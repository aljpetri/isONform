//! The per-cluster pipeline, joining the stages `main` runs in order.
//!
//! ```text
//! fastq -> polyA trim -> batch -> [ intervals -> graph -> simplify -> isoforms ]
//! ```
//!
//! Everything in brackets happens per batch. This module is the wiring; each
//! stage is verified separately by its own oracle, and what is new here is the
//! order, the per-batch parameters, and the outputs each batch leaves behind.
//!
//! # What a batch leaves behind, and why the format differs from the reference
//!
//! The reference writes four things per batch: `{id}batch`, `spoa{id}` and
//! `mapping{id}` as **Python pickles**, and `skip{id}.fa` as text.
//! `isONform_parallel` spawns `main` as a subprocess and later reads those
//! pickles back.
//!
//! This port keeps the same *content* and does not write pickles. They are an
//! internal handoff between two isONform processes, never user-facing output, and
//! reimplementing the pickle protocol would add risk without adding anything a
//! user can see. The consequence is stated rather than buried: **the Rust and
//! Python halves cannot be mixed.** Running Python `isONform_parallel` over
//! output from Rust `main`, or the reverse, will fail to read the intermediates.
//! Each implementation is self-consistent end to end.
//!
//! The user-visible outputs --- `cluster*_merged.fa`/`.fq`, `cluster*_mapping.txt`,
//! `support_*.txt` --- are byte-compatible, and that is what the end-to-end check
//! compares.

use rustc_hash::FxHashMap;

use crate::graph::Interval;
use crate::graph_build::{self, BuildInput, BuildOpts};
use crate::isoforms::{self, MergeOpts, SpoaParasailMerge};
use crate::wis::WisOpts;

/// Everything one batch produces.
#[derive(Debug, Default, Clone)]
pub struct BatchOutput {
    pub batch_id: usize,
    /// The reads the batch saw: `(accession, sequence)`, in order. Written as
    /// `{batch_id}batch` by the reference.
    pub reads: Vec<(String, Vec<u8>)>,
    /// Reads with no usable interval, in the order the reference writes them to
    /// `skip{batch_id}.fa`.
    pub skipped: Vec<(String, Vec<u8>)>,
    /// `spoa{batch_id}`: isoform id -> consensus.
    pub isoforms: Vec<(u32, Vec<u8>)>,
    /// `mapping{batch_id}`: isoform id -> supporting read accessions.
    pub mapping: Vec<(u32, Vec<String>)>,
}

/// The knobs `main` passes down, after its own defaulting.
#[derive(Debug, Clone, Copy)]
pub struct RunParams {
    pub k: usize,
    pub w: usize,
    pub xmin: usize,
    pub xmax: usize,
    pub delta_len: i64,
    pub max_seqs: usize,
    pub max_seqs_to_spoa: usize,
    pub set_w_dynamically: bool,
    pub slow: bool,
    pub delta: f64,
    pub delta_iso_len_3: usize,
    pub delta_iso_len_5: usize,
    pub wis: WisOpts,
    pub build: BuildOpts,
}

/// `w` for one batch, reproducing `main:394-404`.
///
/// A batch of exactly one read is not given a window at all --- the reference
/// writes it straight to the skip file and `continue`s, which is why this returns
/// `None` there rather than a number the caller would have to know to ignore.
pub fn batch_window(params: &RunParams, n_reads: usize) -> Option<usize> {
    if !params.set_w_dynamically {
        return Some(params.w);
    }
    if n_reads == 1 {
        return None;
    }
    if n_reads >= 100 {
        // `min(args.w, args.k + (len(reads) // 100 + 4))`
        Some(params.w.min(params.k + (n_reads / 100 + 4)))
    } else {
        // `args.k + 1 + len(reads) // 30`
        Some(params.k + 1 + n_reads / 30)
    }
}

/// One batch, end to end.
///
/// `reads` are `(r_id, accession, sequence)` with `r_id` the 1-based id the
/// reference assigns across the *whole* file, not within the batch.
pub fn run_batch(
    batch_id: usize,
    reads: &[(u32, String, Vec<u8>)],
    all_read_len: &[(u32, u32)],
    params: &RunParams,
) -> BatchOutput {
    let mut out = BatchOutput {
        batch_id,
        reads: reads
            .iter()
            .map(|(_, a, s)| (a.clone(), s.clone()))
            .collect(),
        ..Default::default()
    };

    let Some(w) = batch_window(params, reads.len()) else {
        // "Not enough reads to work on!" --- the single read is skipped whole.
        out.skipped = reads
            .iter()
            .map(|(_, a, s)| (a.clone(), s.clone()))
            .collect();
        return out;
    };

    // -- intervals ---------------------------------------------------------
    let seqs: Vec<(u32, Vec<u8>)> = reads.iter().map(|(r, _, s)| (*r, s.clone())).collect();
    let iv = crate::intervals::build_batch(
        &seqs,
        w,
        params.k,
        params.xmin,
        params.xmax,
        params.delta_len,
        params.wis,
    );

    let acc_of: FxHashMap<u32, &String> = reads.iter().map(|(r, a, _)| (*r, a)).collect();
    for r_id in &iv.skipped {
        if let Some(a) = acc_of.get(r_id) {
            let seq = reads
                .iter()
                .find(|(r, _, _)| r == r_id)
                .map(|(_, _, s)| s.clone())
                .unwrap_or_default();
            out.skipped.push(((*a).clone(), seq));
        }
    }
    if iv.by_graph_id.is_empty() {
        return out;
    }

    // -- graph -------------------------------------------------------------
    // Ordered by graph id, which is the order the reference's dict holds them
    // and therefore the order nodes are created.
    let mut gids: Vec<u32> = iv.by_graph_id.keys().copied().collect();
    gids.sort_unstable();
    let intervals: Vec<(u32, Vec<Interval>)> = gids
        .iter()
        .map(|g| {
            let v = iv.by_graph_id[g]
                .iter()
                .map(|c| Interval {
                    start: c.start as u32,
                    end: c.stop as u32,
                    weight: c.support as u32,
                    occurrences: c
                        .instance
                        .iter()
                        .flat_map(|(r, a, b)| [*r, *a, *b])
                        .collect(),
                })
                .collect();
            (*g, v)
        })
        .collect();

    // `new_all_reads`: the graph's reads, keyed by graph id.
    let graph_reads: FxHashMap<u32, Vec<u8>> = gids
        .iter()
        .filter_map(|g| {
            let r = iv.read_of_graph_id.get(g)?;
            let s = reads.iter().find(|(rr, _, _)| rr == r)?;
            Some((*g, s.2.clone()))
        })
        .collect();

    // `read_len_dict = get_read_lengths(all_reads)` --- built over the **whole
    // file**, once, outside the batch loop, and *not* re-derived per batch. So it
    // is passed in rather than computed here. See `BuildInput::read_len` and the
    // note on this function.

    let input = BuildInput {
        k: params.k,
        delta_len: params.delta_len,
        intervals: &intervals,
        reads: &graph_reads,
        read_len: all_read_len,
    };
    let Ok((mut g, reads_for_isoforms)) =
        graph_build::generate_graph_from_intervals(&input, params.build)
    else {
        // `BuildError` reproduces a point where the reference raises. Leaving the
        // batch empty is the only honest option; the caller reports it.
        return out;
    };

    // -- simplification ----------------------------------------------------
    crate::simplify::simplify_graph(
        &mut g,
        &graph_reads,
        params.k,
        params.delta_len,
        params.slow,
    );

    // -- isoforms ----------------------------------------------------------
    let mut groups = isoforms::compute_equal_reads(&g, &reads_for_isoforms);
    let mut engine = SpoaParasailMerge;
    let opts = MergeOpts {
        delta: params.delta,
        delta_len: params.delta_len as usize,
        delta_iso_len_3: params.delta_iso_len_3,
        delta_iso_len_5: params.delta_iso_len_5,
        max_seqs_to_spoa: params.max_seqs_to_spoa,
    };
    let consensuses = isoforms::merge_consensuses(&mut engine, &mut groups, &graph_reads, opts);

    // `prepare_consensuses` keeps only the ids still in `equal_reads`, so a group
    // merged away contributes nothing.
    let live: Vec<u32> = groups.iter().map(|(k, _)| *k).collect();
    for id in &live {
        if let Some((_, seq)) = consensuses.iter().find(|(k, _)| k == id) {
            out.isoforms.push((*id, seq.clone()));
        }
    }
    // `mapping[key]` is the *accessions* of the group's reads, looked up through
    // the graph id, not the read id.
    for (id, members) in &groups {
        let accs = members
            .iter()
            .filter_map(|gid| {
                let r = iv.read_of_graph_id.get(gid)?;
                acc_of.get(r).map(|a| (*a).clone())
            })
            .collect();
        out.mapping.push((*id, accs));
    }
    out
}

/// The whole of `main`'s per-cluster work: read the file, trim, batch, run.
///
/// `reads` are already parsed, so the caller owns file IO and this stays testable.
pub fn run_cluster(records: &[crate::fastq::Record], params: &RunParams) -> Vec<BatchOutput> {
    // `{i + 1: (acc, remove_read_polyA_ends(seq, 12, 1), qual)}` --- ids are
    // 1-based over the whole file and survive batching.
    let all: Vec<(u32, String, Vec<u8>)> = records
        .iter()
        .enumerate()
        .map(|(i, r)| {
            (
                i as u32 + 1,
                r.acc.clone(),
                crate::reads::remove_polya_ends(&r.seq, 12, 1),
            )
        })
        .collect();

    // `read_len_dict = get_read_lengths(all_reads)`: the whole file, keyed by
    // global read id, computed once and reused for every batch.
    let all_read_len: Vec<(u32, u32)> = all.iter().map(|(r, _, s)| (*r, s.len() as u32)).collect();

    crate::reads::batch(&all, params.max_seqs)
        .into_iter()
        .enumerate()
        .map(|(batch_id, b)| run_batch(batch_id, &b, &all_read_len, params))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> RunParams {
        RunParams {
            k: 9,
            w: 20,
            xmin: 14,
            xmax: 80,
            delta_len: 5,
            max_seqs: 1000,
            max_seqs_to_spoa: 200,
            set_w_dynamically: false,
            slow: false,
            delta: 0.1,
            delta_iso_len_3: 30,
            delta_iso_len_5: 50,
            wis: WisOpts::default(),
            build: BuildOpts::default(),
        }
    }

    #[test]
    fn a_fixed_window_ignores_the_read_count() {
        let p = params();
        for n in [1, 5, 99, 100, 5000] {
            assert_eq!(batch_window(&p, n), Some(20));
        }
    }

    #[test]
    fn a_dynamic_window_has_three_regimes() {
        let mut p = params();
        p.set_w_dynamically = true;
        // One read: no window at all --- the reference skips the batch entirely.
        assert_eq!(batch_window(&p, 1), None);
        // Under 100: `k + 1 + n // 30`.
        assert_eq!(batch_window(&p, 60), Some(9 + 1 + 2));
        // 100 or more: `min(w, k + (n // 100 + 4))`. At 100 reads that is
        // `min(20, 9 + 5)` --- the formula wins.
        assert_eq!(batch_window(&p, 100), Some(14));
        assert_eq!(batch_window(&p, 100_000), Some(20), "capped by --w");
    }

    #[test]
    fn a_batch_of_one_read_is_skipped_whole() {
        let mut p = params();
        p.set_w_dynamically = true;
        let reads = vec![(1u32, "r1".to_string(), b"ACGTACGTAC".to_vec())];
        let read_len: Vec<(u32, u32)> =
            reads.iter().map(|(r, _, s)| (*r, s.len() as u32)).collect();
        let out = run_batch(0, &reads, &read_len, &p);
        assert_eq!(out.skipped.len(), 1);
        assert!(out.isoforms.is_empty());
        assert_eq!(out.reads.len(), 1, "it is still recorded in the batch file");
    }

    #[test]
    fn reads_sharing_a_route_produce_an_isoform_and_a_mapping() {
        // Varied rather than a repeated motif: a tandem repeat gives every
        // window the same minimizer, so no anchor pair spans anything and the
        // interval stage legitimately returns nothing.
        let mut seed = 42u64;
        let seq: String = (0..600)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ["A", "C", "G", "T"][(seed >> 33) as usize % 4]
            })
            .collect();
        let records: Vec<crate::fastq::Record> = (1..=4)
            .map(|i| crate::fastq::Record {
                acc: format!("read{i}"),
                seq: seq.as_bytes().to_vec(),
                qual: None,
            })
            .collect();
        let out = run_cluster(&records, &params());
        assert_eq!(out.len(), 1, "four reads fit in one batch");
        let b = &out[0];
        assert!(!b.isoforms.is_empty(), "identical reads yield an isoform");
        assert_eq!(
            b.mapping.len(),
            b.isoforms.len(),
            "one mapping per surviving isoform"
        );
        // Mapping carries accessions, not ids.
        assert!(b
            .mapping
            .iter()
            .all(|(_, a)| a.iter().all(|x| x.starts_with("read"))));
    }

    #[test]
    fn read_lengths_span_the_whole_file_not_the_batch() {
        // `read_len_dict = get_read_lengths(all_reads)` is computed **once**,
        // outside the batch loop, over every read in the file --- and
        // `generateGraphfromIntervals` then walks it as `range(1,
        // len(read_len_dict) + 1)`, indexing it by *graph* id. Graph ids restart
        // at 1 in every batch, so from the second batch on the sink node is
        // seeded with the lengths of the file's *first* reads rather than the
        // batch's. That is the reference's behaviour and this port reproduces it.
        //
        // The first version of this driver derived the lengths from the batch
        // instead. Batch 0 was unaffected --- its read ids and graph ids overlap
        // --- so every single-batch corpus passed, and only a cluster split into
        // several batches exposed it: against the reference, batch 0 agreed and
        // batches 1-3 collapsed to roughly one isoform per read. This pins the
        // parameter so it cannot quietly go back to being per-batch.
        let mut seed = 7u64;
        let seq: String = (0..600)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ["A", "C", "G", "T"][(seed >> 33) as usize % 4]
            })
            .collect();
        let records: Vec<crate::fastq::Record> = (1..=8)
            .map(|i| crate::fastq::Record {
                acc: format!("read{i}"),
                seq: seq.as_bytes().to_vec(),
                qual: None,
            })
            .collect();
        let mut p = params();
        p.max_seqs = 4;

        let all: Vec<(u32, String, Vec<u8>)> = records
            .iter()
            .enumerate()
            .map(|(i, r)| {
                (
                    i as u32 + 1,
                    r.acc.clone(),
                    crate::reads::remove_polya_ends(&r.seq, 12, 1),
                )
            })
            .collect();
        let whole_file: Vec<(u32, u32)> =
            all.iter().map(|(r, _, s)| (*r, s.len() as u32)).collect();
        assert_eq!(whole_file.len(), 8, "one entry per read in the file");

        let batches = crate::reads::batch(&all, p.max_seqs);
        assert_eq!(batches.len(), 2, "eight reads at four per batch");
        // The second batch's own ids start at 5, so a per-batch table would not
        // cover graph ids 1..4 at all.
        assert_eq!(batches[1][0].0, 5);

        // What `run_cluster` does must be what `run_batch` does with the
        // whole-file table --- not with the batch's own.
        let via_cluster = run_cluster(&records, &p);
        let direct = run_batch(1, &batches[1], &whole_file, &p);
        assert_eq!(via_cluster.len(), 2);
        assert_eq!(
            via_cluster[1].isoforms, direct.isoforms,
            "run_cluster passes the whole file's read lengths to every batch"
        );
        assert_eq!(via_cluster[1].mapping, direct.mapping);
    }

    #[test]
    fn read_ids_are_one_based_over_the_whole_file_not_per_batch() {
        // `enumerate(readfq(...))` runs before `batch`, so the second batch's
        // first read is not id 1. The graph and the interval stage both key on
        // these ids, so restarting them per batch would silently mis-associate.
        let mut p = params();
        p.max_seqs = 2;
        let records: Vec<crate::fastq::Record> = (1..=5)
            .map(|i| crate::fastq::Record {
                acc: format!("r{i}"),
                seq: b"ACGTACGTACGTACGTACGT".to_vec(),
                qual: None,
            })
            .collect();
        let out = run_cluster(&records, &p);
        assert_eq!(out.len(), 3, "5 reads at max_seqs 2");
        assert_eq!(out[0].reads.len(), 2);
        assert_eq!(out[2].reads.len(), 1, "a lone trailing read");
        assert_eq!(
            out[1].reads[0].0, "r3",
            "second batch starts at the third read"
        );
    }

    #[test]
    fn the_polya_trim_runs_before_batching() {
        // `all_reads` is built with the trim applied, so what the batch records
        // is the trimmed sequence --- which is what the graph and the skip file
        // both see.
        let long_tail = format!("ACGTACGTAC{}", "A".repeat(60));
        let records = vec![crate::fastq::Record {
            acc: "r1".into(),
            seq: long_tail.as_bytes().to_vec(),
            qual: None,
        }];
        let out = run_cluster(&records, &params());
        assert!(
            out[0].reads[0].1.len() < long_tail.len(),
            "the recorded read is the trimmed one"
        );
    }
}
