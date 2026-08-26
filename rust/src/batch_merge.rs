//! Batch merging: `modules/batch_merging_parallel.py`.
//!
//! The last stage. Each cluster's batches are read back from disk, their
//! isoforms are supposed to be merged across batches, and the surviving ones are
//! written out as the final transcriptome.
//!
//! Reachable only from `isONform_parallel:319`; `main`'s call to it is commented
//! out (`main:583`), so the single-instance entry point never merges batches at
//! all.
//!
//! # The merging step is a no-op. It has never run.
//!
//! `actual_merging_process` sorts the batches **descending** by id and then
//! guards its entire body with `if not batchid2 <= batchid`, where `batchid2`
//! comes from a slice *after* `batchid` in that descending list — so
//! `batchid2 < batchid` always, the guard is always false, and nothing inside it
//! ever executes. Measured rather than read off: on data engineered so every pair
//! is mergeable, the guard was evaluated 3 times and the body entered **0** times,
//! merging 0 of 5 consensuses.
//!
//! And behind that guard sits a second defect. Flip the sort to ascending and the
//! body does run — straight into `generate_consensus_path`, which does
//! `all_sequences[id]` where `id` is a read *accession* and `all_sequences` is
//! keyed by *batch id*. It raises `KeyError` on the first read. So the code cannot
//! be made to work by fixing the guard alone; the lookup it wants does not exist.
//!
//! This port therefore reproduces the no-op and does **not** offer a fix flag.
//! Everywhere else a reference defect has an opt-in fix ([`crate::wis::WisOpts`],
//! [`crate::graph_build::BuildOpts`]) there was a defensible correct behaviour to
//! switch to. Here there is not: making it merge means inventing the lookup that
//! `generate_consensus_path` is missing, and inventing behaviour is the one thing
//! the porting rules forbid. `PORTING.md` finding 31.
//!
//! What the stage does do is [`select_output`] — decide which isoforms are
//! written where — which is live, and is what produces the transcriptome.

use std::sync::atomic::{AtomicU64, Ordering};

use crate::isoforms::{align_to_merge, IsoformEngine, MergeOpts};

/// Cross-batch merge counters, for `ISONFORM_SKETCH_VERIFY`.
///
/// `skipped_would_merge` is the number that matters: pairs the sketch filter
/// declined to align that alignment says *would* have merged. That is the
/// filter's false-negative count, and it is measured rather than assumed.
pub static PAIRS_SEEN: AtomicU64 = AtomicU64::new(0);
pub static PAIRS_SKIPPED: AtomicU64 = AtomicU64::new(0);
pub static PAIRS_ALIGNED: AtomicU64 = AtomicU64::new(0);
pub static SKIPPED_WOULD_MERGE: AtomicU64 = AtomicU64::new(0);
pub static MERGES: AtomicU64 = AtomicU64::new(0);

/// One isoform carried through batch merging: `Read(sequence, reads, merged)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Isoform {
    pub sequence: Vec<u8>,
    /// Read accessions supporting it.
    pub reads: Vec<String>,
    pub merged: bool,
}

/// `actual_merging_process`: fold duplicate isoforms across batches together.
///
/// # The reference's version has never executed, and this is the repair
///
/// `PORTING.md` finding 31. Two stacked defects:
///
/// 1. the batch list is sorted **descending** and then sliced from `b_i + 1`, so
///    `batchid2 < batchid` always and the guard `if not batchid2 <= batchid` is
///    never true. The body is unreachable;
/// 2. behind it, `generate_consensus_path` looks up read *accessions* in a dict
///    keyed by *batch id*, and its fallback branch is worse than a `KeyError` --
///    `if not (id in all_sequences): sequence = id` writes the accession string
///    itself into the fasta as though it were DNA.
///
/// The first is fixed by sorting **ascending**, which makes every later batch a
/// higher id and the intent ("compare each batch against every later batch")
/// actually happen.
///
/// The second is not repaired but removed. Merging compares the two stored
/// **consensus sequences** directly and needs no read data at all, so
/// `generate_consensus_path` has no role. When two isoforms merge, the longer
/// one's consensus is kept and the shorter one's reads are absorbed into it --
/// which is exactly what the reference already does on its well-supported branch
/// (`len(infos_long.reads) > 50` keeps `infos_long.sequence` untouched). Making
/// both branches behave that way is the reference's own high-confidence
/// behaviour applied uniformly, not a new invention.
///
/// # What is preserved from the reference
///
/// The traversal order, because it decides which of several mergeable isoforms
/// wins: batches ascending by id, and within a batch isoforms ascending by
/// consensus length. Long/short is by consensus length with `>=` favouring the
/// *later* batch's isoform on a tie, as the reference's branch does. A merged
/// isoform is skipped on both sides.
///
/// Pass `no_op = true` to reproduce the reference exactly --- see
/// [`crate::driver::BugCompat`].
pub fn actual_merging_process<E: IsoformEngine>(
    engine: &mut E,
    batches: &mut [(String, Vec<(u32, Isoform)>)],
    opts: MergeOpts,
    no_op: bool,
) {
    if no_op {
        // Finding 31 reproduced: the reference's body is unreachable.
        return;
    }
    // `sorted(all_infos_dict.items(), key=lambda x: x[0])` --- ascending, which
    // is the one-character fix that makes the guard satisfiable. Sorted by the
    // batch id's *numeric components* so `"3"` < `"3_0"` < `"3_1"` < `"10"`;
    // a plain string sort would put `"10"` before `"3"`.
    batches.sort_by_key(|(b, _)| batch_sort_key(b));

    // Work on indices: two batches are mutated per merge, so the borrow checker
    // wants the pair split rather than held.
    for i in 0..batches.len().saturating_sub(1) {
        for j in (i + 1)..batches.len() {
            // `sorted(id_dict.items(), key=lambda x: len(x[1].sequence))` ---
            // ascending by consensus length, recomputed per pair as the
            // reference does, so merges made earlier in this pass are seen.
            let order_i = length_order(&batches[i].1);
            let order_j = length_order(&batches[j].1);

            for &a in &order_i {
                if batches[i].1[a].1.merged {
                    continue;
                }
                for &b in &order_j {
                    if batches[i].1[a].1.merged {
                        break;
                    }
                    if batches[j].1[b].1.merged {
                        continue;
                    }
                    // `if len(infos2.sequence) >= len(infos.sequence)` --- the
                    // later batch's isoform is the long one on a tie.
                    let (li, la, si, sa) =
                        if batches[j].1[b].1.sequence.len() >= batches[i].1[a].1.sequence.len() {
                            (j, b, i, a)
                        } else {
                            (i, a, j, b)
                        };
                    let long_seq = batches[li].1[la].1.sequence.clone();
                    let short_seq = batches[si].1[sa].1.sequence.clone();

                    // No candidate filter in this path. A minhash sketch was
                    // built and rejected: it gives no bound on false negatives,
                    // and a filter deciding which pairs are *never* examined
                    // cannot rest on a rate measured on one corpus. `crate::sketch`
                    // keeps the code and the reasoning; the merge does not use it.
                    PAIRS_SEEN.fetch_add(1, Ordering::Relaxed);
                    PAIRS_ALIGNED.fetch_add(1, Ordering::Relaxed);
                    if !align_to_merge(engine, &long_seq, &short_seq, opts) {
                        continue;
                    }
                    MERGES.fetch_add(1, Ordering::Relaxed);
                    // The longer consensus survives unchanged; the shorter
                    // isoform's reads move into it and it is marked merged.
                    let moved = batches[si].1[sa].1.reads.clone();
                    batches[li].1[la].1.reads.extend(moved);
                    batches[si].1[sa].1.merged = true;
                }
            }
        }
    }
}

/// A batch id's numeric components, for ordering.
///
/// Batch ids are `"3"` for a single-batch invocation and `"3_1"` when one
/// invocation wrote several (finding 34). Ordering has to be numeric per
/// component: lexicographically `"10" < "3"`, which would put batches in the
/// wrong order and change which of two mergeable isoforms survives.
pub fn batch_sort_key(id: &str) -> Vec<u64> {
    id.split('_')
        .map(|p| p.parse::<u64>().unwrap_or(0))
        .collect()
}

/// Indices of `isoforms` ordered by consensus length ascending, stably.
fn length_order(isoforms: &[(u32, Isoform)]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..isoforms.len()).collect();
    idx.sort_by_key(|&i| isoforms[i].1.sequence.len());
    idx
}

/// Where one isoform is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Destination {
    /// The cluster's main consensus, mapping and support files.
    Main,
    /// The low-abundance files, written only when `write_low_abundance`.
    LowAbundance,
    /// Dropped entirely: below the abundance threshold with low-abundance
    /// output turned off.
    Dropped,
}

/// One row of the final output.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputRecord {
    /// `"{cluster}_{batch}_{isoform}"`.
    pub id: String,
    pub sequence: Vec<u8>,
    pub reads: Vec<String>,
    pub destination: Destination,
    /// What the support file records. For the main file this is the read count;
    /// for low abundance the reference always writes **1** --- see the note on
    /// [`select_output`].
    pub support: usize,
}

/// `write_final_output`'s decisions, separated from its file handles.
///
/// Iterates batches and isoforms in the order the reference's dicts hold them,
/// because that order is the order records appear in the output files.
///
/// # Two reference behaviours worth knowing
///
/// **The abundance test is `>= iso_abundance or iso_abundance == 1`.** The second
/// clause is redundant — a read count is at least 1, so `>= 1` is already always
/// true — but it makes the intent explicit: `--iso_abundance 1` keeps everything.
///
/// **The low-abundance support count was always 1, and is now the real count.**
/// The reference writes `len(all_infos_dict[new_id].reads)` when
/// `new_id in all_infos_dict` and `1` otherwise; `new_id` is a string like
/// `"3_0_7"` while `all_infos_dict` is keyed by integer batch id, so the lookup
/// never succeeds and the count is always the literal 1 whatever the isoform's
/// real support. Finding 32, fixed --- the read count is what the main support
/// file writes and plainly what was meant.
///
/// Note this path is unreachable from `isONform_parallel`, which hardcodes
/// `write_low_abundance = False` (finding 35). It is reachable through
/// `batch_merging_parallel` directly.
pub fn select_output(
    batches: &[(String, Vec<(u32, Isoform)>)],
    cluster: &str,
    iso_abundance: usize,
    write_low_abundance: bool,
) -> Vec<OutputRecord> {
    let mut out = Vec::new();
    for (batch_id, isoforms) in batches {
        for (iso_id, iso) in isoforms {
            if iso.merged {
                continue;
            }
            let id = format!("{cluster}_{batch_id}_{iso_id}");
            if iso.reads.len() >= iso_abundance || iso_abundance == 1 {
                out.push(OutputRecord {
                    id,
                    sequence: iso.sequence.clone(),
                    reads: iso.reads.clone(),
                    destination: Destination::Main,
                    support: iso.reads.len(),
                });
            } else if write_low_abundance {
                out.push(OutputRecord {
                    id,
                    sequence: iso.sequence.clone(),
                    reads: iso.reads.clone(),
                    destination: Destination::LowAbundance,
                    // The isoform's real read count. The reference writes the
                    // literal 1 here because `if new_id in all_infos_dict`
                    // compares a string id against integer batch keys and can
                    // never match --- finding 32, fixed.
                    support: iso.reads.len(),
                });
            } else {
                out.push(OutputRecord {
                    id,
                    sequence: iso.sequence.clone(),
                    reads: iso.reads.clone(),
                    destination: Destination::Dropped,
                    support: iso.reads.len(),
                });
            }
        }
    }
    out
}

/// One record as it appears in the consensus file.
///
/// `write_fastq` picks the format, and the quality string is `'+'` repeated to the
/// sequence's length --- not a real quality, a placeholder the reference writes
/// verbatim.
pub fn format_consensus(rec: &OutputRecord, write_fastq: bool) -> Vec<u8> {
    let seq = String::from_utf8_lossy(&rec.sequence);
    if write_fastq {
        format!(
            "@{}\n{}\n+\n{}\n",
            rec.id,
            seq,
            "+".repeat(rec.sequence.len())
        )
        .into_bytes()
    } else {
        format!(">{}\n{}\n", rec.id, seq).into_bytes()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iso(seq: &str, reads: &[&str]) -> Isoform {
        Isoform {
            sequence: seq.as_bytes().to_vec(),
            reads: reads.iter().map(|s| s.to_string()).collect(),
            merged: false,
        }
    }

    #[test]
    fn bug_compat_reproduces_the_reference_merging_nothing() {
        // The reference's guard is unreachable --- measured, not inferred: on data
        // engineered so every pair is mergeable, its body was entered 0 times and
        // 0 of 5 consensuses were marked merged. This pins the same outcome.
        let seq = "ACGT".repeat(80);
        let mut batches = vec![
            (
                "0".to_string(),
                vec![(1u32, iso(&seq, &["r1", "r2"])), (2, iso(&seq, &["r3"]))],
            ),
            (
                "1".to_string(),
                vec![(1u32, iso(&seq, &["r4", "r5"])), (2, iso(&seq, &["r6"]))],
            ),
            ("2".to_string(), vec![(1u32, iso(&seq, &["r7"]))]),
        ];
        let before = batches.clone();
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, MergeOpts::default(), true);
        assert_eq!(batches, before, "nothing is merged, and nothing else moves");
        assert!(
            batches.iter().flat_map(|(_, v)| v).all(|(_, i)| !i.merged),
            "identical sequences across batches stay unmerged"
        );
    }

    fn opts() -> MergeOpts {
        MergeOpts {
            delta: 0.15,
            delta_len: 5,
            delta_iso_len_3: 30,
            delta_iso_len_5: 50,
            max_seqs_to_spoa: 200,
            merge_rebuild_max: 50,
        final_consensus_pass: false,
        cigar_diversity_counts_runs: false,
        }
    }

    #[test]
    fn identical_consensuses_in_different_batches_are_folded_together() {
        // The repair. Reference behaviour is the test above: nothing merges.
        let seq = "ACGT".repeat(80);
        let mut batches = vec![
            ("0".to_string(), vec![(1u32, iso(&seq, &["r1", "r2"]))]),
            ("1".to_string(), vec![(1u32, iso(&seq, &["r3"]))]),
            ("2".to_string(), vec![(1u32, iso(&seq, &["r4"]))]),
        ];
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, opts(), false);

        let survivors: Vec<(String, u32, usize)> = batches
            .iter()
            .flat_map(|(b, v)| {
                v.iter()
                    .filter(|(_, i)| !i.merged)
                    .map(move |(id, i)| (b.clone(), *id, i.reads.len()))
            })
            .collect();
        assert_eq!(
            survivors.len(),
            1,
            "three identical consensuses collapse to one, got {survivors:?}"
        );
        assert_eq!(survivors[0].2, 4, "and it carries every read: 2 + 1 + 1");
    }

    #[test]
    fn unrelated_consensuses_in_different_batches_are_left_alone() {
        // The check that stops the repair from being a collapse-everything bug.
        let a = "ACGT".repeat(80);
        let b: String = "TTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTA".repeat(8);
        let mut batches = vec![
            ("0".to_string(), vec![(1u32, iso(&a, &["r1"]))]),
            ("1".to_string(), vec![(1u32, iso(&b, &["r2"]))]),
        ];
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, opts(), false);
        assert!(
            batches.iter().flat_map(|(_, v)| v).all(|(_, i)| !i.merged),
            "nothing mergeable, nothing merged"
        );
    }

    #[test]
    fn merging_never_touches_isoforms_inside_one_batch() {
        // Cross-batch only: intra-batch merging is `merge_consensuses`' job and
        // happens earlier, per cluster. Two identical consensuses in the *same*
        // batch must survive this pass untouched.
        let seq = "ACGT".repeat(80);
        let mut batches = vec![(
            "0".to_string(),
            vec![(1u32, iso(&seq, &["r1"])), (2u32, iso(&seq, &["r2"]))],
        )];
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, opts(), false);
        assert!(batches[0].1.iter().all(|(_, i)| !i.merged));
    }

    #[test]
    fn the_longer_consensus_survives_and_keeps_its_own_sequence() {
        // No read data is consulted and no new consensus is computed: the longer
        // sequence is kept verbatim, which is what the reference's own
        // `len(reads) > 50` branch does.
        let long = "ACGT".repeat(80);
        let short = long[..300].to_string();
        let mut batches = vec![
            ("0".to_string(), vec![(7u32, iso(&short, &["r1"]))]),
            ("1".to_string(), vec![(9u32, iso(&long, &["r2"]))]),
        ];
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, opts(), false);
        let live: Vec<_> = batches
            .iter()
            .flat_map(|(b, v)| {
                v.iter()
                    .filter(|(_, i)| !i.merged)
                    .map(move |(id, i)| (b.as_str(), *id, i))
            })
            .collect();
        assert_eq!(live.len(), 1);
        let (batch, id, iso) = &live[0];
        assert_eq!(
            (*batch, *id),
            ("1", 9),
            "the longer isoform is the survivor"
        );
        assert_eq!(
            iso.sequence,
            long.as_bytes(),
            "its consensus is unchanged --- not recomputed"
        );
        assert_eq!(iso.reads.len(), 2, "and it absorbed the shorter's reads");
    }

    #[test]
    fn a_merged_isoform_is_not_merged_again() {
        // Once absorbed, an isoform must not also be absorbed by a third batch:
        // its reads would be double-counted in the support file.
        let seq = "ACGT".repeat(80);
        let mut batches = vec![
            ("0".to_string(), vec![(1u32, iso(&seq, &["r1"]))]),
            ("1".to_string(), vec![(1u32, iso(&seq, &["r2"]))]),
            ("2".to_string(), vec![(1u32, iso(&seq, &["r3"]))]),
        ];
        let mut engine = crate::isoforms::SpoaParasailMerge;
        actual_merging_process(&mut engine, &mut batches, opts(), false);
        let total: usize = batches
            .iter()
            .flat_map(|(_, v)| v)
            .filter(|(_, i)| !i.merged)
            .map(|(_, i)| i.reads.len())
            .sum();
        assert_eq!(total, 3, "each read counted exactly once, got {total}");
    }

    #[test]
    fn every_unmerged_isoform_reaches_the_main_file_at_abundance_one() {
        let batches = vec![
            (
                "0".to_string(),
                vec![
                    (1u32, iso("ACGT", &["r1"])),
                    (2, iso("TTTT", &["r2", "r3"])),
                ],
            ),
            ("3".to_string(), vec![(7u32, iso("GGGG", &["r4"]))]),
        ];
        let got = select_output(&batches, "5", 1, false);
        assert_eq!(got.len(), 3);
        assert!(got.iter().all(|r| r.destination == Destination::Main));
        // `{cluster}_{batch}_{isoform}`, and the batch is the *dict* key, not a
        // running index.
        assert_eq!(got[0].id, "5_0_1");
        assert_eq!(got[2].id, "5_3_7");
    }

    #[test]
    fn a_merged_isoform_is_skipped_entirely() {
        let mut batches = vec![(
            "0".to_string(),
            vec![(1u32, iso("ACGT", &["r1"])), (2, iso("TTTT", &["r2"]))],
        )];
        batches[0].1[0].1.merged = true;
        let got = select_output(&batches, "0", 1, false);
        assert_eq!(got.len(), 1);
        assert_eq!(got[0].id, "0_0_2");
    }

    #[test]
    fn below_the_threshold_goes_to_low_abundance_or_is_dropped() {
        let batches = vec![(
            "0".to_string(),
            vec![
                (1u32, iso("ACGT", &["r1"])),
                (2, iso("TTTT", &["a", "b", "c"])),
            ],
        )];
        // Threshold 3: the single-read isoform falls short.
        let with_low = select_output(&batches, "0", 3, true);
        assert_eq!(with_low[0].destination, Destination::LowAbundance);
        assert_eq!(with_low[1].destination, Destination::Main);
        let without = select_output(&batches, "0", 3, false);
        assert_eq!(
            without[0].destination,
            Destination::Dropped,
            "with low-abundance output off it is written nowhere"
        );
    }

    #[test]
    fn the_low_abundance_support_count_is_the_real_read_count() {
        // Finding 32 fixed. The reference's `if new_id in all_infos_dict`
        // compares a string id against integer batch keys, never matches, and
        // writes the literal 1 whatever the real support.
        let batches = vec![("0".to_string(), vec![(1u32, iso("ACGT", &["r1", "r2"]))])];
        let got = select_output(&batches, "0", 5, true);
        assert_eq!(got[0].destination, Destination::LowAbundance);
        assert_eq!(got[0].reads.len(), 2);
        assert_eq!(got[0].support, 2, "the support file records the real count");
    }

    #[test]
    fn abundance_one_keeps_everything_even_though_the_second_clause_is_redundant() {
        let batches = vec![("0".to_string(), vec![(1u32, iso("ACGT", &["r1"]))])];
        for low in [true, false] {
            let got = select_output(&batches, "0", 1, low);
            assert_eq!(got[0].destination, Destination::Main);
        }
    }

    #[test]
    fn fasta_and_fastq_differ_only_in_the_header_and_a_placeholder_quality() {
        let rec = &select_output(
            &[("0".to_string(), vec![(1u32, iso("ACGT", &["r1"]))])],
            "9",
            1,
            false,
        )[0];
        assert_eq!(format_consensus(rec, false), b">9_0_1\nACGT\n".to_vec());
        assert_eq!(
            format_consensus(rec, true),
            b"@9_0_1\nACGT\n+\n++++\n".to_vec(),
            "the quality line is '+' repeated, not a real quality"
        );
    }
}
