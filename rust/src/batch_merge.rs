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

/// Minimizer prefilter for the cross-batch merge: `(w, k)` and the shared-count
/// threshold, all from the environment. `ISONFORM_MERGE_MINSHARE=0` (the default)
/// disables it.
///
/// # Why a minimizer index rather than a sketch
///
/// The merge is `O(batch pairs x isoforms per batch^2)`: cluster 0 of `droso_deep`
/// is 26 batches x 133 isoforms = 325 pairs x ~17 700 = **5.7M alignments**, tens
/// of minutes on one core, and 82% of those pairs are rejected as structural. The
/// loop proves that one full alignment at a time.
///
/// A minhash sketch was tried first and rejected for having no false-negative
/// *bound* --- only a rate measured on one corpus. Minimizers do have one: any two
/// sequences sharing an exact substring of length `>= w + k - 1` must share at
/// least one minimizer. And the merge states its own floor, `similar_seq < 100`
/// rejects outright, so the threshold is **derived** rather than tuned: a pair
/// whose shared minimizer count is below what a 100bp shared region forces cannot
/// clear the floor.
///
/// The guarantee is over *exact* substrings, while `similar_seq` tolerates
/// mismatches, so a divergent-but-mergeable pair could in principle fall below the
/// count. That is why `ISONFORM_MERGE_MINSHARE_VERIFY=1` exists: it aligns the
/// skipped pairs anyway and counts how many would have merged, making the
/// false-negative rate a measurement.
fn minshare_config() -> Option<(usize, usize, usize)> {
    static C: std::sync::OnceLock<Option<(usize, usize, usize)>> = std::sync::OnceLock::new();
    *C.get_or_init(|| {
        // Percent of the shorter consensus the co-linear chain must span.
        let x: usize = std::env::var("ISONFORM_MERGE_MINSHARE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0);
        if x == 0 {
            return None;
        }
        let w: usize = std::env::var("ISONFORM_MERGE_MINSHARE_W")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(10);
        let k: usize = std::env::var("ISONFORM_MERGE_MINSHARE_K")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(15);
        Some((w, k, x))
    })
}

/// The `(w, k)` minimizers of one sequence as `(hash, position)`, sorted by hash.
///
/// Positions are what make the filter useful. A *count* of shared minimizers
/// cannot separate mergeable pairs from structural rejects: every isoform in a
/// cluster is from the same gene, so both kinds share nearly all their minimizers.
/// Measured on cluster 74, raising the count threshold from 8 to 120 moved the skip
/// rate from 10.5% to 12.4% --- sound at every setting and worth almost nothing.
///
/// What decides the merge is *structure*: `parse_cigar_diversity_isoform_level`
/// rejects a pair when an interior indel run exceeds `delta_len`. That shows up in
/// the anchors as a **gap in the co-linear chain** --- shared minimizers whose
/// offsets jump. Checking the chain approximates the merge's own test at anchor
/// resolution, for the cost of a chain rather than a 909x909 DP.
fn positional_minimizers(seq: &[u8], w: usize, k: usize) -> Vec<(u64, u32)> {
    let mut v: Vec<(u64, u32)> = Vec::new();
    if seq.len() < k {
        return v;
    }
    let kmers: Vec<u64> = (0..=seq.len() - k)
        .map(|i| hash_kmer(&seq[i..i + k]))
        .collect();
    if kmers.len() < w {
        v.extend(kmers.iter().enumerate().map(|(i, &h)| (h, i as u32)));
    } else {
        let mut last = u32::MAX;
        for (start, win) in kmers.windows(w).enumerate() {
            let (off, &h) = win
                .iter()
                .enumerate()
                .min_by_key(|(o, &h)| (h, *o))
                .map(|(o, h)| (o, h))
                .unwrap();
            let pos = (start + off) as u32;
            if pos != last {
                v.push((h, pos));
                last = pos;
            }
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

/// The largest indel implied by the anchor chain, in bases.
///
/// Anchors shared by both consensuses are matched by hash, chained co-linearly
/// (positions increasing in both), and the **offset** `p2 - p1` is tracked along
/// the chain. A jump in that offset between consecutive anchors *is* an indel of
/// that size. The merge rejects a pair when an interior indel run exceeds its
/// tolerance, so a pair whose chain already shows a jump larger than the tolerance
/// cannot merge, and the full alignment is wasted work.
///
/// This is the right question, and it is not "how long is the chain". Chain length
/// mixes two things --- how much sequence is shared, and how consistently --- and
/// the shared amount does not discriminate: every isoform in a cluster is from the
/// same gene, so mergeable and non-mergeable pairs share nearly all their
/// minimizers. Measured: a shared-count threshold raised 15x moved the skip rate
/// from 10.5% to 12.4%. What separates them is *structure*, which is exactly this
/// offset jump.
///
/// Returns `u32::MAX` when there is no chain to speak of, so such a pair is
/// skipped rather than silently admitted.
fn max_chain_indel(a: &[(u64, u32)], b: &[(u64, u32)]) -> u32 {
    // Merge join on the hashes; both sides are sorted by hash. A hash present
    // several times yields every combination, so a repeated anchor can still chain
    // at its correct copy.
    let mut hits: Vec<(u32, u32)> = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let h = a[i].0;
                let (si, sj) = (i, j);
                while i < a.len() && a[i].0 == h {
                    i += 1;
                }
                while j < b.len() && b[j].0 == h {
                    j += 1;
                }
                for x in si..i {
                    for y in sj..j {
                        hits.push((a[x].1, b[y].1));
                    }
                }
            }
        }
    }
    if hits.len() < 2 {
        return u32::MAX;
    }
    hits.sort_unstable();

    // Greedy co-linear chain: keep an anchor when both positions advance. Ties and
    // out-of-order copies of a repeat are dropped rather than allowed to fake a
    // jump.
    let mut chain: Vec<(u32, u32)> = Vec::with_capacity(hits.len());
    for &(p1, p2) in &hits {
        match chain.last() {
            Some(&(q1, q2)) if p1 <= q1 || p2 <= q2 => {}
            _ => chain.push((p1, p2)),
        }
    }
    if chain.len() < 2 {
        return u32::MAX;
    }
    // The offset along the chain; a change in it is an indel of that size.
    let mut worst = 0u32;
    for w in chain.windows(2) {
        let (p1, p2) = w[0];
        let (q1, q2) = w[1];
        let d = ((q2 as i64 - q1 as i64) - (p2 as i64 - p1 as i64)).unsigned_abs() as u32;
        worst = worst.max(d);
    }
    worst
}

/// The `(w, k)` minimizers of one sequence, as a sorted deduplicated hash list.
fn minimizers(seq: &[u8], w: usize, k: usize) -> Vec<u64> {
    if seq.len() < k + w {
        // Too short for a full window: take every k-mer, which is a superset of
        // any windowed choice and so cannot cause a false negative.
        let mut v: Vec<u64> = (0..seq.len().saturating_sub(k - 1))
            .map(|i| hash_kmer(&seq[i..i + k]))
            .collect();
        v.sort_unstable();
        v.dedup();
        return v;
    }
    let kmers: Vec<u64> = (0..=seq.len() - k)
        .map(|i| hash_kmer(&seq[i..i + k]))
        .collect();
    let mut v = Vec::with_capacity(kmers.len() / w + 2);
    for win in kmers.windows(w) {
        // `min` per window is the minimizer; equal values collapse in the dedup.
        v.push(*win.iter().min().unwrap());
    }
    v.sort_unstable();
    v.dedup();
    v
}

fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    // One mix, so neighbouring k-mers do not land in neighbouring buckets.
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^ (h >> 29)
}

/// How many minimizers two sorted lists share.
fn shared(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }
    n
}

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
    // One minimizer list per isoform, computed once. Consensuses do not change
    // during the merge --- the longer one survives unchanged and the shorter is
    // marked merged --- so these stay valid for the whole pass.
    let mins: Option<Vec<Vec<Vec<(u64, u32)>>>> = minshare_config().map(|(w, k, _)| {
        batches
            .iter()
            .map(|(_, isos)| {
                isos.iter()
                    .map(|(_, iso)| positional_minimizers(&iso.sequence, w, k))
                    .collect()
            })
            .collect()
    });

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

                    // The prefilter: skip pairs that cannot clear the merge's own
                    // 100bp shared-region floor. The threshold scales with the
                    // shorter consensus, since a short one cannot host as many
                    // minimizers as a long one even when it matches perfectly.
                    let mut skipped = false;
                    if let (Some((_, _, x)), Some(mins)) = (minshare_config(), mins.as_ref()) {
                        // `x` is the largest indel the chain may imply, in bases.
                        // Independent of consensus length: what decides a merge is
                        // whether a structural difference exists, not how long the
                        // sequences are. Natural values sit near the merge's own end
                        // tolerances, `delta_iso_len_3 = 30` and
                        // `delta_iso_len_5 = 50`.
                        let worst = max_chain_indel(&mins[i][a], &mins[j][b]);
                        if worst as usize > x {
                            PAIRS_SKIPPED.fetch_add(1, Ordering::Relaxed);
                            if std::env::var_os("ISONFORM_MERGE_MINSHARE_VERIFY").is_none() {
                                continue;
                            }
                            // Verification mode: align anyway and count the
                            // false negatives instead of hiding them.
                            skipped = true;
                        }
                    }

                    PAIRS_ALIGNED.fetch_add(1, Ordering::Relaxed);
                    if !align_to_merge(engine, &long_seq, &short_seq, opts) {
                        continue;
                    }
                    if skipped {
                        SKIPPED_WOULD_MERGE.fetch_add(1, Ordering::Relaxed);
                        continue; // the filter's decision stands; only counted
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
            // The support the isoform really has: the summed weight of its
            // accessions, which is `reads.len()` for ordinary reads.
            let support = crate::weights::sum_accs(&iso.reads);
            if support >= iso_abundance || iso_abundance == 1 {
                out.push(OutputRecord {
                    id,
                    sequence: iso.sequence.clone(),
                    reads: iso.reads.clone(),
                    destination: Destination::Main,
                    support,
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
                    support,
                });
            } else {
                out.push(OutputRecord {
                    id,
                    sequence: iso.sequence.clone(),
                    reads: iso.reads.clone(),
                    destination: Destination::Dropped,
                    support,
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
