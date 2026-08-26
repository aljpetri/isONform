//! Isoform generation: `modules/IsoformGeneration.py`'s live path.
//!
//! The simplified graph goes in and isoforms come out, in three steps:
//!
//! 1. [`compute_equal_reads`] walks the graph once per read and groups the reads
//!    that take the same route;
//! 2. [`merge_consensuses`] builds a consensus per group with spoa, then merges
//!    pairs whose consensuses are similar enough, re-running spoa on the union;
//! 3. the caller writes the result out.
//!
//! # Nearly a third of the module is unreachable
//!
//! Of 24 functions in `IsoformGeneration.py`, four are referenced **nowhere in
//! the repository** — `search_last_entries`, `search_first_entries`,
//! `parse_cigar_diversity_isoform_level` (the older sibling of the `_new` one
//! that is live) and `write_transcriptome_single` — and two more are reachable
//! only from commented-out call sites: `generate_isoform_using_spoa` and
//! `generate_isoforms_new`. That is roughly 183 of 631 lines. `PORTING.md`
//! finding 29.
//!
//! Only the live path is ported.

use crate::align::CigarOp;
use crate::graph::{Graph, NodeId, NodeKey};
use rustc_hash::{FxHashMap, FxHashSet};

/// The two external tools this stage needs.
///
/// `spoa` is the same invocation the simplification stage uses. `align_merge` is
/// **not** the same alignment: `align_to_merge` overrides the mismatch penalty to
/// −2 ([`crate::parasail::Scoring::MERGE`]) where bubble popping uses −8, and it
/// needs the *aligned strings* as well as the CIGAR, because the significant-match
/// scan below compares them column by column.
pub trait IsoformEngine {
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8>;
    /// Returns `(cigar_ops, s1_aligned, s2_aligned)`.
    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>);
}

/// The real thing: [`crate::poa`] and [`crate::parasail`].
#[derive(Debug, Default)]
pub struct SpoaParasailMerge;

impl IsoformEngine for SpoaParasailMerge {
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8> {
        crate::poa::consensus(seqs)
            .map(|s| s.into_bytes())
            .unwrap_or_default()
    }

    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let sc = crate::parasail::Scoring::MERGE;
        let aln = crate::wfa::enabled()
            .then(|| crate::wfa::semiglobal(s1, s2, sc))
            .flatten()
            .unwrap_or_else(|| crate::parasail::semiglobal(s1, s2, sc));
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        (aln.ops, a, b)
    }
}

// ===========================================================================
// Step 1: grouping reads that take the same route
// ===========================================================================

/// `is_Sublist`: is `s` a **contiguous** run inside `l`?
///
/// Faithful to the reference, including that an empty `s` is a sublist of
/// anything and that `s == l` short-circuits before the length check.
///
/// The reference can index out of bounds here — `l[i + n]` with `i` near the end
/// of `l` — but only if `s` is longer than the tail, which the `len(s) > len(l)`
/// guard plus the walk make unreachable in Python's bounds. Written with an
/// explicit window check rather than relying on that.
pub fn is_sublist<T: PartialEq>(l: &[T], s: &[T]) -> bool {
    if s.is_empty() || s == l {
        return true;
    }
    if s.len() > l.len() {
        return false;
    }
    l.windows(s.len()).any(|w| w == s)
}

/// `merge_isoform_paths`: fold a group whose member list is a contiguous run of
/// another group's into that other group.
///
/// # This is not doing what its names say
///
/// The parameter is called `isoforms` and its values `vis_nodes_set`, and the
/// docstring says "merges subisoforms into larger isoforms". But
/// `compute_equal_reads` passes `equal_reads`, whose values are **lists of read
/// ids**, not paths — and the path it computes (`visited_nodes`) is assembled and
/// then never read. So the containment test runs over read-id lists whose order
/// comes from Python set iteration, and "one group's reads form a contiguous run
/// inside another's" is not a statement about isoform structure at all.
///
/// Reproduced exactly, because it decides which isoforms are emitted.
/// `PORTING.md` finding 28.
///
/// The reference also deletes from `isoforms` while walking a merge list built
/// beforehand, so a group that is a sublist of two others is deleted twice and
/// raises `KeyError`. Guarded here rather than panicking; see the tests.
pub fn merge_isoform_paths(isoforms: &mut Vec<(u32, Vec<u32>)>) {
    let mut merges: Vec<(u32, u32)> = Vec::new();
    for (i, (id1, a)) in isoforms.iter().enumerate() {
        for (id2, b) in isoforms.iter().skip(i + 1) {
            // The reference's `for id1 ... for id2 ... if id1 < id2` over a dict,
            // which visits every ordered pair once with the smaller id first.
            let (lo, hi) = if id1 < id2 { (id1, id2) } else { (id2, id1) };
            let (lo_v, hi_v) = if id1 < id2 { (a, b) } else { (b, a) };
            if is_sublist(hi_v, lo_v) {
                merges.push((*lo, *hi));
            } else if is_sublist(lo_v, hi_v) {
                merges.push((*hi, *lo));
            }
        }
    }
    for (sub, large) in merges {
        // `isoforms[largeiso].extend(isoforms[subiso]); del isoforms[subiso]`.
        // Both lookups raise if a previous merge already removed the entry; the
        // reference would crash, so skipping is the only non-inventing option and
        // it is recorded rather than silent.
        let Some(si) = isoforms.iter().position(|(k, _)| *k == sub) else {
            continue;
        };
        let taken = isoforms[si].1.clone();
        let Some(li) = isoforms.iter().position(|(k, _)| *k == large) else {
            continue;
        };
        isoforms[li].1.extend(taken);
        isoforms.retain(|(k, _)| *k != sub);
    }
}

/// `compute_equal_reads`: group reads by the route they take through the graph.
///
/// Returns `(representative_read_id, reads_in_group)` in the order the reference
/// inserts them.
///
/// Two reference behaviours are load-bearing and neither is obvious:
///
/// * **`list(current_node_support)[0]` and the member order are CPython set order.**
///   That order names the isoform and decides the order sequences reach spoa,
///   which is order-sensitive. This port uses **ascending** order instead — a
///   deliberate divergence, taken because reproducing CPython's would mean
///   modelling its probing, resize rule and insertion-order propagation through
///   `set()`/`.intersection()`/`-=`, to preserve an order that carries no meaning
///   and is not even seed-dependent (read ids are ints, so `hash(int) == int`).
///   Measured: it fires on 28 of 114 recorded cases and leaves 93% of isoform
///   consensuses byte-identical there. `PORTING.md` finding 28.
/// * **the pop is immediately undone.** `current_node_support = node_support_left`
///   aliases the set, and `.add(read)` puts the popped read straight back —
///   exactly the aliasing `find_paths` has. So the walk starts from the full
///   remaining support, not the support minus the seed.
pub fn compute_equal_reads(g: &Graph, support: &[u32]) -> Vec<(u32, Vec<u32>)> {
    let (Some(s), Some(t)) = (g.lookup(&NodeKey::Source), g.lookup(&NodeKey::Sink)) else {
        return Vec::new();
    };
    let mut left: std::collections::BTreeSet<u32> = support.iter().copied().collect();
    let mut out: Vec<(u32, Vec<u32>)> = Vec::new();

    while let Some(&read) = left.iter().next() {
        // `pop()` then `current_node_support = node_support_left` then `.add(read)`
        // --- the alias means the removal never took effect.
        let mut current: FxHashSet<u32> = left.iter().copied().collect();
        let mut node: NodeId = s;
        while node != t {
            let mut next_found = false;
            for &nb in g.successors(node) {
                if let Some(supp) = g.edge_support(node, nb) {
                    if supp.contains(&read) {
                        let keep: FxHashSet<u32> = supp.iter().copied().collect();
                        current.retain(|r| keep.contains(r));
                        node = nb;
                        next_found = true;
                        break;
                    }
                }
            }
            if !next_found {
                break;
            }
        }
        if current.is_empty() {
            // The reference leaves `node_support_left` unchanged here, so the
            // same read is popped again and the loop never ends. It cannot
            // happen: `read` is in `current` at entry and survives every
            // intersection along its own path. Break rather than hang.
            break;
        }
        let mut group: Vec<u32> = current.iter().copied().collect();
        group.sort_unstable();
        // `id = list(current_node_support)[0]` --- a set-order-dependent
        // representative. Smallest, for the reason in the doc comment.
        let id = group[0];
        for r in &group {
            left.remove(r);
        }
        out.push((id, group));
    }
    merge_isoform_paths(&mut out);
    out
}

// ===========================================================================
// Step 2: consensus and merging
// ===========================================================================

/// `get_overall_alignment_len`: total CIGAR length.
fn overall_alignment_len(ops: &[CigarOp]) -> usize {
    ops.iter().map(|(n, _)| *n).sum()
}

/// `find_first_significant_match`: the first window of `windowsize` alignment
/// columns whose identity **exceeds** `threshold`, or `None`.
pub fn find_first_significant_match(
    a: &[u8],
    b: &[u8],
    windowsize: usize,
    threshold: f64,
) -> Option<usize> {
    // `zip` stops at the shorter, as Python's does.
    let n = a.len().min(b.len());
    if windowsize == 0 || n < windowsize {
        return None;
    }
    let matches: Vec<u32> = (0..n).map(|i| u32::from(a[i] == b[i])).collect();
    let mut run: u32 = matches[..windowsize].iter().sum();
    for i in 0..=(n - windowsize) {
        if i > 0 {
            run = run - matches[i - 1] + matches[i + windowsize - 1];
        }
        if run as f64 / windowsize as f64 > threshold {
            return Some(i);
        }
    }
    None
}

/// Tuning constants for the merge decision, as `main` passes them.
#[derive(Debug, Clone, Copy, Default)]
pub struct MergeOpts {
    pub delta: f64,
    pub delta_len: usize,
    pub delta_iso_len_3: usize,
    pub delta_iso_len_5: usize,
    pub max_seqs_to_spoa: usize,
    /// Reproduce finding 30: accumulate `delta_len` per interior mismatch *run*
    /// instead of the run's actual length.
    ///
    /// This is the one bug-compat switch that is a judgement call rather than a
    /// clear defect. The accumulator is named `miss_match_length` and sits under
    /// a comment saying "we want to add up all missmatches", but adds
    /// `delta_len`; counting runs rather than bases is still a coherent
    /// similarity measure, so fixing it changes the *method*, not just
    /// correctness. Fixed by default under the bug-fix policy, and easy to put
    /// back.
    pub cigar_diversity_counts_runs: bool,
}

/// `parse_cigar_diversity_isoform_level_new`.
///
/// Splits the alignment at the first and last significant match and asks three
/// questions: is anything structural in the middle, do the two ends diverge more
/// than `delta_iso_len_5`/`_3` allows, and is the shared middle similar enough.
///
/// # Every interior mismatch counted as exactly `delta_len` (finding 30, fixed)
///
/// The mismatch accumulator is `miss_match_length += delta_len`, **not**
/// `+= cig_len`, immediately under a comment saying "we want to add up all
/// missmatches to compare to sequence length". A one-base mismatch and a
/// `delta_len`-base one therefore contribute identically, and since anything
/// longer than `delta_len` has already returned `False`, the total is just
/// `delta_len × (number of interior mismatch runs)`. `PORTING.md` finding 28.
#[allow(clippy::too_many_arguments)]
pub fn parse_cigar_diversity_isoform_level(
    ops: &[CigarOp],
    opts: MergeOpts,
    overall_len: usize,
    first_match: usize,
    last_match: usize,
) -> bool {
    let mut miss_match_length = 0usize;
    let mut alignment_len = 0usize;
    let (mut after_last_matches, mut after_last_nomatch) = (0usize, 0usize);
    let (mut before_first_matches, mut before_first_nomatch) = (0usize, 0usize);

    for &(cig_len, cig_type) in ops {
        let this_start_pos = alignment_len;
        alignment_len += cig_len;
        if cig_type != b'=' && cig_type != b'M' {
            if this_start_pos < first_match {
                if cig_type != b'D' {
                    before_first_nomatch += cig_len;
                }
            } else if this_start_pos >= last_match || this_start_pos + cig_len >= last_match {
                if cig_type != b'D' {
                    after_last_nomatch += cig_len;
                }
            } else if cig_len > opts.delta_len {
                // Structural difference: not mergeable.
                bump(|s| {
                    s.structural += 1;
                    match cig_len {
                        0..=10 => s.run_6_10 += 1,
                        11..=20 => s.run_11_20 += 1,
                        21..=50 => s.run_21_50 += 1,
                        51..=100 => s.run_51_100 += 1,
                        101..=200 => s.run_101_200 += 1,
                        _ => s.run_201_plus += 1,
                    }
                    if s.diff_struct_min == 0 || s.last_diff < s.diff_struct_min {
                        s.diff_struct_min = s.last_diff;
                    }
                });
                return false;
            } else if opts.cigar_diversity_counts_runs {
                // finding 30 reproduced: `+= delta_len`, not `+= cig_len`.
                miss_match_length += opts.delta_len;
            } else {
                // finding 30 fixed: the run's actual length, which is what the
                // variable is named for and what the comment above it describes.
                miss_match_length += cig_len;
            }
        } else if this_start_pos < first_match {
            before_first_matches += cig_len;
        } else if this_start_pos >= last_match {
            after_last_matches += cig_len;
        }
    }

    let mergeable_start = before_first_matches + before_first_nomatch < opts.delta_iso_len_5;
    let mergeable_end = after_last_matches + after_last_nomatch < opts.delta_iso_len_3;
    if !mergeable_end || !mergeable_start {
        bump(|s| s.ends_too_long += 1);
        return false;
    }

    let similar_seq = last_match.saturating_sub(first_match);
    // "just to make sure that we only merge reads that have at least 100 nt
    // similar" --- a hard floor, independent of every tuning parameter.
    if similar_seq < 100 {
        bump(|s| s.shared_under_100 += 1);
        return false;
    }
    let diversity = miss_match_length as f64 / similar_seq as f64;
    // Short alignments get twice the tolerance.
    let delta = if overall_len < 200 {
        2.0 * opts.delta
    } else {
        opts.delta
    };
    let max_bp_diff = (delta * similar_seq as f64).max(opts.delta_len as f64);
    let ok = diversity <= max_bp_diff / similar_seq as f64;
    bump(|s| {
        if ok {
            s.merged += 1;
            if s.last_diff > s.diff_merged_max {
                s.diff_merged_max = s.last_diff;
            }
        } else {
            s.too_diverse += 1
        }
    });
    ok
}

/// Why a pair was rejected, for `ISONFORM_MERGE_STATS`.
///
/// Every call to [`align_to_merge`] pays a full O(n*m) alignment before any
/// condition is evaluated, and profiling says that alignment is ~100% of
/// isONform's runtime. A prefilter is only worth building for the condition that
/// actually rejects most pairs, so this counts them rather than guessing.
///
/// Set `ISONFORM_MERGE_STATS=1` and the totals are printed to stderr when the
/// process exits.
#[derive(Default, Debug, Clone, Copy)]
pub struct MergeStats {
    pub calls: u64,
    pub no_start_match: u64,
    pub no_end_match: u64,
    pub ends_too_long: u64,
    pub shared_under_100: u64,
    pub too_diverse: u64,
    pub structural: u64,
    pub merged: u64,
    /// Histogram of the internal indel/mismatch run length that triggered the
    /// `structural` rejection, bucketed: 6-10, 11-20, 21-50, 51-100, 101-200,
    /// 201+. This is what decides whether raising `delta_len` would merge
    /// consensus artefacts (short runs) or real isoform differences (long ones).
    pub run_6_10: u64,
    pub run_11_20: u64,
    pub run_21_50: u64,
    pub run_51_100: u64,
    pub run_101_200: u64,
    pub run_201_plus: u64,
    /// Sum of `min(len1, len2)` over all calls, for sizing a length prefilter.
    pub sum_min_len: u64,
    pub min_len_under_100: u64,
    /// `|len1 - len2|` histogram, split by outcome, to test whether a length
    /// prefilter could separate the structural rejects from the merges.
    pub diff_merged_max: u64,
    pub diff_struct_min: u64,
    pub diff_struct_under_merged_max: u64,
    pub last_diff: u64,
}

thread_local! {
    static MERGE_STATS: std::cell::RefCell<MergeStats> =
        const { std::cell::RefCell::new(MergeStats {
            calls: 0, no_start_match: 0, no_end_match: 0, ends_too_long: 0,
            shared_under_100: 0, too_diverse: 0, structural: 0, merged: 0,
            run_6_10: 0, run_11_20: 0, run_21_50: 0, run_51_100: 0,
            run_101_200: 0, run_201_plus: 0,
            sum_min_len: 0, min_len_under_100: 0,
            diff_merged_max: 0, diff_struct_min: 0,
            diff_struct_under_merged_max: 0, last_diff: 0,
        }) };
}

fn bump(f: impl Fn(&mut MergeStats)) {
    MERGE_STATS.with(|s| f(&mut s.borrow_mut()));
}

/// The counts gathered so far on this thread.
pub fn merge_stats() -> MergeStats {
    MERGE_STATS.with(|s| *s.borrow())
}

/// `align_to_merge`: should these two consensuses become one isoform?
pub fn align_to_merge<E: IsoformEngine>(
    engine: &mut E,
    consensus1: &[u8],
    consensus2: &[u8],
    opts: MergeOpts,
) -> bool {
    let min_len = consensus1.len().min(consensus2.len()) as u64;
    let diff = consensus1.len().abs_diff(consensus2.len()) as u64;
    bump(|s| {
        s.calls += 1;
        s.sum_min_len += min_len;
        s.last_diff = diff;
        if min_len < 100 {
            s.min_len_under_100 += 1;
        }
    });
    let (ops, a, b) = engine.align_merge(consensus1, consensus2);
    let overall_len = overall_alignment_len(&ops);
    const WINDOW: usize = 30;
    const THRESHOLD: f64 = 0.7;

    let Some(start_match) = find_first_significant_match(&a, &b, WINDOW, THRESHOLD) else {
        bump(|s| s.no_start_match += 1);
        return false;
    };
    // The same scan over the reversed alignments, which is the last significant
    // match counted from the end.
    let (ra, rb): (Vec<u8>, Vec<u8>) = (
        a.iter().rev().copied().collect(),
        b.iter().rev().copied().collect(),
    );
    let Some(end_match) = find_first_significant_match(&ra, &rb, WINDOW, THRESHOLD) else {
        bump(|s| s.no_end_match += 1);
        return false;
    };
    let end_match_pos = overall_len.saturating_sub(end_match);
    parse_cigar_diversity_isoform_level(&ops, opts, overall_len, start_match, end_match_pos)
}

/// `generate_all_consensuses`: one consensus per group.
///
/// A group of one read uses that read verbatim rather than calling spoa. Larger
/// groups pass the **first** `max_seqs_to_spoa` reads (`if i < max_seqs_to_spoa`,
/// a plain cap — note this is *not* the off-by-one `>` cutoff the correction
/// stage uses).
fn generate_all_consensuses<E: IsoformEngine>(
    engine: &mut E,
    groups: &[(u32, Vec<u32>)],
    reads: &FxHashMap<u32, Vec<u8>>,
    max_seqs_to_spoa: usize,
) -> Vec<(u32, Vec<u8>)> {
    let mut out = Vec::new();
    for (id, members) in groups {
        if members.len() == 1 {
            // Note it looks up `reads[id]`, the group key, not `members[0]`.
            let seq = reads.get(id).cloned().unwrap_or_default();
            out.push((*id, seq));
        } else {
            let seqs: Vec<&[u8]> = members
                .iter()
                .take(max_seqs_to_spoa)
                .filter_map(|q| reads.get(q).map(|v| v.as_slice()))
                .collect();
            out.push((*id, engine.spoa(&seqs)));
        }
    }
    out
}

/// `generate_new_full_consensus`: spoa over both groups' reads, capped.
fn generate_new_full_consensus<E: IsoformEngine>(
    engine: &mut E,
    a: &[u32],
    b: &[u32],
    reads: &FxHashMap<u32, Vec<u8>>,
    max_seqs_to_spoa: usize,
) -> Vec<u8> {
    let mut all: Vec<u32> = a.to_vec();
    all.extend_from_slice(b);
    // `all_consensus_infos[0:max_seqs_to_spoa]` --- a slice, so also a plain cap.
    all.truncate(max_seqs_to_spoa);
    let seqs: Vec<&[u8]> = all
        .iter()
        .filter_map(|q| reads.get(q).map(|v| v.as_slice()))
        .collect();
    engine.spoa(&seqs)
}

/// `merge_consensuses`: build a consensus per group, then fold similar pairs.
///
/// Returns `(id, consensus)` for the surviving groups, and mutates `groups` to
/// match --- a merged group's reads move into the group it merged with, and its
/// entry is removed, exactly as `add_merged_reads` does.
///
/// The scan order matters and is the reference's: consensuses are sorted by
/// **length ascending** (a stable sort, so equal lengths keep group order), each
/// is compared against every longer one, and the first successful merge ends that
/// consensus's scan.
pub fn merge_consensuses<E: IsoformEngine>(
    engine: &mut E,
    groups: &mut Vec<(u32, Vec<u32>)>,
    reads: &FxHashMap<u32, Vec<u8>>,
    opts: MergeOpts,
) -> Vec<(u32, Vec<u8>)> {
    let all_consensuses = generate_all_consensuses(engine, groups, reads, opts.max_seqs_to_spoa);
    let mut alternative = all_consensuses.clone();
    // `sort(key=lambda x: len(x[1]))` --- stable, so ties keep insertion order.
    alternative.sort_by_key(|(_, seq)| seq.len());

    let mut new_consensuses: Vec<(u32, Vec<u8>)> = Vec::new();
    let get =
        |v: &Vec<(u32, Vec<u8>)>, id: u32| v.iter().find(|(k, _)| *k == id).map(|(_, s)| s.clone());

    for i in 0..alternative.len().saturating_sub(1) {
        let id1 = alternative[i].0;
        // If id1 has already been merged into, use the merged consensus.
        let seq1 = get(&new_consensuses, id1).unwrap_or_else(|| alternative[i].1.clone());

        for cand in alternative.iter().skip(i + 1) {
            let id2 = cand.0;
            let seq2 = cand.1.clone();
            // Shorter first, whichever it is.
            let (c1, c2) = if seq2.len() < seq1.len() {
                (seq2.clone(), seq1.clone())
            } else {
                (seq1.clone(), seq2.clone())
            };
            if !align_to_merge(engine, &c1, &c2, opts) {
                continue;
            }
            let Some(pos2) = groups.iter().position(|(k, _)| *k == id2) else {
                break;
            };
            let Some(pos1) = groups.iter().position(|(k, _)| *k == id1) else {
                break;
            };
            let merged = if groups[pos2].1.len() > 50 {
                // Over fifty supporting reads: keep the existing consensus
                // rather than re-running spoa on the union.
                cand.1.clone()
            } else {
                // `generate_new_full_consensus(id1, id2)` concatenates
                // `curr_best_seqs[id1] + curr_best_seqs[id2]` --- id1's reads
                // first. spoa is order-sensitive, so this is not free to swap.
                generate_new_full_consensus(
                    engine,
                    &groups[pos1].1,
                    &groups[pos2].1,
                    reads,
                    opts.max_seqs_to_spoa,
                )
            };
            // `add_merged_reads(curr_best_seqs, id2, id1)`: id1's reads move into
            // id2, and id1 is removed.
            let moved = groups[pos1].1.clone();
            groups[pos2].1.extend(moved);
            groups.retain(|(k, _)| *k != id1);

            match new_consensuses.iter_mut().find(|(k, _)| *k == id2) {
                Some(slot) => slot.1 = merged,
                None => new_consensuses.push((id2, merged)),
            }
            // The reference breaks out of the inner loop on the first merge.
            break;
        }
    }

    // Anything never merged keeps its original consensus.
    for (id, _) in groups.iter() {
        if !new_consensuses.iter().any(|(k, _)| k == id) {
            if let Some(seq) = get(&all_consensuses, *id) {
                new_consensuses.push((*id, seq));
            }
        }
    }
    new_consensuses
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::ReadInfo;

    fn ri() -> ReadInfo {
        ReadInfo {
            start_mini_end: 0,
            end_mini_start: 0,
            original_support: true,
        }
    }

    /// s -> a -> t and s -> b -> t, reads 1,2 on a and read 3 on b.
    fn two_route_graph() -> Graph {
        let mut g = Graph::new();
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(NodeKey::Interval {
            start: 1,
            end: 2,
            r_id: 1,
        });
        let b = g.add_node(NodeKey::Interval {
            start: 3,
            end: 4,
            r_id: 3,
        });
        let t = g.add_node(NodeKey::Sink);
        for (u, v, rs) in [
            (s, a, vec![1u32, 2]),
            (a, t, vec![1, 2]),
            (s, b, vec![3]),
            (b, t, vec![3]),
        ] {
            g.add_edge(u, v, 0);
            g.upsert_edge_support(u, v, rs);
        }
        for n in [s, a, b, t] {
            for r in [1u32, 2, 3] {
                g.set_read(n, r, ri());
            }
        }
        g
    }

    #[test]
    fn reads_taking_the_same_route_are_grouped() {
        let g = two_route_graph();
        let got = compute_equal_reads(&g, &[1, 2, 3]);
        let mut groups: Vec<Vec<u32>> = got.iter().map(|(_, v)| v.clone()).collect();
        groups.sort();
        assert_eq!(groups, vec![vec![1, 2], vec![3]]);
    }

    #[test]
    fn a_group_is_keyed_by_one_of_its_own_members() {
        let g = two_route_graph();
        for (id, members) in compute_equal_reads(&g, &[1, 2, 3]) {
            assert!(members.contains(&id), "{id} is not in its own group");
        }
    }

    #[test]
    fn is_sublist_wants_a_contiguous_run() {
        assert!(is_sublist(&[1, 2, 3, 4], &[2, 3]));
        assert!(!is_sublist(&[1, 2, 3, 4], &[2, 4]), "not contiguous");
        assert!(is_sublist(&[1, 2], &[]), "empty is a sublist of anything");
        assert!(is_sublist(&[1, 2], &[1, 2]), "equal counts");
        assert!(!is_sublist(&[1], &[1, 2]), "longer cannot be contained");
    }

    #[test]
    fn merge_isoform_paths_folds_a_contained_group_into_its_container() {
        // Read-id lists, not paths --- which is the whole oddity. [2, 3] is a
        // contiguous run inside [1, 2, 3], so group 2 is folded into group 1.
        let mut isoforms = vec![(1u32, vec![1u32, 2, 3]), (2u32, vec![2u32, 3])];
        merge_isoform_paths(&mut isoforms);
        assert_eq!(isoforms.len(), 1);
        assert_eq!(isoforms[0].0, 1);
        assert_eq!(isoforms[0].1, vec![1, 2, 3, 2, 3], "extend, not union");
    }

    #[test]
    fn a_group_contained_in_two_others_does_not_panic() {
        // The reference deletes from the dict while walking a merge list built
        // beforehand, so this raises KeyError on the second delete. The port
        // skips the already-removed entry instead.
        let mut isoforms = vec![
            (1u32, vec![9u32, 9]),
            (2u32, vec![9u32, 9]),
            (3u32, vec![9u32, 9]),
        ];
        merge_isoform_paths(&mut isoforms);
        assert!(!isoforms.is_empty(), "survives rather than panicking");
    }

    #[test]
    fn a_significant_match_needs_to_exceed_the_threshold_not_merely_reach_it() {
        // 21 of 30 is 0.7 exactly, and the test is `>`, so it is not significant.
        let mut a = vec![b'A'; 30];
        let b = vec![b'A'; 30];
        for x in a.iter_mut().take(9) {
            *x = b'C';
        }
        assert_eq!(
            find_first_significant_match(&a, &b, 30, 0.7),
            None,
            "exactly 0.7 does not qualify"
        );
        a[0] = b'A'; // 22 of 30
        assert_eq!(find_first_significant_match(&a, &b, 30, 0.7), Some(0));
    }

    #[test]
    fn no_window_fits_when_the_alignment_is_shorter_than_the_window() {
        assert_eq!(find_first_significant_match(b"AAA", b"AAA", 30, 0.7), None);
    }

    #[test]
    fn a_shared_region_under_a_hundred_bases_never_merges() {
        // The hard floor, independent of every tuning parameter.
        let opts = MergeOpts {
            delta: 0.9,
            delta_len: 100,
            delta_iso_len_3: 1000,
            delta_iso_len_5: 1000,
            max_seqs_to_spoa: 200,
            cigar_diversity_counts_runs: false,
        };
        let ops: Vec<CigarOp> = vec![(99, b'=')];
        assert!(
            !parse_cigar_diversity_isoform_level(&ops, opts, 99, 0, 99),
            "99 shared bases is under the floor even with everything else wide open"
        );
        let ops: Vec<CigarOp> = vec![(150, b'=')];
        assert!(parse_cigar_diversity_isoform_level(&ops, opts, 150, 0, 150));
    }

    #[test]
    fn an_interior_mismatch_longer_than_delta_len_is_structural() {
        let opts = MergeOpts {
            delta: 0.9,
            delta_len: 5,
            delta_iso_len_3: 1000,
            delta_iso_len_5: 1000,
            max_seqs_to_spoa: 200,
            cigar_diversity_counts_runs: false,
        };
        // A 6-base interior mismatch: longer than delta_len, so structural.
        let ops: Vec<CigarOp> = vec![(100, b'='), (6, b'X'), (100, b'=')];
        assert!(!parse_cigar_diversity_isoform_level(
            &ops, opts, 206, 0, 206
        ));
        // Five is tolerated.
        let ops: Vec<CigarOp> = vec![(100, b'='), (5, b'X'), (100, b'=')];
        assert!(parse_cigar_diversity_isoform_level(&ops, opts, 205, 0, 205));
    }

    #[test]
    fn bug_compat_counts_every_interior_mismatch_as_delta_len() {
        // finding 30 reproduced: the accumulator is `+= delta_len`, not
        // `+= cig_len`. So one 1-base mismatch and one 5-base mismatch
        // contribute the same, and two 1-base mismatches contribute twice as
        // much as one 5-base one.
        let opts = MergeOpts {
            delta: 0.0,
            delta_len: 5,
            delta_iso_len_3: 1000,
            delta_iso_len_5: 1000,
            max_seqs_to_spoa: 200,
            cigar_diversity_counts_runs: true,
        };
        // delta 0 means max_bp_diff is delta_len (5), so the budget is 5.
        // One mismatch of any interior length contributes exactly 5: at the
        // budget, and the test is `<=`, so it merges.
        let one_short: Vec<CigarOp> = vec![(200, b'='), (1, b'X'), (200, b'=')];
        assert!(parse_cigar_diversity_isoform_level(
            &one_short, opts, 401, 0, 401
        ));
        let one_long: Vec<CigarOp> = vec![(200, b'='), (5, b'X'), (200, b'=')];
        assert!(parse_cigar_diversity_isoform_level(
            &one_long, opts, 405, 0, 405
        ));
        // Two 1-base mismatches contribute 10, over the budget.
        let two_short: Vec<CigarOp> =
            vec![(150, b'='), (1, b'X'), (150, b'='), (1, b'X'), (150, b'=')];
        assert!(!parse_cigar_diversity_isoform_level(
            &two_short, opts, 452, 0, 452
        ));
    }

    #[test]
    fn the_fixed_diversity_measure_counts_bases_not_runs() {
        // finding 30 fixed. Same budget, same two shapes, opposite verdicts from
        // the bug-compat test above: counting bases, two 1-base mismatches cost
        // 2 and one 5-base mismatch costs 5.
        let opts = MergeOpts {
            delta: 0.0,
            delta_len: 5,
            delta_iso_len_3: 1000,
            delta_iso_len_5: 1000,
            max_seqs_to_spoa: 200,
            cigar_diversity_counts_runs: false,
        };
        // Budget is delta_len = 5 bases.
        let two_singles: Vec<CigarOp> =
            vec![(100, b'='), (1, b'X'), (50, b'='), (1, b'X'), (100, b'=')];
        assert!(
            parse_cigar_diversity_isoform_level(&two_singles, opts, 252, 0, 252),
            "two single-base mismatches cost 2 of a 5-base budget"
        );
        // Under bug-compat the same shape costs 2 x delta_len = 10 and fails.
        let runs = MergeOpts {
            cigar_diversity_counts_runs: true,
            ..opts
        };
        assert!(
            !parse_cigar_diversity_isoform_level(&two_singles, runs, 252, 0, 252),
            "and the reference charges 2 x delta_len = 10 for it"
        );
    }

    #[test]
    fn short_alignments_get_twice_the_diversity_tolerance() {
        // `if overall_len < 200: delta = 2 * delta`.
        let opts = MergeOpts {
            delta: 0.02,
            delta_len: 1,
            delta_iso_len_3: 1000,
            delta_iso_len_5: 1000,
            max_seqs_to_spoa: 200,
            cigar_diversity_counts_runs: false,
        };
        // similar_seq 150, four interior mismatches -> miss_match_length 4.
        // budget = max(delta * 150, 1); at delta 0.02 that is 3 -> refuse.
        // At 2 * delta it is 6 -> allow. overall_len decides which.
        let ops: Vec<CigarOp> = vec![
            (30, b'='),
            (1, b'X'),
            (30, b'='),
            (1, b'X'),
            (30, b'='),
            (1, b'X'),
            (30, b'='),
            (1, b'X'),
            (26, b'='),
        ];
        let total: usize = ops.iter().map(|(n, _)| n).sum();
        assert_eq!(total, 150);
        assert!(
            parse_cigar_diversity_isoform_level(&ops, opts, 150, 0, 150),
            "overall_len 150 is under 200, so delta doubles"
        );
        assert!(
            !parse_cigar_diversity_isoform_level(&ops, opts, 250, 0, 150),
            "the same alignment with overall_len 250 gets the single delta"
        );
    }
}
