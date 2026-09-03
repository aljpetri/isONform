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

/// Order the reads handed to spoa by **decreasing length**.
///
/// POA is incremental: the first sequence becomes the graph's backbone and every
/// later one is aligned onto it, so seeding with the longest read gives a
/// full-length backbone. There is a second, separate effect: the cap keeps the
/// **first** `max_seqs_to_spoa` reads, so sorting changes *which* reads a group
/// larger than the cap uses.
///
/// **Measured and rejected as a default** (finding 44): no runtime benefit (~5%
/// slower), neutral on `sirv_sim_gene`, better on droso (FSM 476 -> 481), but
/// clearly worse on **real** ONT data --- `sirv_real` strict F1 0.735 -> 0.700 and
/// lenient 0.884 -> 0.865. The simulated/real split fits the failure mode: the
/// longest real reads are the likeliest to be chimeric, and a chimeric backbone
/// propagates into the consensus. Kept behind `ISONFORM_SPOA_SORT=1`, off.
pub fn spoa_sort_by_length() -> bool {
    static F: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *F.get_or_init(|| {
        matches!(
            std::env::var("ISONFORM_SPOA_SORT").ok().as_deref(),
            Some("1") | Some("on")
        )
    })
}

/// The reads of `members`, longest first when [`spoa_sort_by_length`] is on, then
/// capped. Sorting happens **before** the cap so the cap keeps the longest.
/// How the sequences handed to spoa are ordered. `ISONFORM_SPOA_ORDER`.
///
/// spoa seeds its graph with the **first** sequence, so that one is the backbone
/// every later sequence is aligned into --- the order is not cosmetic.
///
/// `len` (equivalently `ISONFORM_SPOA_SORT=1`) sorts longest-first *before* the
/// `max_seqs_to_spoa` cap, so it changes the backbone **and** which reads survive
/// the cap. Finding 44 measured that and rejected it: clearly worse on real ONT
/// data (`sirv_real` strict F1 0.735 -> 0.700), because the longest real reads are
/// the likeliest to be chimeric and a chimeric backbone propagates.
///
/// The two variants below separate the knobs finding 44 conflated. Both cap
/// first, in the port's ordinary ascending order, so selection is untouched and
/// only the backbone changes:
///
/// * `cap_len` --- longest of the survivors as backbone. The backbone benefit
///   without letting long chimeras displace ordinary reads at the cap.
/// * `median` --- the **median-length** survivor as backbone. Still a full-length
///   read, but not the outlier that chimeras concentrate in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpoaOrder {
    /// Ascending read id, capped. The port's default.
    None,
    /// Longest first, then capped. Finding 44's variant.
    Len,
    /// Capped, then longest of the survivors first.
    CapLen,
    /// Capped, then the median-length survivor first.
    Median,
}

pub fn spoa_order() -> SpoaOrder {
    static O: std::sync::OnceLock<SpoaOrder> = std::sync::OnceLock::new();
    *O.get_or_init(
        || match std::env::var("ISONFORM_SPOA_ORDER").ok().as_deref() {
            Some("len") => SpoaOrder::Len,
            Some("cap_len") => SpoaOrder::CapLen,
            Some("median") => SpoaOrder::Median,
            _ if spoa_sort_by_length() => SpoaOrder::Len,
            _ => SpoaOrder::None,
        },
    )
}

fn ordered_for_spoa(
    members: &[u32],
    reads: &FxHashMap<u32, Vec<u8>>,
    max_seqs_to_spoa: usize,
) -> Vec<u32> {
    let mut v = members.to_vec();
    let len = |q: &u32| reads.get(q).map_or(0, |r| r.len());
    match spoa_order() {
        SpoaOrder::None => {}
        SpoaOrder::Len => {
            // Descending by read length; stable, so equal lengths keep their
            // order. Before the cap, so the cap keeps the longest.
            v.sort_by(|x, y| len(y).cmp(&len(x)));
        }
        SpoaOrder::CapLen => {
            v.truncate(max_seqs_to_spoa);
            v.sort_by(|x, y| len(y).cmp(&len(x)));
        }
        SpoaOrder::Median => {
            v.truncate(max_seqs_to_spoa);
            if v.len() > 2 {
                // Index of the median-length survivor, by a stable ordering so
                // the pick does not depend on sort implementation.
                let mut by_len: Vec<u32> = v.clone();
                by_len.sort_by(|x, y| len(x).cmp(&len(y)).then(x.cmp(y)));
                let med = by_len[by_len.len() / 2];
                if let Some(i) = v.iter().position(|&r| r == med) {
                    v.remove(i);
                    v.insert(0, med);
                }
            }
        }
    }
    v.truncate(max_seqs_to_spoa);
    v
}

/// The two POA-scheduling defaults, overridable from the environment.
///
/// `ISONFORM_MERGE_REBUILD_MAX` (default 0) and `ISONFORM_FINAL_CONSENSUS`
/// (default on, `0`/`off` to disable). Together at their defaults these are
/// finding 44's scheme; `MERGE_REBUILD_MAX=50 FINAL_CONSENSUS=0` restores the
/// reference's. Read once, and only by the binaries --- library callers pass
/// [`MergeOpts`] explicitly so tests and oracles are not at the mercy of the
/// environment.
pub fn poa_schedule_from_env() -> (usize, bool) {
    static S: std::sync::OnceLock<(usize, bool)> = std::sync::OnceLock::new();
    *S.get_or_init(|| {
        // Faithful mode is the reference's schedule: rebuild after every merge
        // (capped at 50) and no final pass. Finding 44's end-only scheme is the
        // optimisation, opted into with `ISONFORM_FAITHFUL=0` or by setting these
        // explicitly.
        let (d_max, d_fin) = if crate::faithful() {
            (50, false)
        } else {
            (0, true)
        };
        let max = std::env::var("ISONFORM_MERGE_REBUILD_MAX")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(d_max);
        let fin = match std::env::var("ISONFORM_FINAL_CONSENSUS").ok().as_deref() {
            Some("0") | Some("off") => false,
            Some(_) => true,
            None => d_fin,
        };
        (max, fin)
    })
}

/// Nanoseconds spent inside consensus building and inside pairwise alignment.
///
/// The `merge` stage timer covers both `generate_all_consensuses` (partial-order
/// alignment) and `merge_consensuses` (the pairwise scan), and on
/// `sirv_sim_gene`'s critical-path batch that stage is 94% of the wall clock. A
/// faster pairwise aligner only helps to the extent the second of the two
/// dominates, so the two are timed apart rather than argued about.
pub static POA_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static ALIGN_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static POA_CALLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static ALIGN_CALLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn add(c: &std::sync::atomic::AtomicU64, n: u64) {
    c.fetch_add(n, std::sync::atomic::Ordering::Relaxed);
}

/// The real thing: [`crate::poa`] and [`crate::parasail`].
#[derive(Debug, Default)]
pub struct SpoaParasailMerge;

impl IsoformEngine for SpoaParasailMerge {
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8> {
        let t = std::time::Instant::now();
        let out = crate::poa::consensus(seqs)
            .map(|s| s.into_bytes())
            .unwrap_or_default();
        add(&POA_NS, t.elapsed().as_nanos() as u64);
        add(&POA_CALLS, 1);
        out
    }

    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let t = std::time::Instant::now();
        let sc = crate::parasail::Scoring::MERGE;
        let aln = crate::wfa::enabled_merge()
            .then(|| crate::wfa::semiglobal(s1, s2, sc))
            .flatten()
            .unwrap_or_else(|| crate::parasail::semiglobal(s1, s2, sc));
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        add(&ALIGN_NS, t.elapsed().as_nanos() as u64);
        add(&ALIGN_CALLS, 1);
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
        // The insertion order the reference's set has, which decides its slot
        // layout and so `list(...)`. Each `intersection(edge_supp)` builds a
        // **fresh** set by iterating `edge_supp` --- a list, so in list order ---
        // and inserting the members it finds, so the order is simply the last
        // traversed edge's support filtered to the survivors. `None` until the
        // first intersection, while `current_node_support` is still the alias.
        let mut insertion: Option<Vec<u32>> = None;
        let mut node: NodeId = s;
        while node != t {
            let mut next_found = false;
            for &nb in g.successors(node) {
                if let Some(supp) = g.edge_support(node, nb) {
                    if supp.contains(&read) {
                        let keep: FxHashSet<u32> = supp.iter().copied().collect();
                        insertion = Some(
                            supp.iter()
                                .copied()
                                .filter(|r| current.contains(r))
                                .collect(),
                        );
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
        // `list(current_node_support)` is CPython set-iteration order, which is
        // slot order and therefore a function of the values and the order they
        // were inserted. `crate::pyset` models that exactly, so faithful mode
        // reproduces the reference's member order and its `[0]` representative.
        // Finding 28's ascending order is the divergence, kept for the optimised
        // configuration so the two remain comparable.
        let group: Vec<u32> = match (&insertion, crate::reference_semantics()) {
            (Some(ins), true) => crate::pyset::order_of(ins.iter().map(|&r| r as u64))
                .into_iter()
                .map(|r| r as u32)
                .collect(),
            _ => {
                let mut v: Vec<u32> = current.iter().copied().collect();
                v.sort_unstable();
                v
            }
        };
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
#[derive(Debug, Clone, Copy)]
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
    /// Above how many supporting reads a merge skips rebuilding the consensus.
    ///
    /// The reference hard-codes 50, testing the *target* group's read count
    /// rather than the size of the union it would POA. Since most groups sit far
    /// below any threshold the check almost never fires, so every merge pays for
    /// a full POA over both groups' raw reads --- about 259 of the 278 POA calls
    /// on `sirv_sim_gene`'s critical-path batch, where POA is 77% of the wall
    /// clock (finding 42).
    ///
    /// **Defaults to 0**, i.e. never rebuild mid-merge, which only makes sense
    /// together with [`MergeOpts::final_consensus_pass`]. Set to 50 for the
    /// reference's behaviour.
    pub merge_rebuild_max: usize,
    /// Rebuild each surviving group's consensus once, after all merging is done.
    ///
    /// The per-merge rebuild serves two unrelated purposes: it supplies the
    /// consensus later merge *decisions* compare against, and it supplies the
    /// base accuracy of the emitted *sequence*. Finding 43 separated them ---
    /// switching it off kept the isoform structure but cost accuracy. Only the
    /// second purpose needs the union of both read sets, and it needs it once.
    ///
    /// So merge with the consensuses already in hand, then POA each surviving
    /// group once over its final membership: 93 POA calls where the reference
    /// makes 278, and better accuracy on all three corpora (finding 44).
    ///
    /// **Defaults to true.** A deviation: the reference concatenates id1's reads
    /// then id2's, while this POAs the member list in its own order, and spoa is
    /// order-sensitive --- so the consensus differs even on an identical read set,
    /// and merge decisions differ because they compare un-rebuilt consensuses.
    pub final_consensus_pass: bool,
}

impl Default for MergeOpts {
    /// Finding 44's POA schedule, matching what the binaries do without any
    /// environment set. Deriving this would silently give
    /// `merge_rebuild_max: 0` **and** `final_consensus_pass: false` --- rebuild
    /// off with nothing replacing it, which is the one configuration measured to
    /// lose accuracy (finding 43). Anything comparing against the reference must
    /// set `50` / `false` explicitly.
    fn default() -> Self {
        Self {
            delta: 0.0,
            delta_len: 0,
            delta_iso_len_3: 0,
            delta_iso_len_5: 0,
            max_seqs_to_spoa: 0,
            merge_rebuild_max: 0,
            final_consensus_pass: true,
            cigar_diversity_counts_runs: false,
        }
    }
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
/// `align_to_merge`, with WFA2's **positive** verdicts confirmed by parasail.
///
/// WFA2 merges more than parasail does --- 9% of parasail's merges flip --- and
/// the extra merges rebuild the consensus over a larger read union, so the ends
/// creep out and the isoform lands just under the scorer's identity cutoff. The
/// asymmetry is what makes this cheap: merges are rare (1 650 of 51 447 recorded
/// calls, 3.2%), so re-deciding only the positives with the exact aligner costs a
/// few percent of the alignment work and makes every merge exactly parasail's.
///
/// It does **not** catch pairs WFA2 wrongly refuses; that would need confirming
/// the negatives too, which is every pair and no saving. Whether those matter is
/// a measurement, not an argument.
///
/// `ISONFORM_WFA2_CONFIRM=1`.
pub fn wfa_confirm() -> bool {
    static C: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *C.get_or_init(|| {
        matches!(
            std::env::var("ISONFORM_WFA2_CONFIRM").ok().as_deref(),
            Some("1") | Some("on")
        )
    })
}

/// The exact-aligner verdict for a pair, for [`wfa_confirm`].
struct ParasailMergeOnly;

impl IsoformEngine for ParasailMergeOnly {
    fn spoa(&mut self, _seqs: &[&[u8]]) -> Vec<u8> {
        unreachable!("confirmation never builds a consensus")
    }
    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let aln = crate::parasail::semiglobal(s1, s2, crate::parasail::Scoring::MERGE);
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        (aln.ops, a, b)
    }
}

pub fn align_to_merge<E: IsoformEngine>(
    engine: &mut E,
    consensus1: &[u8],
    consensus2: &[u8],
    opts: MergeOpts,
) -> bool {
    let verdict = align_to_merge_inner(engine, consensus1, consensus2, opts);
    // Only a positive needs confirming, and only when the fast aligner produced
    // it --- with WFA2 off the two are the same aligner and this is a no-op.
    if verdict && wfa_confirm() && crate::wfa::enabled_merge() {
        return align_to_merge_inner(&mut ParasailMergeOnly, consensus1, consensus2, opts);
    }
    verdict
}

fn align_to_merge_inner<E: IsoformEngine>(
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
            let picked = ordered_for_spoa(members, reads, max_seqs_to_spoa);
            let seqs: Vec<&[u8]> = picked
                .iter()
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
    // `ordered_for_spoa` applies that same cap, after sorting when asked to.
    let all = ordered_for_spoa(&all, reads, max_seqs_to_spoa);
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
            let merged = if groups[pos2].1.len() > opts.merge_rebuild_max {
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

    if opts.final_consensus_pass {
        // One POA per surviving group, over the union its members now hold. The
        // membership in `groups` is final at this point, so this is exactly
        // `generate_all_consensuses` on the merged groups --- and it is what makes
        // skipping the per-merge rebuild safe for the emitted sequence.
        let rebuilt = generate_all_consensuses(engine, groups, reads, opts.max_seqs_to_spoa);
        for (id, seq) in rebuilt {
            if seq.is_empty() {
                continue;
            }
            match new_consensuses.iter_mut().find(|(k, _)| *k == id) {
                Some(slot) => slot.1 = seq,
                None => new_consensuses.push((id, seq)),
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
            merge_rebuild_max: 50,
            final_consensus_pass: false,
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
