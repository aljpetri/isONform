//! Construction of the interval graph.
//!
//! Ports `modules/GraphGeneration.py::generateGraphfromIntervals` and the five
//! helpers it dispatches to. Structure follows the reference branch for branch,
//! because that is what makes a stage-level diff meaningful; the tidying comes
//! after the port is proven faithful, in its own commit.
//!
//! # The two reference defects this reproduces on purpose
//!
//! Both are `PORTING.md` findings 8 and 9. Reproducing rather than fixing is the
//! working agreement: a behaviour change does not land in the same commit as the
//! port. Each is behind a flag in [`BuildOpts`] so the fixed behaviour can be
//! measured against the reference without editing this file.
//!
//! * **Finding 8, latent.** `known_old_node_action` calls
//!   `DG.add_edge(previous_node, name, this_len)` — a third *positional*
//!   argument, which networkx rejects with `TypeError: add_edge() takes 3
//!   positional arguments but 4 were given`. Verified against networkx 2.8.4. It
//!   is a crash, not a wrong answer, and it needs the same `old_node` to be hit
//!   a third time with a matching predecessor and a close length; measured over
//!   `sirv_real`, `droso` and `sirv_sim_gene`, the filter that guards it produced
//!   4 541 candidate hits and never once matched, so the crash never fires. This
//!   port returns [`BuildError::ReferenceWouldCrash`] there rather than
//!   inventing a behaviour the reference does not have.
//! * **Finding 11, live.** `cycle_added` creates a node to route around a cycle
//!   but the caller keeps the *old* node as `previous_node`, because the rebinding
//!   is local to the helper. The new node ends up with an in-edge and no out-edge
//!   and the read's path is severed. Found by the oracle on real Drosophila data,
//!   not by reading.
//! * **Finding 9, live.** One branch reads `seq` without binding it, so it uses
//!   whatever `seq` held from a previous loop iteration — across reads, that is
//!   *another read's sequence*, and the k-mer taken from it becomes the node's
//!   `end_mini_seq`. Measured: 252 of 255 occurrences on `sirv_real` and 913 of
//!   941 on `droso` used a stale sequence.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::{occurrence_hash, Graph, Interval, NodeId, NodeKey, ReadInfo};

/// Inputs to graph construction, mirroring the reference's five parameters.
pub struct BuildInput<'a> {
    /// Minimizer length.
    pub k: usize,
    /// `--delta_len`.
    pub delta_len: i64,
    /// `all_intervals_for_graph`, as an ordered list rather than a map: the
    /// reference iterates `.items()`, so insertion order decides which node is
    /// created first and it must be explicit rather than left to a hash map.
    pub intervals: &'a [(u32, Vec<Interval>)],
    /// `all_reads` — sequences for the reads the graph sees (`new_all_reads` in
    /// `main`), keyed by the same ids as `intervals`.
    pub reads: &'a FxHashMap<u32, Vec<u8>>,
    /// `read_len_dict`, ordered by read id.
    ///
    /// Built in `main` over **`all_reads`**, not `new_all_reads`, so it covers
    /// every input read including those the interval filter skipped. The sink
    /// node's read set is therefore larger than the graph: measured on
    /// `sirv_small`, 100 entries against 84 reads that have a path. Faithful; see
    /// `PORTING.md`.
    pub read_len: &'a [(u32, u32)],
}

/// Bug fixes. **On by default** since the bug-fix policy changed; use
/// [`BuildOpts::reference`] to reproduce the reference's behaviour, which is what
/// the oracles replay against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BuildOpts {
    /// Bind `seq` from the read being processed, as every other branch does.
    /// Fixes finding 9.
    pub fix_stale_seq: bool,
    /// Continue the read's path from the node `cycle_added` just created, rather
    /// than from the one whose incoming edge it removed. Fixes finding 11.
    pub fix_cycle_continuation: bool,
    /// Seed the sink with the **graph's** reads and their own lengths, as the
    /// source loop two lines above already does, rather than with one entry per
    /// read in the whole *file* indexed as though those were graph ids.
    /// Fixes finding 33.
    ///
    /// The reference's two loops are adjacent and disagree:
    ///
    /// ```python
    /// for i in range(1, len(all_intervals_for_graph) + 1):   # source: graph ids
    /// for i in range(1, len(read_len_dict) + 1):             # sink: file read ids
    /// ```
    ///
    /// Every read has one path, source to sink, so these are the same set. The
    /// source's range is the correct one --- `reads_for_isoforms` is built from
    /// it. Two consequences follow from the sink's being wrong: phantom entries
    /// for graph ids that do not exist, and, for the ids that do, some other
    /// read's length.
    pub fix_sink_read_len: bool,
}

impl Default for BuildOpts {
    /// Correct behaviour: bugs fixed.
    fn default() -> Self {
        Self {
            fix_stale_seq: true,
            fix_cycle_continuation: true,
            fix_sink_read_len: true,
        }
    }
}

impl BuildOpts {
    /// Reproduce the reference's bugs exactly. What the oracles replay against.
    pub fn reference() -> Self {
        Self {
            fix_stale_seq: false,
            fix_cycle_continuation: false,
            fix_sink_read_len: false,
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum BuildError {
    /// The reference raises `TypeError` here. See finding 8 in the module docs.
    ReferenceWouldCrash {
        detail: &'static str,
    },
    /// The reference would raise `NameError`/`UnboundLocalError`: it reads `seq`
    /// on a path that never bound it, and no earlier iteration bound it either.
    UnboundSeq,
    MissingRead(u32),
}

/// The `k`-mer the reference takes as a node's `end_mini_seq`:
/// `seq[inter[1] : inter[1] + k]`.
///
/// Python slicing clamps, so an interval end within `k` of the sequence end
/// yields a short slice rather than an error. Reproduced.
fn end_mini(seq: &[u8], end: u32, k: usize) -> &[u8] {
    let lo = (end as usize).min(seq.len());
    let hi = (lo + k).min(seq.len());
    &seq[lo..hi]
}

/// How often a read finds a node another read predicted for it.
///
/// Registration is *speculative*: `add_prior_read_infos` records where an
/// interval will sit in every later read that shares it. The prediction only pays
/// off if that read's own WIS solution selects an interval with exactly those
/// coordinates. So the hit rate measures **WIS agreement across reads**, which is
/// what decides whether the graph collapses or fragments.
pub static LOOKUPS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
/// `(read, start bucket) -> [(start, end, node)]`, the tolerant lookup's index.
type TolIndex = FxHashMap<(u32, u32), Vec<(u32, u32, NodeId)>>;

/// Assign this read's intervals to existing nodes by **co-linear chaining**
/// instead of an exact-coordinate dict hit.
///
/// # Why
///
/// Finding 45 measured that on the degenerate instance **every** one of 20 376
/// nodes comes from the first-sight branch --- no cycles, no within-read repeats
/// (2 repetitive intervals in 29 403). Reads simply never agree on interval
/// coordinates, and 94% of the misses have a registered coordinate within 20
/// bases. The exact lookup asks "did somebody register this precise coordinate";
/// what it should ask is "does this read's sequence of anchors follow a path
/// already in the graph".
///
/// Finding 39 tried the obvious weakening --- accept the nearest registered
/// coordinate within a tolerance --- and rejected it: nodes collapsed 20 421 ->
/// 1 754 but exact hits *fell*, 9 027 -> 8 106, because a greedy nearest match
/// attaches a read to whichever copy happens to be closest with no regard for
/// order. Chaining is that idea done properly: the assignment must be monotone in
/// the read's own coordinates, so order and spacing disambiguate.
///
/// # The formulation
///
/// Registered coordinates for read `r` are positions **in read `r`**, and so are
/// the read's own interval starts, which makes this a one-dimensional chaining
/// problem rather than anything to do with graph topology. Candidates for
/// interval `i` are nodes whose registered coordinate for `r` is within `tol`;
/// a chain picks at most one per interval such that registered coordinates
/// strictly increase with interval order. Maximise the number of chained
/// intervals, and among those minimise total cost.
///
/// # Cost has two halves: offset and spacing
///
/// Monotonicity alone is weak. It admits a chain whose nodes are in the right
/// order but at the wrong *distances* --- exactly how a read gets pulled onto a
/// subtly wrong node, and the likeliest source of the `sirv_real` strict loss in
/// finding 46. So each transition also pays for **spacing disagreement**: the gap
/// between two consecutive chained nodes' registered coordinates should match the
/// gap between the read's own two interval starts. Since selected intervals tile
/// the read (coverage 0.974, ~60% abutting), that gap is a strong signal ---
/// consecutive tiles are ~40 bases apart, so a mismatched spacing is a real
/// structural disagreement rather than drift.
///
/// It is a **hard bound**, not a tie-breaker. As a soft cost it was measured to
/// be exactly inert: identical chains at weights 0, 1 and 8 for every tolerance
/// tried. Two reasons, both arithmetic. Within a tolerance the skew
/// `|d_i - d_j|` is bounded by the offset cost `|d_i| + |d_j|` the DP already
/// minimises, so it is near-redundant; and the objective maximises *count* first,
/// so cost only breaks ties among chains of equal length. A bound instead removes
/// candidate transitions, which changes which chains exist at all.
///
/// `ISONFORM_CHAIN_SPACING` is the largest tolerated skew in bases, default 0
/// meaning unbounded --- the monotonicity-only chain finding 46 measured.
///
/// Exact hits are zero-offset candidates, so they are never lost by chaining;
/// the caller consults the chain only where the exact lookup missed.
fn chain_assign(
    starts: &[(u32, u32)],
    r_id: u32,
    tol_index: &TolIndex,
    tol: u32,
    bucket_span: u32,
    max_skew: u32,
) -> Vec<Option<NodeId>> {
    let n = starts.len();
    let mut cands: Vec<Vec<(u32, NodeId, u32)>> = vec![Vec::new(); n];
    for (i, &(st_i, en_i)) in starts.iter().enumerate() {
        let b = st_i / bucket_span;
        for bb in b.saturating_sub(1)..=b + 1 {
            if let Some(v) = tol_index.get(&(r_id, bb)) {
                for &(st, en, nd) in v {
                    if st.abs_diff(st_i) <= tol && en.abs_diff(en_i) <= tol {
                        cands[i].push((st, nd, st.abs_diff(st_i) + en.abs_diff(en_i)));
                    }
                }
            }
        }
        // Cheapest first, so ties in chain length resolve toward the closest.
        cands[i].sort_by_key(|&(st, _, cost)| (cost, st));
        cands[i].dedup_by_key(|&mut (st, nd, _)| (st, nd));
    }

    // dp[i][a] = best (count, -cost) for a chain ending at interval i with
    // candidate a. `n` is ~30 and candidates per interval are few, so the
    // quadratic scan is free.
    let mut dp: Vec<Vec<(u32, u32, Option<(usize, usize)>)>> =
        cands.iter().map(|c| vec![(1, 0, None); c.len()]).collect();
    for i in 0..n {
        for a in 0..cands[i].len() {
            let (st_a, _, cost_a) = cands[i][a];
            let mut best = (1u32, cost_a, None);
            for j in 0..i {
                for b in 0..cands[j].len() {
                    let st_b = cands[j][b].0;
                    if st_b >= st_a {
                        continue; // must strictly increase
                    }
                    // Spacing: the registered coordinates must advance by about
                    // the same amount the read's own interval starts do. Rejecting
                    // the transition is what makes this bite --- as an added cost it
                    // measured completely inert.
                    let reg_gap = st_a - st_b;
                    let own_gap = starts[i].0.saturating_sub(starts[j].0);
                    let skew = reg_gap.abs_diff(own_gap);
                    if max_skew > 0 && skew > max_skew {
                        continue;
                    }
                    let (cnt, cost, _) = dp[j][b];
                    let cand = (cnt + 1, cost + cost_a + skew, Some((j, b)));
                    // More chained intervals wins; then less total offset.
                    if cand.0 > best.0 || (cand.0 == best.0 && cand.1 < best.1) {
                        best = cand;
                    }
                }
            }
            dp[i][a] = best;
        }
    }

    let mut end: Option<(usize, usize)> = None;
    let mut best = (0u32, u32::MAX);
    for i in 0..n {
        for a in 0..cands[i].len() {
            let (cnt, cost, _) = dp[i][a];
            if cnt > best.0 || (cnt == best.0 && cost < best.1) {
                best = (cnt, cost);
                end = Some((i, a));
            }
        }
    }

    let mut out = vec![None; n];
    let mut cur = end;
    while let Some((i, a)) = cur {
        out[i] = Some(cands[i][a].1);
        cur = dp[i][a].2;
    }
    out
}

/// Lookups rescued by `ISONFORM_LOOKUP_TOL` that the exact lookup missed.
pub static TOL_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static REGISTRATIONS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
/// Missed lookups bucketed by distance to the nearest prediction for the same
/// read: exact-but-absent, 1-5, 6-20, 21-100, 101+, and "no prediction at all".
pub static MISS_DIST: [std::sync::atomic::AtomicU64; 6] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];

/// `add_prior_read_infos`: record this interval's occurrences in *later* reads,
/// so that when those reads are processed they find this node.
///
/// The reference's condition is `if not r <= r_id`, i.e. strictly greater, and
/// the recorded coordinates are `(r, pos1 + k, pos2)` — note the `+ k` on the
/// start only. First write wins.
fn add_prior_read_infos(
    prior: &mut FxHashMap<(u32, u32, u32), NodeId>,
    inter: &Interval,
    r_id: u32,
    node: NodeId,
    k: usize,
) {
    for triple in inter.occurrences.chunks_exact(3) {
        let (r, p1, p2) = (triple[0], triple[1], triple[2]);
        if r > r_id {
            REGISTRATIONS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            prior.entry((r, p1 + k as u32, p2)).or_insert(node);
        }
    }
}

/// Reproduces `generateGraphfromIntervals`.
///
/// Returns the graph and `reads_for_isoforms`, which is the reference's second
/// return value — simply `1..=len(all_intervals_for_graph)`, built before any
/// work and never filtered.
pub fn generate_graph_from_intervals(
    input: &BuildInput<'_>,
    opts: BuildOpts,
) -> Result<(Graph, Vec<u32>), BuildError> {
    let k = input.k;
    let delta_len = input.delta_len;
    let mut g = Graph::new();

    // Source and sink, with their read sets, exactly as the reference seeds them.
    let n_graph_reads = input.intervals.len() as u32;
    let s = g.add_node(NodeKey::Source);
    for i in 1..=n_graph_reads {
        g.set_read(
            s,
            i,
            ReadInfo {
                start_mini_end: 0,
                end_mini_start: 0,
                original_support: true,
            },
        );
    }
    let t = g.add_node(NodeKey::Sink);
    if opts.fix_sink_read_len {
        // finding 33 fixed: the sink carries the *graph's* reads, each with its
        // own length --- the same range the source loop two lines above uses.
        for i in 1..=n_graph_reads {
            let len = input.reads.get(&i).map(|sq| sq.len()).unwrap_or(0);
            g.set_read(
                t,
                i,
                ReadInfo {
                    start_mini_end: len as i64,
                    end_mini_start: len as i64,
                    original_support: true,
                },
            );
        }
    } else {
        // finding 33 reproduced: `range(1, len(read_len_dict) + 1)` over the
        // *whole file*, indexed as though it were graph ids.
        for &(r, len) in input.read_len {
            g.set_read(
                t,
                r,
                ReadInfo {
                    start_mini_end: len as i64,
                    end_mini_start: len as i64,
                    original_support: true,
                },
            );
        }
    }
    let reads_for_isoforms: Vec<u32> = (1..=n_graph_reads).collect();

    let mut prior_read_infos: FxHashMap<(u32, u32, u32), NodeId> = FxHashMap::default();
    // Optional tolerant index for `ISONFORM_LOOKUP_TOL`. The reference's lookup
    // is an exact dict hit; this additionally allows a prediction whose start and
    // end are within `tol` bases. Bucketed by start so a lookup scans a handful
    // of candidates rather than the whole map.
    let lookup_tol: u32 = std::env::var("ISONFORM_LOOKUP_TOL")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    let bucket_of = |st: u32| st / (2 * lookup_tol.max(1) + 1);
    // `ISONFORM_CHAIN=<tol>`: co-linear chaining for intervals the exact lookup
    // misses. See `chain_assign`.
    let chain_tol: u32 = std::env::var("ISONFORM_CHAIN")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    // Largest tolerated spacing disagreement, in bases. 0 = unbounded, i.e. the
    // monotonicity-only chain finding 46 measured.
    let chain_spacing: u32 = std::env::var("ISONFORM_CHAIN_SPACING")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1);
    let index_span = 2 * lookup_tol.max(chain_tol).max(1) + 1;
    let want_index = lookup_tol > 0 || chain_tol > 0;
    let mut tol_index: TolIndex = FxHashMap::default();
    let mut chained_hits = 0usize;
    // old_node -> [(alternative node, its predecessor, the predecessor edge's length)]
    let mut alternative_nodes: FxHashMap<NodeId, Vec<(NodeId, NodeId, i64)>> = FxHashMap::default();
    // The reference only ever reports `len(alt_cyc_nodes.keys())`.
    let mut alt_cyc_keys: FxHashSet<NodeId> = FxHashSet::default();
    let (mut n_cycle, mut n_repeat, mut n_first) = (0usize, 0usize, 0usize);
    // `is_repetitive` on its own, independent of which branch the interval then
    // takes. The `repeat_within_read` branch fires only when the interval is ALSO
    // already known, so a repetitive-but-unknown interval lands in `first_sight`
    // and would not otherwise be counted as a repeat at all.
    let (mut n_rep_seen, mut n_rep_reads, mut n_intervals) = (0usize, 0usize, 0usize);
    // Which ANCHOR PAIR each selected interval came from. `occurrence_hash` is
    // the pair's occurrence list across all reads, so two reads that selected the
    // same pair hash identically --- and a node hit requires exactly that. This
    // is the quantity the coordinate dumps cannot show.
    let mut pair_selectors: FxHashMap<u64, usize> = FxHashMap::default();

    // `seq` in the reference is function-scoped, so it survives both loops. That
    // is finding 9; `None` models the state before any branch has bound it, where
    // the reference would raise.
    let mut leaked_seq: Option<&[u8]> = None;

    for (r_id, intervals_for_read) in input.intervals {
        let r_id = *r_id;
        let read_seq: &[u8] = input
            .reads
            .get(&r_id)
            .map(|v| v.as_slice())
            .ok_or(BuildError::MissingRead(r_id))?;

        let mut previous_node = s;
        let mut previous_end: i64 = 0;
        // The reference stores `{hash: [name]}` but only ever tests membership.
        let mut read_hashes: FxHashSet<u64> = FxHashSet::default();
        // `name` must survive the interval loop: the sink edge is drawn from
        // whatever it held last. A read with no intervals would leave it unbound
        // in the reference too, so an empty interval list is skipped rather than
        // guessed at.
        let mut name: Option<NodeId> = None;
        let mut read_had_repeat = false;
        // One chain per read, computed against everything earlier reads
        // registered for it. Indexed by the read's interval position.
        let chain: Vec<Option<NodeId>> = if chain_tol > 0 {
            let starts: Vec<(u32, u32)> = intervals_for_read
                .iter()
                .map(|iv| (iv.start, iv.end))
                .collect();
            chain_assign(
                &starts,
                r_id,
                &tol_index,
                chain_tol,
                index_span,
                chain_spacing,
            )
        } else {
            Vec::new()
        };

        // `enumerate` rather than a hand-rolled counter: the chain is indexed by
        // the read's interval position, and a manual increment desyncs silently
        // if a `continue` is ever added to this loop.
        for (iv_idx, inter) in intervals_for_read.iter().enumerate() {
            let info_tuple = (r_id, inter.start, inter.end);
            let curr_hash = occurrence_hash(&inter.occurrences);
            let is_repetitive = read_hashes.contains(&curr_hash);
            n_intervals += 1;
            *pair_selectors.entry(curr_hash).or_insert(0) += 1;
            if is_repetitive {
                n_rep_seen += 1;
                read_had_repeat = true;
            }
            let this_len = inter.start as i64 - previous_end;

            if std::env::var_os("ISONFORM_SHARE_STATS").is_some() {
                use std::sync::atomic::Ordering::Relaxed;
                LOOKUPS.fetch_add(1, Relaxed);
                if prior_read_infos.contains_key(&info_tuple) {
                    HITS.fetch_add(1, Relaxed);
                } else {
                    // How *near* was the nearest prediction for this read? A miss
                    // by a few bases means WIS chose almost the same interval and
                    // a tolerance would recover it; a miss by hundreds means WIS
                    // chose a structurally different set and no tolerance helps.
                    let mut best = u32::MAX;
                    for &(r, st, en) in prior_read_infos.keys() {
                        if r != r_id {
                            continue;
                        }
                        let d = st.abs_diff(info_tuple.1) + en.abs_diff(info_tuple.2);
                        best = best.min(d);
                    }
                    let bucket = match best {
                        u32::MAX => 5,
                        0 => 0,
                        1..=5 => 1,
                        6..=20 => 2,
                        21..=100 => 3,
                        _ => 4,
                    };
                    MISS_DIST[bucket].fetch_add(1, Relaxed);
                }
            }
            let mut found = prior_read_infos.get(&info_tuple).copied();
            if found.is_none() && lookup_tol > 0 {
                let b = bucket_of(info_tuple.1);
                'search: for bb in b.saturating_sub(1)..=b + 1 {
                    if let Some(v) = tol_index.get(&(r_id, bb)) {
                        for &(st, en, nd) in v {
                            if st.abs_diff(info_tuple.1) <= lookup_tol
                                && en.abs_diff(info_tuple.2) <= lookup_tol
                            {
                                found = Some(nd);
                                TOL_HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                break 'search;
                            }
                        }
                    }
                }
            }
            // `ISONFORM_TRACE_LOOKUP=<r_id>`: for one read, every interval's
            // lookup key, whether it hit, and the closest coordinate anyone
            // actually registered for this read. This is the direct observation of
            // the rejection --- finding 45's missing link.
            if std::env::var("ISONFORM_TRACE_LOOKUP")
                .ok()
                .and_then(|v| v.parse::<u32>().ok())
                == Some(r_id)
            {
                let mut best = (u32::MAX, 0u32, 0u32);
                for &(r, st, en) in prior_read_infos.keys() {
                    if r != r_id {
                        continue;
                    }
                    let d = st.abs_diff(inter.start) + en.abs_diff(inter.end);
                    if d < best.0 {
                        best = (d, st, en);
                    }
                }
                eprintln!(
                    "lookup r={r_id} iv={iv_idx} want=({},{}) hit={} \
                     nearest_registered=({},{}) dist={}",
                    inter.start,
                    inter.end,
                    found.is_some(),
                    best.1,
                    best.2,
                    if best.0 == u32::MAX { -1 } else { best.0 as i64 }
                );
            }
            if found.is_none() && chain_tol > 0 {
                if let Some(nd) = chain.get(iv_idx).copied().flatten() {
                    found = Some(nd);
                    chained_hits += 1;
                }
            }
            let node = match found {
                Some(prior) if !is_repetitive => {
                    // The interval is known from an earlier read and does not
                    // repeat within this one.
                    if !g.has_edge(previous_node, prior) {
                        leaked_seq = Some(read_seq);
                        g.set_read(
                            prior,
                            r_id,
                            ReadInfo {
                                start_mini_end: inter.start as i64,
                                end_mini_start: inter.end as i64,
                                original_support: true,
                            },
                        );
                        g.set_end_mini_seq(prior, end_mini(read_seq, inter.end, k));
                        g.add_edge(previous_node, prior, this_len);
                        if g.reaches_itself(previous_node) {
                            // `cycle_added`: drop the edge that closed the cycle
                            // and route the incoming edge through a read-private
                            // node instead.
                            alt_cyc_keys.insert(prior);
                            n_cycle += 1;
                            g.remove_edge(previous_node, prior);
                            let alt = g.add_node(NodeKey::Interval {
                                start: inter.start,
                                end: inter.end,
                                r_id,
                            });
                            g.set_end_mini_seq(alt, end_mini(read_seq, inter.end, k));
                            g.replace_reads(
                                alt,
                                r_id,
                                ReadInfo {
                                    start_mini_end: inter.start as i64,
                                    end_mini_start: inter.end as i64,
                                    original_support: true,
                                },
                            );
                            g.add_edge(previous_node, alt, this_len);
                            g.set_edge_support(previous_node, alt, r_id);

                            // FINDING 11, and this is the whole of it: the
                            // reference returns `prior`, NOT `alt`.
                            //
                            // `cycle_added` rebinds its own local `name` to the
                            // new node, and Python strings are passed by value,
                            // so the caller's `name` still holds the old node.
                            // `previous_node = name` at
                            // `GraphGeneration.py:458` therefore continues the
                            // read's path from the node whose incoming edge was
                            // just removed, leaving `alt` with an in-edge and no
                            // out-edge --- a dead end --- and severing the read's
                            // path in two.
                            //
                            // The port returned `alt` at first, which is the
                            // sensible behaviour, and the oracle caught it on
                            // real Drosophila data. Faithful now; the fix is
                            // `BuildOpts::fix_cycle_continuation`.
                            if opts.fix_cycle_continuation {
                                alt
                            } else {
                                prior
                            }
                        } else {
                            g.set_edge_support(previous_node, prior, r_id);
                            prior
                        }
                    } else {
                        let prev_len = g
                            .edge_length(previous_node, prior)
                            .expect("edge exists, has_edge just said so");
                        if (this_len - prev_len).abs() < delta_len {
                            leaked_seq = Some(read_seq);
                            g.set_read(
                                prior,
                                r_id,
                                ReadInfo {
                                    start_mini_end: inter.start as i64,
                                    end_mini_start: inter.end as i64,
                                    original_support: true,
                                },
                            );
                            g.set_end_mini_seq(prior, end_mini(read_seq, inter.end, k));
                            g.push_edge_support(previous_node, prior, r_id);
                            prior
                        } else {
                            // Length disagrees by too much: split off an
                            // alternative node for this read.
                            let old_node = prior;
                            let alt_key = NodeKey::Interval {
                                start: inter.start,
                                end: inter.end,
                                r_id,
                            };
                            match alternative_nodes.get(&old_node) {
                                None => {
                                    // FINDING 9 lives here: the reference reads
                                    // `seq` on this path without binding it.
                                    let seq_here = if opts.fix_stale_seq {
                                        read_seq
                                    } else {
                                        leaked_seq.ok_or(BuildError::UnboundSeq)?
                                    };
                                    // Note the empty list: `alt_info_tuple` is
                                    // NOT appended here, only on the next hit,
                                    // which is why the filter below needs a
                                    // *third* visit before it can match.
                                    alternative_nodes.insert(old_node, Vec::new());
                                    let alt = g.add_node(alt_key);
                                    g.set_end_mini_seq(alt, end_mini(seq_here, inter.end, k));
                                    g.replace_reads(
                                        alt,
                                        r_id,
                                        ReadInfo {
                                            start_mini_end: inter.start as i64,
                                            end_mini_start: inter.end as i64,
                                            original_support: true,
                                        },
                                    );
                                    g.add_edge(previous_node, alt, this_len);
                                    g.set_edge_support(previous_node, alt, r_id);
                                    alt
                                }
                                Some(alts) => {
                                    let matched = alts.iter().any(|&(_, prev, plen)| {
                                        prev == previous_node && (this_len - plen).abs() < delta_len
                                    });
                                    if matched {
                                        // FINDING 8: the reference raises here.
                                        return Err(BuildError::ReferenceWouldCrash {
                                            detail: "known_old_node_action passes a third \
                                                     positional argument to DG.add_edge, which \
                                                     networkx rejects with TypeError",
                                        });
                                    }
                                    // `no_node_to_add_to_action`: record the
                                    // alternative and build a read-private node.
                                    let alt = g.add_node(alt_key);
                                    let seq_here = if opts.fix_stale_seq {
                                        read_seq
                                    } else {
                                        leaked_seq.ok_or(BuildError::UnboundSeq)?
                                    };
                                    alternative_nodes
                                        .get_mut(&old_node)
                                        .expect("present")
                                        .push((alt, previous_node, prev_len));
                                    g.set_end_mini_seq(alt, end_mini(seq_here, inter.end, k));
                                    g.replace_reads(
                                        alt,
                                        r_id,
                                        ReadInfo {
                                            start_mini_end: inter.start as i64,
                                            end_mini_start: inter.end as i64,
                                            original_support: true,
                                        },
                                    );
                                    g.add_edge(previous_node, alt, this_len);
                                    g.set_edge_support(previous_node, alt, r_id);
                                    alt
                                }
                            }
                        }
                    }
                }
                Some(_) => {
                    // Known interval that *does* repeat within this read: give it
                    // a read-private node rather than closing a loop.
                    n_repeat += 1;
                    leaked_seq = Some(read_seq);
                    let n = g.add_node(NodeKey::Interval {
                        start: inter.start,
                        end: inter.end,
                        r_id,
                    });
                    g.set_end_mini_seq(n, end_mini(read_seq, inter.end, k));
                    g.replace_reads(
                        n,
                        r_id,
                        ReadInfo {
                            start_mini_end: inter.start as i64,
                            end_mini_start: inter.end as i64,
                            original_support: true,
                        },
                    );
                    g.add_edge(previous_node, n, this_len);
                    g.set_edge_support(previous_node, n, r_id);
                    n
                }
                None => {
                    // `new_interval_action`: first time this interval is seen.
                    n_first += 1;
                    leaked_seq = Some(read_seq);
                    let n = g.add_node(NodeKey::Interval {
                        start: inter.start,
                        end: inter.end,
                        r_id,
                    });
                    g.set_end_mini_seq(n, end_mini(read_seq, inter.end, k));
                    g.replace_reads(
                        n,
                        r_id,
                        ReadInfo {
                            start_mini_end: inter.start as i64,
                            end_mini_start: inter.end as i64,
                            original_support: true,
                        },
                    );
                    add_prior_read_infos(&mut prior_read_infos, inter, r_id, n, k);
                    if want_index {
                        for triple in inter.occurrences.chunks_exact(3) {
                            let (r, p1, p2) = (triple[0], triple[1], triple[2]);
                            if r > r_id {
                                let st = p1 + k as u32;
                                tol_index
                                    .entry((r, st / index_span))
                                    .or_default()
                                    .push((st, p2, n));
                            }
                        }
                    }
                    g.add_edge(previous_node, n, this_len);
                    g.set_edge_support(previous_node, n, r_id);
                    n
                }
            };

            name = Some(node);
            previous_node = node;
            previous_end = inter.end as i64;
            if !is_repetitive {
                read_hashes.insert(curr_hash);
            }
        }

        if read_had_repeat {
            n_rep_reads += 1;
        }

        // Close the read's path at the sink.
        if let Some(last) = name {
            if !g.has_edge(last, t) {
                g.add_edge(last, t, 0);
                g.set_edge_support(last, t, r_id);
            } else {
                g.push_edge_support(last, t, r_id);
            }
        }
    }

    if chain_tol > 0 {
        eprintln!(
            "chain-stats tol={chain_tol} max_skew={chain_spacing} \
             chained_hits={chained_hits}"
        );
    }
    if std::env::var_os("ISONFORM_NODE_ORIGIN").is_some() {
        eprintln!(
            "repeat-scan intervals={n_intervals} repetitive={n_rep_seen} \
             ({:.2}%) reads_with_a_repeat={n_rep_reads}",
            100.0 * n_rep_seen as f64 / n_intervals.max(1) as f64
        );
        let mut mult: Vec<usize> = pair_selectors.values().copied().collect();
        mult.sort_unstable();
        let m = mult.len().max(1);
        let sum: usize = mult.iter().sum();
        let once = mult.iter().filter(|&&c| c == 1).count();
        let many = mult.iter().filter(|&&c| c >= 100).count();
        eprintln!(
            "pair-agreement distinct_pairs={} selections={sum} \
             median_reads_per_pair={} p90={} max={} \
             pairs_chosen_by_one_read={once} ({:.1}%) pairs_chosen_by_100plus={many}",
            mult.len(),
            mult[m / 2],
            mult[m * 9 / 10],
            mult[m - 1],
            100.0 * once as f64 / m as f64
        );
        eprintln!(
            "node-origin cycle_added={n_cycle} repeat_within_read={n_repeat} \
             first_sight={n_first}"
        );
    }
    let _ = alt_cyc_keys; // the reference only prints its length
    Ok((g, reads_for_isoforms))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: u32, end: u32, occ: &[u32]) -> Interval {
        Interval {
            start,
            end,
            weight: (occ.len() / 3) as u32,
            occurrences: occ.to_vec(),
        }
    }

    /// `(intervals, read sequences)` --- what a `BuildInput` needs the caller to
    /// own.
    type Fixture = (Vec<(u32, Vec<Interval>)>, FxHashMap<u32, Vec<u8>>);

    fn reads(pairs: &[(u32, &str)]) -> FxHashMap<u32, Vec<u8>> {
        pairs
            .iter()
            .map(|(r, s)| (*r, s.as_bytes().to_vec()))
            .collect()
    }

    #[test]
    fn a_single_read_becomes_a_path_from_source_to_sink() {
        let seq = "A".repeat(60);
        let r = reads(&[(1, &seq)]);
        let intervals = vec![(
            1u32,
            vec![iv(10, 20, &[1, 0, 20]), iv(30, 40, &[1, 20, 40])],
        )];
        let input = BuildInput {
            k: 4,
            delta_len: 5,
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 60)],
        };
        let (g, rfi) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();
        assert_eq!(rfi, vec![1]);
        // s, t, and one node per interval.
        assert_eq!(g.node_count(), 4);
        let s = g.lookup(&NodeKey::Source).unwrap();
        let t = g.lookup(&NodeKey::Sink).unwrap();
        let a = g
            .lookup(&NodeKey::Interval {
                start: 10,
                end: 20,
                r_id: 1,
            })
            .unwrap();
        let b = g
            .lookup(&NodeKey::Interval {
                start: 30,
                end: 40,
                r_id: 1,
            })
            .unwrap();
        assert!(g.has_edge(s, a) && g.has_edge(a, b) && g.has_edge(b, t));
        // length = inter.start - previous_end
        assert_eq!(g.edge_length(s, a), Some(10));
        assert_eq!(g.edge_length(a, b), Some(10));
        assert_eq!(g.edge_length(b, t), Some(0));
    }

    #[test]
    fn the_fixed_sink_carries_the_graphs_reads_with_their_own_lengths() {
        // finding 33 fixed. Same input as the bug-compat test below: the file has
        // three reads of lengths 60, 55, 51, but only read 1 has intervals, so
        // the graph has exactly one read. Its own length is 60.
        let seq = "A".repeat(60);
        let r = reads(&[(1, &seq)]);
        let intervals = vec![(1u32, vec![iv(10, 20, &[1, 0, 20])])];
        let input = BuildInput {
            k: 4,
            delta_len: 5,
            intervals: &intervals,
            reads: &r,
            // Deliberately a *different* length for id 1 than the read really
            // has, so a test that passed by coincidence cannot.
            read_len: &[(1, 999), (2, 55), (3, 51)],
        };
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
        let t = g.lookup(&NodeKey::Sink).unwrap();
        let ids: Vec<u32> = g.reads(t).iter().map(|(r, _)| *r).collect();
        assert_eq!(ids, vec![1], "no phantom entries for reads with no path");
        let s = g.lookup(&NodeKey::Source).unwrap();
        assert_eq!(
            g.reads(s).len(),
            g.reads(t).len(),
            "source and sink carry the same reads --- every path runs between them"
        );
        let (_, info) = g.reads(t)[0];
        assert_eq!(
            info.end_mini_start, 60,
            "the graph read's own length, not read_len_dict's 999"
        );
        assert_eq!(info.start_mini_end, 60);
    }

    #[test]
    fn bug_compat_seeds_the_sink_from_read_len_not_from_the_graph() {
        // read_len_dict is built over all reads in `main`, including those the
        // interval filter skipped, so the sink claims reads with no path. This
        // asserts the faithful behaviour, not the desirable one; see PORTING.md.
        let seq = "A".repeat(60);
        let r = reads(&[(1, &seq)]);
        let intervals = vec![(1u32, vec![iv(10, 20, &[1, 0, 20])])];
        let input = BuildInput {
            k: 4,
            delta_len: 5,
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 60), (2, 55), (3, 51)],
        };
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();
        let t = g.lookup(&NodeKey::Sink).unwrap();
        let ids: Vec<u32> = g.reads(t).iter().map(|(r, _)| *r).collect();
        assert_eq!(ids, vec![1, 2, 3]);
        let s = g.lookup(&NodeKey::Source).unwrap();
        assert_eq!(
            g.reads(s).len(),
            1,
            "the source only knows the graph's reads"
        );
    }

    #[test]
    fn a_shared_interval_is_reused_and_the_edge_gains_support() {
        // Read 1 declares an occurrence in read 2 at (p1, p2), so read 2 finds
        // read 1's node via prior_read_infos. Recorded key is (r, p1 + k, p2).
        let seq = "ACGT".repeat(30);
        let r = reads(&[(1, &seq), (2, &seq)]);
        let k = 4usize;
        // read 2's interval must be (p1 + k, p2) = (14, 40) to match.
        let intervals = vec![
            (1u32, vec![iv(10, 20, &[1, 0, 20, 2, 10, 40])]),
            (2u32, vec![iv(14, 40, &[2, 14, 40])]),
        ];
        let input = BuildInput {
            k,
            delta_len: 100, // large, so the length check takes the merge branch
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 120), (2, 120)],
        };
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();
        let shared = g
            .lookup(&NodeKey::Interval {
                start: 10,
                end: 20,
                r_id: 1,
            })
            .unwrap();
        // Read 2 joined read 1's node rather than creating its own.
        assert!(g
            .lookup(&NodeKey::Interval {
                start: 14,
                end: 40,
                r_id: 2
            })
            .is_none());
        let ids: Vec<u32> = g.reads(shared).iter().map(|(r, _)| *r).collect();
        assert_eq!(
            ids,
            vec![1, 2],
            "both reads support the node, in read order"
        );
        let s = g.lookup(&NodeKey::Source).unwrap();
        assert_eq!(g.edge_support(s, shared), Some(&[1u32, 2][..]));
    }

    #[test]
    fn prior_read_infos_only_points_forward() {
        // `if not r <= r_id` — an occurrence in an *earlier* read is not recorded,
        // so it cannot make a later read attach backwards.
        let mut prior = FxHashMap::default();
        let inter = iv(10, 20, &[5, 0, 20, 3, 1, 21, 9, 2, 22]);
        add_prior_read_infos(&mut prior, &inter, 5, 42, 4);
        let mut keys: Vec<_> = prior.keys().copied().collect();
        keys.sort();
        assert_eq!(
            keys,
            vec![(9, 6, 22)],
            "only read 9 is after read 5; start gains k"
        );
    }

    #[test]
    fn end_mini_seq_clamps_like_python_slicing() {
        let seq = b"ACGTACGT";
        assert_eq!(end_mini(seq, 4, 4), b"ACGT");
        assert_eq!(end_mini(seq, 6, 4), b"GT", "short slice, not an error");
        assert_eq!(end_mini(seq, 20, 4), b"", "beyond the end is empty");
    }

    #[test]
    fn the_unbound_seq_path_is_reported_rather_than_guessed() {
        // Reaching the alternative-node branch before any branch has bound `seq`
        // is a NameError in the reference. The port says so instead of inventing
        // a sequence. Constructing it takes a shared interval whose length
        // disagrees by more than delta_len on the very first interval of a read.
        let seq = "ACGT".repeat(30);
        let r = reads(&[(1, &seq), (2, &seq)]);
        let k = 4usize;
        let intervals = vec![
            (1u32, vec![iv(10, 20, &[1, 0, 20, 2, 10, 40])]),
            // matches read 1's node via (10 + 4, 40); s->node already exists
            // with length 10, this_len is 14, so |14 - 10| = 4.
            (2u32, vec![iv(14, 40, &[2, 14, 40])]),
        ];
        let input = BuildInput {
            k,
            delta_len: 2, // 4 >= 2, so the length check splits
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 120), (2, 120)],
        };
        // Read 1 binds `seq` on its own new-interval branch, so by the time read
        // 2 reaches the split there *is* a leaked value --- and it happens to be
        // the right one here, since both reads share a sequence. The point of
        // this test is that the path is exercised and does not panic.
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();
        assert!(
            g.lookup(&NodeKey::Interval {
                start: 14,
                end: 40,
                r_id: 2
            })
            .is_some(),
            "the length disagreement splits off a read-private node"
        );

        // With the fix the graph is the same shape; only which sequence the
        // end-minimizer came from differs, which this corpus cannot show.
        let (g2, _) = generate_graph_from_intervals(
            &input,
            BuildOpts {
                fix_stale_seq: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(g.node_count(), g2.node_count());
        assert_eq!(g.edge_count(), g2.edge_count());
    }

    /// A fixture that actually reaches `cycle_added`, which is finding 11.
    ///
    /// Read 1 builds `s -> n1 -> n2 -> t` and declares occurrences in read 2 for
    /// both intervals, so read 2 finds them via `prior_read_infos` — but in the
    /// *reverse* order. Read 2 therefore tries to add `n2 -> n1`, closing the
    /// cycle `n1 -> n2 -> n1`, and the reference routes around it.
    fn cycle_fixture() -> Fixture {
        // Distinct sequences so a wrong `seq` would be visible too.
        let s1 = "ACGT".repeat(30);
        let s2 = "TTGA".repeat(30);
        let intervals = vec![
            (
                1u32,
                vec![
                    // registers prior (2, 50 + 4, 60) -> "10, 20, 1"
                    iv(10, 20, &[1, 0, 20, 2, 50, 60]),
                    // registers prior (2, 10 + 4, 20) -> "30, 40, 1"
                    iv(30, 40, &[1, 20, 40, 2, 10, 20]),
                ],
            ),
            (
                2u32,
                // Reverse order: n2 first, then n1.
                //
                // The trailing (9, ...) occurrences exist only to make the two
                // intervals' occurrence *tails* differ. Without them both tails
                // are empty, both hash alike, and the second interval is judged
                // `is_repetitive` --- which sends it down a different branch and
                // never reaches the cycle path at all. That is real reference
                // behaviour, and it cost a debugging round to notice, so it is
                // worth the comment: `convert_array_to_hash` drops the
                // occurrence's own triple, so a lone occurrence hashes as the
                // empty tuple and every such interval in a read collides.
                vec![
                    iv(14, 20, &[2, 14, 20, 9, 1, 2]),
                    iv(54, 60, &[2, 54, 60, 9, 3, 4]),
                ],
            ),
        ];
        (intervals, reads(&[(1, &s1), (2, &s2)]))
    }

    #[test]
    fn the_cycle_path_leaves_a_dead_end_and_continues_from_the_old_node() {
        // Finding 11, reproduced. The reference creates a node to route around
        // the cycle but keeps the OLD node as previous_node, because
        // cycle_added's rebinding of `name` is local to the helper.
        let (intervals, r) = cycle_fixture();
        let input = BuildInput {
            k: 4,
            delta_len: 5,
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 120), (2, 120)],
        };
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();

        let n1 = g
            .lookup(&NodeKey::Interval {
                start: 10,
                end: 20,
                r_id: 1,
            })
            .unwrap();
        let n2 = g
            .lookup(&NodeKey::Interval {
                start: 30,
                end: 40,
                r_id: 1,
            })
            .unwrap();
        let alt = g
            .lookup(&NodeKey::Interval {
                start: 54,
                end: 60,
                r_id: 2,
            })
            .expect("the cycle path created a read-private node");
        let t = g.lookup(&NodeKey::Sink).unwrap();

        // The edge that would have closed the cycle is gone, replaced by one
        // into the new node.
        assert!(!g.has_edge(n2, n1), "the cycle-closing edge was removed");
        assert!(
            g.has_edge(n2, alt),
            "and replaced by an edge into the new node"
        );

        // And here is the defect: the new node is a dead end, while the read's
        // path continues out of the node whose in-edge was just removed.
        assert_eq!(
            g.out_degree(alt),
            0,
            "the routed-to node has no outgoing edge"
        );
        assert!(
            g.has_edge(n1, t),
            "read 2's path continues from the OLD node"
        );
        assert!(!g.has_edge(alt, t));
    }

    #[test]
    fn fixing_the_cycle_continuation_moves_the_outgoing_edge() {
        let (intervals, r) = cycle_fixture();
        let input = BuildInput {
            k: 4,
            delta_len: 5,
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 120), (2, 120)],
        };
        let (g, _) = generate_graph_from_intervals(
            &input,
            BuildOpts {
                fix_cycle_continuation: true,
                ..Default::default()
            },
        )
        .unwrap();
        let n1 = g
            .lookup(&NodeKey::Interval {
                start: 10,
                end: 20,
                r_id: 1,
            })
            .unwrap();
        let alt = g
            .lookup(&NodeKey::Interval {
                start: 54,
                end: 60,
                r_id: 2,
            })
            .unwrap();
        let t = g.lookup(&NodeKey::Sink).unwrap();
        assert!(
            g.has_edge(alt, t),
            "the path now continues from the new node"
        );
        assert!(!g.has_edge(n1, t));
        assert_eq!(g.out_degree(alt), 1);
    }

    #[test]
    fn the_stale_sequence_is_observable_when_reads_differ() {
        // Same shape as above, but the two reads have *different* sequences, so
        // the leaked `seq` produces a different end_mini_seq from the fixed one.
        // This is finding 9 made visible in a unit test.
        let s1 = "AAAACCCC".repeat(15);
        let s2 = "GGGGTTTT".repeat(15);
        let r = reads(&[(1, &s1), (2, &s2)]);
        let k = 4usize;
        let intervals = vec![
            (1u32, vec![iv(10, 20, &[1, 0, 20, 2, 10, 40])]),
            (2u32, vec![iv(14, 40, &[2, 14, 40])]),
        ];
        let input = BuildInput {
            k,
            delta_len: 2,
            intervals: &intervals,
            reads: &r,
            read_len: &[(1, 120), (2, 120)],
        };
        let (g_ref, _) = generate_graph_from_intervals(&input, BuildOpts::reference()).unwrap();
        let (g_fix, _) = generate_graph_from_intervals(
            &input,
            BuildOpts {
                fix_stale_seq: true,
                ..Default::default()
            },
        )
        .unwrap();

        let key = NodeKey::Interval {
            start: 14,
            end: 40,
            r_id: 2,
        };
        let a = g_ref.lookup(&key).unwrap();
        let b = g_fix.lookup(&key).unwrap();
        // The reference takes the k-mer from read 1; the fix takes it from read 2.
        assert_eq!(g_ref.end_mini_seq(a), &s1.as_bytes()[40..44]);
        assert_eq!(g_fix.end_mini_seq(b), &s2.as_bytes()[40..44]);
        assert_ne!(g_ref.end_mini_seq(a), g_fix.end_mini_seq(b));
    }
}
