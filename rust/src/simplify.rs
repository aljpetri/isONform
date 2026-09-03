//! Graph simplification: bubble detection.
//!
//! Ports the front of `modules/SimplifyGraph.py` — the part that decides *which*
//! bubbles exist and *what paths* run through them. The mutation that pops them
//! (`prepare_adding_edges`, `remove_edges`, `linearize_bubble`) comes next; it is
//! kept separate because everything here is pure and testable in isolation, and
//! because the popping side calls spoa while none of this does.
//!
//! # Bubble detection, and where the order comes from
//!
//! A candidate bubble is a pair `(start, end)` where `start` has out-degree > 1,
//! `end` has in-degree > 1, `start` precedes `end` in the topological order, and
//! at least two reads support both. The reference builds these with a nested loop
//! over every start × every end.
//!
//! Every list here is in a deterministic order, and each one traces back to
//! something the port has to reproduce:
//!
//! * starts and ends are collected by iterating `topo_nodes_dict.keys()`, i.e.
//!   **in topological order** — so they inherit
//!   [`crate::graph::Graph::topological_sort`]'s networkx-faithful tie-breaking;
//! * the candidate list is starts × ends in that order;
//! * a candidate's read set is `sorted(set(start) & set(end))`, so it *is* sorted
//!   and the reads maps' own order does not survive into it.
//!
//! That last point is worth stating because it is easy to assume otherwise: the
//! support tuples are built from insertion-ordered dicts, but
//! `generate_combinations` immediately sorts the intersection.

use std::collections::BTreeSet;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::graph::{Graph, NodeId, ReadInfo};

/// A node that could start or end a bubble, with the reads supporting it.
///
/// The reference's `(node, tuple(DG.nodes[node]['reads']))`. Support is kept in
/// the node's own order here, matching the reference; it is sorted later, in
/// [`generate_combinations`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Endpoint {
    pub node: NodeId,
    pub support: Vec<u32>,
}

/// A candidate bubble: `(start, end, shared reads)`, the reads sorted.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Combination {
    pub start: NodeId,
    pub end: NodeId,
    pub support: Vec<u32>,
}

/// `find_possible_starts`: nodes with out-degree > 1, in topological order.
pub fn find_possible_starts(g: &Graph, topo: &[NodeId]) -> Vec<Endpoint> {
    topo.iter()
        .filter(|&&n| g.out_degree(n) > 1)
        .map(|&n| Endpoint {
            node: n,
            support: g.reads(n).iter().map(|(r, _)| *r).collect(),
        })
        .collect()
}

/// `find_possible_ends`: nodes with in-degree > 1, in topological order.
///
/// Note the reference's docstring for this one is a copy of the previous
/// function's and says "at least 2 out nodes"; the code uses `in_degree`. The
/// code is right.
pub fn find_possible_ends(g: &Graph, topo: &[NodeId]) -> Vec<Endpoint> {
    topo.iter()
        .filter(|&&n| g.in_degree(n) > 1)
        .map(|&n| Endpoint {
            node: n,
            support: g.reads(n).iter().map(|(r, _)| *r).collect(),
        })
        .collect()
}

/// `generate_combinations`: every start × end pair that could be a bubble.
///
/// `topo_index` is [`crate::graph::Graph::topological_index`] — position in the
/// topological order, which is what the reference compares.
///
/// Quadratic in the number of branch points, exactly as the reference is. Left
/// that way for now: this is a port, and PORTING.md's rule is to measure before
/// optimising.
///
/// It is an obvious candidate — the pairs could be restricted by walking forward
/// from each start instead of considering every end — but **this order does reach
/// results, narrowly.** `SimplifyGraph.py:946` sorts the candidates by bubble
/// span (`topo[end] - topo[start]`), and Python's `sorted` is stable, so the
/// generation order survives as the tie-break among equal-span bubbles, deciding
/// which of them is popped first. So a cheaper enumeration is a behaviour change
/// needing its own commit and a measurement, not a free win. (Checked, rather
/// than assumed: an earlier draft of this comment named the wrong mechanism.)
pub fn generate_combinations(
    starts: &[Endpoint],
    ends: &[Endpoint],
    topo_index: &[u32],
) -> Vec<Combination> {
    let mut out = Vec::new();
    for s in starts {
        // A set per start rather than per pair: the reference rebuilds
        // `set(startnode[1])` inside the inner loop, which is pure waste and not
        // observable.
        let s_set: FxHashSet<u32> = s.support.iter().copied().collect();
        for e in ends {
            if topo_index[s.node as usize] >= topo_index[e.node as usize] {
                continue;
            }
            let mut inter: Vec<u32> = e
                .support
                .iter()
                .copied()
                .filter(|r| s_set.contains(r))
                .collect();
            inter.sort_unstable();
            inter.dedup();
            if inter.len() >= 2 {
                out.push(Combination {
                    start: s.node,
                    end: e.node,
                    support: inter,
                });
            }
        }
    }
    out
}

/// `filter_combinations`: drop candidates already known to be unpoppable.
///
/// **Dead in the reference.** `new_bubble_popping_routine` inlines the same
/// filter as a comprehension at `SimplifyGraph.py:939` and never calls this
/// function; its only appearance is its own definition. Ported anyway because it
/// is three lines and the inlined version needs exactly this behaviour, but
/// nothing calls it here either yet — it is the shape the routine will use.
pub fn filter_combinations(
    combinations: &[Combination],
    not_viable: &FxHashSet<Combination>,
) -> Vec<Combination> {
    combinations
        .iter()
        .filter(|c| !not_viable.contains(c))
        .cloned()
        .collect()
}

/// `get_dist_to_prev`: mean shift in `end_mini_start` over reads supporting both
/// nodes.
///
/// The reference recomputes `avg_dist = summation / (i + 1)` inside the loop and
/// returns the last value, which is the mean — so despite iterating a Python
/// `set`, the result does not depend on iteration order. Returns 0.0 for a empty
/// intersection, as the reference does (its `avg_dist` starts at 0).
pub fn get_dist_to_prev(g: &Graph, prev_node: NodeId, curr_node: NodeId) -> f64 {
    let prev = g.reads(prev_node);
    let curr = g.reads(curr_node);
    let mut sum: i64 = 0;
    let mut n: i64 = 0;
    for (r, ci) in curr {
        if let Some((_, pi)) = prev.iter().find(|(pr, _)| pr == r) {
            sum += ci.end_mini_start - pi.end_mini_start;
            n += 1;
        }
    }
    if n == 0 {
        0.0
    } else {
        sum as f64 / n as f64
    }
}

/// `get_avg_interval_length`: mean `end_mini_start - start_mini_end` over a
/// node's reads, truncated toward zero.
///
/// The reference's `int(summation / i)` on a Python float truncates toward zero,
/// which differs from integer floor division for negative values — and the
/// difference is reachable, because nothing guarantees `end_mini_start >=
/// start_mini_end`. Reproduced with `trunc`, not with `/`.
pub fn get_avg_interval_length(g: &Graph, node: NodeId) -> i64 {
    let reads = g.reads(node);
    if reads.is_empty() {
        return 0;
    }
    let sum: i64 = reads
        .iter()
        .map(|(_, i)| i.end_mini_start - i.start_mini_end)
        .sum();
    (sum as f64 / reads.len() as f64).trunc() as i64
}

/// One path through a bubble: the nodes walked, and the reads that walk it.
///
/// The reference's tuple is `(visited_nodes, tuple(sorted(support)),
/// final_add_support)`. The third element is **dead** — computed on every path and
/// read nowhere — so it is not carried here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BubblePath {
    /// `visited_nodes`: the bubble start, then each interior node. The bubble
    /// *end* is not included, because the reference appends before testing.
    pub nodes: Vec<NodeId>,
    /// Reads whose edge support survives the whole walk, sorted.
    pub support: Vec<u32>,
}

/// `find_paths`: enumerate the paths through a bubble, one per group of reads.
///
/// Takes the lowest unallocated read, walks the graph following edges that read
/// supports, and attributes the whole intersected support of that walk to one
/// path. Repeats until the support is exhausted or a walk fails.
///
/// # Order, and one aliasing subtlety worth recording
///
/// Two orders reach results and both are reproduced:
///
/// * **which read is taken next.** The reference used `set.pop()`; it now takes
///   the minimum (`PORTING.md` finding 12), which is why a `BTreeSet` is the
///   natural container here rather than a hash set.
/// * **which successor is followed.** The walk takes the *first* successor whose
///   `edge_supp` contains the read and stops looking, so adjacency insertion
///   order decides the path when a read supports more than one outgoing edge.
///
/// The subtlety: the reference writes `current_node_support = node_support_left`,
/// which **aliases** rather than copies, and then `.add(read)` — so the read it
/// just removed goes straight back into `node_support_left`. The pop is
/// effectively undone. That is reproduced by simply not removing it. The alias
/// itself cannot leak further, because the first `intersection()` rebinds
/// `current_node_support` to a fresh set, and `next_found` can only be true if at
/// least one intersection happened — so the later `node_support_left -=
/// current_node_support` is never a set minus itself. Checked rather than
/// assumed; it would be a silent infinite loop if it were wrong.
pub fn find_paths(
    g: &Graph,
    start: NodeId,
    end: NodeId,
    support: &[u32],
    marked: &FxHashSet<NodeId>,
) -> Vec<BubblePath> {
    // The reference would raise `UnboundLocalError` on `next_found` if the walk
    // loop never ran. `generate_combinations` guarantees `start != end` (strictly
    // increasing topological index), so no caller can produce it.
    debug_assert_ne!(
        start, end,
        "the reference raises NameError for start == end"
    );
    if start == end {
        return Vec::new();
    }

    let mut node_support_left: BTreeSet<u32> = support.iter().copied().collect();
    let mut out: Vec<BubblePath> = Vec::new();

    while let Some(&read) = node_support_left.iter().next() {
        // Not removed: see the aliasing note above.
        let mut current: FxHashSet<u32> = node_support_left.iter().copied().collect();
        let mut node = start;
        let mut visited: Vec<NodeId> = Vec::new();
        let mut next_found = false;

        while node != end {
            visited.push(node);
            next_found = false;
            if marked.contains(&node) {
                break;
            }
            let mut advanced = None;
            for &nb in g.successors(node) {
                if let Some(supp) = g.edge_support(node, nb) {
                    if supp.contains(&read) {
                        advanced = Some((nb, supp.to_vec()));
                        break;
                    }
                }
            }
            match advanced {
                Some((nb, supp)) => {
                    let keep: FxHashSet<u32> = supp.into_iter().collect();
                    current.retain(|r| keep.contains(r));
                    node = nb;
                    next_found = true;
                }
                None => break,
            }
        }

        if !current.is_empty() && next_found {
            for r in &current {
                node_support_left.remove(r);
            }
            let mut supp: Vec<u32> = current.into_iter().collect();
            supp.sort_unstable();
            out.push(BubblePath {
                nodes: visited,
                support: supp,
            });
        } else {
            // The reference breaks out entirely rather than trying the next read.
            break;
        }
    }
    out
}

/// An edge lifted out of the graph, with the attributes needed to re-add it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LiftedEdge {
    pub u: NodeId,
    pub v: NodeId,
    /// `None` when the edge exists but carries no `length` --- see
    /// `Graph::upsert_edge_support`. The edge must still be lifted and removed
    /// in that case; only a genuinely *absent* edge is skipped.
    pub length: Option<i64>,
    pub support: Vec<u32>,
}

/// What [`remove_edges`] hands to the edge-adding step.
#[derive(Debug, Default)]
pub struct Lifted {
    /// In insertion order, keyed by `(u, v)` — a later write to the same pair
    /// overwrites, as the reference's dict does.
    pub edges: Vec<LiftedEdge>,
    /// `node_distances`: interior node -> distance from the bubble start.
    pub node_distances: Vec<(NodeId, f64)>,
}

impl Lifted {
    fn put_edge(&mut self, e: LiftedEdge) {
        match self.edges.iter_mut().find(|x| x.u == e.u && x.v == e.v) {
            Some(slot) => *slot = e,
            None => self.edges.push(e),
        }
    }

    fn put_dist(&mut self, n: NodeId, d: f64) {
        match self.node_distances.iter_mut().find(|(x, _)| *x == n) {
            Some(slot) => slot.1 = d,
            None => self.node_distances.push((n, d)),
        }
    }

    pub fn dist(&self, n: NodeId) -> Option<f64> {
        self.node_distances
            .iter()
            .find(|(x, _)| *x == n)
            .map(|(_, d)| *d)
    }
}

/// `new_distance_to_start`: where a node's end-minimizer sits in the path's
/// consensus.
///
/// `consensus.find(node_seq)`, i.e. the first byte offset of `node_seq` in
/// `consensus`, or `-1`. The `-1` is not rare — 16% of calls on `sirv_real` —
/// and the caller has a fallback for it.
pub fn new_distance_to_start(consensus: &[u8], node_seq: &[u8]) -> i64 {
    if node_seq.is_empty() {
        // Python's `"abc".find("")` is 0, not -1.
        return 0;
    }
    if node_seq.len() > consensus.len() {
        return -1;
    }
    for i in 0..=(consensus.len() - node_seq.len()) {
        if &consensus[i..i + node_seq.len()] == node_seq {
            return i as i64;
        }
    }
    -1
}

/// `remove_edges`: lift every edge inside the bubble out of the graph, recording
/// each node's distance from the bubble start on the way.
///
/// `paths` is **mutated**: the reference does `path_node_list.pop(0)` on the
/// caller's own list, dropping the bubble start from each path. Reproduced,
/// because the caller reads the shortened lists afterwards.
///
/// `consensus` maps a path's *first* node to that path's consensus sequence —
/// the reference's `consensus_log`.
///
/// # Two reference quirks reproduced here
///
/// * **`prev_node` never advances** (`SimplifyGraph.py:279`). It is set to the
///   bubble start once per path and never reassigned, so `get_dist_to_prev` always
///   measures from the bubble start and `prev_to_start_dist` is always 0 — which
///   makes the `node_distances[prev_node]` branch beside it unreachable. The
///   comment there ("we are still missing the distance of the previous node to
///   s") describes chaining that does not happen. `PORTING.md` finding 13.
/// * **`pnl_start` is computed twice**, identically, either side of the first
///   `edges_to_delete` write. The second computation is dead.
pub fn remove_edges(
    g: &mut Graph,
    bubble_start: NodeId,
    bubble_end: NodeId,
    paths: &mut [BubblePath],
    consensus: &FxHashMap<NodeId, Vec<u8>>,
) -> Lifted {
    let mut lifted = Lifted::default();

    for path in paths.iter_mut() {
        // Destructive, and the caller depends on it.
        if !path.nodes.is_empty() {
            path.nodes.remove(0);
        }
        let pnl_start = path.nodes.first().copied().unwrap_or(bubble_end);

        // The edge from the bubble start into the path.
        //
        // Gated on existence, not on having a `length` --- a re-added edge from
        // an earlier pop carries no `length` (`Graph::upsert_edge_support`) and
        // still has to be lifted and removed here, exactly as `DG.get_edge_data`
        // returns an attribute dict either way and the reference deletes it
        // unconditionally. Checking `edge_length` here instead of `has_edge`
        // silently skipped exactly the edges a prior pop had re-added, leaving
        // them in the graph and making bubble popping fail to converge.
        let prev_node = bubble_start;
        if g.has_edge(prev_node, pnl_start) {
            lifted.put_edge(LiftedEdge {
                u: prev_node,
                v: pnl_start,
                length: g.edge_length(prev_node, pnl_start),
                support: g.edge_support(prev_node, pnl_start).unwrap_or(&[]).to_vec(),
            });
        } else {
            // `DG.get_edge_data` returns None for a missing edge and the
            // reference stores that None, then calls `DG.remove_edge` on it
            // anyway --- which raises NetworkXError. Recording the pair with no
            // attributes would be inventing behaviour, so it is skipped and
            // noted; if a corpus ever reaches it the oracle will disagree and
            // that is the right outcome.
        }

        let n = path.nodes.len();
        for index in 0..n {
            let path_node = path.nodes[index];

            let inter_dist = consensus
                .get(&pnl_start)
                .map(|c| new_distance_to_start(c, g.end_mini_seq(path_node)))
                .unwrap_or(-1);

            let dist = if inter_dist == -1 {
                // prev_node is the bubble start, always --- see the note above,
                // which is why there is no accumulator here.
                get_dist_to_prev(g, prev_node, path_node)
            } else {
                inter_dist as f64
            };
            lifted.put_dist(path_node, dist);

            let target = if index + 1 < n {
                path.nodes[index + 1]
            } else {
                bubble_end
            };
            if g.has_edge(path_node, target) {
                lifted.put_edge(LiftedEdge {
                    u: path_node,
                    v: target,
                    length: g.edge_length(path_node, target),
                    support: g.edge_support(path_node, target).unwrap_or(&[]).to_vec(),
                });
            }
        }
    }

    for e in &lifted.edges {
        g.remove_edge(e.u, e.v);
    }
    lifted
}

/// `merge_two_dicts`: `dict1` wins on shared keys, `dict2`'s extras are appended.
///
/// Not a symmetric merge and not a `dict2.update(dict1)`: the reference copies
/// every entry of `dict1` first, then adds only the keys of `dict2` it has not
/// already seen. So `dict1`'s ordering leads and `dict1`'s values win.
pub fn merge_two_dicts(d1: &[(u32, ReadInfo)], d2: &[(u32, ReadInfo)]) -> Vec<(u32, ReadInfo)> {
    let mut out: Vec<(u32, ReadInfo)> = d1.to_vec();
    for (k, v) in d2 {
        if !out.iter().any(|(x, _)| x == k) {
            out.push((*k, *v));
        }
    }
    out
}

/// `get_next_node`: a path's next node, or the bubble end once it is down to one.
fn get_next_node(path: &[NodeId], bubble_end: NodeId) -> NodeId {
    if path.len() < 2 {
        bubble_end
    } else {
        path[1]
    }
}

/// `compare_by_length`: pick between two candidate next nodes.
///
/// The bubble end always loses to a real node — which is what stops the walk
/// terminating early while one path still has nodes. Otherwise the earlier
/// topological index wins, and a tie goes to `nextnode2` (the reference's
/// `nextnode1 if nn1 < nn2 else nextnode2`).
// The arms deliberately mirror the reference's four-way structure rather than
// being collapsed: two of them return `nextnode1` for different reasons, and
// merging them would hide which reason applies when reading this against
// `SimplifyGraph.py:384`.
#[allow(clippy::if_same_then_else)]
fn compare_by_length(
    nextnode1: NodeId,
    nextnode2: NodeId,
    bubble_end: NodeId,
    nn1: u32,
    nn2: u32,
) -> NodeId {
    if nextnode1 == bubble_end {
        nextnode2
    } else if nextnode2 == bubble_end {
        nextnode1
    } else if nn1 < nn2 {
        nextnode1
    } else {
        nextnode2
    }
}

/// `find_real_nextnode`: which of the two paths advances next.
///
/// An existing edge between the two candidates settles it — that node must come
/// first — and otherwise the topological order does.
fn find_real_nextnode(
    g: &Graph,
    nextnode1: NodeId,
    nextnode2: NodeId,
    bubble_end: NodeId,
    topo_index: &[u32],
) -> NodeId {
    if g.has_edge(nextnode1, nextnode2) {
        nextnode1
    } else if g.has_edge(nextnode2, nextnode1) {
        nextnode2
    } else {
        compare_by_length(
            nextnode1,
            nextnode2,
            bubble_end,
            topo_index[nextnode1 as usize],
            topo_index[nextnode2 as usize],
        )
    }
}

/// `find_connecting_edges`: edges running between the two bubble paths.
///
/// The reference builds a `set` of `(u, v)` tuples of **node name strings**, so
/// its iteration order is `PYTHONHASHSEED`-dependent — and `test_conn_end` then
/// takes `conn_list[0]`. That is `PORTING.md` finding 14: latent, because
/// measured over 36 027 calls on both real corpora `conn_list` never held more
/// than one entry. A `Vec` in a deterministic order is therefore faithful on all
/// real data and strictly better than reproducing a hash-seeded pick.
fn find_connecting_edges(g: &Graph, path1: &[NodeId], path2: &[NodeId]) -> Vec<(NodeId, NodeId)> {
    let mut out = Vec::new();
    for &node in path1 {
        for &onode in path2 {
            let pair = if g.has_edge(node, onode) {
                Some((node, onode))
            } else if g.has_edge(onode, node) {
                Some((onode, node))
            } else {
                None
            };
            if let Some(p) = pair {
                if !out.contains(&p) {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// `prepare_adding_edges`: relink the bubble as a single chain.
///
/// Walks the two paths together, at each step taking whichever node comes first,
/// and emits one edge per step whose support is the union of both paths' support
/// at that point. The result is a linear chain through every node the bubble
/// contained.
///
/// # What the reference does that is worth knowing
///
/// * **New edges carry no `length`.** `DG.add_edge(u, v, edge_supp=value)` names
///   only the support, so a re-added edge has no `length` attribute at all —
///   measured at 97 of 231 edges on one recorded graph. Harmless, since
///   `GraphGeneration.py:402` is the only reader of `length` and it runs earlier,
///   but observable. See [`Graph::upsert_edge_support`].
/// * **Support accumulates rather than replaces** for an edge that survived:
///   old support values not already present are appended.
/// * `nx.set_node_attributes` is called **inside** the loop, re-applying every
///   accumulated entry each time round. Wasteful, not observable; done once here.
/// * `additional_supp` and `curr_node` are assigned only inside the two branches,
///   so a step that matched neither would read stale values from the previous
///   iteration. That cannot happen: the only node outside both paths is the
///   bubble end, and `compare_by_length` makes the bubble end lose to any real
///   node, while `find_real_nextnode`'s edge tests cannot pick it either — an edge
///   from the bubble end back into the bubble would be a cycle. Checked, because
///   the same pattern *is* a live defect elsewhere (finding 9).
pub fn prepare_adding_edges(
    g: &mut Graph,
    lifted: &Lifted,
    bubble_start: NodeId,
    bubble_end: NodeId,
    paths: &[BubblePath],
    topo_index: &[u32],
) -> Vec<NodeId> {
    let mut path1: Vec<NodeId> = paths.first().map(|p| p.nodes.clone()).unwrap_or_default();
    let mut path2: Vec<NodeId> = paths.get(1).map(|p| p.nodes.clone()).unwrap_or_default();
    let conn_edges = find_connecting_edges(g, &path1, &path2);

    // Emitted in order, keyed by (u, v) with a later write overwriting.
    let mut edge_params: Vec<((NodeId, NodeId), Vec<u32>)> = Vec::new();
    let put = |edge_params: &mut Vec<((NodeId, NodeId), Vec<u32>)>,
               key: (NodeId, NodeId),
               val: Vec<u32>| {
        match edge_params.iter_mut().find(|(k, _)| *k == key) {
            Some(slot) => slot.1 = val,
            None => edge_params.push((key, val)),
        }
    };

    // The reference's `while path1 or path2` has no bound. If a step ever fails
    // to shorten either path the loop spins forever, which is how a `retain`
    // where the reference has `remove` first showed up. Bounded here so that
    // failure mode is a message instead of a hang; the bound is generous enough
    // that no correct run can reach it.
    let max_steps = path1.len() + path2.len() + 2;
    let mut steps = 0usize;

    let mut linearization_order = vec![bubble_start];
    let mut prevnode1 = bubble_start;
    let mut prevnode2 = bubble_start;
    let mut nextnode1 = path1.first().copied().unwrap_or(bubble_end);
    let mut nextnode2 = path2.first().copied().unwrap_or(bubble_end);
    let mut prevnode = bubble_start;

    // The reference indexes `edges_to_delete[prevnode1, nextnode1]` directly, so
    // a missing pair is a `KeyError`. It should not be missing --- `remove_edges`
    // records exactly the pairs this walk asks for --- so an empty default here
    // stands in for "cannot happen", and if it ever does the simplification
    // oracle will disagree, which is the right way to find out.
    let lifted_supp = |u: NodeId, v: NodeId| -> Vec<u32> {
        lifted
            .edges
            .iter()
            .find(|e| e.u == u && e.v == v)
            .map(|e| e.support.clone())
            .unwrap_or_default()
    };

    while !path1.is_empty() || !path2.is_empty() {
        steps += 1;
        assert!(
            steps <= max_steps,
            "prepare_adding_edges did not converge: a step left both paths unchanged"
        );
        let overall_nextnode = find_real_nextnode(g, nextnode1, nextnode2, bubble_end, topo_index);

        // `test_conn_end`: connecting edges that end at the node being added.
        //
        // **Deliberate divergence, and the reference has no answer to match.**
        // `find_connecting_edges` returns a Python `set` of `(name, name)` string
        // tuples, so its iteration order is `PYTHONHASHSEED`-dependent, and
        // `prepare_adding_edges` takes `conn_list[0]`. When two connecting edges
        // end at the same node the pick is therefore a coin flip: replaying one
        // recorded graph under 8 seeds produces 2 distinct outputs. Finding 14
        // called this latent on the evidence then available; a holdout corpus
        // reached it, so it is live --- rare (1 of 140 non-empty results across
        // 19 831 calls) but real.
        //
        // Ordering lexicographically by the source node's *name*, which is the
        // same choice already approved for minimizer selection (finding 1) and is
        // what the reference's own node names sort by. On the one observed
        // multi-candidate case this agrees with `PYTHONHASHSEED=0`, i.e. with
        // every recorded dump --- but that is a sample of one and is not evidence
        // that lexicographic equals seed 0 in general.
        //
        // No env var restores "the reference path" here, because there is no
        // reference path to restore: the behaviour being diverged from is a
        // coin flip, not a defined order.
        let mut conn: Vec<&(NodeId, NodeId)> = conn_edges
            .iter()
            .filter(|(_, v)| *v == overall_nextnode)
            .collect();
        if conn.len() > 1 {
            conn.sort_by_cached_key(|(u, _)| g.key(*u).to_string());
        }

        let new_edge_supp1 = lifted_supp(prevnode1, nextnode1);
        let new_edge_supp2 = lifted_supp(prevnode2, nextnode2);
        let mut full_edge_supp = new_edge_supp1.clone();
        full_edge_supp.extend_from_slice(&new_edge_supp2);
        if let Some(&&(cu, _)) = conn.first() {
            if let Some(extra) = g.edge_support(cu, overall_nextnode) {
                full_edge_supp.extend_from_slice(extra);
            }
        }

        let curr_node;
        if path1.contains(&overall_nextnode) {
            nextnode1 = get_next_node(&path1, bubble_end);
            let additional = additional_node_support(
                g,
                &new_edge_supp2,
                lifted,
                overall_nextnode,
                prevnode2,
                bubble_start,
            );
            remove_first(&mut path1, overall_nextnode);
            prevnode1 = overall_nextnode;
            curr_node = prevnode1;
            let merged = merge_two_dicts(&additional, g.reads(overall_nextnode));
            set_reads(g, overall_nextnode, merged);
        } else if path2.contains(&overall_nextnode) {
            nextnode2 = get_next_node(&path2, bubble_end);
            let additional = additional_node_support(
                g,
                &new_edge_supp1,
                lifted,
                overall_nextnode,
                prevnode1,
                bubble_start,
            );
            remove_first(&mut path2, overall_nextnode);
            prevnode2 = overall_nextnode;
            curr_node = prevnode2;
            let merged = merge_two_dicts(&additional, g.reads(overall_nextnode));
            set_reads(g, overall_nextnode, merged);
        } else {
            // Unreachable; see the note on this function.
            debug_assert!(false, "overall_nextnode belongs to neither path");
            break;
        }

        linearization_order.push(overall_nextnode);
        put(
            &mut edge_params,
            (prevnode, overall_nextnode),
            full_edge_supp,
        );
        prevnode = curr_node;
    }

    // The final hop into the bubble end.
    let mut tail = lifted_supp(prevnode1, bubble_end);
    tail.extend_from_slice(&lifted_supp(prevnode2, bubble_end));
    put(&mut edge_params, (prevnode, bubble_end), tail);
    linearization_order.push(bubble_end);

    for ((u, v), mut support) in edge_params {
        if let Some(old) = g.edge_support(u, v) {
            // The edge survived: keep its support values that the new list lacks.
            let old: Vec<u32> = old.to_vec();
            for o in old {
                if !support.contains(&o) {
                    support.push(o);
                }
            }
        }
        g.upsert_edge_support(u, v, support);
    }

    linearization_order
}

/// Python's `list.remove(x)`: drop the **first** occurrence only.
///
/// `Vec::retain` would drop every occurrence, which is not the same thing and is
/// reachable — a path can revisit a node. Getting this wrong made
/// `prepare_adding_edges` fail to converge, which showed up as a hanging test
/// rather than a wrong answer.
fn remove_first(v: &mut Vec<NodeId>, x: NodeId) {
    if let Some(i) = v.iter().position(|&n| n == x) {
        v.remove(i);
    }
}

/// `additional_node_support`: invent positions for reads that reach a node from
/// the *other* bubble path.
///
/// Only reads absent from the node get an entry, and only when the other path's
/// previous node knows them. `original_support` is `False` on every entry, which
/// is how a read that was never really there is marked.
fn additional_node_support(
    g: &Graph,
    new_support: &[u32],
    lifted: &Lifted,
    node: NodeId,
    other_prevnode: NodeId,
    bubble_start: NodeId,
) -> Vec<(u32, ReadInfo)> {
    let other_dist = if other_prevnode == bubble_start {
        0.0
    } else {
        lifted.dist(other_prevnode).unwrap_or(0.0)
    };
    let this_dist = lifted.dist(node).unwrap_or(0.0);
    let avg_len = get_avg_interval_length(g, node);

    let mut out = Vec::new();
    for &r_id in new_support {
        if g.reads(node).iter().any(|(r, _)| *r == r_id) {
            continue;
        }
        if let Some((_, info)) = g.reads(other_prevnode).iter().find(|(r, _)| *r == r_id) {
            let relative = (this_dist - other_dist).trunc() as i64;
            let newend = info.end_mini_start + relative;
            let newstart = newend - avg_len;
            // NOT clamped. `int()` in the reference truncates toward zero but does
            // not floor at zero, and both of these genuinely come out negative on
            // real data --- 57 nodes across 16 recorded Drosophila
            // simplifications. Clamping here was a divergence the oracle caught.
            out.push((
                r_id,
                ReadInfo {
                    start_mini_end: newstart,
                    end_mini_start: newend,
                    original_support: false,
                },
            ));
        }
    }
    out
}

fn set_reads(g: &mut Graph, node: NodeId, reads: Vec<(u32, ReadInfo)>) {
    // `nx.set_node_attributes(DG, {node: dict}, "reads")` replaces the whole map.
    let mut it = reads.into_iter();
    match it.next() {
        Some((r, i)) => {
            g.replace_reads(node, r, i);
            for (r, i) in it {
                g.set_read(node, r, i);
            }
        }
        None => g.clear_reads(node),
    }
}

/// `linearize_bubble`: lift the bubble's edges out, then relink it as a chain.
///
/// The reference's wrapper, in full: `remove_edges` followed by
/// `prepare_adding_edges`. Returns the linearisation order, which the reference
/// computes and discards.
///
/// `paths` is mutated — `remove_edges` strips the bubble start from each path.
pub fn linearize_bubble(
    g: &mut Graph,
    bubble_start: NodeId,
    bubble_end: NodeId,
    paths: &mut [BubblePath],
    consensus: &FxHashMap<NodeId, Vec<u8>>,
    topo_index: &[u32],
) -> Vec<NodeId> {
    let lifted = remove_edges(g, bubble_start, bubble_end, paths, consensus);
    prepare_adding_edges(g, &lifted, bubble_start, bubble_end, paths, topo_index)
}

// ===========================================================================
// Consensus and the poppability decision
// ===========================================================================

/// The two external tools the poppability decision needs.
///
/// Kept behind a trait because they are the last pieces of this stage that are
/// not pure Rust: spoa is a subprocess in the reference (`spoa <reads> -l 0 -r 0
/// -g -2`, the invocation isONcorrect's `poa.rs` already reproduces via spoars),
/// and the alignment is `parasail.sg_trace_scan_16`, which isONcorrect's
/// `parasail.rs` reproduces exactly. Wiring those in is mechanical; separating
/// them means everything else here is testable now.
pub trait Consensus {
    /// spoa over the given sequences, returning the consensus.
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8>;

    /// `parasail.sg_trace_scan_16(s1, s2, match=2, mismatch=-8, open=12, ext=1)`,
    /// returning the CIGAR as `(length, op)` pairs.
    ///
    /// Note the mismatch penalty: −8 here, where `modules/consensus.py`'s own
    /// default is −2 and isONcorrect uses −8. The call site overrides it.
    fn align(&mut self, s1: &[u8], s2: &[u8]) -> Vec<(u32, u8)>;
}

/// The real thing: spoa via [`crate::poa`] and parasail via [`crate::parasail`].
///
/// Both are pure Rust. [`crate::poa::consensus`] wraps `spoars`, a
/// reimplementation of spoa, and [`crate::parasail::semiglobal`] reproduces
/// `parasail.sg_trace_scan_16` exactly — both carried across from the isONcorrect
/// port with their verification.
///
/// **One of those verifications has not been repeated here.** `poa.rs` matched the
/// `spoa` binary on 505 of 505 real isONcorrect correction intervals; isONform
/// calls spoa on different sequences from a different stage (bubble path
/// consensus, not correction intervals), so that number does not transfer. Until
/// it is re-measured, a simplification oracle that disagrees could be either
/// side. The parasail side needs no such caveat: it is exact by construction and
/// checked at the CIGAR level.
#[derive(Debug, Default)]
pub struct SpoaParasail;

impl Consensus for SpoaParasail {
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8> {
        // The reference writes a fasta and shells out to
        // `spoa <file> -l 0 -r 0 -g -2`; `poa::consensus` is that invocation.
        crate::poa::consensus(seqs)
            .map(|s| s.into_bytes())
            .unwrap_or_default()
    }

    fn align(&mut self, s1: &[u8], s2: &[u8]) -> Vec<(u32, u8)> {
        // `Scoring::BUBBLE` is the poppability scoring: match 2, mismatch -8,
        // open 12, ext 1. Not MERGE (-2), which belongs to isoform merging.
        // A bubble's two paths are delimited by its shared start and end nodes,
        // so they have common endpoints --- a *global* alignment problem. The
        // reference nevertheless calls `parasail.sg_trace_scan_16`, making both
        // end gaps free, which lets the aligner shed exactly the end
        // disagreement that `mergeable_start`/`mergeable_end` exist to catch.
        // `ISONFORM_BUBBLE_GLOBAL=1` charges the end gaps instead. Off by
        // default: it changes which bubbles pop. See `parasail::global`.
        let a = if bubble_global() {
            crate::parasail::global(s1, s2, crate::parasail::Scoring::BUBBLE)
        } else {
            crate::wfa::enabled()
                .then(|| crate::wfa::semiglobal(s1, s2, crate::parasail::Scoring::BUBBLE))
                .flatten()
                .unwrap_or_else(|| {
                    crate::parasail::semiglobal(s1, s2, crate::parasail::Scoring::BUBBLE)
                })
        };
        a.ops.iter().map(|&(len, op)| (len as u32, op)).collect()
    }
}

/// `ISONFORM_BUBBLE_GLOBAL=1`: charge the end gaps when aligning a bubble's two
/// paths. Read once. See [`crate::parasail::global`].
pub fn bubble_global() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("ISONFORM_BUBBLE_GLOBAL").is_some())
}

/// One read's span across a bubble: `(r_id, start, end)`.
///
/// Positions are signed for the same reason [`ReadInfo`]'s are: they come
/// straight from `end_mini_start`, which `additional_node_support` can leave
/// negative.
pub type ConsensusAttr = (u32, i64, i64);

/// Python's `seq[a:b]`, including its negative-index behaviour.
///
/// This is not pedantry. `generate_consensus_path` clamps `pos1` to zero in its
/// spoa branch (`if pos1 < 0: pos1 = 0`) but **not** in its two-read branch, and
/// `align_bubble_nodes` does not clamp its single-read branch either — so a
/// negative position there indexes from the *end* of the read rather than being
/// treated as zero. Since positions really do go negative after a bubble is
/// popped (see [`ReadInfo`]), the difference is reachable, and treating negatives
/// as zero would silently take a different subsequence.
fn py_slice(seq: &[u8], a: i64, b: i64) -> &[u8] {
    let n = seq.len() as i64;
    let norm = |i: i64| -> usize {
        let i = if i < 0 { i + n } else { i };
        i.clamp(0, n) as usize
    };
    let (lo, hi) = (norm(a), norm(b));
    if lo >= hi {
        &[]
    } else {
        &seq[lo..hi]
    }
}

/// `get_consensus_positions`: where each shared read enters and leaves a bubble.
///
/// Both ends come from `end_mini_start`, so a "position" here is the end of the
/// read's minimizer at that node, not the interval start.
///
/// The reference indexes the reads maps directly, so a read in `shared_reads` but
/// absent from either node is a `KeyError`. `generate_combinations` builds the
/// support as the intersection of the two nodes' reads, so it cannot happen from
/// that path; reads missing here are skipped rather than invented.
pub fn get_consensus_positions(
    g: &Graph,
    bubble_start: NodeId,
    bubble_end: NodeId,
    shared_reads: &[u32],
) -> Vec<ConsensusAttr> {
    let mut out = Vec::new();
    for &r_id in shared_reads {
        let start = g.reads(bubble_start).iter().find(|(r, _)| *r == r_id);
        let end = g.reads(bubble_end).iter().find(|(r, _)| *r == r_id);
        if let (Some((_, s)), Some((_, e))) = (start, end) {
            out.push((r_id, s.end_mini_start, e.end_mini_start));
        }
    }
    out
}

/// `generate_consensus_path`: one consensus sequence for a bubble path.
///
/// Two quite different routes depending on how many reads support the path:
///
/// * **more than two** — cut each read's span out, write the ones at least `k`
///   long to a fasta, and run spoa over them. If none reached `k`, return
///   `"X" * max_len` instead, a placeholder the caller then rejects for being
///   too short.
/// * **exactly two** — no spoa at all: take whichever read's span is longer and
///   use that subsequence verbatim. This is the case a port can verify without
///   reproducing spoa, which is why it is worth calling out.
///
/// Called only with two or more attributes (`align_bubble_nodes` handles one
/// separately), so the two-attribute arm cannot index past the end.
pub fn generate_consensus_path<C: Consensus>(
    engine: &mut C,
    reads: &FxHashMap<u32, Vec<u8>>,
    attrs: &[ConsensusAttr],
    k: usize,
) -> Vec<u8> {
    if attrs.len() > 2 {
        let mut to_spoa: Vec<Vec<u8>> = Vec::new();
        let mut max_len = 0usize;
        for &(q_id, pos1, pos2) in attrs {
            let Some(read) = reads.get(&q_id) else {
                continue;
            };
            let mut pos1 = pos1;
            let mut pos2 = pos2;
            if pos2 == 0 {
                // `pos2 == 0` means "no end recorded"; fall back to the read's
                // own end. Can go negative for a read shorter than k, which the
                // slice below then clamps.
                pos2 = read.len() as i64 - k as i64;
            }
            let seq: Vec<u8> = if pos1 < pos2 + k as i64 {
                // This branch *does* clamp, explicitly, so py_slice's
                // negative-index behaviour is deliberately bypassed here.
                if pos1 < 0 {
                    pos1 = 0;
                }
                py_slice(read, pos1, pos2 + k as i64).to_vec()
            } else {
                Vec::new()
            };
            // `endseqlist` in the reference is built here and never read.
            if seq.len() < k {
                if seq.len() > max_len {
                    max_len = seq.len();
                } else if seq.is_empty() {
                    max_len = 1;
                }
            } else {
                to_spoa.push(seq);
            }
        }
        if to_spoa.is_empty() {
            return vec![b'X'; max_len];
        }
        let refs: Vec<&[u8]> = to_spoa.iter().map(|v| v.as_slice()).collect();
        return engine.spoa(&refs);
    }

    // Exactly two: the longer span wins, no spoa.
    let (f_id, fstart, fend) = attrs[0];
    let (e_id, estart, eend) = attrs[1];
    let fdist = fend - fstart;
    let edist = eend - estart;
    let (id, lo, hi) = if fdist > edist {
        (f_id, fstart, fend + k as i64)
    } else {
        (e_id, estart, eend + k as i64)
    };
    // FINDING 21 lives in the line the reference writes next to this one: it
    // files the chosen span under `seq_infos[f_id]` in *both* branches, so when
    // e_id wins, e_id's coordinates are recorded against f_id and e_id gets no
    // entry at all. Harmless, and not reproduced, because `seq_infos` is never
    // read --- see the note on `align_bubble_nodes`.
    match reads.get(&id) {
        // No clamp here: the reference slices directly, so a negative start
        // counts from the end of the read.
        Some(read) => py_slice(read, lo, hi).to_vec(),
        None => Vec::new(),
    }
}

/// `parse_cigar_diversity`: are two consensus sequences close enough to merge?
///
/// Rejects outright if any single non-match CIGAR run is longer than `delta_len`
/// — that is a structural difference, an exon, not an error. Otherwise compares
/// the total mismatch fraction against `max(delta_perc * len, delta_len) / len`.
///
/// `delta_perc` is the hard-coded 0.20 at `SimplifyGraph.py:653`, **not** the
/// `--delta` flag, which never reaches this stage.
pub fn parse_cigar_diversity(cigar: &[(u32, u8)], delta_perc: f64, delta_len: i64) -> bool {
    let mut mismatch = 0i64;
    let mut alignment_len = 0i64;
    for &(len, op) in cigar {
        alignment_len += len as i64;
        if op != b'=' && op != b'M' {
            mismatch += len as i64;
            if len as i64 > delta_len {
                return false;
            }
        }
    }
    if alignment_len == 0 {
        // The reference divides by this unguarded. An empty CIGAR needs an empty
        // input, which the length gate above already excluded.
        return false;
    }
    let diversity = mismatch as f64 / alignment_len as f64;
    let max_bp_diff = (delta_perc * alignment_len as f64).max(delta_len as f64);
    diversity <= max_bp_diff / alignment_len as f64
}

/// The hard-coded diversity threshold at `SimplifyGraph.py:653`.
pub const DIVERSITY_DELTA: f64 = 0.20;

/// `align_bubble_nodes`: build both paths' consensus and decide poppability.
///
/// # Three of the reference's five return values are dead
///
/// It returns `(is_poppable, cigar, seq_infos, consensus_log, spoa_count)`. The
/// caller unpacks all five and uses **two**: `is_poppable` and `consensus_log`.
/// `cigar` and `seq_infos` are never read, and `spoa_count` is only ever
/// incremented and passed back in — never printed, never compared. The cached
/// `multi_consensuses` tuple stores `seq_infos` too, and the cache-hit path reads
/// only elements 0 and 2. So [`AlignVerdict`] carries just the two live values.
///
/// That is what makes finding 21 harmless: `generate_consensus_path` mis-files a
/// `seq_infos` entry, and nothing reads it.
///
/// # The decision
///
/// Purely on the two consensus lengths, then on the alignment:
///
/// * both longer than `delta_len`: reject if their lengths differ by 20% or more,
///   otherwise align and ask [`parse_cigar_diversity`];
/// * both shorter: pop, unconditionally;
/// * one of each: pop when the difference is under `delta_len`.
///
/// Note the strict comparisons leave `shorter_len == delta_len` falling to the
/// third case rather than the second; it gives the same answer, so it is not a
/// defect, only surprising.
pub struct RealAligner<'a, C: Consensus> {
    pub engine: C,
    pub reads: &'a FxHashMap<u32, Vec<u8>>,
    pub k: usize,
    pub delta_len: i64,
    /// `multi_consensuses`: a memo used **only** for megabubbles, keyed by
    /// `(start, end, the path's read set)`. Written only when `is_megabubble`,
    /// read only when `is_megabubble`.
    cache: FxHashMap<(NodeId, NodeId, Vec<u32>), Vec<u8>>,
    graph_snapshot: Option<()>,
    /// `ISONFORM_TRACE_DECIDE=<comma-separated read ids>`: dump the inputs, both
    /// consensus sequences and the verdict for the bubble whose path has exactly
    /// that read set.
    ///
    /// The companion to `ISONFORM_TRACE_POPS`, which says *which* bubbles popped;
    /// this says *why* one did. Both earned a place: this is what showed findings
    /// 24 and 25 --- in each case the two sides computed byte-identical consensus
    /// sequences and disagreed anyway, which is what moved the search out of
    /// `simplify.rs` and into the aligner. Read once at construction, not per
    /// call, because `align` runs thousands of times per graph.
    trace_decide: Option<String>,
}

impl<'a, C: Consensus> RealAligner<'a, C> {
    pub fn new(engine: C, reads: &'a FxHashMap<u32, Vec<u8>>, k: usize, delta_len: i64) -> Self {
        Self {
            engine,
            reads,
            k,
            delta_len,
            cache: FxHashMap::default(),
            graph_snapshot: None,
            trace_decide: std::env::var("ISONFORM_TRACE_DECIDE").ok(),
        }
    }

    /// The consensus for one path, with the megabubble memo applied.
    fn consensus_for(
        &mut self,
        bubble: (NodeId, NodeId),
        attrs: &[ConsensusAttr],
        is_megabubble: bool,
    ) -> Vec<u8> {
        let mut key_reads: Vec<u32> = attrs.iter().map(|a| a.0).collect();
        key_reads.sort_unstable();
        key_reads.dedup();
        let key = (bubble.0, bubble.1, key_reads);

        if is_megabubble {
            if let Some(hit) = self.cache.get(&key) {
                return hit.clone();
            }
        }
        let con = if attrs.len() > 1 {
            let c = generate_consensus_path(&mut self.engine, self.reads, attrs, self.k);
            if is_megabubble {
                self.cache.insert(key, c.clone());
            }
            // A consensus under 3 bases is treated as no consensus at all.
            if c.len() < 3 {
                Vec::new()
            } else {
                c
            }
        } else if attrs.len() == 1 {
            let (q_id, pos1, pos2) = attrs[0];
            let too_short = (pos2 - pos1).abs() < 3;
            if too_short || pos2 < pos1 {
                Vec::new()
            } else {
                match self.reads.get(&q_id) {
                    // Unclamped, as the reference is.
                    Some(read) => py_slice(read, pos1, pos2 + self.k as i64).to_vec(),
                    None => Vec::new(),
                }
            }
        } else {
            Vec::new()
        };
        con
    }

    /// The length-and-alignment gate, given both consensus sequences.
    fn decide(&mut self, c1: &[u8], c2: &[u8]) -> bool {
        let (s1, s2) = (c1.len() as i64, c2.len() as i64);
        let (longer, shorter) = if s1 > s2 { (s1, s2) } else { (s2, s1) };
        if shorter > self.delta_len && longer > self.delta_len {
            if longer == 0 {
                return false;
            }
            if ((longer - shorter) as f64 / longer as f64) < DIVERSITY_DELTA {
                let cigar = self.engine.align(c1, c2);
                parse_cigar_diversity(&cigar, DIVERSITY_DELTA, self.delta_len)
            } else {
                false
            }
        } else if shorter < self.delta_len && longer < self.delta_len {
            true
        } else {
            (longer - shorter) < self.delta_len
        }
    }
}

impl<C: Consensus> BubbleAligner for RealAligner<'_, C> {
    fn align(&mut self, req: &AlignRequest) -> AlignVerdict {
        let _ = self.graph_snapshot;
        // `consensus_infos` is a dict keyed by the path-defining node, so two
        // paths sharing that node would collapse to one entry and the reference
        // would then IndexError on `consensus_list[1]`. It cannot happen: equal
        // path nodes mean the paths share their first interior node, which the
        // caller's intersection check rejects before getting here.
        debug_assert_ne!(req.path_nodes[0], req.path_nodes[1]);

        let bubble = (req.bubble_start, req.bubble_end);
        let mut consensus = FxHashMap::default();
        let mut list: [Vec<u8>; 2] = [Vec::new(), Vec::new()];
        for (slot, (attrs, &node)) in list
            .iter_mut()
            .zip(req.attrs.iter().zip(req.path_nodes.iter()))
        {
            let con = self.consensus_for(bubble, attrs, req.is_megabubble);
            consensus.insert(node, con.clone());
            *slot = con;
        }
        let poppable = self.decide(&list[0], &list[1]);
        if let Some(watch) = self.trace_decide.as_deref() {
            let sig = |a: &[ConsensusAttr]| {
                let mut v: Vec<u32> = a.iter().map(|x| x.0).collect();
                v.sort_unstable();
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            let (s0, s1) = (sig(&req.attrs[0]), sig(&req.attrs[1]));
            // `ISONFORM_TRACE_DECIDE=*` dumps *every* verdict in one line, for
            // diffing against the reference's own sequence of decisions
            // (`dump_reference.py --record-decide`). A read signature watches one
            // bubble instead, which is what findings 24/25 used.
            if watch == "*" {
                let mut reads: Vec<u32> = req.path_support[0]
                    .iter()
                    .chain(req.path_support[1].iter())
                    .copied()
                    .collect();
                reads.sort_unstable();
                reads.dedup();
                eprintln!(
                    "DECIDE\t{}\t{}\t{}\t{}",
                    req.bubble_start,
                    req.bubble_end,
                    reads
                        .iter()
                        .map(|r| r.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    poppable
                );
            }
            if s0 == watch || s1 == watch {
                eprintln!(
                    "DECIDE mega={} p0=[{s0}] p1=[{s1}] len0={} len1={} poppable={}",
                    req.is_megabubble,
                    list[0].len(),
                    list[1].len(),
                    poppable
                );
                eprintln!("  attrs0={:?}", req.attrs[0]);
                eprintln!("  attrs1={:?}", req.attrs[1]);
                eprintln!("  c0={}", String::from_utf8_lossy(&list[0]));
                eprintln!("  c1={}", String::from_utf8_lossy(&list[1]));
            }
        }
        AlignVerdict {
            poppable,
            consensus,
        }
    }
}

// ===========================================================================
// The popping driver
// ===========================================================================

/// What the driver asks of the consensus/alignment step.
///
/// The reference's `align_bubble_nodes` takes the two paths' read positions,
/// builds a consensus for each with spoa, aligns them, and decides from the CIGAR
/// whether their differences are small enough to merge. All of that is behind
/// this trait so the driver can be ported and tested without spoa — see
/// [`NeverPop`] and [`AlwaysPop`].
#[derive(Debug, Clone)]
pub struct AlignRequest {
    pub bubble_start: NodeId,
    pub bubble_end: NodeId,
    /// The node that defines each path (`pathnode1`, `pathnode2`), and the reads
    /// supporting that path.
    pub path_nodes: [NodeId; 2],
    pub path_support: [Vec<u32>; 2],
    /// `get_consensus_positions` for each path: where each supporting read enters
    /// and leaves the bubble. This is what the consensus is built from.
    pub attrs: [Vec<ConsensusAttr>; 2],
    /// True when the bubble had more than two paths before filtering.
    pub is_megabubble: bool,
}

/// The verdict, plus the per-path consensus `linearize_bubble` needs.
#[derive(Debug, Clone, Default)]
pub struct AlignVerdict {
    pub poppable: bool,
    /// `consensus_info_log`: path-defining node -> consensus sequence. Keyed the
    /// way `remove_edges` looks it up.
    pub consensus: FxHashMap<NodeId, Vec<u8>>,
}

pub trait BubbleAligner {
    fn align(&mut self, req: &AlignRequest) -> AlignVerdict;
}

/// Refuses every bubble. Lets the driver's state machine be exercised on its own:
/// with nothing poppable, what is left is exactly the iteration, marking and
/// not-viable bookkeeping.
#[derive(Debug, Default)]
pub struct NeverPop;

impl BubbleAligner for NeverPop {
    fn align(&mut self, _req: &AlignRequest) -> AlignVerdict {
        AlignVerdict::default()
    }
}

/// Pops every bubble it is offered, with an empty consensus.
///
/// An empty consensus means `new_distance_to_start` always misses, so every node
/// distance takes the read-position fallback. Useful for exercising the mutation
/// path without spoa; **not** a stand-in for real behaviour, because which
/// bubbles are poppable is most of what the stage decides.
///
/// It can also drive the loop to a fixed point that still looks poppable —
/// linearisation stops changing the graph while every candidate is still
/// accepted, so nothing is ever recorded as not-viable and the iteration never
/// ends. See [`pop_bubbles`]'s iteration cap: the reference has no such bound.
#[derive(Debug, Default)]
pub struct AlwaysPop;

impl BubbleAligner for AlwaysPop {
    fn align(&mut self, _req: &AlignRequest) -> AlignVerdict {
        AlignVerdict {
            poppable: true,
            consensus: FxHashMap::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct PopOpts {
    /// `--slow`: pop every bubble, i.e. a pop threshold of 1 rather than one
    /// percent of the initial edge count.
    pub slow: bool,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct PopStats {
    pub iterations: usize,
    /// True when the iteration cap stopped the loop rather than the reference's
    /// own exit conditions. Always false on real input; see [`pop_bubbles`].
    pub hit_iteration_cap: bool,
    /// Pops actually performed.
    pub pops: usize,
    /// What the reference *prints* as "Overall number of bubbles popped", which
    /// is not the same number. See `PORTING.md` finding 17.
    pub reported_pops: usize,
}

/// `simplifyGraph`: the stage's entry point.
///
/// The reference's wrapper is one line — `new_bubble_popping_routine(DG,
/// all_reads, work_dir, k_size, delta_len, mode)` — and mutates `DG` in place,
/// returning nothing. Same here, except that the statistics come back rather than
/// being printed.
///
/// `work_dir` has no analogue: it exists only so the reference can write a fasta
/// for spoa to read. [`crate::poa`] takes sequences directly.
pub fn simplify_graph(
    g: &mut Graph,
    reads: &FxHashMap<u32, Vec<u8>>,
    k: usize,
    delta_len: i64,
    slow: bool,
) -> PopStats {
    let mut aligner = RealAligner::new(SpoaParasail, reads, k, delta_len);
    pop_bubbles(g, &mut aligner, PopOpts { slow })
}

/// Wraps a [`Consensus`] and counts what it was asked for.
///
/// The point is attribution. spoa's consensus is the one part of this stage whose
/// equivalence to the reference has **not** been re-verified on isONform's inputs
/// (see `PORTING.md`), so when the simplification oracle disagrees on a case, the
/// first question is whether spoa was involved at all. A case with zero spoa calls
/// that still disagrees is a bug in this port; one with many is not yet
/// attributable.
#[derive(Debug, Default)]
pub struct Counting<C> {
    pub inner: C,
    pub spoa_calls: usize,
    pub align_calls: usize,
}

impl<C: Consensus> Counting<C> {
    pub fn new(inner: C) -> Self {
        Self {
            inner,
            spoa_calls: 0,
            align_calls: 0,
        }
    }
}

impl<C: Consensus> Consensus for Counting<C> {
    fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8> {
        self.spoa_calls += 1;
        self.inner.spoa(seqs)
    }
    fn align(&mut self, s1: &[u8], s2: &[u8]) -> Vec<(u32, u8)> {
        self.align_calls += 1;
        self.inner.align(s1, s2)
    }
}

/// `filter_out_if_marked`: drop paths that touch a marked node, or that retrace a
/// start-to-end pair already popped directly this iteration.
fn filter_out_if_marked(
    all_paths: &mut Vec<BubblePath>,
    marked: &FxHashSet<NodeId>,
    direct_combis: &[(NodeId, NodeId)],
    endnode: NodeId,
) {
    all_paths.retain(|path| {
        if path.nodes.iter().any(|n| marked.contains(n)) {
            return false;
        }
        for &(dstart, dend) in direct_combis {
            for (ident, &node) in path.nodes.iter().enumerate() {
                if ident + 1 < path.nodes.len() {
                    if node == dstart && path.nodes[ident + 1] == dend {
                        return false;
                    }
                } else if node == dstart && endnode == dend {
                    return false;
                }
            }
        }
        true
    });
}

fn union_sorted(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut v: Vec<u32> = a.iter().chain(b).copied().collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// `new_bubble_popping_routine`: iterate bubble detection and popping until it
/// stops paying.
///
/// Ported faithfully, which includes three things worth naming:
///
/// * **`overall_pops` undercounts.** The reference adds the previous iteration's
///   pops at the *top* of the loop, so the final iteration's are never added to
///   the number it prints. `PopStats` reports both the true count and the
///   reference's. `PORTING.md` finding 17.
/// * **An operator-precedence bug** in the multi-path branch, reproduced.
///   `PORTING.md` finding 18.
/// * The loop stops when an iteration pops fewer than `pop_threshold`, which is
///   one percent of the *initial* edge count (at least 1), or 1 under `--slow`.
///   So a graph with few edges pops everything it can, and a large one stops
///   while poppable bubbles may remain — deliberate, and not a defect.
pub fn pop_bubbles<A: BubbleAligner>(g: &mut Graph, aligner: &mut A, opts: PopOpts) -> PopStats {
    let initial_edge_nr = g.edge_count();
    let pop_threshold = if opts.slow {
        1
    } else {
        std::cmp::max(initial_edge_nr / 100, 1)
    };

    let mut not_viable_global: FxHashSet<Combination> = FxHashSet::default();
    let mut not_viable_multibubble: FxHashSet<(NodeId, NodeId, Vec<u32>)> = FxHashSet::default();
    let mut stats = PopStats::default();
    let mut this_it_pops = 0usize;

    // A safety net the reference does not have.
    //
    // `while has_combinations` exits only when no candidate survives filtering or
    // when an iteration pops fewer than the threshold. Neither is guaranteed:
    // if linearisation reaches a fixed point while every candidate is still
    // accepted, the loop pops "successfully" forever without changing anything.
    // The real aligner refuses most bubbles, which records them as not-viable and
    // makes that unreachable in practice --- but nothing in the structure
    // enforces it, and a stub aligner reaches it immediately. Bounded here so the
    // failure is a reported flag rather than a hang; generous enough that no
    // converging run can reach it.
    let iteration_cap = g.node_count() + g.edge_count() + 16;

    // `ISONFORM_TRACE_POPS=1` logs one line per pop: iteration, branch, bubble
    // endpoints, both path supports.
    //
    // Permanent rather than scaffolding, because the reference prints its own
    // per-iteration pop counts and the useful comparison is the *sequence*, not
    // the total. Twice now a simplification-oracle failure has been localised by
    // diffing this against the reference's equivalent and finding the first
    // surplus or missing pop --- once for Finding 24, where the local signal
    // (only synthetic reads differ) pointed at the wrong function entirely.
    // Read once here rather than per pop, so it costs nothing when unset.
    let trace_pops = std::env::var_os("ISONFORM_TRACE_POPS").is_some();

    loop {
        if stats.iterations >= iteration_cap {
            stats.hit_iteration_cap = true;
            break;
        }
        stats.reported_pops += this_it_pops; // top of the loop --- see the note
        this_it_pops = 0;
        stats.iterations += 1;

        let mut marked: FxHashSet<NodeId> = FxHashSet::default();
        let mut direct_combis: Vec<(NodeId, NodeId)> = Vec::new();

        let Some(topo) = g.topological_sort() else {
            // The reference lets nx.topological_sort raise here. A cyclic graph
            // at this point is a construction bug, and stopping is the only
            // honest option since every downstream step needs the order.
            break;
        };
        let topo_index = g.topological_index().expect("acyclic, just checked");

        let starts = find_possible_starts(g, &topo);
        let ends = find_possible_ends(g, &topo);
        let combinations = generate_combinations(&starts, &ends, &topo_index);
        let mut filtered: Vec<Combination> = combinations
            .into_iter()
            .filter(|c| !not_viable_global.contains(c))
            .collect();
        if filtered.is_empty() {
            break;
        }

        // Shortest bubbles first. Python's `sorted` is stable, so the generation
        // order survives as the tie-break among equal spans.
        filtered.sort_by_key(|c| {
            topo_index[c.end as usize] as i64 - topo_index[c.start as usize] as i64
        });

        for combination in &filtered {
            let mut all_paths = find_paths(
                g,
                combination.start,
                combination.end,
                &combination.support,
                &marked,
            );
            let initial_all_paths = all_paths.len();
            if initial_all_paths == 1 {
                not_viable_global.insert(combination.clone());
            }
            filter_out_if_marked(&mut all_paths, &marked, &direct_combis, combination.end);

            if all_paths.len() < 2 {
                continue;
            }
            if all_paths.len() == 2 {
                let is_megabubble = initial_all_paths > 2;
                let this_combi_reads = union_sorted(&all_paths[0].support, &all_paths[1].support);
                let this_combi = (combination.start, combination.end, this_combi_reads);
                if not_viable_multibubble.contains(&this_combi) {
                    continue;
                }

                let pathnode1 = get_next_node(&all_paths[0].nodes, combination.end);
                let pathnode2 = get_next_node(&all_paths[1].nodes, combination.end);
                // Interior nodes only: the bubble start is `nodes[0]`.
                let p_set1: FxHashSet<NodeId> =
                    all_paths[0].nodes.iter().skip(1).copied().collect();
                let intersects = all_paths[1]
                    .nodes
                    .iter()
                    .skip(1)
                    .any(|n| p_set1.contains(n));

                if intersects {
                    if initial_all_paths == 2 {
                        not_viable_global.insert(combination.clone());
                    } else {
                        not_viable_multibubble.insert(this_combi);
                    }
                    continue;
                }

                let req = AlignRequest {
                    bubble_start: combination.start,
                    bubble_end: combination.end,
                    path_nodes: [pathnode1, pathnode2],
                    path_support: [all_paths[0].support.clone(), all_paths[1].support.clone()],
                    attrs: [
                        get_consensus_positions(
                            g,
                            combination.start,
                            combination.end,
                            &all_paths[0].support,
                        ),
                        get_consensus_positions(
                            g,
                            combination.start,
                            combination.end,
                            &all_paths[1].support,
                        ),
                    ],
                    is_megabubble,
                };
                let verdict = aligner.align(&req);
                if verdict.poppable {
                    if trace_pops {
                        eprintln!(
                            "POP it={} 2path {} -> {} p1={:?} p2={:?}",
                            stats.iterations,
                            g.key(combination.start),
                            g.key(combination.end),
                            all_paths[0].support,
                            all_paths[1].support
                        );
                    }
                    linearize_bubble(
                        g,
                        combination.start,
                        combination.end,
                        &mut all_paths,
                        &verdict.consensus,
                        &topo_index,
                    );
                    this_it_pops += 1;
                    stats.pops += 1;
                    for p in &all_paths {
                        for &n in &p.nodes {
                            marked.insert(n);
                        }
                    }
                    if all_paths[0].nodes.is_empty() || all_paths[1].nodes.is_empty() {
                        let pair = (combination.start, combination.end);
                        if !direct_combis.contains(&pair) {
                            direct_combis.push(pair);
                        }
                    }
                } else if initial_all_paths == 2 {
                    not_viable_global.insert(combination.clone());
                } else {
                    not_viable_multibubble.insert(this_combi);
                }
            } else {
                // More than two surviving paths: try each pair.
                let mut directpath_marked = false;
                let pairs: Vec<(usize, usize)> = (0..all_paths.len())
                    .flat_map(|i| ((i + 1)..all_paths.len()).map(move |j| (i, j)))
                    .filter(|&(i, j)| {
                        let reads = union_sorted(&all_paths[i].support, &all_paths[j].support);
                        !not_viable_multibubble.contains(&(
                            combination.start,
                            combination.end,
                            reads,
                        ))
                    })
                    .collect();
                if pairs.is_empty() {
                    not_viable_global.insert(combination.clone());
                    continue;
                }

                for (i, j) in pairs {
                    let reads = union_sorted(&all_paths[i].support, &all_paths[j].support);
                    let this_combi = (combination.start, combination.end, reads);
                    let p_set1: FxHashSet<NodeId> =
                        all_paths[i].nodes.iter().skip(1).copied().collect();
                    if all_paths[j]
                        .nodes
                        .iter()
                        .skip(1)
                        .any(|n| p_set1.contains(n))
                    {
                        not_viable_multibubble.insert(this_combi);
                        continue;
                    }
                    // FINDING 18, reproduced: `and` binds tighter than `or`, so
                    // this reads `(!p1) || (!p2 && directpath_marked)` --- an
                    // empty first path skips unconditionally.
                    if all_paths[i].nodes.is_empty()
                        || (all_paths[j].nodes.is_empty() && directpath_marked)
                    {
                        continue;
                    }
                    if all_paths[i].nodes.iter().any(|n| marked.contains(n))
                        || all_paths[j].nodes.iter().any(|n| marked.contains(n))
                    {
                        continue;
                    }

                    let req = AlignRequest {
                        bubble_start: combination.start,
                        bubble_end: combination.end,
                        path_nodes: [
                            get_next_node(&all_paths[i].nodes, combination.end),
                            get_next_node(&all_paths[j].nodes, combination.end),
                        ],
                        path_support: [all_paths[i].support.clone(), all_paths[j].support.clone()],
                        attrs: [
                            get_consensus_positions(
                                g,
                                combination.start,
                                combination.end,
                                &all_paths[i].support,
                            ),
                            get_consensus_positions(
                                g,
                                combination.start,
                                combination.end,
                                &all_paths[j].support,
                            ),
                        ],
                        // The reference hard-codes True here.
                        is_megabubble: true,
                    };
                    let verdict = aligner.align(&req);
                    if !verdict.poppable {
                        not_viable_multibubble.insert(this_combi);
                        continue;
                    }
                    // The reference does `all_paths_filtered = [p1, p2]`, and
                    // those are the *same tuples* that stay in `all_paths` --- so
                    // when `remove_edges` strips the bubble start from each with
                    // `pop(0)`, later pairs in this same iteration see the
                    // shortened lists. Cloning here instead left them unstripped,
                    // which changed which pairs the emptiness and marked checks
                    // skipped and made the loop pop indefinitely. Mutate in place.
                    if trace_pops {
                        eprintln!(
                            "POP it={} multi {} -> {} p1={:?} p2={:?}",
                            stats.iterations,
                            g.key(combination.start),
                            g.key(combination.end),
                            all_paths[i].support,
                            all_paths[j].support
                        );
                    }
                    let mut pair_paths = vec![
                        std::mem::replace(
                            &mut all_paths[i],
                            BubblePath {
                                nodes: Vec::new(),
                                support: Vec::new(),
                            },
                        ),
                        std::mem::replace(
                            &mut all_paths[j],
                            BubblePath {
                                nodes: Vec::new(),
                                support: Vec::new(),
                            },
                        ),
                    ];
                    linearize_bubble(
                        g,
                        combination.start,
                        combination.end,
                        &mut pair_paths,
                        &verdict.consensus,
                        &topo_index,
                    );
                    this_it_pops += 1;
                    // The reference does NOT increment nr_popped here --- it is
                    // commented out --- but does count this_it_pops.
                    stats.pops += 1;
                    for p in &pair_paths {
                        for &n in &p.nodes {
                            marked.insert(n);
                        }
                    }
                    if pair_paths[0].nodes.is_empty() || pair_paths[1].nodes.is_empty() {
                        directpath_marked = true;
                        let pair = (combination.start, combination.end);
                        if !direct_combis.contains(&pair) {
                            direct_combis.push(pair);
                        }
                    }
                    // Put the mutated paths back, so later pairs see them as the
                    // reference does.
                    all_paths[j] = pair_paths.pop().expect("two paths");
                    all_paths[i] = pair_paths.pop().expect("two paths");
                }
            }
        }

        // `ISONFORM_TRACE_POPS=1`: this iteration's pop count, to line up against
        // the reference's own "This iterations pops" line. Comparing totals hides
        // *where* two runs diverge; comparing the per-iteration series says which
        // iteration to look at.
        if std::env::var_os("ISONFORM_TRACE_POPS").is_some() {
            eprintln!(
                "trace-pops iteration {} pops {this_it_pops}",
                stats.iterations
            );
        }
        if this_it_pops < pop_threshold {
            break;
        }
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{NodeKey, ReadInfo};

    fn ri(a: i64, b: i64) -> ReadInfo {
        ReadInfo {
            start_mini_end: a,
            end_mini_start: b,
            original_support: true,
        }
    }

    /// s -> a -> t and s -> b -> t: `s` starts a bubble, `t` ends one.
    fn diamond() -> (Graph, NodeId, NodeId, NodeId, NodeId) {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let t = g.add_node(NodeKey::Sink);
        g.add_edge(s, a, 0);
        g.add_edge(s, b, 0);
        g.add_edge(a, t, 0);
        g.add_edge(b, t, 0);
        for (n, reads) in [
            (s, &[1u32, 2, 3][..]),
            (a, &[1, 2]),
            (b, &[3]),
            (t, &[1, 2, 3]),
        ] {
            for &r in reads {
                g.set_read(n, r, ri(0, (r * 10) as i64));
            }
        }
        (g, s, a, b, t)
    }

    #[test]
    fn starts_are_out_degree_gt_one_and_ends_in_degree_gt_one() {
        let (g, s, _a, _b, t) = diamond();
        let topo = g.topological_sort().unwrap();
        let starts = find_possible_starts(&g, &topo);
        let ends = find_possible_ends(&g, &topo);
        assert_eq!(starts.iter().map(|e| e.node).collect::<Vec<_>>(), vec![s]);
        assert_eq!(ends.iter().map(|e| e.node).collect::<Vec<_>>(), vec![t]);
    }

    #[test]
    fn endpoints_come_out_in_topological_order() {
        // Two branch points in sequence: the earlier one must be listed first,
        // regardless of the order the nodes were added to the graph.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        // Add the *later* branch point first, so insertion order disagrees with
        // topological order.
        let second = g.add_node(k(9));
        let first = g.add_node(k(1));
        let x = g.add_node(k(2));
        let y = g.add_node(k(3));
        let p = g.add_node(k(4));
        let q = g.add_node(k(5));
        g.add_edge(first, x, 0);
        g.add_edge(first, y, 0);
        g.add_edge(x, second, 0);
        g.add_edge(y, second, 0);
        g.add_edge(second, p, 0);
        g.add_edge(second, q, 0);
        let topo = g.topological_sort().unwrap();
        let starts = find_possible_starts(&g, &topo);
        assert_eq!(
            starts.iter().map(|e| e.node).collect::<Vec<_>>(),
            vec![first, second],
            "collected in topological order, not insertion order"
        );
    }

    #[test]
    fn a_combination_needs_two_shared_reads_and_the_right_direction() {
        let (g, s, _a, _b, t) = diamond();
        let topo = g.topological_sort().unwrap();
        let idx = g.topological_index().unwrap();
        let combis = generate_combinations(
            &find_possible_starts(&g, &topo),
            &find_possible_ends(&g, &topo),
            &idx,
        );
        assert_eq!(combis.len(), 1);
        assert_eq!(combis[0].start, s);
        assert_eq!(combis[0].end, t);
        assert_eq!(combis[0].support, vec![1, 2, 3], "sorted intersection");
    }

    #[test]
    fn a_single_shared_read_is_not_a_bubble() {
        let (mut g, s, _a, _b, t) = diamond();
        // Strip the sink down to one shared read.
        g.replace_reads(t, 1, ri(0, 10));
        let topo = g.topological_sort().unwrap();
        let idx = g.topological_index().unwrap();
        let combis = generate_combinations(
            &find_possible_starts(&g, &topo),
            &find_possible_ends(&g, &topo),
            &idx,
        );
        assert!(combis.is_empty(), "len(inter) >= 2 is required");
        let _ = s;
    }

    #[test]
    fn the_intersection_is_sorted_not_in_map_order() {
        // Build supports whose insertion order is descending, so a port that
        // forgot to sort would be caught.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let t = g.add_node(NodeKey::Sink);
        g.add_edge(s, a, 0);
        g.add_edge(s, b, 0);
        g.add_edge(a, t, 0);
        g.add_edge(b, t, 0);
        for n in [s, t] {
            for &r in &[9u32, 4, 7, 1] {
                g.set_read(n, r, ri(0, 0));
            }
        }
        let topo = g.topological_sort().unwrap();
        let idx = g.topological_index().unwrap();
        let combis = generate_combinations(
            &find_possible_starts(&g, &topo),
            &find_possible_ends(&g, &topo),
            &idx,
        );
        assert_eq!(combis[0].support, vec![1, 4, 7, 9]);
    }

    #[test]
    fn filter_combinations_drops_known_bad_ones() {
        let c1 = Combination {
            start: 0,
            end: 1,
            support: vec![1, 2],
        };
        let c2 = Combination {
            start: 0,
            end: 2,
            support: vec![1, 3],
        };
        let mut not_viable = FxHashSet::default();
        not_viable.insert(c1.clone());
        let kept = filter_combinations(&[c1, c2.clone()], &not_viable);
        assert_eq!(kept, vec![c2]);
    }

    /// `s -> a -> t` carrying reads 1,2 and `s -> b -> t` carrying read 3.
    fn two_path_bubble() -> (Graph, NodeId, NodeId, NodeId, NodeId) {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let t = g.add_node(NodeKey::Sink);
        g.add_edge(s, a, 0);
        g.set_edge_support(s, a, 1);
        g.push_edge_support(s, a, 2);
        g.add_edge(a, t, 0);
        g.set_edge_support(a, t, 1);
        g.push_edge_support(a, t, 2);
        g.add_edge(s, b, 0);
        g.set_edge_support(s, b, 3);
        g.add_edge(b, t, 0);
        g.set_edge_support(b, t, 3);
        // Node reads matter: `find_possible_starts`/`_ends` read them, and
        // `generate_combinations` needs at least two reads shared between the
        // start and the end before it will call this a bubble at all. A first
        // draft of this fixture set only edge support, so `pop_bubbles` found no
        // candidates and every driver test silently measured nothing.
        for (n, rs) in [
            (s, &[1u32, 2, 3][..]),
            (a, &[1, 2][..]),
            (b, &[3][..]),
            (t, &[1, 2, 3][..]),
        ] {
            for &r in rs {
                g.set_read(n, r, ri(0, (r * 10) as i64));
            }
        }
        (g, s, a, b, t)
    }

    #[test]
    fn find_paths_splits_reads_by_the_route_they_take() {
        let (g, s, a, b, t) = two_path_bubble();
        let paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        assert_eq!(paths.len(), 2);
        // Lowest read first, so the 1,2 path comes out before the 3 path.
        assert_eq!(paths[0].nodes, vec![s, a]);
        assert_eq!(paths[0].support, vec![1, 2]);
        assert_eq!(paths[1].nodes, vec![s, b]);
        assert_eq!(paths[1].support, vec![3]);
    }

    #[test]
    fn the_bubble_end_is_not_part_of_a_path() {
        // `visited_nodes.append(node)` happens before the `node != endnode` test,
        // so the end node is never appended.
        let (g, s, a, _b, t) = two_path_bubble();
        let paths = find_paths(&g, s, t, &[1, 2], &FxHashSet::default());
        assert_eq!(paths[0].nodes, vec![s, a]);
        assert!(!paths[0].nodes.contains(&t));
    }

    #[test]
    fn path_order_follows_the_lowest_read_not_the_support_order() {
        // Support given in descending order; the 3-path must still come second,
        // because allocation takes the minimum read id. This is finding 12: the
        // reference used set.pop() here and the order decided which path became
        // path1.
        let (g, s, a, b, t) = two_path_bubble();
        let paths = find_paths(&g, s, t, &[3, 2, 1], &FxHashSet::default());
        assert_eq!(paths[0].nodes, vec![s, a], "reads 1,2 first");
        assert_eq!(paths[1].nodes, vec![s, b]);
    }

    #[test]
    fn a_marked_node_stops_the_walk_and_ends_enumeration() {
        // Hitting a marked node sets next_found = false, which makes the
        // reference break out of the *outer* loop too --- so no further paths are
        // enumerated, not merely this one skipped.
        let (g, s, a, b, t) = two_path_bubble();
        let mut marked = FxHashSet::default();
        marked.insert(a);
        let paths = find_paths(&g, s, t, &[1, 2, 3], &marked);
        assert!(
            paths.is_empty(),
            "read 1 walks into the marked node, which abandons the whole bubble"
        );
        let _ = (b, t);
    }

    #[test]
    fn a_read_with_no_supporting_edge_ends_enumeration() {
        let (g, s, _a, _b, t) = two_path_bubble();
        // Read 9 supports nothing, and it is the highest, so the two real paths
        // are found first and then the walk for 9 fails.
        let paths = find_paths(&g, s, t, &[1, 2, 3, 9], &FxHashSet::default());
        assert_eq!(paths.len(), 2, "the two real paths survive");
        // And when the failing read is the lowest, nothing is found at all.
        let paths = find_paths(&g, s, t, &[0, 1, 2, 3], &FxHashSet::default());
        assert!(paths.is_empty(), "read 0 fails first and aborts the bubble");
    }

    #[test]
    fn the_first_matching_successor_wins() {
        // A read supported by two outgoing edges follows whichever edge was added
        // first, because the reference breaks out of the successor loop.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let first = g.add_node(k(1));
        let second = g.add_node(k(2));
        let t = g.add_node(NodeKey::Sink);
        // Both edges out of s carry read 1; `first` was added first.
        g.add_edge(s, first, 0);
        g.set_edge_support(s, first, 1);
        g.add_edge(s, second, 0);
        g.set_edge_support(s, second, 1);
        g.add_edge(first, t, 0);
        g.set_edge_support(first, t, 1);
        g.add_edge(second, t, 0);
        g.set_edge_support(second, t, 1);
        let paths = find_paths(&g, s, t, &[1, 2], &FxHashSet::default());
        assert_eq!(paths[0].nodes, vec![s, first], "adjacency order decides");
    }

    #[test]
    fn support_is_intersected_along_the_whole_walk() {
        // Read 1 and 2 leave `s` together but only 1 continues past `a`, so the
        // path's support is the intersection, not the starting set.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let t = g.add_node(NodeKey::Sink);
        g.add_edge(s, a, 0);
        g.set_edge_support(s, a, 1);
        g.push_edge_support(s, a, 2);
        g.add_edge(a, t, 0);
        g.set_edge_support(a, t, 1);
        let paths = find_paths(&g, s, t, &[1, 2], &FxHashSet::default());
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].support, vec![1], "2 is dropped at the second edge");
    }

    #[test]
    fn new_distance_to_start_is_pythons_str_find() {
        assert_eq!(new_distance_to_start(b"ACGTACGT", b"GTA"), 2);
        assert_eq!(
            new_distance_to_start(b"ACGTACGT", b"ACGT"),
            0,
            "first match wins"
        );
        assert_eq!(new_distance_to_start(b"ACGTACGT", b"TTT"), -1);
        assert_eq!(
            new_distance_to_start(b"AC", b"ACGT"),
            -1,
            "needle longer than haystack"
        );
        // Python: "abc".find("") == 0
        assert_eq!(new_distance_to_start(b"ACGT", b""), 0);
        assert_eq!(new_distance_to_start(b"", b""), 0);
    }

    #[test]
    fn remove_edges_lifts_the_whole_bubble_and_drops_the_start_from_each_path() {
        let (mut g, s, a, b, t) = two_path_bubble();
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        assert_eq!(paths[0].nodes, vec![s, a], "before: the start is present");

        let lifted = remove_edges(&mut g, s, t, &mut paths, &FxHashMap::default());

        // Destructive: the bubble start is gone from every path.
        assert_eq!(paths[0].nodes, vec![a]);
        assert_eq!(paths[1].nodes, vec![b]);

        // All four bubble edges lifted, and gone from the graph.
        assert_eq!(lifted.edges.len(), 4);
        for (u, v) in [(s, a), (a, t), (s, b), (b, t)] {
            assert!(!g.has_edge(u, v), "edge should have been removed");
            assert!(
                lifted.edges.iter().any(|e| e.u == u && e.v == v),
                "edge should have been recorded"
            );
        }
        // Node set untouched --- the invariant the dumps show holds for the whole
        // stage.
        assert_eq!(g.node_count(), 4);
    }

    #[test]
    fn lifted_edges_keep_their_length_and_support() {
        let (mut g, s, a, _b, t) = two_path_bubble();
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let lifted = remove_edges(&mut g, s, t, &mut paths, &FxHashMap::default());
        let e = lifted.edges.iter().find(|e| e.u == s && e.v == a).unwrap();
        assert_eq!(e.length, Some(0));
        assert_eq!(e.support, vec![1, 2]);
    }

    #[test]
    fn a_lengthless_edge_inside_the_bubble_is_still_lifted_and_removed() {
        // An edge `prepare_adding_edges` re-added on an earlier pop carries no
        // `length` (`Graph::upsert_edge_support`). `remove_edges` used to gate
        // lifting on `edge_length` returning `Some`, which treated "exists with
        // no length" the same as "does not exist" and left the edge in the
        // graph forever --- so a later pop that reused the same bubble start
        // kept finding it, and `pop_bubbles` never converged. Regression for
        // that: the edge must be lifted (and actually removed from the graph)
        // even though it has no length.
        let (mut g, s, a, _b, t) = two_path_bubble();
        // Simulate what an earlier pop leaves behind: the edge was removed and
        // then re-added via `upsert_edge_support`, which --- like the
        // reference's `DG.add_edge(u, v, edge_supp=...)` on a fresh pair ---
        // creates it with no `length`.
        let supp = g.edge_support(s, a).unwrap().to_vec();
        g.remove_edge(s, a);
        g.upsert_edge_support(s, a, supp);
        assert!(g.has_edge(s, a));
        assert_eq!(g.edge_length(s, a), None);

        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let lifted = remove_edges(&mut g, s, t, &mut paths, &FxHashMap::default());

        let e = lifted.edges.iter().find(|e| e.u == s && e.v == a).unwrap();
        assert_eq!(e.length, None);
        assert_eq!(e.support, vec![1, 2]);
        assert!(
            !g.has_edge(s, a),
            "a lengthless edge inside the bubble must still be removed"
        );
    }

    #[test]
    fn a_consensus_hit_gives_the_offset_and_a_miss_falls_back_to_the_read_distance() {
        let (mut g, s, a, b, t) = two_path_bubble();
        // Give the nodes end-minimizers and read positions so both branches can
        // be told apart.
        g.set_end_mini_seq(a, b"GGGG");
        g.set_end_mini_seq(b, b"CCCC");
        g.replace_reads(s, 1, ri(0, 100));
        g.set_read(s, 2, ri(0, 100));
        g.set_read(s, 3, ri(0, 100));
        g.replace_reads(a, 1, ri(0, 140));
        g.set_read(a, 2, ri(0, 140));

        let mut consensus = FxHashMap::default();
        // Path 1 starts at `a`: its consensus contains a's minimizer at offset 7.
        consensus.insert(a, b"TTTTTTTGGGG".to_vec());
        // Path 2 starts at `b`: consensus lacks b's minimizer, forcing -1.
        consensus.insert(b, b"TTTTTTT".to_vec());

        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let lifted = remove_edges(&mut g, s, t, &mut paths, &consensus);

        assert_eq!(
            lifted.dist(a),
            Some(7.0),
            "found in the consensus -> the offset"
        );
        // b was not found in its consensus, so the distance comes from
        // get_dist_to_prev(s, b): they share only read 3, at 100 on s and 30 on b
        // (the fixture's `r * 10`), so the mean shift is 30 - 100 = -70. Negative
        // distances are reachable, which is why get_avg_interval_length has to
        // truncate toward zero rather than floor.
        assert_eq!(
            lifted.dist(b),
            Some(-70.0),
            "-1 -> measured from the bubble start"
        );
    }

    #[test]
    fn the_fallback_always_measures_from_the_bubble_start() {
        // Finding 13. On a three-node path the second interior node's fallback
        // distance is measured from the bubble START, not from the node before
        // it, because `prev_node` never advances. If it chained, the value would
        // be the sum of the two hops.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let x = g.add_node(k(1));
        let y = g.add_node(k(2));
        let other = g.add_node(k(3));
        let t = g.add_node(NodeKey::Sink);
        for (u, v, supp) in [
            (s, x, &[1u32, 2][..]),
            (x, y, &[1, 2]),
            (y, t, &[1, 2]),
            (s, other, &[3]),
            (other, t, &[3]),
        ] {
            g.add_edge(u, v, 0);
            let mut first = true;
            for &r in supp {
                if first {
                    g.set_edge_support(u, v, r);
                    first = false;
                } else {
                    g.push_edge_support(u, v, r);
                }
            }
        }
        // Positions: s at 0, x at 10, y at 30. Chaining would give y = 10 + 20;
        // measuring from the start gives y = 30. Same here --- so make the read
        // sets differ, which is where the two disagree: read 2 is absent from x.
        g.set_read(s, 1, ri(0, 0));
        g.set_read(s, 2, ri(0, 0));
        g.set_read(x, 1, ri(0, 10));
        g.set_read(y, 1, ri(0, 30));
        g.set_read(y, 2, ri(0, 50));

        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        // No consensus at all, so every distance takes the -1 fallback.
        let lifted = remove_edges(&mut g, s, t, &mut paths, &FxHashMap::default());

        // get_dist_to_prev(s, y) over reads 1 and 2: (30-0) and (50-0) -> mean 40.
        // Chaining would have given dist(x) + get_dist_to_prev(x, y) = 10 + 20 = 30.
        assert_eq!(lifted.dist(x), Some(10.0));
        assert_eq!(
            lifted.dist(y),
            Some(40.0),
            "measured from the bubble start, not chained from x"
        );
    }

    #[test]
    fn a_one_edge_path_records_the_start_to_end_edge() {
        // A path with no interior node: `pnl_start` falls back to the bubble end.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let t = g.add_node(NodeKey::Sink);
        g.add_edge(s, a, 0);
        g.set_edge_support(s, a, 1);
        g.add_edge(a, t, 0);
        g.set_edge_support(a, t, 1);
        g.add_edge(s, t, 5); // the direct path, read 2
        g.set_edge_support(s, t, 2);

        let mut paths = find_paths(&g, s, t, &[1, 2], &FxHashSet::default());
        assert_eq!(paths.len(), 2);
        let lifted = remove_edges(&mut g, s, t, &mut paths, &FxHashMap::default());
        assert!(
            lifted
                .edges
                .iter()
                .any(|e| e.u == s && e.v == t && e.length == Some(5)),
            "the direct start->end edge is lifted too"
        );
        assert!(
            paths[1].nodes.is_empty(),
            "the direct path has no interior node"
        );
    }

    #[test]
    fn merge_two_dicts_lets_the_first_win_and_keeps_its_order() {
        let a = [(3u32, ri(0, 3)), (1, ri(0, 1))];
        let b = [(1u32, ri(9, 9)), (2, ri(0, 2))];
        let m = merge_two_dicts(&a, &b);
        let keys: Vec<u32> = m.iter().map(|(k, _)| *k).collect();
        assert_eq!(
            keys,
            vec![3, 1, 2],
            "dict1's order leads, dict2's extras follow"
        );
        assert_eq!(m[1].1, ri(0, 1), "dict1's value wins on a shared key");
    }

    #[test]
    fn the_bubble_end_never_wins_the_next_node_race() {
        // This is what keeps the walk going while one path still has nodes, and
        // what makes the "belongs to neither path" branch unreachable.
        assert_eq!(compare_by_length(9, 5, 9, 0, 100), 5);
        assert_eq!(compare_by_length(5, 9, 9, 100, 0), 5);
        // Otherwise the earlier topological index wins, ties to nextnode2.
        assert_eq!(compare_by_length(1, 2, 9, 3, 7), 1);
        assert_eq!(compare_by_length(1, 2, 9, 7, 3), 2);
        assert_eq!(
            compare_by_length(1, 2, 9, 5, 5),
            2,
            "a tie goes to nextnode2"
        );
    }

    #[test]
    fn linearizing_a_two_path_bubble_produces_one_chain() {
        let (mut g, s, a, b, t) = two_path_bubble();
        let topo = g.topological_index().unwrap();
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let order = linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);

        // Every node still there --- the invariant the recorded dumps show for the
        // whole stage.
        assert_eq!(g.node_count(), 4);
        // A single chain start -> ... -> end, visiting both former branches.
        assert_eq!(order.len(), 4);
        assert_eq!(order[0], s);
        assert_eq!(*order.last().unwrap(), t);
        assert!(order.contains(&a) && order.contains(&b));
        // Consecutive nodes in the order are joined, and nothing else is.
        for w in order.windows(2) {
            assert!(g.has_edge(w[0], w[1]), "chain edge missing");
        }
        assert_eq!(g.edge_count(), 3, "a 4-node chain has exactly 3 edges");
    }

    #[test]
    fn relinked_edges_carry_no_length_and_the_union_of_both_paths_support() {
        let (mut g, s, _a, _b, t) = two_path_bubble();
        let topo = g.topological_index().unwrap();
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let order = linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);

        // Finding: re-added edges have no `length` attribute at all.
        assert_eq!(
            g.edge_length(order[0], order[1]),
            None,
            "prepare_adding_edges never names a length"
        );
        assert!(g.has_edge(order[0], order[1]), "but the edge does exist");

        // The first hop carries both paths' support, since it replaces both
        // branches.
        let mut supp = g.edge_support(order[0], order[1]).unwrap().to_vec();
        supp.sort_unstable();
        assert_eq!(supp, vec![1, 2, 3]);
    }

    #[test]
    fn the_topological_order_decides_which_branch_comes_first() {
        // The same bubble as `two_path_bubble` but with the branches wired in
        // the opposite order, so `b` takes the earlier topological slot and the
        // chain visits b before a.
        //
        // Note what actually decides that: **edge** insertion order, not node
        // insertion order. Only the first generation is seeded from the node
        // list; after that a generation is ordered by parent, then by successor.
        // A first draft of this test added `b` as a node before `a` and asserted
        // b came first --- it does not, because both are reached from `s` and
        // `s`'s successor list is what orders them.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let t = g.add_node(NodeKey::Sink);
        // b's branch is wired first this time.
        g.add_edge(s, b, 0);
        g.set_edge_support(s, b, 3);
        g.add_edge(b, t, 0);
        g.set_edge_support(b, t, 3);
        g.add_edge(s, a, 0);
        g.set_edge_support(s, a, 1);
        g.push_edge_support(s, a, 2);
        g.add_edge(a, t, 0);
        g.set_edge_support(a, t, 1);
        g.push_edge_support(a, t, 2);

        let topo = g.topological_index().unwrap();
        assert!(
            topo[b as usize] < topo[a as usize],
            "b comes first because s->b was added first"
        );
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let order = linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);
        assert_eq!(order, vec![s, b, a, t]);
    }

    #[test]
    fn the_connecting_edge_pick_is_the_lexicographically_smallest_source() {
        // `find_connecting_edges` can return two edges ending at the same node,
        // and `prepare_adding_edges` folds the support of `conn_list[0]` into the
        // new edge. In the reference `conn_list` comes from a Python `set` of
        // string tuples, so which one that is depends on PYTHONHASHSEED --- two
        // distinct outputs across 8 seeds on real data. The port orders by the
        // source node's name instead. PORTING.md finding 14.
        //
        // Bubble s -> t. path1 = [a1, a2], path2 = [c], with BOTH a1 -> c and
        // a2 -> c present, so `test_conn_end` for `c` has two candidates.
        let mut g = Graph::new();
        // Names sort as "100, .." < "99, ..", so lexicographic and numeric order
        // disagree here on purpose: the reference's names are strings.
        let k = |start: u32| NodeKey::Interval {
            start,
            end: start + 10,
            r_id: 1,
        };
        let s = g.add_node(NodeKey::Source);
        let a1 = g.add_node(k(99));
        let a2 = g.add_node(k(100));
        let c = g.add_node(k(200));
        let t = g.add_node(NodeKey::Sink);
        assert!(
            g.key(a2).to_string() < g.key(a1).to_string(),
            "fixture intent: '100, ..' sorts before '99, ..' as a string"
        );

        // path1: s -> a1 -> a2 -> t   path2: s -> c -> t
        for (u, v, r) in [(s, a1, 1u32), (a1, a2, 1), (a2, t, 1), (s, c, 2), (c, t, 2)] {
            g.add_edge(u, v, 0);
            g.set_edge_support(u, v, r);
        }
        // The two connecting edges into `c`, each with its own marker read.
        g.add_edge(a1, c, 0);
        g.set_edge_support(a1, c, 91);
        g.add_edge(a2, c, 0);
        g.set_edge_support(a2, c, 92);

        for n in [s, a1, a2, c, t] {
            g.set_read(n, 1, ri(0, 10));
            g.set_read(n, 2, ri(0, 10));
        }

        let topo = g.topological_index().unwrap();
        let mut paths = vec![
            BubblePath {
                nodes: vec![s, a1, a2],
                support: vec![1],
            },
            BubblePath {
                nodes: vec![s, c],
                support: vec![2],
            },
        ];
        linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);

        // `a2` sorts first as a string, so its edge's support (92) is the one
        // folded in --- not `a1`'s (91), which numeric order would have chosen.
        // The linearised chain runs s -> a1 -> a2 -> c -> t, and the edge into
        // `c` is where `conn_list[0]`'s support gets folded in. `a2` sorts first
        // as a string, so that is 92. Without the ordering the port would take
        // `find_connecting_edges`' insertion order, which visits `a1` first and
        // would fold in 91 --- so this test discriminates.
        let into_c = g.edge_support(a2, c).expect("a2 -> c is in the new chain");
        assert!(
            into_c.contains(&92),
            "the lexicographically smallest source's support is folded in, got {into_c:?}"
        );
        assert!(
            !into_c.contains(&91),
            "and the other candidate's is not, got {into_c:?}"
        );
        // `a1 -> c` was never on a bubble path, so it survives untouched --- which
        // is why the check has to be on one edge rather than on the whole graph.
        assert_eq!(g.edge_support(a1, c).unwrap(), &[91]);
    }

    #[test]
    fn an_existing_edge_keeps_its_length_and_accumulates_support() {
        // A bubble whose start connects directly to its end: that edge is lifted
        // and re-added, so it is created fresh without a length. But an edge that
        // was NOT lifted and gets support added keeps its length --- covered by
        // graph.rs's upsert test; here we check the accumulate half end to end.
        let (mut g, s, a, _b, t) = two_path_bubble();
        let topo = g.topological_index().unwrap();
        // An extra read on the surviving s->a edge that no path carries.
        g.push_edge_support(s, a, 77);
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        let order = linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);
        // 77 rode along on the lifted edge's support, so it is still present.
        let supp = g.edge_support(order[0], order[1]).unwrap();
        assert!(
            supp.contains(&77),
            "lifted support is carried into the new edge"
        );
    }

    #[test]
    fn additional_support_is_marked_as_not_original() {
        // A read reaching a node only via the other branch gets an invented
        // position with original_support = false, which is how the graph records
        // "this read was never really here".
        let (mut g, s, a, b, t) = two_path_bubble();
        let topo = g.topological_index().unwrap();
        g.set_read(s, 3, ri(0, 100));
        g.set_read(b, 3, ri(0, 120));
        let mut paths = find_paths(&g, s, t, &[1, 2, 3], &FxHashSet::default());
        linearize_bubble(&mut g, s, t, &mut paths, &FxHashMap::default(), &topo);
        // Whichever node gained a read from the other path has it flagged.
        let invented: Vec<_> = [a, b]
            .iter()
            .flat_map(|&n| g.reads(n).iter().copied())
            .filter(|(_, i)| !i.original_support)
            .collect();
        assert!(
            !invented.is_empty(),
            "at least one read was carried across the bubble"
        );
    }

    #[test]
    fn never_pop_leaves_the_graph_untouched_and_terminates() {
        let (mut g, _s, _a, _b, _t) = two_path_bubble();
        let before = (g.node_count(), g.edge_count());
        let stats = pop_bubbles(&mut g, &mut NeverPop, PopOpts::default());
        assert_eq!((g.node_count(), g.edge_count()), before);
        assert_eq!(stats.pops, 0);
        assert_eq!(stats.iterations, 1, "one iteration, then nothing pops");
    }

    #[test]
    fn a_refused_two_path_bubble_is_remembered_as_not_viable() {
        // The proof it is remembered: a second call finds the same bubble again,
        // so if it were not cached the iteration count would keep climbing.
        // Instead the routine stops after one iteration both times.
        let (mut g, _s, _a, _b, _t) = two_path_bubble();
        let first = pop_bubbles(&mut g, &mut NeverPop, PopOpts::default());
        let second = pop_bubbles(&mut g, &mut NeverPop, PopOpts::default());
        assert_eq!(first.iterations, 1);
        assert_eq!(second.iterations, 1);
    }

    #[test]
    fn popping_a_two_path_bubble_linearises_it_and_stops() {
        let (mut g, s, a, b, t) = two_path_bubble();
        let stats = pop_bubbles(&mut g, &mut AlwaysPop, PopOpts::default());
        assert_eq!(stats.pops, 1);
        // Same nodes, now a chain.
        assert_eq!(g.node_count(), 4);
        assert_eq!(g.edge_count(), 3);
        for n in [s, a, b, t] {
            assert!(g.lookup(&g.key(n)).is_some());
        }
        // Nothing left to pop: a chain has no branch points.
        let again = pop_bubbles(&mut g, &mut AlwaysPop, PopOpts::default());
        assert_eq!(again.pops, 0);
    }

    #[test]
    fn the_reported_pop_count_agrees_when_the_threshold_is_one() {
        // Finding 17 is real but narrower than it first looks, and this test
        // records the boundary rather than the bug.
        //
        // `overall_pops += this_it_pops` runs at the *top* of the loop, so the
        // final iteration's pops are never added to the number the reference
        // prints. But the loop only exits early when `this_it_pops <
        // pop_threshold`, and the threshold is `max(edges / 100, 1)` --- so below
        // 200 initial edges the threshold is 1, an early exit means zero pops
        // that iteration, and nothing is lost. The undercount needs a graph big
        // enough for a threshold above 1 *and* a final iteration popping between
        // 1 and threshold-1. Measured on real data in PORTING.md.
        let (mut g, _s, _a, _b, _t) = two_path_bubble();
        let stats = pop_bubbles(&mut g, &mut AlwaysPop, PopOpts::default());
        assert_eq!(stats.pops, 1);
        assert_eq!(
            stats.reported_pops, stats.pops,
            "with a threshold of 1 the two agree"
        );
    }

    #[test]
    fn slow_mode_lowers_the_pop_threshold_to_one() {
        // With 4 edges, max(4/100, 1) is 1 either way, so this checks the wiring
        // rather than a behaviour difference: both modes agree here.
        let (mut g1, _, _, _, _) = two_path_bubble();
        let (mut g2, _, _, _, _) = two_path_bubble();
        let fast = pop_bubbles(&mut g1, &mut AlwaysPop, PopOpts { slow: false });
        let slow = pop_bubbles(&mut g2, &mut AlwaysPop, PopOpts { slow: true });
        assert_eq!(fast.pops, slow.pops);
    }

    #[test]
    fn the_aligner_sees_both_path_nodes_and_their_support() {
        // What a real aligner needs, and what the stub proves is wired: the two
        // path-defining nodes and each path's reads.
        struct Recorder(Vec<AlignRequest>);
        impl BubbleAligner for Recorder {
            fn align(&mut self, req: &AlignRequest) -> AlignVerdict {
                self.0.push(req.clone());
                AlignVerdict::default()
            }
        }
        let (mut g, s, a, b, t) = two_path_bubble();
        let mut rec = Recorder(Vec::new());
        pop_bubbles(&mut g, &mut rec, PopOpts::default());
        assert_eq!(rec.0.len(), 1);
        let req = &rec.0[0];
        assert_eq!((req.bubble_start, req.bubble_end), (s, t));
        // Each path is [start, interior], i.e. length 2, so `get_next_node`
        // returns `nodes[1]` --- the interior node --- rather than falling back to
        // the bubble end. The fallback needs a path shorter than 2.
        assert_eq!(req.path_nodes, [a, b]);
        assert_eq!(req.path_support[0], vec![1, 2]);
        assert_eq!(req.path_support[1], vec![3]);
        assert!(!req.is_megabubble);
        let _ = t;
    }

    #[test]
    fn an_aligner_that_never_refuses_still_converges_by_merging_paths_away() {
        // A three-path bubble plus an always-yes aligner: earlier this hit the
        // iteration cap, because `remove_edges` gated lifting an edge on it
        // having a `length` rather than on it existing, so a re-added
        // (lengthless) edge from an earlier pop was never removed on the next
        // one --- see `a_lengthless_edge_inside_the_bubble_is_still_lifted_and_removed`.
        // With that fixed, three parallel paths merge into one linear chain and
        // the loop finds nothing left to pop: no fixed point here after all, the
        // node set stays put (popping relinks edges, never deletes nodes) but
        // the edge count drops from 6 to 4 as the chain forms.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let t = g.add_node(NodeKey::Sink);
        for (i, r) in [(1u32, 1u32), (2, 2), (3, 3)] {
            let m = g.add_node(k(i));
            g.add_edge(s, m, 0);
            g.set_edge_support(s, m, r);
            g.add_edge(m, t, 0);
            g.set_edge_support(m, t, r);
            g.set_read(m, r, ri(0, (r * 10) as i64));
            g.set_read(s, r, ri(0, 0));
            g.set_read(t, r, ri(0, 100));
        }
        let stats = pop_bubbles(&mut g, &mut AlwaysPop, PopOpts::default());
        assert!(!stats.hit_iteration_cap, "this now terminates on its own");
        assert_eq!(stats.pops, 2, "three paths merge pairwise, twice");
        assert_eq!(
            g.node_count(),
            5,
            "popping relinks edges, it never deletes nodes"
        );
        assert_eq!(g.edge_count(), 4, "one linear chain: s -> m -> m -> m -> t");
    }

    #[test]
    fn filter_out_if_marked_drops_paths_touching_marked_nodes() {
        let mut paths = vec![
            BubblePath {
                nodes: vec![1, 2],
                support: vec![1],
            },
            BubblePath {
                nodes: vec![1, 3],
                support: vec![2],
            },
        ];
        let mut marked = FxHashSet::default();
        marked.insert(2u32);
        filter_out_if_marked(&mut paths, &marked, &[], 9);
        assert_eq!(paths.len(), 1);
        assert_eq!(paths[0].nodes, vec![1, 3]);
    }

    #[test]
    fn filter_out_if_marked_drops_a_path_retracing_a_direct_pop() {
        // A path whose last node pairs with the bubble end via a recorded direct
        // combination is dropped, which is what stops the same start-to-end pair
        // being popped twice in one iteration.
        let mut paths = vec![BubblePath {
            nodes: vec![1, 5],
            support: vec![1],
        }];
        filter_out_if_marked(&mut paths, &FxHashSet::default(), &[(5, 9)], 9);
        assert!(paths.is_empty());

        // And an adjacent pair inside the path counts too.
        let mut paths = vec![BubblePath {
            nodes: vec![1, 5, 7],
            support: vec![1],
        }];
        filter_out_if_marked(&mut paths, &FxHashSet::default(), &[(5, 7)], 9);
        assert!(paths.is_empty());
    }

    #[test]
    fn a_three_path_bubble_pops_one_pair_then_gives_up() {
        // s -> {a, b, c} -> t, three disjoint paths.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(NodeKey::Source);
        let t = g.add_node(NodeKey::Sink);
        let mut mids = Vec::new();
        for (i, r) in [(1u32, 1u32), (2, 2), (3, 3)] {
            let m = g.add_node(k(i));
            g.add_edge(s, m, 0);
            g.set_edge_support(s, m, r);
            g.add_edge(m, t, 0);
            g.set_edge_support(m, t, r);
            g.set_read(m, r, ri(0, (r * 10) as i64));
            // Node reads, without which no bubble is detected at all.
            g.set_read(s, r, ri(0, 0));
            g.set_read(t, r, ri(0, 100));
            mids.push(m);
        }
        // An aligner that accepts once and then refuses --- which is what a real
        // one does, and what lets the loop terminate. `AlwaysPop` would spin here:
        // linearisation reaches a fixed point while every candidate is still
        // accepted, so nothing is ever recorded not-viable.
        struct PopOnce(usize);
        impl BubbleAligner for PopOnce {
            fn align(&mut self, _req: &AlignRequest) -> AlignVerdict {
                self.0 += 1;
                AlignVerdict {
                    poppable: self.0 == 1,
                    consensus: FxHashMap::default(),
                }
            }
        }
        let stats = pop_bubbles(&mut g, &mut PopOnce(0), PopOpts::default());
        assert_eq!(stats.pops, 1, "exactly the one pair the aligner accepted");
        assert!(!stats.hit_iteration_cap, "terminated on its own");
        assert_eq!(g.node_count(), 5, "the node set never changes");
        let _ = mids;
    }

    /// A `Consensus` that records what it was asked and answers predictably:
    /// spoa returns the first sequence, the aligner returns a caller-set CIGAR.
    struct StubConsensus {
        cigar: Vec<(u32, u8)>,
        spoa_calls: usize,
        align_calls: usize,
    }

    impl StubConsensus {
        fn new() -> Self {
            Self {
                cigar: vec![(100, b'=')],
                spoa_calls: 0,
                align_calls: 0,
            }
        }
    }

    impl Consensus for StubConsensus {
        fn spoa(&mut self, seqs: &[&[u8]]) -> Vec<u8> {
            self.spoa_calls += 1;
            seqs.first().map(|s| s.to_vec()).unwrap_or_default()
        }
        fn align(&mut self, _s1: &[u8], _s2: &[u8]) -> Vec<(u32, u8)> {
            self.align_calls += 1;
            self.cigar.clone()
        }
    }

    fn reads_map(pairs: &[(u32, &str)]) -> FxHashMap<u32, Vec<u8>> {
        pairs
            .iter()
            .map(|(r, s)| (*r, s.as_bytes().to_vec()))
            .collect()
    }

    #[test]
    fn two_supporting_reads_take_the_longer_span_and_never_call_spoa() {
        // The case a port can verify without reproducing spoa at all.
        let reads = reads_map(&[(1, "AAAACCCCGGGGTTTT"), (2, "TTTTGGGGCCCCAAAA")]);
        let mut engine = StubConsensus::new();
        // read 1 spans 0..4 (dist 4); read 2 spans 0..8 (dist 8) --- read 2 wins.
        let attrs = vec![(1u32, 0, 4), (2u32, 0, 8)];
        let con = generate_consensus_path(&mut engine, &reads, &attrs, 2);
        assert_eq!(con, b"TTTTGGGGCC".to_vec(), "read 2's span plus k");
        assert_eq!(engine.spoa_calls, 0, "two reads never reach spoa");
    }

    #[test]
    fn more_than_two_supporting_reads_go_to_spoa() {
        let reads = reads_map(&[(1, "AAAAAAAAAA"), (2, "CCCCCCCCCC"), (3, "GGGGGGGGGG")]);
        let mut engine = StubConsensus::new();
        let attrs = vec![(1u32, 0, 6), (2, 0, 6), (3, 0, 6)];
        let con = generate_consensus_path(&mut engine, &reads, &attrs, 2);
        assert_eq!(engine.spoa_calls, 1);
        assert_eq!(
            con,
            b"AAAAAAAA".to_vec(),
            "the stub returns the first sequence"
        );
    }

    #[test]
    fn spans_shorter_than_k_are_withheld_from_spoa_and_yield_an_x_placeholder() {
        // A span runs pos1..pos2 + k, so it is only shorter than k when the read
        // itself is --- the slice clamps. (A first draft used 10-base reads with
        // k = 8 and got 9-base spans, i.e. long enough for spoa, so it measured
        // the opposite of what it claimed.) Here the reads are 3 bases and k is 8,
        // so nothing is written and the reference returns "X" * max_len, which the
        // caller then rejects for being under 3 long.
        let reads = reads_map(&[(1, "AAA"), (2, "CCC"), (3, "GGG")]);
        let mut engine = StubConsensus::new();
        let attrs = vec![(1u32, 0, 1), (2, 0, 1), (3, 0, 1)];
        let con = generate_consensus_path(&mut engine, &reads, &attrs, 8);
        assert_eq!(engine.spoa_calls, 0);
        assert!(
            con.iter().all(|&b| b == b'X'),
            "placeholder, not a consensus"
        );
    }

    #[test]
    fn cigar_diversity_rejects_one_long_indel_but_tolerates_scattered_errors() {
        // A single non-match run longer than delta_len is structure, not error.
        assert!(!parse_cigar_diversity(&[(90, b'='), (10, b'I')], 0.20, 5));
        // The same total mismatch spread thinly is tolerated.
        let mut scattered = Vec::new();
        for _ in 0..5 {
            scattered.push((18, b'='));
            scattered.push((2, b'X'));
        }
        assert!(parse_cigar_diversity(&scattered, 0.20, 5));
    }

    #[test]
    fn cigar_diversity_uses_the_larger_of_the_two_thresholds() {
        // max(delta_perc * len, delta_len): on a short alignment delta_len wins,
        // so a proportionally large mismatch still passes.
        // len 20, mismatch 4 -> diversity 0.20; threshold max(0.20*20, 5)=5 -> 0.25.
        assert!(parse_cigar_diversity(&[(16, b'='), (4, b'X')], 0.20, 5));
        // Same shape, longer alignment: the percentage wins and 0.30 > 0.20.
        assert!(!parse_cigar_diversity(
            &[(70, b'='), (3, b'X'), (27, b'M')],
            0.20,
            2
        ));
        // An empty CIGAR: the reference divides by zero here; we refuse instead.
        assert!(!parse_cigar_diversity(&[], 0.20, 5));
    }

    #[test]
    fn m_counts_as_a_match_alongside_eq() {
        // Both '=' and 'M' are matches; only other ops count as mismatch.
        assert!(parse_cigar_diversity(&[(50, b'M'), (50, b'=')], 0.20, 5));
    }

    #[test]
    fn two_short_consensuses_pop_unconditionally() {
        let reads = reads_map(&[(1, "AAAA")]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 2, 10);
        // Both under delta_len = 10.
        assert!(a.decide(b"AAA", b"AAAAA"), "both short -> pop");
        assert_eq!(a.engine.align_calls, 0, "no alignment needed");
    }

    #[test]
    fn a_big_length_difference_is_rejected_without_aligning() {
        let reads = reads_map(&[(1, "A")]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 2, 5);
        // 100 vs 50: 50% difference, well over the 20% gate.
        let long = vec![b'A'; 100];
        let short = vec![b'A'; 50];
        assert!(!a.decide(&long, &short));
        assert_eq!(a.engine.align_calls, 0, "rejected on length alone");
    }

    #[test]
    fn similar_lengths_are_settled_by_the_alignment() {
        let reads = reads_map(&[(1, "A")]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 2, 5);
        let s1 = vec![b'A'; 100];
        let s2 = vec![b'A'; 95];
        assert!(a.decide(&s1, &s2), "a clean CIGAR pops");
        assert_eq!(a.engine.align_calls, 1);
        // Now make the aligner report one long indel.
        a.engine.cigar = vec![(90, b'='), (10, b'D')];
        assert!(!a.decide(&s1, &s2), "a long indel does not");
    }

    #[test]
    fn one_short_one_long_pops_only_within_delta_len() {
        let reads = reads_map(&[(1, "A")]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 2, 10);
        // shorter under delta_len, longer over: the third branch, comparing the
        // raw difference.
        assert!(a.decide(&[b'A'; 9], &[b'A'; 12]), "difference 3 < 10");
        assert!(!a.decide(&[b'A'; 5], &[b'A'; 40]), "difference 35 >= 10");
    }

    #[test]
    fn get_consensus_positions_reads_both_ends_and_skips_absent_reads() {
        let (mut g, s, _a, _b, t) = two_path_bubble();
        g.replace_reads(s, 1, ri(0, 5));
        g.set_read(s, 2, ri(0, 7));
        g.replace_reads(t, 1, ri(0, 90));
        g.set_read(t, 2, ri(0, 95));
        // Read 3 is asked for but absent from both ends.
        let attrs = get_consensus_positions(&g, s, t, &[1, 2, 3]);
        assert_eq!(attrs, vec![(1, 5, 90), (2, 7, 95)]);
    }

    #[test]
    fn the_megabubble_memo_only_applies_to_megabubbles() {
        let reads = reads_map(&[
            (1, "AAAAAAAAAAAAAAAAAAAA"),
            (2, "CCCCCCCCCCCCCCCCCCCC"),
            (3, "GGGGGGGGGGGGGGGGGGGG"),
        ]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 2, 5);
        let attrs = vec![(1u32, 0, 6), (2, 0, 6), (3, 0, 6)];
        // Not a megabubble: no caching, so spoa runs each time.
        a.consensus_for((0, 1), &attrs, false);
        a.consensus_for((0, 1), &attrs, false);
        assert_eq!(a.engine.spoa_calls, 2);
        // A megabubble: the second call is served from the memo.
        a.consensus_for((0, 1), &attrs, true);
        let after_first = a.engine.spoa_calls;
        a.consensus_for((0, 1), &attrs, true);
        assert_eq!(a.engine.spoa_calls, after_first, "memoised");
    }

    #[test]
    fn a_consensus_under_three_bases_counts_as_none() {
        let reads = reads_map(&[(1, "AAAAAAAAAA"), (2, "CC"), (3, "GG")]);
        let mut a = RealAligner::new(StubConsensus::new(), &reads, 1, 5);
        // Spans of length 1 with k = 1 give 2-base sequences: under k, so the
        // placeholder path, and the result is shorter than 3.
        let attrs = vec![(1u32, 0, 1), (2, 0, 1), (3, 0, 1)];
        let con = a.consensus_for((0, 1), &attrs, false);
        assert!(con.is_empty(), "too short to be a consensus");
    }

    #[test]
    fn the_real_engine_produces_a_consensus_and_a_cigar() {
        // Smoke test that spoa and parasail are wired up and behave. Not a
        // verification of either --- parasail's own tests do that, and the spoa
        // side still owes a differential check against the binary on isONform's
        // inputs.
        let mut e = SpoaParasail;
        let a = b"ACGTACGTACGTACGTACGT".to_vec();
        let b = b"ACGTACGTACGTACGTACGT".to_vec();
        let con = e.spoa(&[&a, &b]);
        assert!(!con.is_empty(), "spoa returned something");

        let cigar = e.align(&a, &b);
        assert!(!cigar.is_empty());
        // Identical sequences: one run of matches covering both.
        assert_eq!(cigar.len(), 1);
        assert_eq!(cigar[0], (a.len() as u32, b'='));
    }

    #[test]
    fn spoa_on_an_empty_input_gives_an_empty_consensus() {
        // `poa::consensus` returns None for no sequences, which the wrapper turns
        // into an empty consensus --- and the caller then treats anything under
        // three bases as no consensus at all.
        let mut e = SpoaParasail;
        assert!(e.spoa(&[]).is_empty());
    }

    #[test]
    fn the_real_engine_decides_a_bubble_end_to_end() {
        // A whole bubble through the real consensus and alignment path: two
        // near-identical branches should merge, two unrelated ones should not.
        let reads = reads_map(&[
            (1, "ACGTACGTACGTACGTACGTACGTACGTACGT"),
            (2, "ACGTACGTACGTACGTACGTACGTACGTACGT"),
            (3, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
        ]);
        let mut a = RealAligner::new(SpoaParasail, &reads, 4, 5);
        // Reads 1 and 2 are identical over the same span.
        let same = a.decide(&reads[&1][0..24], &reads[&2][0..24]);
        assert!(same, "identical branches merge");
        // Read 3 is unrelated over the same length.
        let different = a.decide(&reads[&1][0..24], &reads[&3][0..24]);
        assert!(!different, "unrelated branches do not");
    }

    #[test]
    fn dist_to_prev_is_the_mean_shift_and_ignores_order() {
        let (mut g, s, a, _b, _t) = diamond();
        // s: read 1 at 10, read 2 at 20; a: read 1 at 40, read 2 at 60.
        g.replace_reads(s, 1, ri(0, 10));
        g.set_read(s, 2, ri(0, 20));
        g.replace_reads(a, 1, ri(0, 40));
        g.set_read(a, 2, ri(0, 60));
        // shifts are 30 and 40 -> mean 35
        assert_eq!(get_dist_to_prev(&g, s, a), 35.0);
    }

    #[test]
    fn dist_to_prev_of_a_disjoint_pair_is_zero() {
        let (g, _s, a, b, _t) = diamond();
        // a supports reads 1,2 and b supports read 3 --- no overlap.
        assert_eq!(get_dist_to_prev(&g, a, b), 0.0);
    }

    #[test]
    fn avg_interval_length_truncates_toward_zero() {
        let mut g = Graph::new();
        let n = g.add_node(NodeKey::Source);
        // lengths 5 and 8 -> 13/2 = 6.5 -> 6
        g.set_read(n, 1, ri(10, 15));
        g.set_read(n, 2, ri(10, 18));
        assert_eq!(get_avg_interval_length(&g, n), 6);

        // Negative mean: Python's int() truncates toward zero (-3), where Rust's
        // integer division would floor to -4 if written as sum/len.
        let m = g.add_node(NodeKey::Sink);
        g.set_read(m, 1, ri(10, 5)); // -5
        g.set_read(m, 2, ri(10, 8)); // -2
        assert_eq!(
            get_avg_interval_length(&g, m),
            -3,
            "int() truncates, not floors"
        );
        assert_eq!((-7i64) / 2, -3, "Rust integer division also truncates here");
    }

    #[test]
    fn avg_interval_length_of_an_empty_node_is_zero() {
        let mut g = Graph::new();
        let n = g.add_node(NodeKey::Source);
        assert_eq!(get_avg_interval_length(&g, n), 0);
    }
}
