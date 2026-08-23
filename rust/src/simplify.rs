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

use crate::graph::{Graph, NodeId};

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
            sum += ci.end_mini_start as i64 - pi.end_mini_start as i64;
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
        .map(|(_, i)| i.end_mini_start as i64 - i.start_mini_end as i64)
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
    pub length: i64,
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
        let prev_node = bubble_start;
        if let Some(length) = g.edge_length(prev_node, pnl_start) {
            lifted.put_edge(LiftedEdge {
                u: prev_node,
                v: pnl_start,
                length,
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
            if let Some(length) = g.edge_length(path_node, target) {
                lifted.put_edge(LiftedEdge {
                    u: path_node,
                    v: target,
                    length,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{NodeKey, ReadInfo};

    fn ri(a: u32, b: u32) -> ReadInfo {
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
                g.set_read(n, r, ri(0, r * 10));
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
        assert_eq!(e.length, 0);
        assert_eq!(e.support, vec![1, 2]);
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
        // b was not found, so the distance comes from get_dist_to_prev(s, b) ---
        // and b shares no read position with s beyond read 3, whose positions are
        // both 100, so the mean shift is 0.
        assert_eq!(
            lifted.dist(b),
            Some(0.0),
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
                .any(|e| e.u == s && e.v == t && e.length == 5),
            "the direct start->end edge is lifted too"
        );
        assert!(
            paths[1].nodes.is_empty(),
            "the direct path has no interior node"
        );
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
