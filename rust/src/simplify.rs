//! Graph simplification: bubble detection.
//!
//! Ports the front of `modules/SimplifyGraph.py` — the part that decides *which*
//! bubbles exist. Popping them (`new_bubble_popping_routine`, ~270 lines of
//! in-place mutation plus spoa calls) is not ported yet, deliberately: the
//! detection side is pure, order-sensitive and testable in isolation, which makes
//! it the right thing to nail down first.
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

use rustc_hash::FxHashSet;

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
