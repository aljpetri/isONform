//! The interval graph, and its construction from WIS-selected intervals.
//!
//! Ports `modules/GraphGeneration.py::generateGraphfromIntervals`. The reference
//! is normative; where this file and it disagree, it is right — with the two
//! exceptions recorded under *Findings* in `PORTING.md`, which are reproduced
//! here on purpose and called out in comments.
//!
//! # Representation, and why not `petgraph`
//!
//! The reference builds a `networkx.DiGraph` whose node identifiers are
//! **formatted strings**, `"{start}, {end}, {r_id}"`, plus `"s"` and `"t"`. Every
//! node carries a `reads` dict and an `end_mini_seq` string; every edge carries
//! `length` and an `edge_supp` list. That is a Python object per node, a dict per
//! node and a dict per edge, and it is a large part of why the reference is slow.
//!
//! Here: nodes are `u32` indices. A node's identity is the tuple the name is
//! formatted *from* ([`NodeKey`]), so the map is keyed on three integers rather
//! than on a heap string — exactly equivalent, since the formatting is injective.
//! Attributes live in parallel arrays indexed by node id.
//!
//! Two places deliberately keep a `Vec` where a `HashMap` looks natural:
//!
//! * a node's read list, and an edge's support list, are `Vec`s searched
//!   linearly. They hold a handful to a few tens of entries, so a linear scan
//!   beats hashing, and it preserves insertion order for free.
//!
//!   **Is that order observable? Checked, and as the code stands: no.** All five
//!   consumers of `DG.nodes[node]['reads']` in `SimplifyGraph.py` are
//!   order-independent — `:84` and `:103` convert to a `set` and then sort;
//!   `get_dist_to_prev` (`:218`) and `get_avg_interval_length` (`:295`) sum and
//!   divide, which commutes; `additional_node_support` (`:311`) writes one
//!   independent value per read id; `:491` merges dicts whose result only feeds
//!   the above. So order is preserved because it is free and faithful, **not**
//!   because it is known to reach results. An earlier version of this comment
//!   claimed the latter and was wrong.
//! * adjacency is `Vec<Vec<NodeId>>` rather than CSR, because construction is
//!   incremental and one branch *removes* an edge. CSR is the right shape for the
//!   simplification stage, which only reads; freezing to CSR is where that
//!   conversion belongs when it is written.

use std::collections::VecDeque;
use std::fmt;

use rustc_hash::FxHashMap;

pub type NodeId = u32;
pub type EdgeId = u32;

/// A node's identity.
///
/// The reference names nodes `"{start}, {end}, {r_id}"`, so two distinct triples
/// always give distinct names and vice versa: keying on the triple is the same
/// graph with none of the string formatting. [`fmt::Display`] reproduces the
/// reference's spelling, which is what lets a dump be diffed against it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NodeKey {
    Source,
    Sink,
    Interval { start: u32, end: u32, r_id: u32 },
    /// Identity by **content**: a hash of the interval's flanking minimizer pair
    /// `(m1, m2)`, used when `ISONFORM_PAIR_NODES=1`.
    ///
    /// The `Interval` key is per-read coordinates, which makes node identity
    /// depend on coordinate drift, on read order, and on file position --- moving
    /// two reads from the front of a file to the end changes the graph from 26 507
    /// nodes to 500 (finding 48). A pair hash depends on none of those.
    Pair { h: u64 },
}

impl fmt::Display for NodeKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NodeKey::Source => f.write_str("s"),
            NodeKey::Sink => f.write_str("t"),
            // `str(a) + ", " + str(b) + ", " + str(c)` in the reference.
            NodeKey::Interval { start, end, r_id } => write!(f, "{start}, {end}, {r_id}"),
            NodeKey::Pair { h } => write!(f, "p{h:016x}"),
        }
    }
}

/// The reference's `Read_infos` namedtuple: `(start_mini_end, end_mini_start,
/// original_support)`.
///
/// The positions are **signed**, and that is not defensive typing. Graph
/// construction only ever stores non-negative coordinates, but
/// `additional_node_support` invents positions for reads that reach a node via
/// the other branch of a bubble — `newend = prev_end + relative_dist`, then
/// `newstart = newend - avg_len` — and both can come out negative. Measured on 16
/// recorded Drosophila simplifications: **57 nodes carry a negative position
/// afterwards, and none did before**. So a position here is a virtual coordinate
/// that may precede the read's own start, and anything consuming it has to
/// tolerate that.
///
/// An earlier version of this port clamped them to zero. That was a divergence
/// introduced without measuring, and the simplification oracle caught it on its
/// first run against real data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadInfo {
    pub start_mini_end: i64,
    pub end_mini_start: i64,
    pub original_support: bool,
}

/// One interval from the WIS solution: `(start, end, weight, occurrences)`.
///
/// `occurrences` is the reference's `inter[3]`, a flat `array("I")` of
/// `(r_id, pos1, pos2)` triples. **The reference mutates it**:
/// `convert_array_to_hash` pops the first three elements before hashing
/// (`GraphGeneration.py:129-131`), so the hash is over the *other* reads'
/// occurrences only. That is reproduced in [`occurrence_hash`], which reads from
/// index 3 rather than mutating.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub start: u32,
    pub end: u32,
    pub weight: u32,
    pub occurrences: Vec<u32>,
}

/// The graph.
#[derive(Debug, Default)]
pub struct Graph {
    keys: Vec<NodeKey>,
    index: FxHashMap<NodeKey, NodeId>,

    /// Per-node, parallel to `keys`.
    end_mini_seq: Vec<Vec<u8>>,
    reads: Vec<Vec<(u32, ReadInfo)>>,
    succ: Vec<Vec<NodeId>>,
    /// In-degree, maintained incrementally. networkx's topological sort reads
    /// `G.in_degree()`, and reproducing its order needs this cheaply.
    in_deg: Vec<u32>,

    /// Edges. `edge_index` maps `(u, v)` to a slot in the three parallel arrays.
    edge_index: FxHashMap<(NodeId, NodeId), EdgeId>,
    edge_uv: Vec<(NodeId, NodeId)>,
    /// `None` means the edge carries no `length` attribute at all.
    ///
    /// Not hypothetical: `prepare_adding_edges` re-adds bubble edges with
    /// `DG.add_edge(u, v, edge_supp=value)` and no `length`, so after
    /// simplification a large share of edges simply lack it — measured at 97 of
    /// 231 on one recorded Drosophila graph. Nothing reads `length` after
    /// construction (`GraphGeneration.py:402` is its only consumer), so this is
    /// harmless, but it *is* observable and the dump records it as `NA`.
    edge_len: Vec<Option<i64>>,
    edge_supp: Vec<Vec<u32>>,
    /// Slots freed by `remove_edge`, reused so ids stay dense.
    free_edges: Vec<EdgeId>,
}

impl Graph {
    pub fn new() -> Self {
        Self::default()
    }

    // -- nodes ------------------------------------------------------------

    pub fn node_count(&self) -> usize {
        self.keys.len()
    }

    pub fn edge_count(&self) -> usize {
        self.edge_index.len()
    }

    pub fn key(&self, n: NodeId) -> NodeKey {
        self.keys[n as usize]
    }

    pub fn lookup(&self, key: &NodeKey) -> Option<NodeId> {
        self.index.get(key).copied()
    }

    /// Add a node, or return the existing one. `DG.add_node` is idempotent in
    /// networkx and leaves existing attributes alone, which this matches.
    pub fn add_node(&mut self, key: NodeKey) -> NodeId {
        if let Some(&n) = self.index.get(&key) {
            return n;
        }
        let n = self.keys.len() as NodeId;
        self.keys.push(key);
        self.index.insert(key, n);
        self.end_mini_seq.push(Vec::new());
        self.reads.push(Vec::new());
        self.succ.push(Vec::new());
        self.in_deg.push(0);
        n
    }

    pub fn set_end_mini_seq(&mut self, n: NodeId, seq: &[u8]) {
        self.end_mini_seq[n as usize].clear();
        self.end_mini_seq[n as usize].extend_from_slice(seq);
    }

    pub fn end_mini_seq(&self, n: NodeId) -> &[u8] {
        &self.end_mini_seq[n as usize]
    }

    /// `nodes_for_graph[name][r_id] = r_infos`.
    ///
    /// Python dict assignment: an existing key keeps its position and takes the
    /// new value; a new key is appended. Reproduced exactly — see the module
    /// docs for whether that order is observable (as the code stands, it is
    /// not; it is kept because it is free and faithful).
    pub fn set_read(&mut self, n: NodeId, r_id: u32, info: ReadInfo) {
        let list = &mut self.reads[n as usize];
        if let Some(slot) = list.iter_mut().find(|(r, _)| *r == r_id) {
            slot.1 = info;
        } else {
            list.push((r_id, info));
        }
    }

    /// `nodes_for_graph[name] = {r_id: r_infos}` — a *fresh* dict.
    ///
    /// Distinct from [`Graph::set_read`], and the distinction is load-bearing.
    /// Several branches of the reference build a new local `nodelist` and assign
    /// it wholesale, which **discards** whatever read map the node already had.
    /// Mostly those branches create a new node, so the map was empty anyway —
    /// but `cycle_added` and the repetitive-interval branch both reuse a name
    /// that may already exist (`if not DG.has_node(name)`), so the discard is
    /// reachable. Reproduced rather than softened to a merge.
    pub fn replace_reads(&mut self, n: NodeId, r_id: u32, info: ReadInfo) {
        let list = &mut self.reads[n as usize];
        list.clear();
        list.push((r_id, info));
    }

    /// Drop every read from a node, for `set_node_attributes` with an empty map.
    pub fn clear_reads(&mut self, n: NodeId) {
        self.reads[n as usize].clear();
    }

    pub fn reads(&self, n: NodeId) -> &[(u32, ReadInfo)] {
        &self.reads[n as usize]
    }

    pub fn successors(&self, n: NodeId) -> &[NodeId] {
        &self.succ[n as usize]
    }

    // -- edges ------------------------------------------------------------

    pub fn has_edge(&self, u: NodeId, v: NodeId) -> bool {
        self.edge_index.contains_key(&(u, v))
    }

    /// The edge's `length`, or `None` if the edge is absent **or** carries no
    /// `length`. Use [`Graph::has_edge`] to tell those apart.
    pub fn edge_length(&self, u: NodeId, v: NodeId) -> Option<i64> {
        self.edge_index
            .get(&(u, v))
            .and_then(|&e| self.edge_len[e as usize])
    }

    pub fn edge_support(&self, u: NodeId, v: NodeId) -> Option<&[u32]> {
        self.edge_index
            .get(&(u, v))
            .map(|&e| self.edge_supp[e as usize].as_slice())
    }

    /// `DG.add_edge(u, v, length=length)`. Idempotent on `(u, v)`, overwriting
    /// `length`, as networkx is.
    pub fn add_edge(&mut self, u: NodeId, v: NodeId, length: i64) -> EdgeId {
        self.add_edge_inner(u, v, Some(length))
    }

    /// `DG.add_edge(u, v, edge_supp=...)` — no `length` argument.
    ///
    /// networkx merges attribute dicts, so this **creates** an edge with no
    /// `length` when the pair is absent and **preserves** the existing `length`
    /// when it is present. `prepare_adding_edges` relies on both halves.
    /// `support` replaces whatever support the edge had.
    pub fn upsert_edge_support(&mut self, u: NodeId, v: NodeId, support: Vec<u32>) -> EdgeId {
        let e = match self.edge_index.get(&(u, v)) {
            Some(&e) => e, // existing edge: length untouched
            None => self.add_edge_inner(u, v, None),
        };
        self.edge_supp[e as usize] = support;
        e
    }

    fn add_edge_inner(&mut self, u: NodeId, v: NodeId, length: Option<i64>) -> EdgeId {
        if let Some(&e) = self.edge_index.get(&(u, v)) {
            self.edge_len[e as usize] = length;
            return e;
        }
        let e = match self.free_edges.pop() {
            Some(e) => {
                self.edge_uv[e as usize] = (u, v);
                self.edge_len[e as usize] = length;
                self.edge_supp[e as usize].clear();
                e
            }
            None => {
                self.edge_uv.push((u, v));
                self.edge_len.push(length);
                self.edge_supp.push(Vec::new());
                (self.edge_uv.len() - 1) as EdgeId
            }
        };
        self.edge_index.insert((u, v), e);
        self.succ[u as usize].push(v);
        self.in_deg[v as usize] += 1;
        e
    }

    /// `DG.remove_edge(u, v)`. Only `cycle_added` needs this.
    pub fn remove_edge(&mut self, u: NodeId, v: NodeId) -> bool {
        match self.edge_index.remove(&(u, v)) {
            Some(e) => {
                // Remove one occurrence; adjacency never holds duplicates
                // because add_edge is idempotent on (u, v).
                if let Some(pos) = self.succ[u as usize].iter().position(|&w| w == v) {
                    self.succ[u as usize].remove(pos);
                }
                self.in_deg[v as usize] -= 1;
                self.edge_supp[e as usize].clear();
                self.free_edges.push(e);
                true
            }
            None => false,
        }
    }

    /// `edge_support[u, v] = [r_id]` — replaces any existing list.
    ///
    /// The reference's `add_edge_support` assigns a fresh single-element list
    /// rather than appending, so an earlier support list is discarded. Faithful.
    pub fn set_edge_support(&mut self, u: NodeId, v: NodeId, r_id: u32) {
        if let Some(&e) = self.edge_index.get(&(u, v)) {
            let s = &mut self.edge_supp[e as usize];
            s.clear();
            s.push(r_id);
        }
    }

    /// Append `r_id` to an edge's support if not already present, preserving
    /// order.
    pub fn push_edge_support(&mut self, u: NodeId, v: NodeId, r_id: u32) {
        if let Some(&e) = self.edge_index.get(&(u, v)) {
            let s = &mut self.edge_supp[e as usize];
            if !s.contains(&r_id) {
                s.push(r_id);
            }
        }
    }

    pub fn out_degree(&self, n: NodeId) -> usize {
        self.succ[n as usize].len()
    }

    pub fn in_degree(&self, n: NodeId) -> usize {
        self.in_deg[n as usize] as usize
    }

    // -- topological order ------------------------------------------------

    /// Reproduces `networkx.topological_sort` as of 2.8.4 — which matters,
    /// because a topological order is **not unique** and the reference compares
    /// the resulting *indices* to decide which node pairs are candidate bubbles
    /// (`SimplifyGraph.py:115`). A different valid order finds different bubbles.
    ///
    /// In 2.8.4 `topological_sort` simply yields from `topological_generations`,
    /// so it is generation-by-generation Kahn, **not** the LIFO variant:
    ///
    /// ```text
    /// zero_indegree = [v for v, d in G.in_degree() if d == 0]
    /// while zero_indegree:
    ///     this_generation, zero_indegree = zero_indegree, []
    ///     for node in this_generation:
    ///         for child in G.neighbors(node):
    ///             indegree_map[child] -= 1
    ///             if indegree_map[child] == 0:
    ///                 zero_indegree.append(child)
    ///     yield this_generation
    /// ```
    ///
    /// The order within a generation therefore comes from **node insertion
    /// order** for the first generation and from parent-order × successor-order
    /// after that. Both are preserved here: `keys` is in insertion order and
    /// `succ[u]` is in edge-insertion order. This is the place where keeping
    /// insertion order genuinely reaches results.
    ///
    /// Returns `None` if the graph has a cycle, as networkx raises
    /// `NetworkXUnfeasible`.
    /// A degree profile, for measuring what simplification actually did.
    ///
    /// Node and edge counts are the wrong instrument. `SimplifyGraph.py` calls
    /// only `remove_edge`/`add_edge` and **never** touches nodes, so popping
    /// cannot change the node count at all --- measuring nodes alone made a graph
    /// that had been popped heavily look untouched. Edge count moves, but it
    /// cannot distinguish "a bubble was linearised" from "an edge was dropped".
    ///
    /// Degrees can. A bubble is a node with out-degree > 1 reconverging at one
    /// with in-degree > 1; linearising it turns interior nodes into in-1/out-1;
    /// and because the abandoned path's nodes are left behind rather than
    /// removed, they fall to degree 0. So:
    ///
    /// * `branch` and `merge` falling means bubbles were resolved;
    /// * `linear` rising means paths were straightened;
    /// * `orphan` rising counts what popping left behind --- nodes the reference
    ///   never deletes.
    ///
    /// Returns `(branch, merge, linear, orphan)`. A node can be both a branch and
    /// a merge, so those two overlap; `linear` and `orphan` are disjoint from
    /// both.
    pub fn degree_profile(&self) -> (usize, usize, usize, usize) {
        let (mut branch, mut merge, mut linear, mut orphan) = (0, 0, 0, 0);
        for n in 0..self.node_count() as u32 {
            let (i, o) = (self.in_degree(n), self.out_degree(n));
            if o > 1 {
                branch += 1;
            }
            if i > 1 {
                merge += 1;
            }
            if i == 1 && o == 1 {
                linear += 1;
            }
            if i == 0 && o == 0 {
                orphan += 1;
            }
        }
        (branch, merge, linear, orphan)
    }

    pub fn topological_sort(&self) -> Option<Vec<NodeId>> {
        let n = self.keys.len();
        let mut indegree: Vec<u32> = self.in_deg.clone();
        let mut this_generation: Vec<NodeId> = (0..n as NodeId)
            .filter(|&v| indegree[v as usize] == 0)
            .collect();
        let mut out: Vec<NodeId> = Vec::with_capacity(n);

        while !this_generation.is_empty() {
            let mut next: Vec<NodeId> = Vec::new();
            for &node in &this_generation {
                for &child in &self.succ[node as usize] {
                    let d = &mut indegree[child as usize];
                    *d -= 1;
                    if *d == 0 {
                        next.push(child);
                    }
                }
            }
            out.extend_from_slice(&this_generation);
            this_generation = next;
        }

        if out.len() == n {
            Some(out)
        } else {
            None // a cycle: some nodes never reached in-degree zero
        }
    }

    /// `{node: index}` over [`Graph::topological_sort`], which is what the
    /// reference actually uses (`topo_nodes_dict`).
    pub fn topological_index(&self) -> Option<Vec<u32>> {
        let order = self.topological_sort()?;
        let mut idx = vec![u32::MAX; self.keys.len()];
        for (i, &n) in order.iter().enumerate() {
            idx[n as usize] = i as u32;
        }
        Some(idx)
    }

    // -- cycle detection --------------------------------------------------

    /// Reproduces `GraphGeneration.bfs(DG, startnode)`.
    ///
    /// Forward BFS from `start`, returning true as soon as `start` is found
    /// among some node's successors — i.e. "is `start` on a cycle". Note it does
    /// **not** detect cycles that avoid `start`, which is all the reference asks
    /// of it, and the `neighbour == startnode` test happens before the visited
    /// check, so a self-loop reports true.
    pub fn reaches_itself(&self, start: NodeId) -> bool {
        let mut visited = vec![false; self.keys.len()];
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(m) = queue.pop_front() {
            for &nb in &self.succ[m as usize] {
                if !visited[nb as usize] {
                    visited[nb as usize] = true;
                    queue.push_back(nb);
                }
                if nb == start {
                    return true;
                }
            }
        }
        false
    }

    // -- output -----------------------------------------------------------

    /// Node ids sorted by the reference's node *name*, so a dump matches
    /// `sorted(dg.nodes())` in `bench/dump_reference.py`.
    pub fn nodes_sorted_by_name(&self) -> Vec<NodeId> {
        let mut ids: Vec<NodeId> = (0..self.keys.len() as NodeId).collect();
        ids.sort_by_cached_key(|&n| self.keys[n as usize].to_string());
        ids
    }

    /// Edges sorted by `(name(u), name(v))`, matching `sorted(dg.edges())`.
    pub fn edges_sorted_by_name(&self) -> Vec<(NodeId, NodeId)> {
        let mut es: Vec<(NodeId, NodeId)> = self.edge_index.keys().copied().collect();
        es.sort_by_cached_key(|&(u, v)| {
            (
                self.keys[u as usize].to_string(),
                self.keys[v as usize].to_string(),
            )
        });
        es
    }
}

/// The hash the reference computes over an interval's occurrence array.
///
/// `convert_array_to_hash` pops the first three elements — this occurrence's own
/// `(r_id, pos1, pos2)` — and hashes `tuple(rest)`. It is used only for equality
/// ("has this interval been seen before in this read"), never for ordering or
/// bucketing that reaches output, so any injective-enough function will do and
/// CPython's tuple hash need not be reproduced. That is worth stating because
/// `hash()` on a *string* was a genuine determinism defect here (`PORTING.md`
/// finding 1) while this one is not: `hash()` over a tuple of ints is not
/// seed-randomised, verified across four `PYTHONHASHSEED` values.
///
/// Reads from index 3 rather than mutating the input, which is the one place this
/// port deliberately does not copy the reference's *mechanism* — the observable
/// result is identical, and a destructive read of a caller's buffer is not
/// something to carry forward.
pub fn occurrence_hash(occurrences: &[u32]) -> u64 {
    // FNV-1a. Cheap, and collisions only cost a spurious "repetitive" verdict,
    // which the reference's own hash is equally exposed to.
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &x in occurrences.iter().skip(3) {
        for b in x.to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100_0000_01b3);
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_names_match_the_reference_spelling() {
        assert_eq!(NodeKey::Source.to_string(), "s");
        assert_eq!(NodeKey::Sink.to_string(), "t");
        assert_eq!(
            NodeKey::Interval {
                start: 40,
                end: 97,
                r_id: 3
            }
            .to_string(),
            "40, 97, 3"
        );
    }

    #[test]
    fn add_node_is_idempotent_and_keeps_attributes() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Interval {
            start: 1,
            end: 2,
            r_id: 3,
        });
        g.set_end_mini_seq(a, b"ACGT");
        let again = g.add_node(NodeKey::Interval {
            start: 1,
            end: 2,
            r_id: 3,
        });
        assert_eq!(a, again);
        assert_eq!(g.end_mini_seq(a), b"ACGT");
        assert_eq!(g.node_count(), 1);
    }

    #[test]
    fn set_read_overwrites_in_place_and_appends_new_in_order() {
        // Python dict semantics, and the order is contract: SimplifyGraph does
        // tuple(DG.nodes[node]['reads']).
        let mut g = Graph::new();
        let n = g.add_node(NodeKey::Source);
        let ri = |a: i64, b: i64| ReadInfo {
            start_mini_end: a,
            end_mini_start: b,
            original_support: true,
        };
        g.set_read(n, 7, ri(1, 2));
        g.set_read(n, 3, ri(3, 4));
        g.set_read(n, 7, ri(9, 9)); // overwrite, must not move to the back
        let got: Vec<u32> = g.reads(n).iter().map(|(r, _)| *r).collect();
        assert_eq!(got, vec![7, 3]);
        assert_eq!(g.reads(n)[0].1, ri(9, 9));
    }

    #[test]
    fn replace_reads_discards_the_existing_map() {
        let mut g = Graph::new();
        let n = g.add_node(NodeKey::Source);
        let ri = |a: i64, b: i64| ReadInfo {
            start_mini_end: a,
            end_mini_start: b,
            original_support: true,
        };
        g.set_read(n, 1, ri(1, 1));
        g.set_read(n, 2, ri(2, 2));
        g.replace_reads(n, 9, ri(9, 9));
        let got: Vec<u32> = g.reads(n).iter().map(|(r, _)| *r).collect();
        assert_eq!(
            got,
            vec![9],
            "a fresh nodelist assignment drops prior reads"
        );
    }

    #[test]
    fn add_edge_is_idempotent_and_overwrites_length() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        let b = g.add_node(NodeKey::Sink);
        g.add_edge(a, b, 5);
        g.add_edge(a, b, 9);
        assert_eq!(g.edge_count(), 1);
        assert_eq!(g.edge_length(a, b), Some(9));
        assert_eq!(g.successors(a), &[b]); // adjacency must not gain a duplicate
    }

    #[test]
    fn upsert_edge_support_creates_without_a_length_and_preserves_an_existing_one() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        let b = g.add_node(NodeKey::Sink);
        // Absent: created with no length at all --- what `prepare_adding_edges`
        // does to every bubble edge it re-adds.
        g.upsert_edge_support(a, b, vec![1, 2]);
        assert!(g.has_edge(a, b));
        assert_eq!(g.edge_length(a, b), None, "no length attribute");
        assert_eq!(g.edge_support(a, b), Some(&[1u32, 2][..]));

        // Present with a length: networkx merges attribute dicts, so the length
        // survives and only the support is replaced.
        let mut g2 = Graph::new();
        let a = g2.add_node(NodeKey::Source);
        let b = g2.add_node(NodeKey::Sink);
        g2.add_edge(a, b, 42);
        g2.upsert_edge_support(a, b, vec![7]);
        assert_eq!(g2.edge_length(a, b), Some(42), "existing length preserved");
        assert_eq!(g2.edge_support(a, b), Some(&[7u32][..]));
    }

    #[test]
    fn edge_length_none_means_absent_or_unlabelled_and_has_edge_tells_them_apart() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        let b = g.add_node(NodeKey::Sink);
        assert_eq!(g.edge_length(a, b), None);
        assert!(!g.has_edge(a, b), "absent");
        g.upsert_edge_support(a, b, vec![]);
        assert_eq!(g.edge_length(a, b), None);
        assert!(g.has_edge(a, b), "present but unlabelled");
    }

    #[test]
    fn set_edge_support_replaces_and_push_appends_once() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        let b = g.add_node(NodeKey::Sink);
        g.add_edge(a, b, 0);
        g.set_edge_support(a, b, 1);
        g.push_edge_support(a, b, 2);
        g.push_edge_support(a, b, 2); // already there
        assert_eq!(g.edge_support(a, b), Some(&[1u32, 2][..]));
        // add_edge_support in the reference *replaces* the list.
        g.set_edge_support(a, b, 9);
        assert_eq!(g.edge_support(a, b), Some(&[9u32][..]));
    }

    #[test]
    fn remove_edge_clears_adjacency_and_index() {
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        let b = g.add_node(NodeKey::Sink);
        g.add_edge(a, b, 3);
        assert!(g.remove_edge(a, b));
        assert!(!g.has_edge(a, b));
        assert!(g.successors(a).is_empty());
        assert_eq!(g.edge_count(), 0);
        assert!(!g.remove_edge(a, b)); // second removal is a no-op
                                       // The freed slot is reused rather than leaked.
        g.add_edge(a, b, 4);
        assert_eq!(g.edge_length(a, b), Some(4));
    }

    #[test]
    fn reaches_itself_detects_only_cycles_through_the_start() {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let c = g.add_node(k(3));
        g.add_edge(a, b, 0);
        g.add_edge(b, c, 0);
        assert!(!g.reaches_itself(a));
        g.add_edge(c, b, 0); // cycle b->c->b, which does NOT contain a
        assert!(
            !g.reaches_itself(a),
            "bfs() only reports cycles through start"
        );
        assert!(g.reaches_itself(b));
        assert!(g.reaches_itself(c));
    }

    #[test]
    fn reaches_itself_reports_a_self_loop() {
        // The reference tests `neighbour == startnode` before the visited check,
        // so a self-loop counts.
        let mut g = Graph::new();
        let a = g.add_node(NodeKey::Source);
        g.add_edge(a, a, 0);
        assert!(g.reaches_itself(a));
    }

    #[test]
    fn occurrence_hash_skips_this_occurrence() {
        // The reference pops the first three entries before hashing, so two
        // intervals differing only in their own (r_id, p1, p2) hash the same.
        let a = occurrence_hash(&[1, 10, 20, 5, 50, 60]);
        let b = occurrence_hash(&[9, 90, 99, 5, 50, 60]);
        assert_eq!(a, b);
        let c = occurrence_hash(&[1, 10, 20, 6, 50, 60]);
        assert_ne!(a, c);
    }

    #[test]
    fn occurrence_hash_of_a_lone_occurrence_is_the_empty_hash() {
        assert_eq!(occurrence_hash(&[1, 2, 3]), occurrence_hash(&[]));
    }

    #[test]
    fn degrees_track_edges() {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let (a, b, c) = (g.add_node(k(1)), g.add_node(k(2)), g.add_node(k(3)));
        g.add_edge(a, c, 0);
        g.add_edge(b, c, 0);
        assert_eq!((g.out_degree(a), g.in_degree(a)), (1, 0));
        assert_eq!((g.out_degree(c), g.in_degree(c)), (0, 2));
        g.add_edge(a, c, 7); // idempotent: must not double-count
        assert_eq!(g.in_degree(c), 2);
        g.remove_edge(b, c);
        assert_eq!(g.in_degree(c), 1);
    }

    #[test]
    fn topological_sort_is_generation_based_not_lifo() {
        // The distinction is the whole point. Graph: a->c, b->c, c->d, plus an
        // isolated-source e added last.
        //   generation 0 = [a, b, e]  (node insertion order among in-degree 0)
        //   generation 1 = [c]
        //   generation 2 = [d]
        // A LIFO Kahn would emit e, b, a, ... instead.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let a = g.add_node(k(1));
        let b = g.add_node(k(2));
        let c = g.add_node(k(3));
        let d = g.add_node(k(4));
        let e = g.add_node(k(5));
        g.add_edge(a, c, 0);
        g.add_edge(b, c, 0);
        g.add_edge(c, d, 0);
        let order = g.topological_sort().expect("acyclic");
        assert_eq!(order, vec![a, b, e, c, d]);
    }

    #[test]
    fn within_a_generation_order_follows_parent_then_successor_order() {
        // s has successors x then y; both become ready in that order.
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let s = g.add_node(k(0));
        // Add y as a node BEFORE x, so node insertion order and successor order
        // disagree --- successor order must win for generation 1.
        let y = g.add_node(k(2));
        let x = g.add_node(k(1));
        g.add_edge(s, x, 0);
        g.add_edge(s, y, 0);
        let order = g.topological_sort().expect("acyclic");
        assert_eq!(order, vec![s, x, y], "successor insertion order decides");
    }

    #[test]
    fn topological_sort_reports_a_cycle_as_none() {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let (a, b) = (g.add_node(k(1)), g.add_node(k(2)));
        g.add_edge(a, b, 0);
        g.add_edge(b, a, 0);
        assert!(g.topological_sort().is_none());
    }

    #[test]
    fn topological_index_maps_node_to_position() {
        let mut g = Graph::new();
        let k = |i| NodeKey::Interval {
            start: i,
            end: i,
            r_id: 0,
        };
        let (a, b) = (g.add_node(k(1)), g.add_node(k(2)));
        g.add_edge(a, b, 0);
        let idx = g.topological_index().unwrap();
        assert!(idx[a as usize] < idx[b as usize]);
    }

    #[test]
    fn sorting_matches_pythons_string_sort() {
        // The dump sorts by name, and these are *strings*: "100, ..." sorts
        // before "40, ...". Getting this wrong would make every dump diff noise.
        let mut g = Graph::new();
        g.add_node(NodeKey::Interval {
            start: 40,
            end: 1,
            r_id: 1,
        });
        g.add_node(NodeKey::Interval {
            start: 100,
            end: 1,
            r_id: 1,
        });
        g.add_node(NodeKey::Source);
        let names: Vec<String> = g
            .nodes_sorted_by_name()
            .into_iter()
            .map(|n| g.key(n).to_string())
            .collect();
        assert_eq!(names, vec!["100, 1, 1", "40, 1, 1", "s"]);
    }
}
