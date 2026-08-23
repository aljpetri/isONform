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

/// Deliberate divergences, off by default.
#[derive(Debug, Clone, Copy, Default)]
pub struct BuildOpts {
    /// Bind `seq` from the read being processed, as every other branch does.
    /// Fixes finding 9.
    pub fix_stale_seq: bool,
    /// Continue the read's path from the node `cycle_added` just created, rather
    /// than from the one whose incoming edge it removed. Fixes finding 11.
    pub fix_cycle_continuation: bool,
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
    for &(r, len) in input.read_len {
        g.set_read(
            t,
            r,
            ReadInfo {
                start_mini_end: len,
                end_mini_start: len,
                original_support: true,
            },
        );
    }
    let reads_for_isoforms: Vec<u32> = (1..=n_graph_reads).collect();

    let mut prior_read_infos: FxHashMap<(u32, u32, u32), NodeId> = FxHashMap::default();
    // old_node -> [(alternative node, its predecessor, the predecessor edge's length)]
    let mut alternative_nodes: FxHashMap<NodeId, Vec<(NodeId, NodeId, i64)>> = FxHashMap::default();
    // The reference only ever reports `len(alt_cyc_nodes.keys())`.
    let mut alt_cyc_keys: FxHashSet<NodeId> = FxHashSet::default();

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

        for inter in intervals_for_read {
            let info_tuple = (r_id, inter.start, inter.end);
            let curr_hash = occurrence_hash(&inter.occurrences);
            let is_repetitive = read_hashes.contains(&curr_hash);
            let this_len = inter.start as i64 - previous_end;

            let node = match prior_read_infos.get(&info_tuple).copied() {
                Some(prior) if !is_repetitive => {
                    // The interval is known from an earlier read and does not
                    // repeat within this one.
                    if !g.has_edge(previous_node, prior) {
                        leaked_seq = Some(read_seq);
                        g.set_read(
                            prior,
                            r_id,
                            ReadInfo {
                                start_mini_end: inter.start,
                                end_mini_start: inter.end,
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
                                    start_mini_end: inter.start,
                                    end_mini_start: inter.end,
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
                                    start_mini_end: inter.start,
                                    end_mini_start: inter.end,
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
                                            start_mini_end: inter.start,
                                            end_mini_start: inter.end,
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
                                            start_mini_end: inter.start,
                                            end_mini_start: inter.end,
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
                            start_mini_end: inter.start,
                            end_mini_start: inter.end,
                            original_support: true,
                        },
                    );
                    g.add_edge(previous_node, n, this_len);
                    g.set_edge_support(previous_node, n, r_id);
                    n
                }
                None => {
                    // `new_interval_action`: first time this interval is seen.
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
                            start_mini_end: inter.start,
                            end_mini_start: inter.end,
                            original_support: true,
                        },
                    );
                    add_prior_read_infos(&mut prior_read_infos, inter, r_id, n, k);
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
        let (g, rfi) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
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
    fn the_sink_read_set_comes_from_read_len_not_from_the_graph() {
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
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
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
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
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
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
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
        let (g, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();

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
        let (g_ref, _) = generate_graph_from_intervals(&input, BuildOpts::default()).unwrap();
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
