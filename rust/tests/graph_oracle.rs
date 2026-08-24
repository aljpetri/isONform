//! Differential oracle for graph construction.
//!
//! Replays the inputs recorded by `bench/dump_reference.py` through
//! [`isonform::graph_build::generate_graph_from_intervals`] and diffs the result
//! against the graph the reference produced. This is method point 2: an
//! end-to-end comparison tells you *that* something is wrong, never *where*.
//!
//! ```text
//! bench/dump_reference.py --fastq-folder bench/corpus/sirv_small --outdir /tmp/d
//! ISONFORM_GRAPH_DUMPS=/tmp/d cargo test --test graph_oracle -- --nocapture
//! ```
//!
//! With the variable unset the test **skips loudly** rather than passing
//! silently: a green oracle that never ran is worse than a red one. CI generates
//! the dumps and sets it, so it is a real gate there.
//!
//! # What is compared, and what would be missed if it were sloppier
//!
//! Everything the reference's graph carries, structurally *and in order*:
//!
//! * the node set, by the reference's own node names;
//! * each node's `end_mini_seq`;
//! * each node's `reads` map **in insertion order**. Stricter than the reference
//!   strictly requires — that order turns out not to be observable, see the
//!   audit in `graph.rs` — but comparing it costs nothing and would catch a
//!   divergence the moment a consumer started depending on it. The report says
//!   `ORDER ONLY` when order is the only difference, so a failure of this kind
//!   is not mistaken for a structural one;
//! * the edge set, each edge's `length`, and each edge's `edge_supp` in
//!   insertion order;
//! * `reads_for_isoforms`;
//! * and the **exact `nx.topological_sort` order**, not merely a valid one. A
//!   topological order is not unique, and the reference compares the resulting
//!   indices to decide which node pairs are candidate bubbles
//!   (`SimplifyGraph.py:115`), so a port emitting a different valid order would
//!   find different bubbles and this comparison is what catches that.
//!
//! Node ids differ between the two — the reference uses formatted strings, the
//! port uses `u32` — so the comparison maps through the name, which is what makes
//! the representation change invisible to the diff.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isonform::graph::{Interval, ReadInfo};
use isonform::graph_build::{generate_graph_from_intervals, BuildInput, BuildOpts};
use rustc_hash::FxHashMap;

/// One recorded call: the inputs, and the graph the reference built from them.
#[derive(Default)]
struct Case {
    path: PathBuf,
    k: usize,
    delta_len: i64,
    reads: FxHashMap<u32, Vec<u8>>,
    read_len: Vec<(u32, u32)>,
    intervals: Vec<(u32, Vec<Interval>)>,
    /// name -> (end_mini_seq, reads in insertion order)
    exp_nodes: BTreeMap<String, (String, Vec<(u32, ReadInfo)>)>,
    /// (u, v) -> (length, support in insertion order)
    exp_edges: BTreeMap<(String, String), (i64, Vec<u32>)>,
    exp_reads_for_isoforms: Vec<u32>,
    exp_counts: (usize, usize),
    /// `nx.topological_sort` order, by node name.
    ///
    /// Three states, and conflating any two of them is a bug this file has
    /// already had: `None` means the dump carries no `T` record at all (it
    /// predates the record, so the order is simply unchecked), `Some(None)` means
    /// the reference reported a cycle, `Some(Some(order))` is an order to
    /// compare. Reading "absent" as "cyclic" made 30 of 30 real cases fail with
    /// a message that pointed nowhere.
    exp_topo: Option<Option<Vec<String>>>,
}

fn parse_reads_attr(field: &str) -> Vec<(u32, ReadInfo)> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .map(|item| {
            let mut p = item.split(':');
            let r: u32 = p.next().unwrap().parse().unwrap();
            // Signed: positions really do go negative after a bubble is popped.
            let a: i64 = p.next().unwrap().parse().unwrap();
            let b: i64 = p.next().unwrap().parse().unwrap();
            let o: u8 = p.next().unwrap().parse().unwrap();
            (
                r,
                ReadInfo {
                    start_mini_end: a,
                    end_mini_start: b,
                    original_support: o != 0,
                },
            )
        })
        .collect()
}

fn parse_case(path: &Path) -> Case {
    let text = std::fs::read_to_string(path).expect("dump readable");
    let mut c = Case {
        path: path.to_path_buf(),
        ..Default::default()
    };
    // r_id -> intervals, kept in first-seen order so the replay iterates the
    // reads in the same order the reference did.
    let mut order: Vec<u32> = Vec::new();
    let mut by_read: FxHashMap<u32, Vec<Interval>> = FxHashMap::default();

    for line in text.lines() {
        if line.starts_with("# params") {
            for tok in line.split_whitespace() {
                if let Some(v) = tok.strip_prefix("k=") {
                    c.k = v.parse().unwrap();
                } else if let Some(v) = tok.strip_prefix("delta_len=") {
                    c.delta_len = v.parse().unwrap();
                }
            }
            continue;
        }
        let mut f = line.splitn(2, ' ');
        let tag = f.next().unwrap_or("");
        let rest = f.next().unwrap_or("");
        match tag {
            "R" => {
                let mut p = rest.splitn(2, ' ');
                let r: u32 = p.next().unwrap().parse().unwrap();
                c.reads
                    .insert(r, p.next().unwrap_or("").as_bytes().to_vec());
            }
            "L" => {
                let mut p = rest.split_whitespace();
                let r: u32 = p.next().unwrap().parse().unwrap();
                let l: u32 = p.next().unwrap().parse().unwrap();
                c.read_len.push((r, l));
            }
            "I" => {
                let p: Vec<&str> = rest.split_whitespace().collect();
                let r: u32 = p[0].parse().unwrap();
                // p[1] is the position within the read, implied by order here.
                let start: u32 = p[2].parse().unwrap();
                let end: u32 = p[3].parse().unwrap();
                let weight: u32 = p[4].parse().unwrap();
                let occ: Vec<u32> = if p.len() > 5 && !p[5].is_empty() {
                    p[5].split(',').map(|x| x.parse().unwrap()).collect()
                } else {
                    Vec::new()
                };
                if !by_read.contains_key(&r) {
                    order.push(r);
                }
                by_read.entry(r).or_default().push(Interval {
                    start,
                    end,
                    weight,
                    occurrences: occ,
                });
            }
            "N" => {
                // `N <name...> <end_mini_seq> <reads>` --- the name contains
                // spaces ("40, 97, 3"), so split from the right.
                let parts: Vec<&str> = rest.rsplitn(3, ' ').collect();
                let reads_field = parts[0];
                let ems = parts[1];
                let name = parts[2];
                c.exp_nodes.insert(
                    name.to_string(),
                    (
                        if ems == "-" {
                            String::new()
                        } else {
                            ems.to_string()
                        },
                        parse_reads_attr(reads_field),
                    ),
                );
            }
            "E" => {
                // `E <u...> -> <v...> <length> <supp>`
                let arrow = rest.find(" -> ").expect("edge line has an arrow");
                let u = rest[..arrow].to_string();
                let tail = &rest[arrow + 4..];
                let parts: Vec<&str> = tail.rsplitn(3, ' ').collect();
                let supp_field = parts[0];
                let length: i64 = parts[1].parse().unwrap();
                let v = parts[2].to_string();
                let supp: Vec<u32> = if supp_field.is_empty() {
                    Vec::new()
                } else {
                    supp_field.split(',').map(|x| x.parse().unwrap()).collect()
                };
                c.exp_edges.insert((u, v), (length, supp));
            }
            "T" => {
                c.exp_topo = Some(if rest.starts_with("CYCLIC") {
                    None
                } else {
                    Some(rest.split('|').map(|s| s.to_string()).collect())
                });
            }
            "F" => {
                c.exp_reads_for_isoforms = if rest.trim().is_empty() {
                    Vec::new()
                } else {
                    rest.trim().split(',').map(|x| x.parse().unwrap()).collect()
                };
            }
            "S" => {
                let mut p = rest.split_whitespace();
                c.exp_counts = (
                    p.next().unwrap().parse().unwrap(),
                    p.next().unwrap().parse().unwrap(),
                );
            }
            _ => {}
        }
    }

    c.intervals = order
        .into_iter()
        .map(|r| (r, by_read.remove(&r).unwrap()))
        .collect();
    c
}

/// Compare one case, returning a human-readable report of every disagreement.
fn check(c: &Case) -> Vec<String> {
    let mut problems = Vec::new();

    let input = BuildInput {
        k: c.k,
        delta_len: c.delta_len,
        intervals: &c.intervals,
        reads: &c.reads,
        read_len: &c.read_len,
    };
    let (g, rfi) = match generate_graph_from_intervals(&input, BuildOpts::default()) {
        Ok(v) => v,
        Err(e) => {
            problems.push(format!("port returned an error: {e:?}"));
            return problems;
        }
    };

    // -- nodes -------------------------------------------------------------
    let mut got_nodes: BTreeMap<String, (String, Vec<(u32, ReadInfo)>)> = BTreeMap::new();
    for n in g.nodes_sorted_by_name() {
        got_nodes.insert(
            g.key(n).to_string(),
            (
                String::from_utf8_lossy(g.end_mini_seq(n)).to_string(),
                g.reads(n).to_vec(),
            ),
        );
    }

    for name in c.exp_nodes.keys() {
        if !got_nodes.contains_key(name) {
            problems.push(format!("missing node {name:?}"));
        }
    }
    for name in got_nodes.keys() {
        if !c.exp_nodes.contains_key(name) {
            problems.push(format!("extra node {name:?}"));
        }
    }
    for (name, (exp_ems, exp_reads)) in &c.exp_nodes {
        if let Some((got_ems, got_reads)) = got_nodes.get(name) {
            if exp_ems != got_ems {
                problems.push(format!(
                    "node {name:?} end_mini_seq: reference {exp_ems:?}, port {got_ems:?}"
                ));
            }
            if exp_reads != got_reads {
                // Order is part of the comparison, so say whether it is only
                // the order that differs --- that distinction points at
                // different code.
                let mut a = exp_reads.clone();
                let mut b = got_reads.clone();
                a.sort_by_key(|(r, _)| *r);
                b.sort_by_key(|(r, _)| *r);
                let kind = if a == b { "ORDER ONLY" } else { "contents" };
                problems.push(format!(
                    "node {name:?} reads differ ({kind}): reference {} entries, port {}",
                    exp_reads.len(),
                    got_reads.len()
                ));
            }
        }
    }

    // -- edges -------------------------------------------------------------
    let mut got_edges: BTreeMap<(String, String), (i64, Vec<u32>)> = BTreeMap::new();
    for (u, v) in g.edges_sorted_by_name() {
        got_edges.insert(
            (g.key(u).to_string(), g.key(v).to_string()),
            (
                g.edge_length(u, v).unwrap(),
                g.edge_support(u, v).unwrap_or(&[]).to_vec(),
            ),
        );
    }
    for e in c.exp_edges.keys() {
        if !got_edges.contains_key(e) {
            problems.push(format!("missing edge {:?} -> {:?}", e.0, e.1));
        }
    }
    for e in got_edges.keys() {
        if !c.exp_edges.contains_key(e) {
            problems.push(format!("extra edge {:?} -> {:?}", e.0, e.1));
        }
    }
    for (e, (exp_len, exp_supp)) in &c.exp_edges {
        if let Some((got_len, got_supp)) = got_edges.get(e) {
            if exp_len != got_len {
                problems.push(format!(
                    "edge {:?}->{:?} length: reference {exp_len}, port {got_len}",
                    e.0, e.1
                ));
            }
            if exp_supp != got_supp {
                let mut a = exp_supp.clone();
                let mut b = got_supp.clone();
                a.sort_unstable();
                b.sort_unstable();
                let kind = if a == b { "ORDER ONLY" } else { "contents" };
                problems.push(format!(
                    "edge {:?}->{:?} edge_supp differs ({kind}): reference {exp_supp:?}, port {got_supp:?}",
                    e.0, e.1
                ));
            }
        }
    }

    // -- the topological order --------------------------------------------
    //
    // Compared exactly, not just checked for validity. A topological order is
    // not unique and the reference compares the resulting indices to pick
    // candidate bubbles, so producing *a* valid order is not enough.
    let got_topo = g.topological_sort();
    match (&c.exp_topo, &got_topo) {
        // No `T` record: the dump predates it, so the order is simply unchecked.
        // Not silently ignored --- the caller counts these and reports how many.
        (None, _) => {}
        (Some(None), None) => {
            // Both agree the graph is cyclic. Worth knowing, but not a mismatch.
        }
        (Some(None), Some(_)) => {
            problems.push("reference reported a cycle, port found a topological order".to_string())
        }
        (Some(Some(_)), None) => {
            problems.push("port reported a cycle, reference found a topological order".to_string())
        }
        (Some(Some(exp)), Some(got)) => {
            let got_names: Vec<String> = got.iter().map(|&n| g.key(n).to_string()).collect();
            if *exp != got_names {
                let first = exp
                    .iter()
                    .zip(got_names.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(exp.len().min(got_names.len()));
                problems.push(format!(
                    "topological order diverges at index {first}: reference {:?}, port {:?} \
                     (lengths {} vs {})",
                    exp.get(first),
                    got_names.get(first),
                    exp.len(),
                    got_names.len()
                ));
            }
        }
    }

    // -- the rest ----------------------------------------------------------
    if rfi != c.exp_reads_for_isoforms {
        problems.push(format!(
            "reads_for_isoforms: reference {} entries, port {}",
            c.exp_reads_for_isoforms.len(),
            rfi.len()
        ));
    }
    if (g.node_count(), g.edge_count()) != c.exp_counts {
        problems.push(format!(
            "counts: reference {:?}, port {:?}",
            c.exp_counts,
            (g.node_count(), g.edge_count())
        ));
    }

    problems
}

#[test]
fn graph_construction_matches_the_reference() {
    let Some(dir) = std::env::var_os("ISONFORM_GRAPH_DUMPS") else {
        // Loud skip. A silently-green oracle is worse than a red one.
        eprintln!(
            "SKIPPED: set ISONFORM_GRAPH_DUMPS to a directory produced by\n\
             bench/dump_reference.py to run the graph oracle. Without it this\n\
             test asserts nothing."
        );
        return;
    };
    let dir = PathBuf::from(dir);

    let mut cases: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display()))
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|s| s.to_str())
                .is_some_and(|s| s.starts_with("graph_") && s.ends_with(".txt"))
        })
        .collect();
    cases.sort();

    assert!(
        !cases.is_empty(),
        "no graph_*.txt under {} --- did bench/dump_reference.py run?",
        dir.display()
    );

    let mut report = String::new();
    let mut failed = 0usize;
    let mut n_nodes = 0usize;
    let mut n_edges = 0usize;
    let mut n_topo_unchecked = 0usize;

    for path in &cases {
        let c = parse_case(path);
        n_nodes += c.exp_counts.0;
        n_edges += c.exp_counts.1;
        if c.exp_topo.is_none() {
            n_topo_unchecked += 1;
        }
        let problems = check(&c);
        if !problems.is_empty() {
            failed += 1;
            let _ = writeln!(report, "\n=== {} ===", c.path.display());
            // Cap the per-case output: one wrong branch produces thousands of
            // identical-looking lines and the first few are what localise it.
            for p in problems.iter().take(20) {
                let _ = writeln!(report, "  {p}");
            }
            if problems.len() > 20 {
                let _ = writeln!(report, "  ... and {} more", problems.len() - 20);
            }
        }
    }

    eprintln!(
        "graph oracle: {} case(s), {} reference nodes and {} edges in total",
        cases.len(),
        n_nodes,
        n_edges
    );
    if n_topo_unchecked > 0 {
        eprintln!(
            "  NOTE: {n_topo_unchecked} case(s) carry no T record, so their \
             topological order was NOT compared. Re-run bench/dump_reference.py \
             to record it."
        );
    }
    assert!(
        failed == 0,
        "{failed} of {} case(s) disagree with the reference:{report}",
        cases.len()
    );
}
