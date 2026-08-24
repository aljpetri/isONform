//! Differential oracle for graph simplification.
//!
//! Rebuilds the graph the reference had on entry to `simplifyGraph`, runs the
//! port's [`isonform::simplify::simplify_graph`] over it, and diffs the result
//! against the graph the reference produced.
//!
//! ```text
//! bench/dump_reference.py --fastq-folder bench/corpus/sirv_small \
//!     --outdir /tmp/d --stage simplify
//! ISONFORM_SIMPLIFY_DUMPS=/tmp/d cargo test --test simplify_oracle -- --nocapture
//! ```
//!
//! Skips loudly with the variable unset, like the graph oracle.
//!
//! # What a disagreement here means
//!
//! Graph construction is pure; this stage calls **spoa**, whose consensus decides
//! which bubbles are poppable. So this oracle used to split its verdict — cases
//! that never touched spoa failed the build, cases that did were only reported —
//! because `crate::poa`'s equivalence to the `spoa` binary had been measured on
//! isONcorrect's correction intervals and nothing said it transferred to
//! isONform's different sequences from a different stage.
//!
//! **It transfers, and that is measured**: 3 277 of 3 277 recorded isONform spoa
//! calls give an identical consensus to the binary, 3 074 of them from
//! `SimplifyGraph.py:570`, the site these verdicts depend on. See `crate::poa`
//! and `bench/dump_reference.py --record-spoa`.
//!
//! So every disagreement now fails the build. One boundary is worth keeping in
//! view: that result licenses "same inputs ⇒ same consensus", not "the port
//! computes the same inputs". If span extraction diverges, spoa is handed
//! different sequences and faithfully returns a different consensus — still a bug
//! on this side, which is why the gate is unconditional either way.
//!
//! `spoa_calls` is still reported per case, now as a diagnostic rather than a
//! verdict: it says whether a failure went through the consensus path at all,
//! which is the first thing worth knowing when localising one.
//!
//! # Rebuilding the "before" graph
//!
//! From the `BN`/`BE` records, replayed **in order**. That is load-bearing:
//! networkx emits them in insertion order and `nx.topological_sort`'s output
//! depends on node insertion order and per-node adjacency order, which in turn
//! decides which node pairs are candidate bubbles. Sorted records would rebuild a
//! graph with a different (still valid) topological order and every case would
//! disagree for a reason that has nothing to do with the algorithm.
//!
//! The rebuild is checked rather than trusted: each case asserts the port's
//! topological order equals the recorded `BT` **before** running the stage. A
//! failure there is a reconstruction problem, reported separately so it is not
//! mistaken for an algorithmic one.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isonform::graph::{Graph, NodeKey, ReadInfo};
use isonform::simplify::{pop_bubbles, Counting, PopOpts, RealAligner, SpoaParasail};
use rustc_hash::FxHashMap;

/// A node's `reads` map and `end_mini_seq`, as recorded.
type NodeRec = (String, Vec<(u32, ReadInfo)>);
/// An edge's `length` (absent means the attribute is missing) and support.
type EdgeRec = (Option<i64>, Vec<u32>);

#[derive(Default)]
struct Case {
    path: PathBuf,
    k: usize,
    delta_len: i64,
    slow: bool,
    reads: FxHashMap<u32, Vec<u8>>,
    /// In record order, which is insertion order.
    before_nodes: Vec<(String, NodeRec)>,
    before_edges: Vec<((String, String), EdgeRec)>,
    before_topo: Option<Vec<String>>,
    after_nodes: BTreeMap<String, NodeRec>,
    after_edges: BTreeMap<(String, String), EdgeRec>,
}

fn parse_reads_attr(field: &str) -> Vec<(u32, ReadInfo)> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .filter(|s| !s.is_empty())
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

/// `<name...> <end_mini_seq> <reads>` — the name contains spaces, so split from
/// the right.
fn parse_node(rest: &str) -> (String, NodeRec) {
    let parts: Vec<&str> = rest.rsplitn(3, ' ').collect();
    let reads = parse_reads_attr(parts[0]);
    let ems = if parts[1] == "-" {
        String::new()
    } else {
        parts[1].to_string()
    };
    (parts[2].to_string(), (ems, reads))
}

/// `<u...> -> <v...> <length|NA> <supp>`
fn parse_edge(rest: &str) -> ((String, String), EdgeRec) {
    let arrow = rest.find(" -> ").expect("edge record has an arrow");
    let u = rest[..arrow].to_string();
    let tail = &rest[arrow + 4..];
    let parts: Vec<&str> = tail.rsplitn(3, ' ').collect();
    let supp: Vec<u32> = parts[0]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|x| x.parse().unwrap())
        .collect();
    // `NA` is not a placeholder for "unknown": it is how the dump records an edge
    // that genuinely has no `length` attribute, which is what
    // `prepare_adding_edges` produces for every edge it re-adds.
    let length = if parts[1] == "NA" {
        None
    } else {
        Some(parts[1].parse().unwrap())
    };
    ((u, parts[2].to_string()), (length, supp))
}

fn parse_case(path: &Path) -> Case {
    let text = std::fs::read_to_string(path).expect("dump readable");
    let mut c = Case {
        path: path.to_path_buf(),
        ..Default::default()
    };
    for line in text.lines() {
        if line.starts_with("# params") {
            for tok in line.split_whitespace() {
                if let Some(v) = tok.strip_prefix("k=") {
                    c.k = v.parse().unwrap();
                } else if let Some(v) = tok.strip_prefix("delta_len=") {
                    c.delta_len = v.parse().unwrap();
                } else if let Some(v) = tok.strip_prefix("mode=") {
                    // argparse's `--slow` reaches this as a Python bool repr.
                    c.slow = v == "True";
                }
            }
            continue;
        }
        let mut f = line.splitn(2, ' ');
        let (tag, rest) = (f.next().unwrap_or(""), f.next().unwrap_or(""));
        match tag {
            "R" => {
                let mut p = rest.splitn(2, ' ');
                let r: u32 = p.next().unwrap().parse().unwrap();
                c.reads
                    .insert(r, p.next().unwrap_or("").as_bytes().to_vec());
            }
            "BN" => c.before_nodes.push(parse_node(rest)),
            "BE" => c.before_edges.push(parse_edge(rest)),
            "BT" => {
                c.before_topo = if rest.starts_with("CYCLIC") {
                    None
                } else {
                    Some(rest.split('|').map(|s| s.to_string()).collect())
                }
            }
            "AN" => {
                let (name, rec) = parse_node(rest);
                c.after_nodes.insert(name, rec);
            }
            "AE" => {
                let (key, rec) = parse_edge(rest);
                c.after_edges.insert(key, rec);
            }
            _ => {}
        }
    }
    c
}

/// The reference's node name back into a [`NodeKey`].
fn key_from_name(name: &str) -> NodeKey {
    match name {
        "s" => NodeKey::Source,
        "t" => NodeKey::Sink,
        _ => {
            let f: Vec<&str> = name.split(", ").collect();
            assert_eq!(f.len(), 3, "unexpected node name {name:?}");
            NodeKey::Interval {
                start: f[0].parse().unwrap(),
                end: f[1].parse().unwrap(),
                r_id: f[2].parse().unwrap(),
            }
        }
    }
}

fn rebuild(c: &Case) -> Graph {
    let mut g = Graph::new();
    // Nodes first, in record order, so insertion order matches the reference's.
    for (name, (ems, reads)) in &c.before_nodes {
        let n = g.add_node(key_from_name(name));
        g.set_end_mini_seq(n, ems.as_bytes());
        for &(r, info) in reads {
            g.set_read(n, r, info);
        }
    }
    // Then edges, in record order, so each node's adjacency order matches.
    for ((u, v), (length, supp)) in &c.before_edges {
        let (u, v) = (
            g.lookup(&key_from_name(u)).expect("edge endpoint exists"),
            g.lookup(&key_from_name(v)).expect("edge endpoint exists"),
        );
        match length {
            Some(l) => {
                g.add_edge(u, v, *l);
            }
            None => {
                g.upsert_edge_support(u, v, Vec::new());
            }
        }
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
    g
}

struct Outcome {
    problems: Vec<String>,
    spoa_calls: usize,
    /// Set when the rebuild itself is wrong, which is not an algorithm failure.
    rebuild_failed: bool,
    /// Set when the iteration cap stopped the loop.
    ///
    /// This is attributable **regardless of spoa**: the reference terminates on
    /// these graphs and the port does not, and no consensus difference can turn a
    /// terminating loop into a non-terminating one --- it can only change which
    /// bubbles get popped. So a cap hit is a port bug and is treated as one, even
    /// on a case that called spoa thousands of times.
    hit_cap: bool,
}

fn check(c: &Case) -> Outcome {
    let mut problems = Vec::new();
    let mut g = rebuild(c);

    // Precondition: the rebuild has to reproduce the reference's topological
    // order, or nothing downstream means anything.
    if let Some(exp) = &c.before_topo {
        match g.topological_sort() {
            Some(got) => {
                let got: Vec<String> = got.iter().map(|&n| g.key(n).to_string()).collect();
                if *exp != got {
                    let at = exp
                        .iter()
                        .zip(got.iter())
                        .position(|(a, b)| a != b)
                        .unwrap_or(0);
                    return Outcome {
                        problems: vec![format!(
                            "REBUILD: topological order diverges at index {at}: \
                             reference {:?}, rebuilt {:?}",
                            exp.get(at),
                            got.get(at)
                        )],
                        spoa_calls: 0,
                        rebuild_failed: true,
                        hit_cap: false,
                    };
                }
            }
            None => {
                return Outcome {
                    problems: vec!["REBUILD: rebuilt graph is cyclic".to_string()],
                    spoa_calls: 0,
                    rebuild_failed: true,
                    hit_cap: false,
                }
            }
        }
    }

    let mut aligner = RealAligner::new(Counting::new(SpoaParasail), &c.reads, c.k, c.delta_len);
    let stats = pop_bubbles(&mut g, &mut aligner, PopOpts { slow: c.slow });
    let spoa_calls = aligner.engine.spoa_calls;
    let hit_cap = stats.hit_iteration_cap;
    if hit_cap {
        problems.push(
            "PORT BUG: the iteration cap stopped the loop. The reference \
             terminates on this graph and the port does not, which no consensus \
             difference can explain."
                .to_string(),
        );
    }

    // -- nodes -------------------------------------------------------------
    let mut got_nodes: BTreeMap<String, NodeRec> = BTreeMap::new();
    for n in g.nodes_sorted_by_name() {
        got_nodes.insert(
            g.key(n).to_string(),
            (
                String::from_utf8_lossy(g.end_mini_seq(n)).to_string(),
                g.reads(n).to_vec(),
            ),
        );
    }
    for name in c.after_nodes.keys() {
        if !got_nodes.contains_key(name) {
            problems.push(format!("missing node {name:?}"));
        }
    }
    for name in got_nodes.keys() {
        if !c.after_nodes.contains_key(name) {
            problems.push(format!("extra node {name:?}"));
        }
    }
    for (name, (exp_ems, exp_reads)) in &c.after_nodes {
        if let Some((got_ems, got_reads)) = got_nodes.get(name) {
            if exp_ems != got_ems {
                problems.push(format!("node {name:?} end_mini_seq differs"));
            }
            if exp_reads != got_reads {
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
    let mut got_edges: BTreeMap<(String, String), EdgeRec> = BTreeMap::new();
    for (u, v) in g.edges_sorted_by_name() {
        got_edges.insert(
            (g.key(u).to_string(), g.key(v).to_string()),
            (
                g.edge_length(u, v),
                g.edge_support(u, v).unwrap_or(&[]).to_vec(),
            ),
        );
    }
    for e in c.after_edges.keys() {
        if !got_edges.contains_key(e) {
            problems.push(format!("missing edge {:?} -> {:?}", e.0, e.1));
        }
    }
    for e in got_edges.keys() {
        if !c.after_edges.contains_key(e) {
            problems.push(format!("extra edge {:?} -> {:?}", e.0, e.1));
        }
    }
    for (e, (exp_len, exp_supp)) in &c.after_edges {
        if let Some((got_len, got_supp)) = got_edges.get(e) {
            if exp_len != got_len {
                problems.push(format!(
                    "edge {:?}->{:?} length: reference {exp_len:?}, port {got_len:?}",
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
                    "edge {:?}->{:?} support differs ({kind}): {} vs {} entries",
                    e.0,
                    e.1,
                    exp_supp.len(),
                    got_supp.len()
                ));
            }
        }
    }

    Outcome {
        problems,
        spoa_calls,
        rebuild_failed: false,
        hit_cap,
    }
}

#[test]
fn simplification_matches_the_reference() {
    let Some(dir) = std::env::var_os("ISONFORM_SIMPLIFY_DUMPS") else {
        eprintln!(
            "SKIPPED: set ISONFORM_SIMPLIFY_DUMPS to a directory produced by\n\
             bench/dump_reference.py --stage simplify to run the simplification\n\
             oracle. Without it this test asserts nothing."
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
                .is_some_and(|s| s.starts_with("simplify_") && s.ends_with(".txt"))
        })
        .collect();
    cases.sort();
    assert!(
        !cases.is_empty(),
        "no simplify_*.txt under {} --- did bench/dump_reference.py run with \
         --stage simplify?",
        dir.display()
    );

    let mut report = String::new();
    let (mut ok_no_spoa, mut ok_spoa) = (0usize, 0usize);
    let (mut bad_no_spoa, mut bad_spoa) = (0usize, 0usize);
    let mut rebuild_failures = 0usize;

    for path in &cases {
        let c = parse_case(path);
        let out = check(&c);
        let name = c
            .path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        if out.rebuild_failed {
            rebuild_failures += 1;
        } else if out.problems.is_empty() {
            if out.spoa_calls == 0 {
                ok_no_spoa += 1;
            } else {
                ok_spoa += 1;
            }
            continue;
        } else if out.spoa_calls == 0 || out.hit_cap {
            // A cap hit is attributable whatever spoa did --- see `Outcome`.
            bad_no_spoa += 1;
        } else {
            bad_spoa += 1;
        }

        let _ = writeln!(
            report,
            "\n=== {name} === ({} spoa call(s){})",
            out.spoa_calls,
            if out.hit_cap { ", and the cap hit" } else { "" }
        );
        for p in out.problems.iter().take(12) {
            let _ = writeln!(report, "  {p}");
        }
        if out.problems.len() > 12 {
            let _ = writeln!(report, "  ... and {} more", out.problems.len() - 12);
        }
    }

    eprintln!(
        "simplify oracle: {} case(s)\n  \
         agree: {ok_no_spoa} without spoa, {ok_spoa} with\n  \
         disagree: {bad_no_spoa} without spoa, {bad_spoa} with\n  \
         rebuild failures: {rebuild_failures}",
        cases.len()
    );

    // The gate used to be split: cases that called spoa were reported but did not
    // fail, because `crate::poa`'s equivalence had only ever been measured on
    // isONcorrect's correction intervals and could not be assumed to carry over.
    //
    // It does carry over, and that is now measured rather than assumed: 3 277 of
    // 3 277 recorded isONform spoa calls produce an identical consensus to the
    // binary, 3 074 of them from `SimplifyGraph.py:570`, the site these verdicts
    // depend on (`crate::poa`, and `--record-spoa` in `bench/dump_reference.py`).
    //
    // So the escape hatch is gone and every disagreement fails. Note what the
    // spoa result does and does not license: it says that *given the same input
    // sequences in the same order* the consensus agrees. It does not say the port
    // computes the same inputs — if span extraction diverges, spoa is handed
    // different sequences and correctly returns a different consensus. Either way
    // the bug is on this side, which is exactly why these cases now fail.
    assert_eq!(
        rebuild_failures, 0,
        "the before-graph rebuild is wrong, so nothing else here means anything:{report}"
    );
    assert_eq!(
        bad_no_spoa + bad_spoa,
        0,
        "{} case(s) disagree ({bad_no_spoa} without spoa, {bad_spoa} with). \
         spoa equivalence is measured, so all of these are attributable:{report}",
        bad_no_spoa + bad_spoa
    );
}
