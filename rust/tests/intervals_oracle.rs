//! Differential oracle for the front half of `main`: reads in, intervals out.
//!
//! Replays a batch's reads through [`isonform::intervals::build_batch`] and diffs
//! the result against what the reference chose --- the `all_intervals_for_graph`
//! that `generateGraphfromIntervals` receives, which is the boundary between the
//! shared isONcorrect machinery and isONform's own algorithm.
//!
//! ```text
//! bench/dump_reference.py --fastq-folder bench/corpus/sirv_small \
//!     --outdir /tmp/d --stage intervals
//! ISONFORM_INTERVAL_DUMPS=/tmp/d cargo test --test intervals_oracle -- --nocapture
//! ```
//!
//! Skips loudly with the variable unset, like the other oracles.
//!
//! # What this covers that the graph oracle does not
//!
//! The graph oracle starts *from* these intervals: it replays them and checks the
//! graph. So everything upstream --- minimizer selection, the anchor database,
//! span support, weighted interval scheduling --- was unverified until now, and a
//! wrong minimizer would have produced a graph the graph oracle happily agreed
//! with, because it would have been fed the reference's intervals rather than the
//! port's.
//!
//! It is also the first stage whose input is a *read* rather than a recorded
//! intermediate, which is what makes the port runnable end to end.
//!
//! # Three things are compared, not one
//!
//! * the **intervals** per `graph_id`: start, stop, support, and the full
//!   `(r_id, p1, p2)` instance list *in order*, because the order reaches the
//!   graph;
//! * the **`graph_id` assignment**: which read got which id, since the reference
//!   numbers them 1-up in read order and the graph is keyed by it;
//! * which reads were **skipped**, by implication --- a read absent from the
//!   assignment produced no intervals.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isonform::intervals::{build_batch, ChosenInterval};
use isonform::wis::WisOpts;

#[derive(Default)]
struct Case {
    path: PathBuf,
    k: usize,
    w: usize,
    x_low: usize,
    x_high: usize,
    delta_len: i64,
    /// Every read the batch saw, in `r_id` order.
    reads: Vec<(u32, Vec<u8>)>,
    /// `graph_id -> r_id`, as the reference assigned them.
    graph_ids: BTreeMap<u32, u32>,
    /// `graph_id -> intervals`, in record order.
    intervals: BTreeMap<u32, Vec<ChosenInterval>>,
}

fn parse_case(path: &Path) -> Case {
    let text = std::fs::read_to_string(path).expect("dump readable");
    let mut c = Case {
        path: path.to_path_buf(),
        ..Default::default()
    };
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("# params") {
            for tok in rest.split_whitespace() {
                let (key, val) = tok.split_once('=').expect("params are key=value");
                match key {
                    "k" => c.k = val.parse().unwrap(),
                    "w" => c.w = val.parse().unwrap(),
                    "xmin" => c.x_low = val.parse().unwrap(),
                    "xmax" => c.x_high = val.parse().unwrap(),
                    "delta_len" => c.delta_len = val.parse().unwrap(),
                    _ => {}
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
                    .push((r, p.next().unwrap_or("").as_bytes().to_vec()));
            }
            "G" => {
                let mut p = rest.split_whitespace();
                let gid: u32 = p.next().unwrap().parse().unwrap();
                let rid: i64 = p.next().unwrap().parse().unwrap();
                assert!(
                    rid >= 0,
                    "{}: unmatched read for graph id {gid}",
                    path.display()
                );
                c.graph_ids.insert(gid, rid as u32);
            }
            "I" => {
                let f: Vec<&str> = rest.split_whitespace().collect();
                let gid: u32 = f[0].parse().unwrap();
                // f[1] is the position within the read's list; records are in
                // order, so pushing preserves it.
                let start: usize = f[2].parse().unwrap();
                let stop: usize = f[3].parse().unwrap();
                let support: usize = f[4].parse().unwrap();
                let flat: Vec<u32> = f
                    .get(5)
                    .map(|s| {
                        s.split(',')
                            .filter(|x| !x.is_empty())
                            .map(|x| x.parse().unwrap())
                            .collect()
                    })
                    .unwrap_or_default();
                let instance = flat.chunks(3).map(|t| (t[0], t[1], t[2])).collect();
                c.intervals.entry(gid).or_default().push(ChosenInterval {
                    start,
                    stop,
                    support,
                    instance,
                });
            }
            _ => {}
        }
    }
    c.reads.sort_by_key(|(r, _)| *r);
    c
}

fn check(c: &Case) -> Vec<String> {
    let mut problems = Vec::new();
    let got = build_batch(
        &c.reads,
        c.w,
        c.k,
        c.x_low,
        c.x_high,
        c.delta_len,
        // The reference is off by one in `fill_p2`; the goldens contain that.
        WisOpts::default(),
    );

    // -- graph_id assignment ----------------------------------------------
    let got_ids: BTreeMap<u32, u32> = got.read_of_graph_id.iter().map(|(g, r)| (*g, *r)).collect();
    if got_ids != c.graph_ids {
        let only_ref: Vec<_> = c
            .graph_ids
            .iter()
            .filter(|(g, r)| got_ids.get(g) != Some(r))
            .take(4)
            .collect();
        let only_port: Vec<_> = got_ids
            .iter()
            .filter(|(g, r)| c.graph_ids.get(g) != Some(r))
            .take(4)
            .collect();
        problems.push(format!(
            "graph_id assignment differs: reference {} entries, port {} \
             | first ref-only {only_ref:?} | first port-only {only_port:?}",
            c.graph_ids.len(),
            got_ids.len()
        ));
    }

    // -- intervals ---------------------------------------------------------
    for (gid, want) in &c.intervals {
        let Some(have) = got.by_graph_id.get(gid) else {
            problems.push(format!("graph_id {gid}: no intervals from the port"));
            continue;
        };
        if have.len() != want.len() {
            problems.push(format!(
                "graph_id {gid}: reference chose {} interval(s), port {}",
                want.len(),
                have.len()
            ));
            continue;
        }
        for (i, (a, b)) in want.iter().zip(have.iter()).enumerate() {
            if a.start != b.start || a.stop != b.stop {
                problems.push(format!(
                    "graph_id {gid} interval {i}: span {}..{} vs {}..{}",
                    a.start, a.stop, b.start, b.stop
                ));
            } else if a.support != b.support {
                problems.push(format!(
                    "graph_id {gid} interval {i} ({}..{}): support {} vs {}",
                    a.start, a.stop, a.support, b.support
                ));
            } else if a.instance != b.instance {
                let kind = {
                    let (mut x, mut y) = (a.instance.clone(), b.instance.clone());
                    x.sort_unstable();
                    y.sort_unstable();
                    if x == y {
                        "ORDER ONLY"
                    } else {
                        "contents"
                    }
                };
                problems.push(format!(
                    "graph_id {gid} interval {i} ({}..{}): instance differs ({kind}), \
                     {} vs {} entries",
                    a.start,
                    a.stop,
                    a.instance.len(),
                    b.instance.len()
                ));
            }
        }
    }
    for gid in got.by_graph_id.keys() {
        if !c.intervals.contains_key(gid) {
            problems.push(format!(
                "graph_id {gid}: port produced intervals, reference none"
            ));
        }
    }
    problems
}

#[test]
fn interval_selection_matches_the_reference() {
    let Some(dir) = std::env::var_os("ISONFORM_INTERVAL_DUMPS") else {
        eprintln!(
            "SKIPPED: set ISONFORM_INTERVAL_DUMPS to a directory produced by\n\
             bench/dump_reference.py --stage intervals to run the interval\n\
             oracle. Without it this test asserts nothing."
        );
        return;
    };
    let dir = PathBuf::from(dir);
    let mut cases: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display()))
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("intervals_") && n.ends_with(".txt"))
        })
        .collect();
    cases.sort();
    assert!(!cases.is_empty(), "no intervals_*.txt in {}", dir.display());

    let (mut ok, mut bad) = (0usize, 0usize);
    let (mut reads, mut kept) = (0usize, 0usize);
    let mut report = String::new();
    for path in &cases {
        let c = parse_case(path);
        reads += c.reads.len();
        kept += c.graph_ids.len();
        let problems = check(&c);
        if problems.is_empty() {
            ok += 1;
            continue;
        }
        bad += 1;
        let name = c.path.file_name().unwrap().to_string_lossy();
        let _ = writeln!(
            report,
            "\n=== {name} === ({} reads in, {} kept by the reference)",
            c.reads.len(),
            c.graph_ids.len()
        );
        for p in problems.iter().take(10) {
            let _ = writeln!(report, "  {p}");
        }
        if problems.len() > 10 {
            let _ = writeln!(report, "  ... and {} more", problems.len() - 10);
        }
    }
    eprintln!(
        "interval oracle: {} case(s), {ok} agree, {bad} disagree \
         ({reads} reads in, {kept} kept)",
        cases.len()
    );
    assert_eq!(bad, 0, "{bad} case(s) disagree:{report}");
}
