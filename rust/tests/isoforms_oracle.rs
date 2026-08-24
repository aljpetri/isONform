//! Differential oracle for isoform generation.
//!
//! Rebuilds the simplified graph the reference had on entry to
//! `generate_isoforms`, replays both halves of the stage, and diffs each against
//! what the reference produced:
//!
//! * **grouping** — `compute_equal_reads`, including the sub-isoform merge, against
//!   the recorded `Q` records;
//! * **merging** — `merge_consensuses`, which runs spoa and parasail, against the
//!   recorded `C` records.
//!
//! ```text
//! bench/dump_reference.py --fastq-folder bench/corpus/sirv_small \
//!     --outdir /tmp/d --stage isoforms
//! ISONFORM_ISOFORM_DUMPS=/tmp/d cargo test --test isoforms_oracle -- --nocapture
//! ```
//!
//! Skips loudly with the variable unset, like the others.
//!
//! # The two halves are reported separately, because they fail for different reasons
//!
//! Grouping is a pure graph walk: a disagreement there is unambiguously this
//! port's. Merging runs spoa and parasail, both of which are now verified exactly
//! on isONform's own inputs (`PORTING.md` findings 24 and 25), so a disagreement
//! there is also this port's — but it is worth knowing which half moved, because
//! a wrong grouping feeds wrong sequences into the merge and would show up as
//! both.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isonform::graph::{Graph, NodeKey, ReadInfo};
use isonform::isoforms::{compute_equal_reads, merge_consensuses, MergeOpts, SpoaParasailMerge};
use rustc_hash::FxHashMap;

#[derive(Default)]
struct Case {
    path: PathBuf,
    opts: MergeOpts,
    reads: FxHashMap<u32, Vec<u8>>,
    support: Vec<u32>,
    nodes: Vec<String>,
    edges: Vec<(String, String, Vec<u32>)>,
    /// `equal_reads` after grouping and sub-isoform merging.
    groups: BTreeMap<u32, Vec<u32>>,
    /// `new_consensuses`.
    consensuses: BTreeMap<u32, Vec<u8>>,
}

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

fn parse_case(path: &Path) -> Case {
    let text = std::fs::read_to_string(path).expect("dump readable");
    let mut c = Case {
        path: path.to_path_buf(),
        ..Default::default()
    };
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("# params") {
            for tok in rest.split_whitespace() {
                let (k, v) = tok.split_once('=').expect("key=value");
                match k {
                    "delta" => c.opts.delta = v.parse().unwrap(),
                    "delta_len" => c.opts.delta_len = v.parse().unwrap(),
                    "d3" => c.opts.delta_iso_len_3 = v.parse().unwrap(),
                    "d5" => c.opts.delta_iso_len_5 = v.parse().unwrap(),
                    "max_seqs" => c.opts.max_seqs_to_spoa = v.parse().unwrap(),
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
                    .insert(r, p.next().unwrap_or("").as_bytes().to_vec());
            }
            "S" => {
                c.support = rest
                    .split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect()
            }
            "N" => c.nodes.push(rest.to_string()),
            "E" => {
                let (u, tail) = rest.split_once(" -> ").expect("edge arrow");
                let (v, supp) = tail.rsplit_once(' ').unwrap_or((tail, ""));
                let supp = supp
                    .split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect();
                c.edges.push((u.to_string(), v.to_string(), supp));
            }
            "Q" => {
                let (k, v) = rest.split_once(' ').unwrap_or((rest, ""));
                let members = v
                    .split(',')
                    .filter(|x| !x.is_empty())
                    .map(|x| x.parse().unwrap())
                    .collect();
                c.groups.insert(k.parse().unwrap(), members);
            }
            "C" => {
                let (k, v) = rest.split_once(' ').unwrap_or((rest, ""));
                c.consensuses
                    .insert(k.parse().unwrap(), v.as_bytes().to_vec());
            }
            _ => {}
        }
    }
    c
}

/// Nodes then edges, both in record order, so adjacency order matches the
/// reference's --- `compute_equal_reads` takes the **first** out-edge carrying the
/// read, so that order decides the walk.
fn rebuild(c: &Case) -> Graph {
    let mut g = Graph::new();
    let blank = ReadInfo {
        start_mini_end: 0,
        end_mini_start: 0,
        original_support: true,
    };
    for name in &c.nodes {
        let n = g.add_node(key_from_name(name));
        // `compute_equal_reads` never reads node attributes, only edge support.
        let _ = (n, blank);
    }
    for (u, v, supp) in &c.edges {
        let (u, v) = (
            g.lookup(&key_from_name(u)).expect("edge endpoint exists"),
            g.lookup(&key_from_name(v)).expect("edge endpoint exists"),
        );
        g.upsert_edge_support(u, v, supp.clone());
    }
    g
}

#[derive(Default)]
struct Outcome {
    /// The groups themselves differ --- a read is in the wrong group, or the
    /// number of groups differs. Unambiguously this port's bug.
    grouping: Vec<String>,
    /// The groups partition the reads identically, but the *order* within a group
    /// or the chosen representative id differs. That is CPython set-iteration
    /// order, which the port does not model --- see `PORTING.md` finding 28.
    set_order: Vec<String>,
    merging: Vec<String>,
}

fn check(c: &Case) -> Outcome {
    let mut out = Outcome::default();
    let g = rebuild(c);

    // -- grouping ----------------------------------------------------------
    let mut groups = compute_equal_reads(&g, &c.support);
    let got: BTreeMap<u32, Vec<u32>> = groups
        .iter()
        .map(|(k, v)| {
            let mut v = v.clone();
            v.sort_unstable();
            (*k, v)
        })
        .collect();
    let want: BTreeMap<u32, Vec<u32>> = c
        .groups
        .iter()
        .map(|(k, v)| {
            let mut v = v.clone();
            v.sort_unstable();
            (*k, v)
        })
        .collect();
    // Compare the *partition* first: same reads in the same groups, ignoring
    // both the representative id and the order within a group. That is the part
    // this port is responsible for.
    let partition = |m: &BTreeMap<u32, Vec<u32>>| {
        let mut v: Vec<Vec<u32>> = m.values().cloned().collect();
        v.sort();
        v
    };
    if partition(&got) != partition(&want) {
        out.grouping.push(format!(
            "the partition differs: reference {} group(s), port {}",
            want.len(),
            got.len()
        ));
        for (k, w) in want.iter().take(4) {
            if let Some(h) = got.get(k) {
                if h != w {
                    out.grouping.push(format!(
                        "group {k}: only-ref {:?} | only-port {:?}",
                        w.iter()
                            .filter(|x| !h.contains(x))
                            .take(6)
                            .collect::<Vec<_>>(),
                        h.iter()
                            .filter(|x| !w.contains(x))
                            .take(6)
                            .collect::<Vec<_>>()
                    ));
                }
            }
        }
    } else {
        // Same partition. Do the labels and the within-group order agree?
        let raw_got: BTreeMap<u32, Vec<u32>> = groups.iter().cloned().collect();
        if raw_got.keys().collect::<Vec<_>>() != c.groups.keys().collect::<Vec<_>>() {
            let d: Vec<u32> = c
                .groups
                .keys()
                .filter(|k| !raw_got.contains_key(k))
                .copied()
                .take(4)
                .collect();
            out.set_order.push(format!(
                "same partition, different representative id(s): reference uses {d:?} where the port does not"
            ));
        }
        for (k, w) in &c.groups {
            if let Some(h) = raw_got.get(k) {
                if h != w {
                    out.set_order.push(format!(
                        "group {k}: same members, different order --- reference {:?}, port {:?}",
                        &w[..w.len().min(6)],
                        &h[..h.len().min(6)]
                    ));
                    break;
                }
            }
        }
    }

    // -- merging -----------------------------------------------------------
    // Only meaningful if the grouping agreed; otherwise the merge is being fed
    // different input and would report a difference that is not its own.
    if !out.grouping.is_empty() || !out.set_order.is_empty() {
        out.merging.push(
            "not evaluated: the grouping or its order differs, so the merge sees \
             different input --- and spoa is order-sensitive"
                .into(),
        );
        return out;
    }
    let mut engine = SpoaParasailMerge;
    let made = merge_consensuses(&mut engine, &mut groups, &c.reads, c.opts);
    let got: BTreeMap<u32, Vec<u8>> = made.into_iter().collect();
    if got.len() != c.consensuses.len() {
        out.merging.push(format!(
            "reference produced {} isoform(s), port {}",
            c.consensuses.len(),
            got.len()
        ));
    }
    let mut differing = 0usize;
    for (k, want) in &c.consensuses {
        match got.get(k) {
            None => out
                .merging
                .push(format!("isoform {k}: missing from the port")),
            Some(have) if have != want => {
                differing += 1;
                if differing <= 3 {
                    out.merging.push(format!(
                        "isoform {k}: consensus differs, {} vs {} bases",
                        want.len(),
                        have.len()
                    ));
                }
            }
            _ => {}
        }
    }
    if differing > 3 {
        out.merging.push(format!(
            "... and {} more differing consensuses",
            differing - 3
        ));
    }
    for k in got.keys() {
        if !c.consensuses.contains_key(k) {
            out.merging
                .push(format!("isoform {k}: extra, not in the reference"));
        }
    }
    out
}

#[test]
fn isoform_generation_matches_the_reference() {
    let Some(dir) = std::env::var_os("ISONFORM_ISOFORM_DUMPS") else {
        eprintln!(
            "SKIPPED: set ISONFORM_ISOFORM_DUMPS to a directory produced by\n\
             bench/dump_reference.py --stage isoforms to run the isoform oracle.\n\
             Without it this test asserts nothing."
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
                .is_some_and(|n| n.starts_with("isoforms_") && n.ends_with(".txt"))
        })
        .collect();
    cases.sort();
    assert!(!cases.is_empty(), "no isoforms_*.txt in {}", dir.display());

    let (mut ok, mut bad_group, mut bad_order, mut bad_merge) = (0usize, 0usize, 0usize, 0usize);
    let mut report = String::new();
    for path in &cases {
        let c = parse_case(path);
        let o = check(&c);
        if o.grouping.is_empty() && o.set_order.is_empty() && o.merging.is_empty() {
            ok += 1;
            continue;
        }
        if !o.grouping.is_empty() {
            bad_group += 1;
        }
        if !o.set_order.is_empty() {
            bad_order += 1;
        }
        if !o.merging.is_empty() && o.grouping.is_empty() && o.set_order.is_empty() {
            bad_merge += 1;
        }
        let name = c.path.file_name().unwrap().to_string_lossy();
        let _ = writeln!(
            report,
            "\n=== {name} === ({} reads, {} reference group(s), {} isoform(s))",
            c.reads.len(),
            c.groups.len(),
            c.consensuses.len()
        );
        for p in o.grouping.iter().take(4) {
            let _ = writeln!(report, "  GROUPING:  {p}");
        }
        for p in o.set_order.iter().take(2) {
            let _ = writeln!(report, "  SET-ORDER: {p}");
        }
        if o.grouping.is_empty() && o.set_order.is_empty() {
            for p in o.merging.iter().take(4) {
                let _ = writeln!(report, "  MERGING:   {p}");
            }
        }
    }
    eprintln!(
        "isoform oracle: {} case(s), {ok} fully agree | {bad_group} wrong partition, \
         {bad_order} set-order only, {bad_merge} merging",
        cases.len()
    );

    // The gate. A wrong *partition* is this port's bug and fails. A difference
    // that is only CPython's set-iteration order --- same reads in the same
    // groups, different representative or different order within the group --- is
    // reported and counted but does not fail: reproducing it means modelling
    // CPython's set internals, which is an open decision recorded as `PORTING.md`
    // finding 28, not something to pretend is already resolved.
    assert_eq!(
        bad_group, 0,
        "{bad_group} case(s) partition the reads differently:{report}"
    );
    assert_eq!(
        bad_merge, 0,
        "{bad_merge} case(s) agree on grouping and order but differ in merging:{report}"
    );
    if bad_order > 0 {
        eprintln!(
            "\nNOTE: {bad_order} case(s) differ only in CPython set-iteration order \
             (same reads, same groups). Not failing the build --- see PORTING.md \
             finding 28.{report}"
        );
    }
}
