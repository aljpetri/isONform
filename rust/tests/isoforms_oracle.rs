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
    /// The same keys in record order.
    group_order: Vec<u32>,
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
                let key: u32 = k.parse().unwrap();
                c.group_order.push(key);
                c.groups.insert(key, members);
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
    /// What the ordering divergence actually costs, for cases where it fires:
    /// `(reference isoforms, port isoforms, byte-identical consensuses)`.
    cost: Option<(usize, usize, usize)>,
}

fn check(c: &Case) -> Outcome {
    let mut out = Outcome::default();
    let g = rebuild(c);

    // -- grouping ----------------------------------------------------------
    let groups = compute_equal_reads(&g, &c.support);
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
        //
        // Compared group-by-group on *membership*, not by key: when the labels
        // differ, key `k` names a different group on each side, and diffing those
        // reports members as "missing" that are simply in the group next door.
        let raw_got: Vec<(u32, Vec<u32>)> = groups.clone();
        let mut label_diffs = 0usize;
        let mut order_diffs = 0usize;
        let mut first_order: Option<String> = None;
        for (want_k, want_v) in c
            .group_order
            .iter()
            .filter_map(|k| c.groups.get(k).map(|v| (*k, v)))
        {
            // The port's group holding the same reads, whatever it called it.
            let mut sorted_want = want_v.clone();
            sorted_want.sort_unstable();
            let Some((got_k, got_v)) = raw_got.iter().find(|(_, v)| {
                let mut sv = (*v).clone();
                sv.sort_unstable();
                sv == sorted_want
            }) else {
                continue;
            };
            if *got_k != want_k {
                label_diffs += 1;
            }
            if got_v != want_v {
                order_diffs += 1;
                if first_order.is_none() {
                    first_order = Some(format!(
                        "e.g. reference {:?}, port {:?}",
                        &want_v[..want_v.len().min(6)],
                        &got_v[..got_v.len().min(6)]
                    ));
                }
            }
        }
        if label_diffs > 0 {
            out.set_order.push(format!(
                "{label_diffs} group(s) carry a different representative id \
                 (`list(set)[0]` --- finding 28)"
            ));
        }
        if order_diffs > 0 {
            out.set_order.push(format!(
                "{order_diffs} group(s) hold the same reads in a different order. {}",
                first_order.unwrap_or_default()
            ));
        }
    }

    // -- merging -----------------------------------------------------------
    //
    // Seeded with the **reference's own** grouping, read straight from the `Q`
    // records, rather than with the port's. That is deliberate and it is what
    // makes this a real gate: the port uses ascending order where the reference
    // uses CPython's set order (finding 28), so on 28 of 114 cases the port's
    // groups reach `merge_consensuses` in a different order. Feeding the port's
    // groups here would mean those cases could never be checked at all --- the
    // merge would be judged on input it was never meant to see.
    //
    // With the reference's order in, the merge must reproduce the reference's
    // consensuses exactly, on every case. Anything else is this port's bug.
    let mut seeded: Vec<(u32, Vec<u32>)> = c.groups.iter().map(|(k, v)| (*k, v.clone())).collect();
    // `equal_reads` is a dict; `merge_consensuses` iterates it in insertion
    // order, which for the recorded dump is the order `compute_equal_reads`
    // inserted. The `Q` records are written sorted by key, so restore insertion
    // order from the recorded sequence instead.
    seeded.sort_by_key(|(k, _)| {
        c.group_order
            .iter()
            .position(|x| x == k)
            .unwrap_or(usize::MAX)
    });
    let mut engine = SpoaParasailMerge;
    let made = merge_consensuses(&mut engine, &mut seeded, &c.reads, c.opts);
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
    // -- what the ordering divergence costs ---------------------------------
    //
    // Not a gate. The port uses ascending order where the reference uses
    // CPython's set order (finding 28, decided rather than pending), so on the
    // cases where those differ the output differs too. Running the port's *own*
    // grouping through the merge and comparing measures that, instead of leaving
    // "it diverges" as an unquantified word.
    if !out.set_order.is_empty() {
        let mut own = groups.clone();
        let mut engine = SpoaParasailMerge;
        let made = merge_consensuses(&mut engine, &mut own, &c.reads, c.opts);
        let identical = made
            .iter()
            .filter(|(_, seq)| c.consensuses.values().any(|w| w == seq))
            .count();
        out.cost = Some((c.consensuses.len(), made.len(), identical));
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
    let (mut cost_ref, mut cost_port, mut cost_same) = (0usize, 0usize, 0usize);
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
        if !o.merging.is_empty() {
            bad_merge += 1;
        }
        if let Some((r, p, same)) = o.cost {
            cost_ref += r;
            cost_port += p;
            cost_same += same;
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
        for p in o.merging.iter().take(4) {
            let _ = writeln!(report, "  MERGING:   {p}");
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
        "{bad_merge} case(s) differ in merging, given the reference's own grouping:{report}"
    );
    if bad_order > 0 && cost_ref > 0 {
        eprintln!(
            "  ordering divergence, on the {bad_order} case(s) where it fires: \
             reference emitted {cost_ref} isoform(s), the port {cost_port}, \
             {cost_same} of them byte-identical"
        );
    }
    if bad_order > 0 {
        eprintln!(
            "\nNOTE: {bad_order} case(s) differ only in CPython set-iteration order \
             (same reads, same groups). Not failing the build --- see PORTING.md \
             finding 28.{report}"
        );
    }
}
