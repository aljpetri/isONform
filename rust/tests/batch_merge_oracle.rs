//! Differential oracle for batch merging, the last stage.
//!
//! Replays the isoforms the reference had on entry through
//! [`isonform::batch_merge`] and diffs two things:
//!
//! * the state **after** `actual_merging_process` (the `A` records) --- which is
//!   the state before it, because that function is a no-op;
//! * every record `write_final_output` emitted (the `F` records): destination,
//!   id, support and sequence, in order.
//!
//! ```text
//! bench/dump_reference.py --fastq-folder bench/corpus/sirv_small \
//!     --outdir /tmp/d --stage batchmerge
//! ISONFORM_BATCHMERGE_DUMPS=/tmp/d cargo test --test batch_merge_oracle -- --nocapture
//! ```
//!
//! Skips loudly with the variable unset, like the others.
//!
//! # The dump comes from `isONform_parallel`, not `main`
//!
//! This is the one stage `main` never reaches --- its call is commented out at
//! `main:583`. `isONform_parallel` forks the per-cluster work but calls
//! `join_back_via_batch_merging` in the parent, after `pool.join()`, so a
//! recording wrapper survives to see it.
//!
//! # It asserts the *reference's* merging does nothing, on purpose
//!
//! The reference's `actual_merging_process` body is unreachable (`PORTING.md`
//! finding 31), so its `A` records equal its `B` records. This replays with
//! `no_op = true` and checks that, which keeps it a real gate on the recorded
//! data: if a future reference ever fixes the guard, this fails loudly instead of
//! the dump silently disagreeing.
//!
//! The port itself **does** merge across batches by default; that repair is
//! covered by the unit tests in `batch_merge.rs` and measured in `PORTING.md`.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use isonform::batch_merge::{actual_merging_process, select_output, Destination, Isoform};
use isonform::isoforms::{MergeOpts, SpoaParasailMerge};

#[derive(Default)]
struct Case {
    path: PathBuf,
    iso_abundance: usize,
    write_low_abundance: bool,
    cluster: String,
    /// `(batch_id, [(isoform_id, isoform)])`, in record order.
    before: Vec<(u32, Vec<(u32, Isoform)>)>,
    /// `(batch_id, isoform_id, merged, n_reads)` after merging.
    after: Vec<(u32, u32, bool, usize)>,
    /// `(destination, id, support, sequence)`.
    written: Vec<(String, String, usize, Vec<u8>)>,
}

fn push_isoform(v: &mut Vec<(u32, Vec<(u32, Isoform)>)>, b: u32, i: u32, iso: Isoform) {
    match v.iter_mut().find(|(k, _)| *k == b) {
        Some(slot) => slot.1.push((i, iso)),
        None => v.push((b, vec![(i, iso)])),
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
                    "iso_abundance" => c.iso_abundance = v.parse().unwrap(),
                    "low" => c.write_low_abundance = v == "1",
                    "cluster" => c.cluster = v.to_string(),
                    _ => {}
                }
            }
            continue;
        }
        let f: Vec<&str> = line.splitn(2, ' ').collect();
        if f.len() != 2 {
            continue;
        }
        match f[0] {
            "B" => {
                let p: Vec<&str> = f[1].splitn(5, ' ').collect();
                let n_reads: usize = p[3].parse().unwrap();
                push_isoform(
                    &mut c.before,
                    p[0].parse().unwrap(),
                    p[1].parse().unwrap(),
                    Isoform {
                        sequence: p.get(4).map(|s| s.as_bytes().to_vec()).unwrap_or_default(),
                        // Accessions are not recorded, only the count --- the
                        // stage uses `len(reads)` and never the names, except in
                        // the mapping file which this oracle does not compare.
                        reads: (0..n_reads).map(|i| i.to_string()).collect(),
                        merged: p[2] == "1",
                    },
                );
            }
            "A" => {
                let p: Vec<&str> = f[1].split_whitespace().collect();
                c.after.push((
                    p[0].parse().unwrap(),
                    p[1].parse().unwrap(),
                    p[2] == "1",
                    p[3].parse().unwrap(),
                ));
            }
            "F" => {
                let p: Vec<&str> = f[1].splitn(4, ' ').collect();
                c.written.push((
                    p[0].to_string(),
                    p[1].to_string(),
                    p[2].parse().unwrap(),
                    p.get(3).map(|s| s.as_bytes().to_vec()).unwrap_or_default(),
                ));
            }
            _ => {}
        }
    }
    c
}

fn dest_name(d: Destination) -> &'static str {
    match d {
        Destination::Main => "main",
        Destination::LowAbundance => "low",
        Destination::Dropped => "dropped",
    }
}

fn check(c: &Case) -> Vec<String> {
    let mut problems = Vec::new();
    let mut batches = c.before.clone();

    // The oracle replays the *reference*, whose merging step never executes
    // (finding 31). The port fixes that by default, so this asks for the bug
    // back explicitly --- exactly as the interval and graph oracles ask for
    // `WisOpts::reference()`.
    let mut engine = SpoaParasailMerge;
    actual_merging_process(
        &mut engine,
        &mut batches,
        // Unused with `no_op = true`; the reference's dump records no merge
        // thresholds because its merge never ran to need them.
        MergeOpts::default(),
        true,
    );

    // The `A` records: what the reference had after merging. Since merging is a
    // no-op this equals `B`, and the point of checking is to notice if that ever
    // stops being true upstream.
    let got_after: Vec<(u32, u32, bool, usize)> = batches
        .iter()
        .flat_map(|(b, v)| {
            v.iter()
                .map(move |(i, iso)| (*b, *i, iso.merged, iso.reads.len()))
        })
        .collect();
    let mut want_after = c.after.clone();
    let mut got_sorted = got_after.clone();
    want_after.sort();
    got_sorted.sort();
    if got_sorted != want_after {
        let changed = c
            .after
            .iter()
            .filter(|(b, i, m, n)| {
                !c.before.iter().any(|(bb, v)| {
                    bb == b
                        && v.iter()
                            .any(|(ii, iso)| ii == i && iso.merged == *m && iso.reads.len() == *n)
                })
            })
            .count();
        problems.push(format!(
            "post-merge state differs: reference {} entries, port {} \
             ({changed} of the reference's changed during merging --- if that is \
             non-zero the reference no longer no-ops, and finding 31 needs revisiting)",
            want_after.len(),
            got_sorted.len()
        ));
    }

    // The written records.
    let got = select_output(&batches, &c.cluster, c.iso_abundance, c.write_low_abundance);
    if got.len() != c.written.len() {
        problems.push(format!(
            "reference wrote {} record(s), port {}",
            c.written.len(),
            got.len()
        ));
    }
    for (i, (want, have)) in c.written.iter().zip(got.iter()).enumerate() {
        let (wdest, wid, wsupp, wseq) = want;
        if *wid != have.id {
            problems.push(format!("record {i}: id {wid:?} vs {:?}", have.id));
        } else if wdest != dest_name(have.destination) {
            problems.push(format!(
                "record {i} ({wid}): destination {wdest} vs {}",
                dest_name(have.destination)
            ));
        } else if *wsupp != have.support {
            problems.push(format!(
                "record {i} ({wid}): support {wsupp} vs {}",
                have.support
            ));
        } else if *wseq != have.sequence {
            problems.push(format!(
                "record {i} ({wid}): sequence differs, {} vs {} bases",
                wseq.len(),
                have.sequence.len()
            ));
        }
    }
    problems
}

#[test]
fn batch_merging_matches_the_reference() {
    let Some(dir) = std::env::var_os("ISONFORM_BATCHMERGE_DUMPS") else {
        eprintln!(
            "SKIPPED: set ISONFORM_BATCHMERGE_DUMPS to a directory produced by\n\
             bench/dump_reference.py --stage batchmerge to run the batch-merge\n\
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
                .is_some_and(|n| n.starts_with("batchmerge_") && n.ends_with(".txt"))
        })
        .collect();
    cases.sort();
    assert!(
        !cases.is_empty(),
        "no batchmerge_*.txt in {}",
        dir.display()
    );

    let (mut ok, mut bad, mut records) = (0usize, 0usize, 0usize);
    let mut report = String::new();
    for path in &cases {
        let c = parse_case(path);
        records += c.written.len();
        let problems = check(&c);
        if problems.is_empty() {
            ok += 1;
            continue;
        }
        bad += 1;
        let name = c.path.file_name().unwrap().to_string_lossy();
        let _ = writeln!(
            report,
            "\n=== {name} === (cluster {}, {} record(s))",
            c.cluster,
            c.written.len()
        );
        for p in problems.iter().take(8) {
            let _ = writeln!(report, "  {p}");
        }
    }
    eprintln!(
        "batch-merge oracle: {} case(s), {ok} agree, {bad} disagree ({records} output records)",
        cases.len()
    );
    assert_eq!(bad, 0, "{bad} case(s) disagree:{report}");
}
