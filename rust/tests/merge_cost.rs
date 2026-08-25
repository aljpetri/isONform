//! Cross-batch merge cost on one real cluster, with and without the filter.
//!
//! Repairing finding 31 traded a correctness bug for a performance cliff: the
//! merge is O(batch pairs × isoforms per batch²), and on `droso_deep`'s deepest
//! cluster that measured **26 batches, ~133 isoforms each, ~5.7M alignments, 93
//! minutes at 99% CPU without finishing**. This is the harness for fixing that
//! without a multi-hour pipeline run.
//!
//! ```text
//! ISONFORM_MERGE_COST_DIR=<a cluster's output folder> \
//!   cargo test --release --test merge_cost -- --nocapture
//! ```
//!
//! The folder is one produced by `main` under `--parallel`: `{id}_batch`,
//! `spoa{id}`, `mapping{id}` per batch. A killed pipeline run leaves them behind,
//! which is how this gets a real deep cluster to work on.
//!
//! It reports, for the filter off and on: wall time, pairs seen, pairs actually
//! aligned, merges, and — the number that decides whether the filter is
//! acceptable — how many skipped pairs *would* have merged.

use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::time::Instant;

use isonform::batch_merge::{
    actual_merging_process, MERGES, PAIRS_ALIGNED, PAIRS_SEEN, PAIRS_SKIPPED, SKIPPED_WOULD_MERGE,
};
use isonform::isoforms::{MergeOpts, SpoaParasailMerge};
use isonform::parallel::read_cluster;

fn reset() {
    for c in [
        &PAIRS_SEEN,
        &PAIRS_SKIPPED,
        &PAIRS_ALIGNED,
        &SKIPPED_WOULD_MERGE,
        &MERGES,
    ] {
        c.store(0, Ordering::Relaxed);
    }
}

#[test]
fn cross_batch_merge_cost() {
    let Some(dir) = std::env::var_os("ISONFORM_MERGE_COST_DIR") else {
        eprintln!(
            "SKIPPED: set ISONFORM_MERGE_COST_DIR to a cluster output folder\n\
             (one holding `{{id}}_batch`, `spoa{{id}}`, `mapping{{id}}` files).\n\
             Without it this test measures nothing."
        );
        return;
    };
    let dir = PathBuf::from(dir);
    let cluster = read_cluster(&dir).expect("cluster readable");
    let n_batches = cluster.batches.len();
    let n_iso: usize = cluster.batches.iter().map(|(_, v)| v.len()).sum();
    assert!(n_batches > 0, "no batch files in {}", dir.display());
    eprintln!(
        "cluster {}: {n_batches} batches, {n_iso} isoforms, {} batch pairs",
        dir.display(),
        n_batches * (n_batches - 1) / 2
    );

    // `--delta 0.15` is what `main` hard-codes; the parallel driver's default is
    // 0.1. Using main's keeps this comparable to the pipeline run it came from.
    let opts = MergeOpts {
        delta: 0.15,
        delta_len: 5,
        delta_iso_len_3: 30,
        delta_iso_len_5: 50,
        max_seqs_to_spoa: 200,
        cigar_diversity_counts_runs: false,
    };

    let mut batches = cluster.batches.clone();
    let mut engine = SpoaParasailMerge;
    reset();
    let t0 = Instant::now();
    actual_merging_process(&mut engine, &mut batches, opts, false);
    let elapsed = t0.elapsed().as_secs_f64();

    let seen = PAIRS_SEEN.load(Ordering::Relaxed);
    let skipped = PAIRS_SKIPPED.load(Ordering::Relaxed);
    let aligned = PAIRS_ALIGNED.load(Ordering::Relaxed);
    let fns = SKIPPED_WOULD_MERGE.load(Ordering::Relaxed);
    let merges = MERGES.load(Ordering::Relaxed);
    let surviving = batches
        .iter()
        .flat_map(|(_, v)| v)
        .filter(|(_, i)| !i.merged)
        .count();

    eprintln!("  wall:              {elapsed:.1}s");
    eprintln!("  pairs seen:        {seen}");
    eprintln!(
        "  pairs aligned:     {aligned}  ({:.1}% of seen)",
        if seen > 0 {
            100.0 * aligned as f64 / seen as f64
        } else {
            0.0
        }
    );
    eprintln!(
        "  pairs skipped:     {skipped}  ({:.1}% of seen)",
        if seen > 0 {
            100.0 * skipped as f64 / seen as f64
        } else {
            0.0
        }
    );
    eprintln!("  merges:            {merges}");
    eprintln!("  isoforms surviving: {surviving} of {n_iso}");
    if skipped > 0 {
        eprintln!(
            "  skipped-but-would-merge: {fns}   <-- the filter's false negatives \
             (0 unless ISONFORM_SKETCH_VERIFY=1, which pays for the skipped alignments)"
        );
    }
}
