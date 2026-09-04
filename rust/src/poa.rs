//! Consensus by partial-order alignment, matching `IsoformGeneration.run_spoa`
//! (`modules/IsoformGeneration.py:143`).
//!
//! That citation is deliberate: this module arrived from the isONcorrect port
//! naming `create_augmented_reference.run_spoa`, which in isONform is dead code
//! in a module nothing imports (PORTING.md, reconnaissance correction 2). The
//! invocation is character-for-character the same, so the code was right and only
//! the attribution was wrong — but the live one is what a reader needs to find.
//!
//! The reference shells out per call:
//!
//! ```text
//! spoa <reads.fa> -l 0 -r 0 -g -2
//! ```
//!
//! Resolved against spoa's CLI defaults (`m=5, n=-4, g=-8, e=-6, q=-10, c=-4`)
//! and `AlignmentEngine::Create`'s subtype dispatch (`g >= e` ⇒ linear, then
//! `e := g`), the effective configuration is **local (kSW) alignment, linear gap
//! −2, match +5, mismatch −4, consensus only**. No affine, no convex, and the
//! MSA output path is unused.
//!
//! This wraps [`spoars`], a pure-Rust reimplementation of spoa, rather than
//! binding the C++ library. That choice is backed by a differential oracle, not
//! by the crate's own claim, and it has now been measured twice — on two
//! different tools' inputs, which is the part that matters:
//!
//! * isONcorrect: 505 of 505 real correction intervals.
//! * **isONform: 3 277 of 3 277 recorded calls**, across `sirv_small`, 56
//!   Drosophila clusters and 7 real SIRV clusters — 2 to 96 sequences per call,
//!   consensus 20 to 1 742 bases. All three live call sites are covered, 3 074 of
//!   the cases coming from `SimplifyGraph.py:570` (bubble-path consensus), which
//!   is the one the simplification oracle's verdicts depend on.
//!
//! The isONcorrect number does **not** imply the isONform one: different
//! sequences from a different stage, so it had to be re-measured rather than
//! inherited. Recording is `bench/dump_reference.py --record-spoa`; the `oracle`
//! test below replays it from `SPOA_CASES`.
//!
//! # Measured and rejected: reusing the engine across intervals
//!
//! `SimdEngine` owns grow-only DP scratch that it reuses across `align` calls, so
//! building one per interval — 11 775 per ten real clusters — looks like pure
//! waste. Hoisting it into a thread-local is **worth nothing**: `run_spoa` went
//! 21.28 s to 21.21 s, inside the noise. The per-call allocation is not the cost;
//! the POA dynamic programming is.
//!
//! That also bears on whether to write a native POA (see *Deferred improvements*).
//! Allocation reuse was one of the two arguments for it, and this measures that
//! argument at zero. A native implementation constrained to reproduce spoa's
//! output exactly would do the same asymptotic work.
//!
//! **Sequence insertion order changes the consensus.** POA is order-sensitive,
//! so the order `get_best_corrections` writes segments — which is `to_add`'s
//! insertion order from `find_most_supported_span` — must be preserved,
//! including the `i > max_seqs_to_spoa` cutoff (note `>`, not `>=`, so it admits
//! `max_seqs_to_spoa + 1` sequences).

use spoars::align::{AlignmentEngine, AlignmentType, Scoring, SimdEngine};
use spoars::graph::Graph;

/// spoa's CLI defaults for the parameters isONcorrect does not override.
const MATCH: i8 = 5;
const MISMATCH: i8 = -4;
/// `-g -2`, and linear means the extension penalty equals it.
const GAP: i8 = -2;

/// Consensus over `seqs`, in the order given.
///
/// Returns `None` for an empty input, matching the fact that the reference
/// never calls spoa with no sequences.
pub fn consensus<S: AsRef<[u8]>>(seqs: &[S]) -> Option<String> {
    if seqs.is_empty() {
        return None;
    }
    // Linear gaps: spoa collapses e, q and c onto g for this subtype.
    let scoring = Scoring::new(MATCH, MISMATCH, GAP, GAP, GAP, GAP).ok()?;
    let mut engine = SimdEngine::new(AlignmentType::Local, scoring);
    let mut graph = Graph::new();

    for s in seqs {
        let s = s.as_ref();
        if s.is_empty() {
            continue;
        }
        let (aln, _score) = engine.align(s, &graph);
        let weights = vec![1u32; s.len()];
        graph.add_alignment(&aln, s, &weights).ok()?;
    }
    Some(graph.generate_consensus())
}

/// Segments to hand to the POA, applying the reference's cutoff.
///
/// `get_best_corrections` breaks when `i > max_seqs_to_spoa`, so it admits
/// `max_seqs_to_spoa + 1` sequences — an off-by-one worth preserving.
pub fn apply_spoa_cutoff<S: Clone>(segments: &[S], max_seqs_to_spoa: usize) -> Vec<S> {
    segments
        .iter()
        .take(max_seqs_to_spoa + 1)
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_has_no_consensus() {
        let empty: [&[u8]; 0] = [];
        assert!(consensus(&empty).is_none());
    }

    #[test]
    fn a_single_sequence_is_its_own_consensus() {
        assert_eq!(consensus(&[b"ACGTACGT".as_slice()]).unwrap(), "ACGTACGT");
    }

    #[test]
    fn identical_sequences_yield_that_sequence() {
        let s = b"ACGTACGTACGT".as_slice();
        assert_eq!(consensus(&[s, s, s]).unwrap(), "ACGTACGTACGT");
    }

    #[test]
    fn the_majority_base_wins() {
        // Three reads agree on position 4, one disagrees.
        let seqs: Vec<&[u8]> = vec![b"ACGTAACGT", b"ACGTAACGT", b"ACGTAACGT", b"ACGTCACGT"];
        assert_eq!(consensus(&seqs).unwrap(), "ACGTAACGT");
    }

    #[test]
    fn cutoff_admits_one_more_than_the_limit() {
        // `if i > max_seqs_to_spoa: break` runs the body for i == limit.
        let segs: Vec<usize> = (0..10).collect();
        assert_eq!(apply_spoa_cutoff(&segs, 3).len(), 4);
        assert_eq!(apply_spoa_cutoff(&segs, 0).len(), 1);
        // Fewer segments than the limit is left alone.
        assert_eq!(apply_spoa_cutoff(&segs, 100).len(), 10);
    }
}

/// Differential oracle against the real `spoa` binary.
///
/// Generate the cases inside the reference environment, then run:
///
/// ```bash
/// bench/dump_reference.py --fastq test_data/isoncorrect/0.fastq --outdir /tmp/poa
/// SPOA_CASES=/tmp/poa/spoa_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml poa::oracle -- --nocapture
/// ```
///
/// This is the check that justifies using `spoars` at all, so it should be
/// re-run whenever the crate is upgraded or a new corpus becomes available.
#[cfg(test)]
mod oracle {
    use super::*;

    #[test]
    fn matches_the_spoa_binary_on_recorded_intervals() {
        let Ok(path) = std::env::var("SPOA_CASES") else {
            return; // no recorded cases in this checkout
        };
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("SPOA_CASES={path} unreadable: {e}"));

        let (mut n, mut bad) = (0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let mut f = line.split('\t');
            let Some(want) = f.next() else { continue };
            let seqs: Vec<&[u8]> = f.map(str::as_bytes).collect();
            if seqs.is_empty() {
                continue;
            }
            let got = consensus(&seqs).expect("non-empty input");
            if got != want {
                if bad < 3 {
                    eprintln!(
                        "mismatch on {} seqs:\n  spoa  : {want}\n  spoars: {got}",
                        seqs.len()
                    );
                }
                bad += 1;
            }
            n += 1;
        }
        assert!(n > 0, "no cases in {path}");
        assert_eq!(
            bad, 0,
            "{bad} of {n} consensus sequences differed from spoa"
        );
        eprintln!("consensus matched the spoa binary on all {n} recorded intervals");
    }
}
