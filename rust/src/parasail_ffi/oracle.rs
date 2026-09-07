//! Does linking parasail change any **verdict**?
//!
//! Finding 41: score equality is the wrong gate, and so is CIGAR equality. What
//! ships is two booleans --- whether two consensuses merge
//! (`IsoformGeneration.py:381`) and whether a bubble pops
//! (`SimplifyGraph.py:657`) --- and those are what have to agree. An aligner can
//! return a different co-optimal CIGAR and decide identically, or return the
//! same score and decide differently.
//!
//! The FFI backend *should* agree on all three of score, CIGAR and verdict,
//! because it is the same C library the reference calls rather than a second
//! implementation of it. This is the check that says so rather than assuming
//! it, on isONform's own recorded calls:
//!
//! ```text
//! PARASAIL_CASES=<recorded calls> cargo test --release --lib parasail_ffi::oracle -- --nocapture
//! ```

use crate::align::CigarOp;
use crate::isoforms::{align_to_merge, IsoformEngine, MergeOpts};
use crate::parasail::Scoring;

/// `align_to_merge` driven by the linked C library.
struct FfiMerge;
/// `align_to_merge` driven by the scalar reimplementation.
struct ScalarMerge;

impl IsoformEngine for FfiMerge {
    fn spoa(&mut self, _seqs: &[&[u8]]) -> Vec<u8> {
        unreachable!("the verdict gate never builds a consensus")
    }
    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let aln = super::semiglobal(s1, s2, Scoring::MERGE);
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        (aln.ops, a, b)
    }
}

impl IsoformEngine for ScalarMerge {
    fn spoa(&mut self, _seqs: &[&[u8]]) -> Vec<u8> {
        unreachable!()
    }
    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let aln = crate::parasail::semiglobal(s1, s2, Scoring::MERGE);
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        (aln.ops, a, b)
    }
}

fn cases() -> Option<Vec<(Vec<u8>, Vec<u8>)>> {
    let path = std::env::var("PARASAIL_CASES").ok()?;
    let raw = std::fs::read_to_string(path).ok()?;
    Some(
        raw.lines()
            .filter(|l| !l.starts_with('#'))
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (f.len() >= 2).then(|| (f[0].as_bytes().to_vec(), f[1].as_bytes().to_vec()))
            })
            .filter(|(a, b)| !a.is_empty() && !b.is_empty())
            .collect(),
    )
}

/// The bubble verdict, as `SimplifyGraph.py:657` reaches it: align at the bubble
/// scoring, then `parse_cigar_diversity` at the hard-coded 0.20 / `delta_len`.
fn bubble_verdict(cigar: &[(u32, u8)]) -> bool {
    crate::simplify::parse_cigar_diversity(cigar, 0.20, 5)
}

fn merge_opts() -> MergeOpts {
    MergeOpts {
        delta: 0.15,
        delta_len: 5,
        delta_iso_len_3: 30,
        delta_iso_len_5: 50,
        max_seqs_to_spoa: 200,
        merge_rebuild_max: 50,
        final_consensus_pass: false,
        cigar_diversity_counts_runs: false,
    }
}

#[test]
fn verdicts_match_the_scalar_implementation() {
    let Some(cases) = cases() else {
        eprintln!("SKIPPED: set PARASAIL_CASES to the recorded calls");
        return;
    };
    let o = merge_opts();
    let (mut m_agree, mut m_disagree, mut m_merges) = (0usize, 0usize, 0usize);
    let (mut b_agree, mut b_disagree, mut b_pops) = (0usize, 0usize, 0usize);
    // Score and CIGAR are reported too --- not as the gate, but because for this
    // backend they are expected to be exact, and a drop below 100% would say the
    // FFI is not calling what we think it is.
    let (mut score_eq, mut cigar_eq) = (0usize, 0usize);
    let mut examples: Vec<String> = Vec::new();

    for (s1, s2) in &cases {
        // --- the merge verdict
        let want = align_to_merge(&mut ScalarMerge, s1, s2, o);
        let got = align_to_merge(&mut FfiMerge, s1, s2, o);
        if want {
            m_merges += 1;
        }
        if want == got {
            m_agree += 1;
        } else {
            m_disagree += 1;
            if examples.len() < 10 {
                examples.push(format!(
                    "    merge: len1={} len2={} scalar={want} ffi={got}",
                    s1.len(),
                    s2.len()
                ));
            }
        }

        // --- the bubble verdict, at the bubble scoring
        let ps = crate::parasail::semiglobal(s1, s2, Scoring::BUBBLE);
        let fs = super::semiglobal(s1, s2, Scoring::BUBBLE);
        let pc: Vec<(u32, u8)> = ps.ops.iter().map(|&(n, t)| (n as u32, t)).collect();
        let fc: Vec<(u32, u8)> = fs.ops.iter().map(|&(n, t)| (n as u32, t)).collect();
        let (bw, bg) = (bubble_verdict(&pc), bubble_verdict(&fc));
        if bw {
            b_pops += 1;
        }
        if bw == bg {
            b_agree += 1;
        } else {
            b_disagree += 1;
            if examples.len() < 10 {
                examples.push(format!(
                    "    bubble: len1={} len2={} scalar={bw} ffi={bg}",
                    s1.len(),
                    s2.len()
                ));
            }
        }
        if ps.score == fs.score {
            score_eq += 1;
        }
        if ps.cigar == fs.cigar {
            cigar_eq += 1;
        }
    }

    let n = cases.len();
    let pct = |k: usize| 100.0 * k as f64 / n as f64;
    println!(
        "parasail-ffi oracle: {n} cases\n  \
         merge verdicts:  {m_agree} agree, {m_disagree} disagree ({m_merges} merge)\n  \
         bubble verdicts: {b_agree} agree, {b_disagree} disagree ({b_pops} pop)\n  \
         score==scalar:   {score_eq} ({:.1}%)\n  \
         cigar==scalar:   {cigar_eq} ({:.1}%)",
        pct(score_eq),
        pct(cigar_eq)
    );
    for e in &examples {
        println!("{e}");
    }
    assert_eq!(m_disagree, 0, "merge verdicts must not change");
    assert_eq!(b_disagree, 0, "bubble verdicts must not change");
}
