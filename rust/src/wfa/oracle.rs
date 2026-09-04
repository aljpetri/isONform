//! Does swapping in WFA2 change any **merge verdict**?
//!
//! Score equality is the wrong gate. `align_to_merge` returns a bool, and that
//! bool is what decides whether two isoforms collapse; two co-optimal alignments
//! with different CIGARs can score identically and still disagree, while two
//! alignments with different scores can agree perfectly. The earlier
//! block-aligner work gated on scores and so spent its effort on 270 cases
//! holding 0.1% of the runtime. This gates on the decision.
//!
//! ```text
//! PARASAIL_CASES=<recorded calls> cargo test --release --lib wfa::oracle -- --nocapture
//! ```
//!
//! The corpus is the recorded real `parasail_alignment` calls: two tab-separated
//! sequences plus the scoring the caller passed.

use crate::align::CigarOp;
use crate::isoforms::{align_to_merge, IsoformEngine, MergeOpts};
use crate::parasail::Scoring;

/// Drives `align_to_merge` with WFA2 in place of parasail, falling back exactly
/// as the real site would so the comparison covers the code that would ship.
struct WfaMerge {
    fell_back: usize,
}

impl IsoformEngine for WfaMerge {
    fn spoa(&mut self, _seqs: &[&[u8]]) -> Vec<u8> {
        unreachable!("the verdict gate never builds a consensus")
    }
    fn align_merge(&mut self, s1: &[u8], s2: &[u8]) -> (Vec<CigarOp>, Vec<u8>, Vec<u8>) {
        let aln = match crate::wfa::semiglobal(s1, s2, Scoring::MERGE) {
            Some(a) => a,
            None => {
                self.fell_back += 1;
                crate::parasail::semiglobal(s1, s2, Scoring::MERGE)
            }
        };
        let (a, b) = crate::align::ops_to_seq(&aln.ops, s1, s2).unwrap_or_default();
        (aln.ops, a, b)
    }
}

struct ParaMerge;

impl IsoformEngine for ParaMerge {
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
            .collect(),
    )
}

/// `MergeOpts` as `main` configures it: `delta 0.15` and the fixed cigar
/// diversity of finding 30.
fn opts() -> MergeOpts {
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
fn verdicts_match_parasail() {
    let Some(cases) = cases() else {
        eprintln!("SKIPPED: set PARASAIL_CASES to the recorded calls");
        return;
    };
    let o = opts();
    let mut wfa = WfaMerge { fell_back: 0 };
    let mut para = ParaMerge;
    let (mut agree, mut disagree, mut merges) = (0usize, 0usize, 0usize);
    let mut examples: Vec<String> = Vec::new();
    let mut disagreeing: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();

    for (s1, s2) in &cases {
        if s1.is_empty() || s2.is_empty() {
            continue;
        }
        let want = align_to_merge(&mut para, s1, s2, o);
        let got = align_to_merge(&mut wfa, s1, s2, o);
        if want {
            merges += 1;
        }
        if want == got {
            agree += 1;
        } else {
            disagree += 1;
            disagreeing.push((s1.clone(), s2.clone()));
            if examples.len() < 10 {
                examples.push(format!(
                    "    len1={} len2={} parasail={want} wfa2={got}",
                    s1.len(),
                    s2.len()
                ));
            }
        }
    }

    // Why do they disagree: a worse alignment, or a co-optimal one with a
    // different shape? Only the first is a defect this layer could fix.
    let (mut worse, mut equal, mut better) = (0usize, 0usize, 0usize);
    let (mut worst_delta, mut sum_delta) = (0i64, 0i64);
    for (s1, s2) in &disagreeing {
        let p = crate::parasail::semiglobal(s1, s2, Scoring::MERGE).score;
        let w = crate::wfa::semiglobal(s1, s2, Scoring::MERGE).map_or(p, |a| a.score);
        let d = (w - p) as i64;
        sum_delta += d;
        worst_delta = worst_delta.min(d);
        match d.cmp(&0) {
            std::cmp::Ordering::Less => worse += 1,
            std::cmp::Ordering::Equal => equal += 1,
            std::cmp::Ordering::Greater => better += 1,
        }
    }

    // Score quality over *every* case, not only the disagreeing ones. Verdicts
    // are a threshold on the score, so a change that improves alignments without
    // crossing the threshold is invisible above --- which is exactly what the
    // dangle DP (`ISONFORM_WFA2_DANGLE_DP`) is meant to do most of the time.
    let (mut sworse, mut sequal, mut sbetter) = (0usize, 0usize, 0usize);
    let (mut ssum, mut sworst) = (0i64, 0i64);
    let mut dangle_hist: Vec<usize> = Vec::new();
    let mut two_sided = 0usize;
    for (s1, s2) in &cases {
        if s1.is_empty() || s2.is_empty() {
            continue;
        }
        let p = crate::parasail::semiglobal(s1, s2, Scoring::MERGE).score;
        let Some(w) = crate::wfa::semiglobal(s1, s2, Scoring::MERGE) else {
            continue; // declined; the call site uses parasail, so nothing is lost
        };
        // How much sequence is left *outside* the aligned core: the leading and
        // trailing pure-gap runs of the returned alignment. This is the material
        // the dangle DP has to work with, so if it is ~0 on the cases WFA2 loses,
        // the loss is in the core and no end repair can reach it.
        let al = |&(_, t): &(usize, u8)| t == b'=' || t == b'X';
        let lead: usize = w.ops.iter().take_while(|o| !al(o)).map(|&(n, _)| n).sum();
        let trail: usize = w
            .ops
            .iter()
            .rev()
            .take_while(|o| !al(o))
            .map(|&(n, _)| n)
            .sum();
        // Two-sided or one-sided? `anchored_dp` can only align a dangle where
        // BOTH sequences contribute bases; a pure overhang (one side only) is a
        // free end for parasail too, so there is nothing there to win.
        let side = |it: &mut dyn Iterator<Item = (usize, u8)>| -> (usize, usize) {
            let (mut i, mut dd) = (0usize, 0usize);
            for (n, t) in it {
                match t {
                    b'I' => i += n,
                    b'D' => dd += n,
                    _ => break,
                }
            }
            (i, dd)
        };
        let (li, ld) = side(&mut w.ops.iter().copied());
        let (ti, td) = side(&mut w.ops.iter().rev().copied());
        let d = (w.score - p) as i64;
        if d < 0 {
            dangle_hist.push(lead + trail);
            if li.min(ld) > 0 || ti.min(td) > 0 {
                two_sided += 1;
            }
        }
        ssum += d;
        sworst = sworst.min(d);
        match d.cmp(&0) {
            std::cmp::Ordering::Less => sworse += 1,
            std::cmp::Ordering::Equal => sequal += 1,
            std::cmp::Ordering::Greater => sbetter += 1,
        }
    }
    let scored = sworse + sequal + sbetter;
    eprintln!("score quality over {scored} non-declined cases");
    eprintln!("  worse:  {sworse}");
    eprintln!("  equal:  {sequal}");
    eprintln!("  better: {sbetter}  (must be 0; parasail is exact)");
    eprintln!(
        "  mean delta {:.4}, worst {sworst}",
        ssum as f64 / scored.max(1) as f64
    );

    if !dangle_hist.is_empty() {
        dangle_hist.sort_unstable();
        let n = dangle_hist.len();
        let zero = dangle_hist.iter().filter(|&&v| v == 0).count();
        eprintln!(
            "  unaligned end bases on the {n} worse cases: zero on {zero} ({:.1}%), \
             median {}, p90 {}, max {}; two-sided (DP can act) on {two_sided}",
            100.0 * zero as f64 / n as f64,
            dangle_hist[n / 2],
            dangle_hist[n * 9 / 10],
            dangle_hist[n - 1]
        );
    }

    eprintln!("verdict gate over {} cases", agree + disagree);
    eprintln!("  agree:      {agree}");
    eprintln!("  disagree:   {disagree}");
    eprintln!("  parasail merged: {merges}  (the verdicts that can actually change output)");
    eprintln!("  wfa2 fell back to parasail: {}", wfa.fell_back);
    eprintln!("  of the disagreements, WFA2 scored:");
    eprintln!("    worse:  {worse}   (a suboptimal alignment --- fixable in principle)");
    eprintln!("    equal:  {equal}   (co-optimal, different shape --- a tie-break difference)");
    eprintln!("    better: {better}  (must be 0; parasail is exact)");
    if !disagreeing.is_empty() {
        eprintln!(
            "    worst delta {worst_delta}, mean {:.2}",
            sum_delta as f64 / disagreeing.len() as f64
        );
    }
    for e in &examples {
        eprintln!("{e}");
    }
    // Left as a report, not a hard gate: the decision is whether the rate is
    // tolerable end to end, which this test cannot answer on its own.
    eprintln!(
        "  verdict flip rate {:.2}% of all calls, {:.2}% of parasail's merges",
        100.0 * disagree as f64 / (agree + disagree) as f64,
        100.0 * disagree as f64 / merges as f64
    );
}
