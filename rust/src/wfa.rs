//! WFA2 at the two hot alignment sites, reconciled to parasail's objective.
//!
//! # Why this is not a drop-in
//!
//! Both hot sites --- [`crate::isoforms::SpoaParasailMerge::align_merge`] and
//! `align_bubble_nodes` --- call `parasail.sg_trace_scan_16`, whose semi-global
//! objective **maximises a score** with a positive match reward. WFA2 instead
//! **minimises a penalty** with `match = 0`. Those are not the same problem, and
//! the difference is not cosmetic: with all four ends free, the empty alignment
//! costs nothing and is therefore always optimal for WFA2. Measured, on two
//! *identical* 20bp sequences:
//!
//! ```text
//! cigar=DDDDDDDDDDDDDDDDDDDDIIIIIIIIIIIIIIIIIIII    (aligns nothing, penalty 0)
//! ```
//!
//! Bounding the free ends removes that degeneracy, but a second gap remains:
//! WFA2 pays nothing either way for terminal *matches*, so it leaves them
//! unaligned, while parasail is paid `match_score` per column to include them.
//! On two identical 200bp sequences with a 50-base budget it shaved 48 bases off
//! the ends --- the budget, essentially in full.
//!
//! # The reconciliation, and why it is exact rather than approximate
//!
//! Under both scoring schemes in use here, parasail never pays penalty to buy
//! coverage: each aligned column is worth `match_score / 2 = 1` per base, while
//! the cheapest way to extend is a mismatch at 4 (MERGE) or 10 (BUBBLE), and a
//! length-`L` gap costs `open + ext*(L-1)`. So a parasail-optimal alignment is
//! always among the **minimum-penalty** alignments, and differs from WFA2's only
//! in preferring maximum coverage within that set. Hence:
//!
//! 1. bound the free ends, so the empty alignment is unreachable;
//! 2. greedily extend the aligned core outward over matching base pairs, which
//!    is free in WFA2's model and paid in parasail's;
//! 3. score the resulting columns with **parasail's own rules**.
//!
//! Step 3 is what makes the penalty transform unnecessary. Rather than mapping
//! WFA2's penalty back into a parasail score --- which is where the free-end
//! bookkeeping went wrong twice --- WFA2 is used only to *find* the column
//! assignment, and [`Scoring`] then scores it. A disagreement with parasail can
//! then only mean the two found different alignments, never that the arithmetic
//! drifted.
//!
//! # What is NOT guaranteed
//!
//! Step 2 extends over exact matches only, so it reaches maximum coverage when
//! the shaved region is pure matches --- the case that motivated it. It is not
//! proven maximal when that region also holds mismatches or indels, and bounding
//! the ends is itself a departure from parasail's unbounded `sg`. Both are
//! reasons this is gated on merge verdicts rather than assumed; see
//! [`oracle::verdicts_match_parasail`].

use biodiff_wfa2_sys as sys;
use std::ffi::c_char;

use crate::align::CigarOp;
use crate::parasail::{Alignment, Scoring};

/// How many bases at each end may go unpenalised.
///
/// Tied to [`crate::isoforms::MergeOpts::delta_iso_len_5`], the largest end
/// divergence the merge tolerates before rejecting a pair outright: a pair whose
/// ends differ by more than this is refused regardless, so a wider budget could
/// not change a verdict. Not unbounded, because unbounded is degenerate (above).
pub const FREE_ENDS: i32 = 50;

/// Is WFA2 enabled at the hot sites? **On by default.**
///
/// 4.5x faster than parasail on the 54 884 recorded calls, and it flips 2.19% of
/// merge verdicts (finding 41). End to end, with the finding-44 POA schedule and
/// finding-48 node identity, that combination is **2.35x--3.74x faster than the
/// coordinate-keyed parasail port** on the three corpora --- droso 24.9s -> 10.6s,
/// sirv_sim_gene 107.8s -> 29.4s, sirv_real 220.1s -> 58.9s --- for one real
/// accuracy cost: `sirv_real` strict F1 0.735 -> 0.711, of which 0.017 is WFA2's
/// and 0.006 finding 48's.
///
/// `ISONFORM_WFA2=0` restores parasail. The known fix for the larger half of that
/// cost is an exact DP over the `<=50`-base dangles (finding 41), not yet written.
pub fn enabled() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| {
        !matches!(std::env::var("ISONFORM_WFA2").ok().as_deref(), Some("0") | Some("off"))
    })
}

/// One WFA2 aligner, reused across calls. Creating one per alignment would cost
/// more than the alignment.
struct Aligner {
    raw: *mut sys::wavefront_aligner_t,
    sc: Scoring,
}

impl Aligner {
    /// WFA2 takes penalties (positive, `match = 0`); parasail takes scores. The
    /// standard affine mapping, given parasail charges `open + ext*(L-1)` per
    /// gap: `mismatch = a + b`, `open = o - e`, `ext = e + a/2`.
    fn new(sc: Scoring) -> Self {
        unsafe {
            let mut attr = sys::wavefront_aligner_attr_default;
            attr.distance_metric = sys::distance_metric_t_gap_affine;
            attr.affine_penalties.match_ = 0;
            attr.affine_penalties.mismatch = sc.match_score - sc.mismatch;
            attr.affine_penalties.gap_opening = sc.open - sc.ext;
            attr.affine_penalties.gap_extension = sc.ext + sc.match_score / 2;
            attr.alignment_scope = sys::alignment_scope_t_compute_alignment;
            Self { raw: sys::wavefront_aligner_new(&mut attr), sc }
        }
    }

    /// The raw column string, in WFA2's own `M`/`X`/`I`/`D` alphabet.
    fn columns(&mut self, s1: &[u8], s2: &[u8]) -> Vec<u8> {
        unsafe {
            // Each bound must not exceed its own sequence, or WFA2 aborts the
            // process with a diagnostic rather than returning an error.
            //
            // It must also leave a core that cannot be skipped. Clamping only to
            // the length lets a short sequence be *entirely* free, and then the
            // empty alignment costs nothing and wins --- finding 40's degeneracy,
            // which a 20bp identical pair reproduced exactly (2 cigar ops instead
            // of 1). Capping each bound at a third of its own sequence forces at
            // least a third of both to align, so the empty path is unreachable at
            // every length.
            let f1 = FREE_ENDS.min(s1.len() as i32 / 3);
            let f2 = FREE_ENDS.min(s2.len() as i32 / 3);
            sys::wavefront_aligner_set_alignment_free_ends(self.raw, f1, f1, f2, f2);
            sys::wavefront_align(
                self.raw,
                s1.as_ptr() as *const c_char,
                s1.len() as i32,
                s2.as_ptr() as *const c_char,
                s2.len() as i32,
            );
            let cig = (*self.raw).cigar;
            let (b, e) = ((*cig).begin_offset, (*cig).end_offset);
            if e <= b {
                return Vec::new();
            }
            std::slice::from_raw_parts(
                ((*cig).operations as *const u8).add(b as usize),
                (e - b) as usize,
            )
            .to_vec()
        }
    }
}

impl Drop for Aligner {
    fn drop(&mut self) {
        unsafe { sys::wavefront_aligner_delete(self.raw) }
    }
}

thread_local! {
    /// One aligner per scoring scheme actually used. Two, in practice.
    static ALIGNERS: std::cell::RefCell<Vec<Aligner>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Semi-global alignment via WFA2, in parasail's convention and scored by
/// parasail's rules. `None` when the inputs are degenerate or the reconciled
/// columns fail to account for both sequences, so the caller can fall back.
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Option<Alignment> {
    if s1.is_empty() || s2.is_empty() {
        return None;
    }
    let raw = ALIGNERS.with(|a| {
        let mut a = a.borrow_mut();
        if !a.iter().any(|x| x.sc == sc) {
            a.push(Aligner::new(sc));
        }
        let al = a.iter_mut().find(|x| x.sc == sc)?;
        Some(al.columns(s1, s2))
    })?;
    if raw.is_empty() {
        return None;
    }

    // WFA2's `I` consumes the *text* (s2) and `D` the *pattern* (s1); ours is the
    // opposite (`crate::align::ops_to_seq`: `I` takes the query, `D` the
    // reference). Swap while classifying, and re-derive `=`/`X` from the bases
    // rather than trusting `M`, which WFA2 emits for both.
    let mut cols: Vec<u8> = Vec::with_capacity(raw.len());
    let (mut i, mut j) = (0usize, 0usize);
    for &c in &raw {
        match c {
            b'M' | b'X' => {
                let (x, y) = (*s1.get(i)?, *s2.get(j)?);
                cols.push(if x == y { b'=' } else { b'X' });
                i += 1;
                j += 1;
            }
            b'D' => {
                cols.push(b'I');
                i += 1;
            }
            b'I' => {
                cols.push(b'D');
                j += 1;
            }
            _ => return None,
        }
    }
    if i != s1.len() || j != s2.len() {
        return None;
    }

    extend_ends_over_matches(&mut cols, s1, s2);
    let ops = run_length_encode(&cols);

    // Guard against the co-optimal offset. On periodic sequence a shift costs
    // nothing, so WFA2 may align fewer columns for the same penalty and greedy
    // extension cannot recover the full diagonal --- the shifted dangle is all `I`
    // with no `D` to pair against (finding 40). Real reads are not periodic, but
    // AT-rich ONT data comes close enough that a 20bp `ACGTACGT...` fixture
    // reproduces it exactly: `4I16=4D` instead of `20=`.
    //
    // Rather than return a worse alignment, decline and let the caller use
    // parasail. The test is coverage of the shorter sequence: a genuine
    // semi-global alignment of two similar reads covers nearly all of it, while a
    // shift covers `len - offset`.
    let covered: usize = ops
        .iter()
        .filter(|&&(_, t)| t == b'=' || t == b'X')
        .map(|&(n, _)| n)
        .sum();
    let min_len = s1.len().min(s2.len());
    if covered * 10 < min_len * 9 {
        return None;
    }

    let score = score_ops(&ops, sc);
    Some(Alignment { score, cigar: crate::align::encode_cigar(&ops), ops })
}

/// Step 2: pull the aligned core outward through matching base pairs.
///
/// WFA2 is indifferent to terminal matches; parasail is paid for them. The
/// terminal region is a block of `I`s (s1 only) and `D`s (s2 only) in some order;
/// wherever the innermost `I` and innermost `D` refer to equal bases, that pair
/// becomes one aligned column. Repeated until they differ.
fn extend_ends_over_matches(cols: &mut Vec<u8>, s1: &[u8], s2: &[u8]) {
    let aligned = |c: u8| c == b'=' || c == b'X';
    for end in [End::Left, End::Right] {
        loop {
            let Some(core) = (match end {
                End::Left => cols.iter().position(|&c| aligned(c)),
                End::Right => cols.iter().rposition(|&c| aligned(c)),
            }) else {
                return; // nothing aligned at all; not a case to repair here
            };
            // The dangle on this side, and the base each side's next column
            // would consume.
            let dangle: Vec<usize> = match end {
                End::Left => (0..core).collect(),
                End::Right => (core + 1..cols.len()).rev().collect(),
            };
            let (mut last_i, mut last_d) = (None, None);
            for &p in &dangle {
                match cols[p] {
                    b'I' => last_i = Some(p),
                    b'D' => last_d = Some(p),
                    _ => {}
                }
            }
            let (Some(pi), Some(pd)) = (last_i, last_d) else { break };
            // Which bases those columns consume: count the consumption of each
            // sequence up to that column.
            let idx = |upto: usize, op: u8| -> usize {
                cols[..upto]
                    .iter()
                    .filter(|&&c| c == op || aligned(c))
                    .count()
            };
            let (bi, bd) = (idx(pi, b'I'), idx(pd, b'D'));
            let (Some(&x), Some(&y)) = (s1.get(bi), s2.get(bd)) else { break };
            if x != y {
                break;
            }
            // Fuse the pair into one aligned column and drop the other.
            let (keep, drop) = (pi.min(pd), pi.max(pd));
            cols[keep] = b'=';
            cols.remove(drop);
        }
    }
}

enum End {
    Left,
    Right,
}

fn run_length_encode(cols: &[u8]) -> Vec<CigarOp> {
    let mut ops: Vec<CigarOp> = Vec::new();
    for &c in cols {
        match ops.last_mut() {
            Some((n, t)) if *t == c => *n += 1,
            _ => ops.push((1, c)),
        }
    }
    ops
}

/// Score the columns exactly as parasail's semi-global does.
///
/// `match_score` per matched column, `mismatch` per mismatched one, and
/// `open + ext*(L-1)` per gap run --- except for the single outermost run at each
/// end, which is free.
///
/// # Only the *outermost* run at each end is free
///
/// parasail's `sg` zeroes row 0 and column 0, so a path may start at `(i, 0)` or
/// `(0, j)`: it skips a prefix of **one** sequence for nothing, not of both. An
/// alignment whose left dangle is `10I 5D` therefore pays for the `5D`, because
/// only one of the two can be the free start. Treating every outside-the-core gap
/// as free made this layer score *above* the exact DP on 51 of the recorded
/// calls, which is how the error surfaced.
fn score_ops(ops: &[CigarOp], sc: Scoring) -> i32 {
    let aligned = |t: u8| t == b'=' || t == b'X';
    let first = ops.iter().position(|&(_, t)| aligned(t));
    let last = ops.iter().rposition(|&(_, t)| aligned(t));
    let (Some(first), Some(last)) = (first, last) else {
        return 0;
    };
    let mut score = 0i32;
    for (k, &(len, t)) in ops.iter().enumerate() {
        let len = len as i32;
        match t {
            b'=' => score += sc.match_score * len,
            b'X' => score += sc.mismatch * len,
            // The free start (k == 0, before the core) or the free finish
            // (k == last op, after the core). Every other gap is charged.
            _ if (k == 0 && k < first) || (k + 1 == ops.len() && k > last) => {}
            _ => score -= sc.open + sc.ext * (len - 1),
        }
    }
    score
}

pub mod oracle;

#[cfg(test)]
mod tests {
    use super::*;

    /// A pseudorandom sequence. Deliberately NOT periodic: a periodic fixture
    /// makes an offset alignment co-optimal (a 48-base shift of a period-12
    /// string is still all matches), which tests the fixture rather than the
    /// aligner. That case is covered on purpose in
    /// `a_periodic_sequence_admits_a_co_optimal_offset`.
    fn seq(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed | 1;
        (0..len)
            .map(|_| {
                x ^= x >> 33;
                x = x.wrapping_mul(0xff51afd7ed558ccd);
                x ^= x >> 29;
                b"ACGT"[(x >> 40) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn identical_sequences_align_end_to_end_and_score_every_column() {
        // The case bounded ends-free gets wrong on its own: WFA2 shaves terminal
        // matches because skipping them is free. Step 2 must restore them.
        let s = seq(1, 200);
        let a = semiglobal(&s, &s, Scoring::BUBBLE).expect("aligns");
        assert_eq!(a.cigar, "200=");
        assert_eq!(a.score, 400);
    }

    #[test]
    fn a_terminal_overhang_is_free_but_still_appears_in_the_cigar() {
        let s = seq(2, 200);
        let mut tail = s.clone();
        tail.extend_from_slice(&[b'G'; 20]);
        let a = semiglobal(&tail, &s, Scoring::BUBBLE).expect("aligns");
        assert_eq!(a.score, 400, "the 20-base overhang must cost nothing: {}", a.cigar);
        // `ops_to_seq` needs every base accounted for, so the dangle stays.
        assert!(a.cigar.ends_with("20I"), "dangle missing: {}", a.cigar);
    }

    #[test]
    fn an_interior_mismatch_costs_exactly_the_mismatch_penalty() {
        let s = seq(3, 200);
        let mut m = s.clone();
        m[100] = if m[100] == b'A' { b'C' } else { b'A' };
        let a = semiglobal(&s, &m, Scoring::BUBBLE).expect("aligns");
        assert_eq!(a.score, 199 * 2 - 8, "cigar was {}", a.cigar);
    }

    #[test]
    fn an_interior_gap_is_charged_the_affine_cost() {
        let s = seq(4, 200);
        let mut d = s[..100].to_vec();
        d.extend_from_slice(&s[103..]); // a 3-base deletion, interior
        let a = semiglobal(&s, &d, Scoring::BUBBLE).expect("aligns");
        // 197 matches at +2, one length-3 gap at open + ext*(L-1) = 12 + 2.
        assert_eq!(a.score, 197 * 2 - (12 + 2), "cigar was {}", a.cigar);
    }

    #[test]
    fn the_ops_account_for_both_sequences_completely() {
        // `align_to_merge` feeds these to `ops_to_seq`, which returns None unless
        // the ops consume both inputs exactly --- so a lost dangle would silently
        // disable merging rather than fail loudly.
        let s1 = seq(5, 300);
        let mut s2 = s1[7..280].to_vec();
        s2[40] = b'A';
        s2[41] = b'A';
        let a = semiglobal(&s1, &s2, Scoring::MERGE).expect("aligns");
        assert!(
            crate::align::ops_to_seq(&a.ops, &s1, &s2).is_some(),
            "ops {} do not consume both sequences",
            a.cigar
        );
    }

    #[test]
    fn a_periodic_sequence_is_declined_or_correct_never_worse() {
        // The limitation, recorded rather than hidden: when a shift costs nothing
        // because the sequence repeats, WFA2 may align fewer columns for the same
        // penalty, and extending the ends cannot recover the full diagonal --- the
        // left dangle is all `I` with no `D` to pair against. The coverage guard
        // catches the bad ones and declines; whatever comes back must never be
        // worse than the exact score.
        let s: Vec<u8> = (0..200u32).map(|i| b"ACGTTGCAAGCT"[(i as usize * 7) % 12]).collect();
        let exact = crate::parasail::semiglobal(&s, &s, Scoring::BUBBLE);
        assert_eq!(exact.score, 400);
        match semiglobal(&s, &s, Scoring::BUBBLE) {
            None => {}
            Some(a) => assert!(
                a.score <= exact.score,
                "never beat the exact score: {} vs {}",
                a.score,
                exact.score
            ),
        }
    }
}

#[cfg(test)]
mod short_seq_tests {
    use super::*;

    /// Non-periodic, so a shift is not co-optimal.
    fn seq(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed | 1;
        (0..len)
            .map(|_| {
                x ^= x >> 33;
                x = x.wrapping_mul(0xff51afd7ed558ccd);
                x ^= x >> 29;
                b"ACGT"[(x >> 40) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn short_pairs_align_completely_when_not_periodic() {
        // Below 2 * FREE_ENDS the bounded ends-free mode can skip everything, so
        // these are the lengths where the degeneracy would show.
        for n in [12usize, 20, 30, 45, 60, 90] {
            let s = seq(n as u64 + 7, n);
            let a = semiglobal(&s, &s, Scoring::BUBBLE)
                .unwrap_or_else(|| panic!("declined at length {n}"));
            assert_eq!(a.cigar, format!("{n}="), "length {n}");
            assert_eq!(a.score, 2 * n as i32);
        }
    }

    #[test]
    fn a_periodic_pair_is_declined_rather_than_answered_badly() {
        // Period 4, so a 4-base shift is all-matches and costs nothing; WFA2
        // returns `4I16=4D` and greedy extension cannot undo a whole-alignment
        // shift. The coverage guard turns that into a decline, and the caller
        // falls back to parasail --- which is why `simplify`'s smoke test, whose
        // fixture is exactly this, still gets `20=`.
        let s = b"ACGTACGTACGTACGTACGT";
        assert!(
            semiglobal(s, s, Scoring::BUBBLE).is_none(),
            "a co-optimal shift must be declined, not returned"
        );
        // And the production path still answers correctly, via the fallback.
        let exact = crate::parasail::semiglobal(s, s, Scoring::BUBBLE);
        assert_eq!(exact.cigar, "20=");
    }
}
