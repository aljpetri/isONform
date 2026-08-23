//! Semi-global affine alignment matching `parasail.sg_trace_scan_16`.
//!
//! This is the alignment behind the structural-overcorrection guard: the
//! corrected read is aligned back against the original, and [`crate::guard`]
//! reverts any indel run longer than 10 columns on the assumption it is real
//! structure (an exon difference) rather than an error.
//!
//! `correct_read` calls it as
//!
//! ```text
//! parasail_alignment(seq, corr, match_score=4, mismatch_penalty=-8,
//!                    opening_penalty=12, gap_ext=1)
//! ```
//!
//! — note 12, not the function's own default of 24.
//!
//! # What "sg" means, measured
//!
//! parasail's plain `sg` leaves gaps at **both ends of both sequences** free.
//! Verified against the library rather than read off documentation:
//!
//! ```text
//! sg("ACGTACGT", "TTTACGTACGTTTT") -> 3D8=3D, score 32 == 8 matches * 4
//! sg("TTTACGTACGTTTT", "ACGTACGT") -> 3I8=3I, score 32
//! ```
//!
//! Both end gaps contribute nothing to the score, and the CIGAR still spans
//! both sequences in full. A gap of length `L` inside the alignment costs
//! `open + (L - 1) * ext`, so a single-base gap costs `open` — also measured:
//! `sg("ACGTACGT", "ACGACGT")` scores 16, which is the free-end-gap alignment
//! tying exactly with `7 * 4 - 12` for the one-base-deletion alignment.
//!
//! That tie is the whole problem. Free end gaps plus affine costs leave many
//! alignments at the optimal score, and parasail returns one of them. Which one
//! changes where indel runs fall, which changes what the guard reverts, which
//! changes the output — so the choice has to be reproduced, not merely matched
//! on score.
//!
//! # CIGAR convention
//!
//! `I` consumes `s1` (the query), `D` consumes `s2` (the database), matching
//! what `cigar_to_seq(cigar, s1, s2)` then expects — the same convention as
//! edlib's, so [`crate::align::cigar_to_seq`] expands these unchanged.

use crate::align::CigarOp;

/// Scoring, in parasail's sign convention: `match_score` positive, `mismatch`
/// negative, and both gap penalties positive numbers that get subtracted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub open: i32,
    pub ext: i32,
}

impl Scoring {
    /// What `align_bubble_nodes` passes when deciding whether a bubble is
    /// poppable (`SimplifyGraph.py:656`): `match=2, mismatch=-8, open=12, ext=1`.
    ///
    /// Note the mismatch penalty overrides `parasail_alignment`'s own default of
    /// −2. isONcorrect's equivalent was `match=4, mismatch=-8`; the constant it
    /// used (`GUARD`, for its overcorrection guard) is gone, because neither the
    /// guard nor that scoring exists here.
    pub const BUBBLE: Scoring = Scoring {
        match_score: 2,
        mismatch: -8,
        open: 12,
        ext: 1,
    };

    /// What `IsoformGeneration.align_to_merge` passes (`IsoformGeneration.py:379`),
    /// which is also `parasail_alignment`'s own default: `match=2, mismatch=-2,
    /// open=12, ext=1`.
    ///
    /// The two differ only in the mismatch penalty, and that difference is the
    /// whole distinction between "are these two bubble paths the same" and "are
    /// these two isoforms the same" --- worth not conflating.
    pub const MERGE: Scoring = Scoring {
        match_score: 2,
        mismatch: -2,
        open: 12,
        ext: 1,
    };
}

/// Which way a tie is broken. Every field was determined by sweeping against
/// recorded parasail output; see the `sweep` test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TieBreak {
    /// Order to try the three predecessors of an `H` cell.
    pub h_order: [Step; 3],
    /// On a tie inside a gap chain, prefer opening a new gap over extending.
    pub prefer_open: bool,
    /// When several cells on the last row/column share the best score, scan the
    /// last column before the last row.
    pub column_first: bool,
    /// Keep the last equally-good end cell rather than the first.
    pub last_max: bool,
}

/// A predecessor of an `H` cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Step {
    /// Diagonal: a match or mismatch.
    Diagonal,
    /// From the gap-in-`s1` chain, emitting `D`.
    Delete,
    /// From the gap-in-`s2` chain, emitting `I`.
    Insert,
}

impl TieBreak {
    /// The ordering that reproduces `sg_trace_scan_16`, measured rather than
    /// guessed — see the `sweep` test for the numbers.
    ///
    /// `column_first` and `last_max` only matter when several cells on the last
    /// row or column tie for the best score, which never happens on the recorded
    /// corpus: the two sequences are a read and its own correction, so a
    /// full-length alignment always beats ending early. Both values are
    /// therefore **unpinned by evidence**, and the sweep reports all four
    /// combinations as equally perfect. If a corpus ever reaches such a tie,
    /// re-run the sweep before trusting these two.
    pub const PARASAIL: TieBreak = TieBreak {
        h_order: [Step::Diagonal, Step::Insert, Step::Delete],
        prefer_open: false,
        column_first: false,
        last_max: false,
    };

    /// Every combination, for the sweep that picks the real one.
    #[cfg(test)]
    fn all() -> Vec<TieBreak> {
        use Step::*;
        let orders = [
            [Diagonal, Delete, Insert],
            [Diagonal, Insert, Delete],
            [Delete, Diagonal, Insert],
            [Delete, Insert, Diagonal],
            [Insert, Diagonal, Delete],
            [Insert, Delete, Diagonal],
        ];
        let mut out = Vec::new();
        for h_order in orders {
            for prefer_open in [true, false] {
                for column_first in [true, false] {
                    for last_max in [true, false] {
                        out.push(TieBreak {
                            h_order,
                            prefer_open,
                            column_first,
                            last_max,
                        });
                    }
                }
            }
        }
        out
    }
}

/// Result of a semi-global alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub score: i32,
    pub cigar: String,
    pub ops: Vec<CigarOp>,
}

const NEG: i32 = i32::MIN / 4;

/// Align `s1` against `s2`, reproducing `parasail.sg_trace_scan_16`.
///
/// `O(n*m)` time and memory. The guard runs this once per read against a
/// sequence of nearly the same length, so the table is the read length squared —
/// the one place in the port where that is a real cost, and the reason
/// *Deferred improvements* notes banding it.
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    semiglobal_with(s1, s2, sc, TieBreak::PARASAIL)
}

/// The four traceback bits per cell, plus the scores the end-cell scan needs.
///
/// # Why the decisions are stored rather than the matrices
///
/// The obvious implementation keeps all three `i32` matrices and re-derives each
/// traceback step by comparing them. That is 12 bytes per cell — 12 MB for a
/// 1 000 bp read against its own correction, once per read — and it was **65% of
/// the port's runtime** once the edit-distance hot spot was fixed.
///
/// Two rows of `H`/`E`/`F` are enough for the forward pass. What the traceback
/// actually needs is not the scores but the **choices**, and there are only four
/// bits of those per cell:
///
/// * which of the three predecessors won the `H` cell (2 bits), decided by
///   [`TieBreak::h_order`];
/// * whether a `D` chain leaves for `H` here or extends (1 bit);
/// * the same for an `I` chain (1 bit).
///
/// One byte per cell that is `n * m` — 1 MB instead of 12 MB — and it is the
/// *same* computation, so the recorded CIGARs must not move.
///
/// Four bits would halve it again, and was tried: the read-modify-write and
/// shifting cost about 11% of runtime while the DP is already streaming
/// row-by-row, so the cache never benefited. One byte per cell is the better
/// trade — measured, not assumed.
struct Table {
    m: usize,
    packed: Vec<u8>,
    /// `H[i][m]` for every row, for the last-column scan.
    last_col: Vec<i32>,
    /// `H[n][j]` for every column, for the last-row scan.
    last_row: Vec<i32>,
}

impl Table {
    #[inline]
    fn get(&self, i: usize, j: usize) -> u8 {
        self.packed[(i - 1) * self.m + (j - 1)]
    }

    #[inline]
    fn h_choice(&self, i: usize, j: usize) -> Step {
        match self.get(i, j) & 0b11 {
            0 => Step::Diagonal,
            1 => Step::Delete,
            _ => Step::Insert,
        }
    }

    #[inline]
    fn delete_leaves(&self, i: usize, j: usize) -> bool {
        self.get(i, j) & 0b100 != 0
    }

    #[inline]
    fn insert_leaves(&self, i: usize, j: usize) -> bool {
        self.get(i, j) & 0b1000 != 0
    }
}

/// Row bounds of the band, clamped to the matrix.
#[inline]
fn band(i: usize, half: usize, m: usize) -> (usize, usize) {
    (i.saturating_sub(half).max(1), (i + half).min(m))
}

/// Forward pass: fill the decision table row by row, keeping only two rows of
/// scores alive.
///
/// `half` restricts each row to `|i - j| <= half`; `None` computes every cell.
/// Cells outside the band keep `NEG`, so they can never win the end-cell scan.
fn forward(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak, half: Option<usize>) -> Table {
    let (n, m) = (s1.len(), s2.len());

    // Free end gaps: an unaligned prefix of either sequence costs nothing, so
    // row 0 and column 0 start at zero — including the gap chains, which is
    // what keeps a *leading* gap free rather than charging `open`.
    let mut h_prev = vec![0i32; m + 1];
    let mut e_prev = vec![0i32; m + 1];
    let mut f_prev = vec![NEG; m + 1];
    f_prev[0] = 0;

    let mut h_cur = vec![NEG; m + 1];
    let mut e_cur = vec![NEG; m + 1];
    let mut f_cur = vec![NEG; m + 1];

    let mut table = Table {
        m,
        packed: vec![0u8; n * m],
        last_col: vec![NEG; n + 1],
        last_row: vec![NEG; m + 1],
    };
    table.last_col[0] = 0;

    // Mask of achieving predecessors -> the winner under `tb.h_order`. Bit 0 is
    // the diagonal, 1 is delete, 2 is insert. Mask 0 cannot happen, since the
    // maximum is achieved by something.
    let choice_lut: [u8; 8] = std::array::from_fn(|mask| {
        for step in tb.h_order {
            let (bit, code) = match step {
                Step::Diagonal => (0, 0u8),
                Step::Delete => (1, 1),
                Step::Insert => (2, 2),
            };
            if mask >> bit & 1 == 1 {
                return code;
            }
        }
        0
    });

    for i in 1..=n {
        let (lo, hi) = match half {
            Some(h) => band(i, h, m),
            None => (1, m),
        };
        if half.is_some() {
            // Outside the band a cell is unreachable, not zero.
            h_cur.fill(NEG);
            e_cur.fill(NEG);
            f_cur.fill(NEG);
        }
        if lo == 1 {
            h_cur[0] = 0;
            e_cur[0] = NEG;
            f_cur[0] = 0;
        }
        let c1 = s1[i - 1];
        // One row of decisions, sliced out once: the index arithmetic and bounds
        // check would otherwise be repeated for every cell.
        let row = &mut table.packed[(i - 1) * m..i * m];

        for j in lo..=hi {
            let open_d = h_cur[j - 1] - sc.open;
            let ext_d = e_cur[j - 1] - sc.ext;
            e_cur[j] = open_d.max(ext_d);
            let d_leaves = if tb.prefer_open {
                open_d >= ext_d
            } else {
                open_d > ext_d
            };

            let open_i = h_prev[j] - sc.open;
            let ext_i = f_prev[j] - sc.ext;
            f_cur[j] = open_i.max(ext_i);
            let i_leaves = if tb.prefer_open {
                open_i >= ext_i
            } else {
                open_i > ext_i
            };

            let diag = h_prev[j - 1]
                + if c1 == s2[j - 1] {
                    sc.match_score
                } else {
                    sc.mismatch
                };
            let best = diag.max(e_cur[j]).max(f_cur[j]);
            h_cur[j] = best;

            // Which predecessors achieve the cell, as a 3-bit mask resolved
            // through the table above. Walking `tb.h_order` with an early break
            // instead is a data-dependent branch chain over a runtime array,
            // executed once per cell, and the predictor cannot learn it.
            let mask = usize::from(diag == best)
                | usize::from(e_cur[j] == best) << 1
                | usize::from(f_cur[j] == best) << 2;
            row[j - 1] = choice_lut[mask] | (u8::from(d_leaves) << 2) | (u8::from(i_leaves) << 3);
        }

        // The last column is only reachable when the band reaches it.
        table.last_col[i] = if hi == m { h_cur[m] } else { NEG };
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut e_prev, &mut e_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }

    // After the last swap `h_prev` holds row n --- or row 0 when the loop never
    // ran, which is all zeros and correct.
    table.last_row.copy_from_slice(&h_prev);
    table
}

/// [`semiglobal`] with an explicit tie-break, for the sweep.
///
/// # Banding: measured, and rejected
///
/// The full DP is `O(n*m)` and is the largest single cost in the port. The two
/// sequences are a read and its own correction, so the optimal path hugs the
/// diagonal — over 2 408 recorded guard alignments the largest excursion is 63
/// cells, and a half-band of 64 reproduces **every** parasail CIGAR exactly.
/// The `banding` tests below are that measurement.
///
/// It is still not used, for two reasons that only became clear once both were
/// measured:
///
/// * **the cheap band is not provably correct.** There is a real bound — every
///   column is diagonal or a gap, so `2c + g = n + m`, and reaching deviation
///   `W` needs `g >= W`, hence any such path scores at most `2(n+m) - 2W`. But
///   the scores here sit far below the all-match maximum (median slack 753), so
///   the band that *proves* optimality is `~slack/2 ≈ 380`, not 64. And even
///   that proves only the score: a cell near the band edge can have an
///   under-estimated `E`/`F`, which could flip a tie-break and change the CIGAR;
/// * **the provable band is worth almost nothing.** Measured on a 200-read gene
///   cluster: full DP 3.4 s, provable band 3.1 s, unproven ±64 band 2.1 s. The
///   rigorous version buys 9%.
///
/// So the choice was 9% with a proof, or 35% on evidence that would have to be
/// re-established for every new kind of input. For a tool whose entire contract
/// is byte-identical output, and whose failure mode here is silent, neither is a
/// good trade. The remaining win is vectorising this loop, which is exact
/// because it does not change the recurrence.
pub fn semiglobal_with(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> Alignment {
    // Every cell, always: banding was implemented, measured and rejected above.
    traceback(s1, s2, tb, &forward(s1, s2, sc, tb, None))
}

// ---------------------------------------------------------------------------
// Removed on the way across from isONcorrect
// ---------------------------------------------------------------------------
//
// Two blocks of research scaffolding were dropped here, both isONcorrect-only:
//
//  * `global_affine`, `global_check` and `global_check_report` with their
//    counters, which existed to answer whether exact affine *global* alignment
//    could replace this semi-global one. That question was about isONcorrect's
//    structural-overcorrection guard --- the comparison was whether
//    `guard::fix_correction` landed on the same corrected sequence --- and
//    isONform has no guard, so there is nothing to compare.
//  * `band_check` and `band_check_report`, which cross-checked a banded
//    alignment against the exact one on arbitrary real data. Banding was
//    measured and rejected in the isONcorrect port; the exact path is what
//    shipped, and this file keeps it.
//
// Both ran only when an `ISONCORRECT_*` environment variable was set, so
// dropping them changes nothing on any normal run. If banding is ever revisited
// here, that scaffolding is worth copying back rather than rewriting.

impl Table {
    /// The end cell the alignment starts its traceback from.
    fn end_cell(&self, tb: TieBreak) -> (i32, usize, usize) {
        let n = self.last_col.len() - 1;
        let m = self.m;
        let mut best: (i32, usize, usize) = (NEG, n, m);
        let row = (0..=m).map(|j| (self.last_row[j], n, j));
        let col = (0..=n).map(|i| (self.last_col[i], i, m));
        let consider = |cand: (i32, usize, usize), best: &mut (i32, usize, usize)| {
            if cand.0 > best.0 || (tb.last_max && cand.0 == best.0) {
                *best = cand;
            }
        };
        if tb.column_first {
            for c in col {
                consider(c, &mut best);
            }
            for r in row {
                consider(r, &mut best);
            }
        } else {
            for r in row {
                consider(r, &mut best);
            }
            for c in col {
                consider(c, &mut best);
            }
        }
        best
    }
}

/// Walk the stored decisions back from the end cell.
fn traceback(s1: &[u8], s2: &[u8], tb: TieBreak, table: &Table) -> Alignment {
    let (score, i, j) = table.end_cell(tb);
    let mut aln = traceback_from(s1, s2, table, i, j);
    aln.score = score;
    aln
}

/// The traceback walk itself, from an explicit end cell.
fn traceback_from(s1: &[u8], s2: &[u8], table: &Table, i0: usize, j0: usize) -> Alignment {
    let (n, m) = (s1.len(), s2.len());
    let (mut i, mut j) = (i0, j0);

    let mut ops: Vec<u8> = Vec::with_capacity(n + m);
    // Free trailing gaps first: whatever is left of either sequence.
    ops.extend(std::iter::repeat_n(b'D', m - j));
    ops.extend(std::iter::repeat_n(b'I', n - i));

    // `Some(Step::Delete)` means we are inside a D chain and must decide
    // open-vs-extend before leaving it; `None` means we are on H.
    let mut chain: Option<Step> = None;
    while i > 0 && j > 0 {
        match chain {
            Some(Step::Delete) => {
                let leaves = table.delete_leaves(i, j);
                ops.push(b'D');
                j -= 1;
                chain = if leaves { None } else { Some(Step::Delete) };
            }
            Some(Step::Insert) => {
                let leaves = table.insert_leaves(i, j);
                ops.push(b'I');
                i -= 1;
                chain = if leaves { None } else { Some(Step::Insert) };
            }
            Some(Step::Diagonal) => unreachable!("the diagonal is not a chain"),
            None => match table.h_choice(i, j) {
                Step::Diagonal => {
                    ops.push(if s1[i - 1] == s2[j - 1] { b'=' } else { b'X' });
                    i -= 1;
                    j -= 1;
                }
                step => chain = Some(step),
            },
        }
    }
    // Free leading gaps.
    ops.extend(std::iter::repeat_n(b'D', j));
    ops.extend(std::iter::repeat_n(b'I', i));

    ops.reverse();
    let (cigar, runs) = encode(&ops);
    // The caller fills in the score; the walk does not know it.
    Alignment {
        score: 0,
        cigar,
        ops: runs,
    }
}

fn encode(ops: &[u8]) -> (String, Vec<CigarOp>) {
    let mut out = String::new();
    let mut runs = Vec::new();
    let mut k = 0;
    while k < ops.len() {
        let c = ops[k];
        let start = k;
        while k < ops.len() && ops[k] == c {
            k += 1;
        }
        out.push_str(&(k - start).to_string());
        out.push(c as char);
        runs.push((k - start, c));
    }
    (out, runs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The scoring these semantics tests were originally verified against.
    ///
    /// `match=4, mismatch=-8, open=12, ext=1` --- isONcorrect's, not isONform's.
    /// Kept deliberately: the expected scores below were read off the parasail
    /// *library* at this scoring, so they are measurements, and re-deriving them
    /// arithmetically under `Scoring::BUBBLE` would turn evidence into
    /// assumption. What they establish --- free end gaps on both sequences, a
    /// one-base gap costing `open`, longer gaps costing `open + (L-1) * ext` ---
    /// is a property of parasail's semi-global mode and holds at any scoring.
    ///
    /// isONform's own scorings are [`Scoring::BUBBLE`] and [`Scoring::MERGE`];
    /// the tie-break and CIGAR tests below use those.
    const VERIFIED: Scoring = Scoring {
        match_score: 4,
        mismatch: -8,
        open: 12,
        ext: 1,
    };

    fn sg(a: &str, b: &str) -> Alignment {
        semiglobal(a.as_bytes(), b.as_bytes(), VERIFIED)
    }

    #[test]
    fn isonforms_two_scorings_differ_only_in_the_mismatch_penalty() {
        // The distinction is load-bearing: -8 for "are these bubble paths the
        // same", -2 for "are these isoforms the same". Conflating them would
        // silently change both decisions.
        assert_eq!(Scoring::BUBBLE.match_score, Scoring::MERGE.match_score);
        assert_eq!(Scoring::BUBBLE.open, Scoring::MERGE.open);
        assert_eq!(Scoring::BUBBLE.ext, Scoring::MERGE.ext);
        assert_eq!(Scoring::BUBBLE.mismatch, -8);
        assert_eq!(Scoring::MERGE.mismatch, -2);
    }

    #[test]
    fn the_mismatch_penalty_changes_where_a_gap_beats_a_mismatch() {
        // Same input, the two scorings, different alignments --- which is why they
        // are separate constants rather than one shared default.
        let a = b"ACGTACGT";
        let b = b"ACGAACGT";
        let strict = semiglobal(a, b, Scoring::BUBBLE);
        let loose = semiglobal(a, b, Scoring::MERGE);
        assert_ne!(
            strict.score, loose.score,
            "a mismatch is priced differently by the two"
        );
    }

    /// Ground truth read off the parasail library itself, not assumed.
    #[test]
    fn end_gaps_are_free_on_both_sequences() {
        let a = sg("ACGTACGT", "TTTACGTACGTTTT");
        assert_eq!(a.cigar, "3D8=3D");
        assert_eq!(a.score, 32);

        let b = sg("TTTACGTACGTTTT", "ACGTACGT");
        assert_eq!(b.cigar, "3I8=3I");
        assert_eq!(b.score, 32);
    }

    #[test]
    fn identical_sequences_score_every_match() {
        let a = sg("ACGTACGT", "ACGTACGT");
        assert_eq!(a.cigar, "8=");
        assert_eq!(a.score, 32);
    }

    #[test]
    fn a_single_base_gap_costs_the_opening_penalty() {
        // Both sides long enough that end gaps cannot win.
        let s1 = "ACGTACGTACGTTGCAACGT";
        let s2 = "ACGTACGTACGTGCAACGT"; // one T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 19 * 4 - 12);
        assert_eq!(a.cigar.matches('I').count(), 1);
    }

    #[test]
    fn a_longer_gap_costs_open_plus_extensions() {
        let s1 = "ACGTACGTACGTTTTTGCAACGTACGT";
        let s2 = "ACGTACGTACGTGCAACGTACGT"; // four T dropped
        let a = sg(s1, s2);
        assert_eq!(a.score, 23 * 4 - (12 + 3));
    }

    #[test]
    fn the_cigar_spans_both_sequences() {
        for (s1, s2) in [
            ("ACGTACGT", "TTTACGTACGTTTT"),
            ("ACGTACGTACGTTGCAACGT", "ACGTACGTACGTGCAACGT"),
            ("ACGT", "ACGT"),
            ("", "ACGT"),
            ("ACGT", ""),
        ] {
            let a = semiglobal(s1.as_bytes(), s2.as_bytes(), Scoring::BUBBLE);
            let (mut q, mut r) = (0usize, 0usize);
            for (len, op) in &a.ops {
                match op {
                    b'=' | b'X' => {
                        q += len;
                        r += len;
                    }
                    b'I' => q += len,
                    b'D' => r += len,
                    _ => panic!("unexpected op"),
                }
            }
            assert_eq!((q, r), (s1.len(), s2.len()), "on {s1:?} vs {s2:?}");
        }
    }

    #[test]
    fn an_empty_sequence_is_all_free_gaps() {
        let a = sg("", "ACGT");
        assert_eq!((a.cigar.as_str(), a.score), ("4D", 0));
        let b = sg("ACGT", "");
        assert_eq!((b.cigar.as_str(), b.score), ("4I", 0));
        let c = sg("", "");
        assert_eq!((c.cigar.as_str(), c.score), ("", 0));
    }
}

/// Differential oracle against the parasail library.
///
/// ```bash
/// bench/dump_reference.py --fastq cluster.fastq --outdir /tmp/d
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --manifest-path rust/Cargo.toml parasail::oracle -- --nocapture
/// ```
///
/// `sweep` is the tool that *chose* [`TieBreak::PARASAIL`]: it runs every
/// combination against the recorded cases and prints the score. Re-run it if a
/// mismatch ever appears — a new corpus reaching a tie the current ordering gets
/// wrong would show up as a different winner.
#[cfg(test)]
mod oracle {
    use super::*;

    struct Case {
        s1: Vec<u8>,
        s2: Vec<u8>,
        cigar: String,
        score: i32,
        scoring: Scoring,
    }

    fn load() -> Option<Vec<Case>> {
        let path = std::env::var("PARASAIL_CASES").ok()?;
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("PARASAIL_CASES={path} unreadable: {e}"));
        let mut cases = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            assert!(f.len() >= 8, "short row");
            cases.push(Case {
                s1: f[0].as_bytes().to_vec(),
                s2: f[1].as_bytes().to_vec(),
                cigar: f[2].to_string(),
                score: f[3].parse().unwrap(),
                scoring: Scoring {
                    match_score: f[4].parse().unwrap(),
                    mismatch: f[5].parse().unwrap(),
                    open: f[6].parse().unwrap(),
                    ext: f[7].parse().unwrap(),
                },
            });
        }
        Some(cases)
    }

    #[test]
    fn matches_parasail_exactly() {
        let Some(cases) = load() else { return };
        assert!(!cases.is_empty(), "no cases");

        let (mut bad_score, mut bad_cigar, mut shown) = (0usize, 0usize, 0usize);
        for c in &cases {
            let got = semiglobal(&c.s1, &c.s2, c.scoring);
            if got.score != c.score {
                bad_score += 1;
            }
            if got.cigar != c.cigar {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "cigar mismatch on {}x{}:\n  parasail: {}\n  rust    : {}",
                        c.s1.len(),
                        c.s2.len(),
                        c.cigar,
                        got.cigar
                    );
                }
                bad_cigar += 1;
            }
        }
        eprintln!(
            "parasail oracle: {} cases, {bad_score} score mismatches, {bad_cigar} cigar mismatches",
            cases.len()
        );
        assert_eq!(bad_score, 0, "scores differed");
        assert_eq!(bad_cigar, 0, "CIGARs differed");
    }

    /// Not an assertion --- a measurement. Prints how each tie-break ordering
    /// scores so the right one can be picked from evidence.
    #[test]
    fn sweep() {
        let Some(mut cases) = load() else { return };
        let Ok(setting) = std::env::var("PARASAIL_SWEEP") else {
            return;
        };
        // 96 configurations times full-length reads is hours of quadratic DP, so
        // the variable doubles as a case cap: PARASAIL_SWEEP=200 sweeps the
        // first 200. Ties are what matter here, not volume.
        if let Ok(cap) = setting.parse::<usize>() {
            if cap > 1 {
                cases.truncate(cap);
            }
        }
        let mut results: Vec<(usize, TieBreak)> = TieBreak::all()
            .into_iter()
            .map(|tb| {
                let hits = cases
                    .iter()
                    .filter(|c| semiglobal_with(&c.s1, &c.s2, c.scoring, tb).cigar == c.cigar)
                    .count();
                (hits, tb)
            })
            .collect();
        results.sort_by_key(|(hits, _)| std::cmp::Reverse(*hits));
        eprintln!("tie-break sweep over {} cases:", cases.len());
        for (hits, tb) in results.iter().take(8) {
            eprintln!(
                "  {hits:>6} / {}  h_order={:?} open={} col_first={} last_max={}",
                cases.len(),
                tb.h_order,
                tb.prefer_open,
                tb.column_first,
                tb.last_max
            );
        }
    }
}

/// Does banding change the answer? A measurement, not a feature.
///
/// The guard is ~65% of runtime and its DP is `O(n*m)` over a read against its
/// own correction. Those two sequences differ very little, so a band around the
/// diagonal should lose nothing --- but `fix_correction` walks the *path*, so an
/// equally-optimal-but-different path changes the output. This computes the
/// banded answer and compares CIGARs against what parasail recorded.
///
/// ```bash
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --release --manifest-path rust/Cargo.toml parasail::banding -- --nocapture
/// ```
#[cfg(test)]
mod banding {
    use super::*;

    /// Same recurrence, restricted to `|i - j - drift| <= half`, where `drift`
    /// centres the band on the diagonal the length difference implies.
    fn banded(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak, half: usize) -> Option<String> {
        let (n, m) = (s1.len(), s2.len());
        if n == 0 || m == 0 {
            return Some(semiglobal_with(s1, s2, sc, tb).cigar);
        }
        let lo_of = |i: usize| -> usize { i.saturating_sub(half).max(1) };
        let hi_of = |i: usize| -> usize { (i + half).min(m) };

        let mut h_prev = vec![0i32; m + 1];
        let mut e_prev = vec![0i32; m + 1];
        let mut f_prev = vec![NEG; m + 1];
        f_prev[0] = 0;
        let mut h_cur = vec![NEG; m + 1];
        let mut e_cur = vec![NEG; m + 1];
        let mut f_cur = vec![NEG; m + 1];

        let mut packed = vec![0u8; n * m];
        let mut last_col = vec![NEG; n + 1];
        last_col[0] = 0;
        let mut touched_edge = false;

        for i in 1..=n {
            let (lo, hi) = (lo_of(i), hi_of(i));
            // Cells outside the band are unreachable, not zero.
            for v in h_cur.iter_mut() {
                *v = NEG;
            }
            for v in e_cur.iter_mut() {
                *v = NEG;
            }
            for v in f_cur.iter_mut() {
                *v = NEG;
            }
            if lo == 1 {
                h_cur[0] = 0;
                f_cur[0] = 0;
            }
            let c1 = s1[i - 1];
            for j in lo..=hi {
                let open_d = h_cur[j - 1] - sc.open;
                let ext_d = e_cur[j - 1] - sc.ext;
                e_cur[j] = open_d.max(ext_d);
                let d_leaves = if tb.prefer_open {
                    open_d >= ext_d
                } else {
                    open_d > ext_d
                };
                let open_i = h_prev[j] - sc.open;
                let ext_i = f_prev[j] - sc.ext;
                f_cur[j] = open_i.max(ext_i);
                let i_leaves = if tb.prefer_open {
                    open_i >= ext_i
                } else {
                    open_i > ext_i
                };
                let diag = h_prev[j - 1]
                    + if c1 == s2[j - 1] {
                        sc.match_score
                    } else {
                        sc.mismatch
                    };
                let best = diag.max(e_cur[j]).max(f_cur[j]);
                h_cur[j] = best;
                let mut choice = 0u8;
                for step in tb.h_order {
                    let (hit, code) = match step {
                        Step::Diagonal => (diag == best, 0u8),
                        Step::Delete => (e_cur[j] == best, 1),
                        Step::Insert => (f_cur[j] == best, 2),
                    };
                    if hit {
                        choice = code;
                        break;
                    }
                }
                packed[(i - 1) * m + (j - 1)] =
                    choice | (u8::from(d_leaves) << 2) | (u8::from(i_leaves) << 3);
            }
            last_col[i] = if hi == m { h_cur[m] } else { NEG };
            std::mem::swap(&mut h_prev, &mut h_cur);
            std::mem::swap(&mut e_prev, &mut e_cur);
            std::mem::swap(&mut f_prev, &mut f_cur);
        }

        let mut best = (NEG, n, m);
        for (j, &v) in h_prev.iter().enumerate() {
            if v > best.0 {
                best = (v, n, j);
            }
        }
        for (i, &v) in last_col.iter().enumerate() {
            if v > best.0 {
                best = (v, i, m);
            }
        }
        let (_score, mut i, mut j) = best;

        let mut ops: Vec<u8> = Vec::with_capacity(n + m);
        ops.extend(std::iter::repeat_n(b'D', m - j));
        ops.extend(std::iter::repeat_n(b'I', n - i));
        let mut chain: Option<Step> = None;
        while i > 0 && j > 0 {
            // A *real* band boundary, not the sequence boundary: `lo`/`hi` are
            // clamped to [1, m], and hitting those is not clipping.
            let (lo, hi) = (lo_of(i), hi_of(i));
            if (j == lo && lo > 1) || (j == hi && hi < m) {
                touched_edge = true;
            }
            let bits = packed[(i - 1) * m + (j - 1)];
            match chain {
                Some(Step::Delete) => {
                    ops.push(b'D');
                    j -= 1;
                    chain = if bits & 0b100 != 0 {
                        None
                    } else {
                        Some(Step::Delete)
                    };
                }
                Some(Step::Insert) => {
                    ops.push(b'I');
                    i -= 1;
                    chain = if bits & 0b1000 != 0 {
                        None
                    } else {
                        Some(Step::Insert)
                    };
                }
                Some(Step::Diagonal) => unreachable!(),
                None => match bits & 0b11 {
                    0 => {
                        ops.push(if s1[i - 1] == s2[j - 1] { b'=' } else { b'X' });
                        i -= 1;
                        j -= 1;
                    }
                    1 => chain = Some(Step::Delete),
                    _ => chain = Some(Step::Insert),
                },
            }
        }
        ops.extend(std::iter::repeat_n(b'D', j));
        ops.extend(std::iter::repeat_n(b'I', i));
        ops.reverse();
        if touched_edge {
            return None; // the band clipped the path; a wider one is needed
        }
        Some(encode(&ops).0)
    }

    /// Would "same score at W and at 2W" have caught the cases a too-narrow
    /// band gets wrong? If so it is a sound-in-practice acceptance test: a
    /// better path within 2W is found by the wider run and the scores diverge.
    #[test]
    fn does_widening_detect_a_band_that_is_too_narrow() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let (mut wrong_caught, mut wrong_missed, mut ok) = (0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let (s1, s2) = (f[0].as_bytes(), f[1].as_bytes());
            let sc = Scoring {
                match_score: f[4].parse().unwrap(),
                mismatch: f[5].parse().unwrap(),
                open: f[6].parse().unwrap(),
                ext: f[7].parse().unwrap(),
            };
            // Deliberately too narrow, so wrong answers actually occur.
            let narrow = banded_score(s1, s2, sc, TieBreak::PARASAIL, 32);
            let wider = banded_score(s1, s2, sc, TieBreak::PARASAIL, 64);
            let truth: i32 = f[3].parse().unwrap();
            match (narrow, wider) {
                (Some(a), Some(b)) => {
                    if a == truth {
                        ok += 1;
                    } else if a != b {
                        wrong_caught += 1;
                    } else {
                        wrong_missed += 1;
                    }
                }
                _ => ok += 1,
            }
        }
        eprintln!(
            "narrow band correct {ok}, wrong-and-caught-by-widening {wrong_caught}, \
             wrong-and-missed {wrong_missed}"
        );
    }

    /// Banded score only, without the traceback or the edge check.
    fn banded_score(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak, half: usize) -> Option<i32> {
        let (n, m) = (s1.len(), s2.len());
        if n == 0 || m == 0 {
            return None;
        }
        let mut h_prev = vec![0i32; m + 1];
        let mut e_prev = vec![0i32; m + 1];
        let mut f_prev = vec![NEG; m + 1];
        f_prev[0] = 0;
        let mut h_cur = vec![NEG; m + 1];
        let mut e_cur = vec![NEG; m + 1];
        let mut f_cur = vec![NEG; m + 1];
        let mut best = NEG;
        for i in 1..=n {
            let lo = i.saturating_sub(half).max(1);
            let hi = (i + half).min(m);
            h_cur.fill(NEG);
            e_cur.fill(NEG);
            f_cur.fill(NEG);
            if lo == 1 {
                h_cur[0] = 0;
                f_cur[0] = 0;
            }
            let c1 = s1[i - 1];
            for j in lo..=hi {
                e_cur[j] = (h_cur[j - 1] - sc.open).max(e_cur[j - 1] - sc.ext);
                f_cur[j] = (h_prev[j] - sc.open).max(f_prev[j] - sc.ext);
                let diag = h_prev[j - 1]
                    + if c1 == s2[j - 1] {
                        sc.match_score
                    } else {
                        sc.mismatch
                    };
                h_cur[j] = diag.max(e_cur[j]).max(f_cur[j]);
            }
            if hi == m {
                best = best.max(h_cur[m]);
            }
            std::mem::swap(&mut h_prev, &mut h_cur);
            std::mem::swap(&mut e_prev, &mut e_cur);
            std::mem::swap(&mut f_prev, &mut f_cur);
            let _ = tb;
        }
        for &v in h_prev.iter() {
            best = best.max(v);
        }
        Some(best)
    }

    #[test]
    fn how_often_does_banding_change_the_cigar() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let mut cases = Vec::new();
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            cases.push((
                f[0].as_bytes().to_vec(),
                f[1].as_bytes().to_vec(),
                f[2].to_string(),
                Scoring {
                    match_score: f[4].parse().unwrap(),
                    mismatch: f[5].parse().unwrap(),
                    open: f[6].parse().unwrap(),
                    ext: f[7].parse().unwrap(),
                },
            ));
        }
        eprintln!("{} recorded alignments", cases.len());
        for half in [32usize, 64, 128, 256, 512] {
            let (mut same, mut differ, mut clipped) = (0usize, 0usize, 0usize);
            for (s1, s2, want, sc) in &cases {
                match banded(s1, s2, *sc, TieBreak::PARASAIL, half) {
                    None => clipped += 1,
                    Some(got) if &got == want => same += 1,
                    Some(_) => differ += 1,
                }
            }
            eprintln!(
                "  half-band {half:>4}: identical {same:>5}  DIFFERENT {differ:>5}  \
                 clipped(fallback) {clipped:>5}"
            );
        }
    }
}

/// Does the guard actually need *semi*-global alignment?
///
/// `correct_read` aligns `seq` against `corr`, and `corr` is built by splicing
/// corrected spans into `seq` — so the two share an **exact prefix and suffix**.
/// Matching those scores `+match` per base while a free end gap scores 0, so
/// free ends should never win. If exact affine *global* alignment agrees with
/// parasail's `sg` on real data, the semi-globalness is an artifact of reusing
/// `parasail_alignment`'s default, and any global aligner is a candidate.
///
/// ```bash
/// PARASAIL_CASES=/tmp/d/parasail_cases.tsv \
///   cargo test --release --manifest-path rust/Cargo.toml parasail::global_vs_sg -- --nocapture
/// ```
#[cfg(test)]
mod global_vs_sg {
    use super::*;

    /// Exact affine global (Needleman-Wunsch-Gotoh), same tie-break as
    /// [`TieBreak::PARASAIL`]. Self-contained so the production path is
    /// untouched by this measurement.
    fn global(s1: &[u8], s2: &[u8], sc: Scoring, tb: TieBreak) -> String {
        let (n, m) = (s1.len(), s2.len());
        let at = |i: usize, j: usize| i * (m + 1) + j;
        let mut h = vec![NEG; (n + 1) * (m + 1)];
        let mut e = vec![NEG; (n + 1) * (m + 1)];
        let mut f = vec![NEG; (n + 1) * (m + 1)];
        let gap = |l: usize| -> i32 {
            if l == 0 {
                0
            } else {
                -(sc.open + (l as i32 - 1) * sc.ext)
            }
        };
        h[at(0, 0)] = 0;
        for i in 1..=n {
            h[at(i, 0)] = gap(i);
            f[at(i, 0)] = gap(i);
        }
        for j in 1..=m {
            h[at(0, j)] = gap(j);
            e[at(0, j)] = gap(j);
        }
        for i in 1..=n {
            for j in 1..=m {
                e[at(i, j)] = (h[at(i, j - 1)] - sc.open).max(e[at(i, j - 1)] - sc.ext);
                f[at(i, j)] = (h[at(i - 1, j)] - sc.open).max(f[at(i - 1, j)] - sc.ext);
                let diag = h[at(i - 1, j - 1)]
                    + if s1[i - 1] == s2[j - 1] {
                        sc.match_score
                    } else {
                        sc.mismatch
                    };
                h[at(i, j)] = diag.max(e[at(i, j)]).max(f[at(i, j)]);
            }
        }

        let (mut i, mut j) = (n, m);
        let mut ops: Vec<u8> = Vec::new();
        let mut chain: Option<Step> = None;
        while i > 0 && j > 0 {
            match chain {
                Some(Step::Delete) => {
                    let open_d = h[at(i, j - 1)] - sc.open;
                    let ext_d = e[at(i, j - 1)] - sc.ext;
                    let leaves = if tb.prefer_open {
                        open_d >= ext_d
                    } else {
                        open_d > ext_d
                    };
                    ops.push(b'D');
                    j -= 1;
                    chain = if leaves { None } else { Some(Step::Delete) };
                }
                Some(Step::Insert) => {
                    let open_i = h[at(i - 1, j)] - sc.open;
                    let ext_i = f[at(i - 1, j)] - sc.ext;
                    let leaves = if tb.prefer_open {
                        open_i >= ext_i
                    } else {
                        open_i > ext_i
                    };
                    ops.push(b'I');
                    i -= 1;
                    chain = if leaves { None } else { Some(Step::Insert) };
                }
                Some(Step::Diagonal) => unreachable!(),
                None => {
                    let cur = h[at(i, j)];
                    let same = s1[i - 1] == s2[j - 1];
                    let diag =
                        h[at(i - 1, j - 1)] + if same { sc.match_score } else { sc.mismatch };
                    let mut taken = Step::Diagonal;
                    for step in tb.h_order {
                        let hit = match step {
                            Step::Diagonal => diag == cur,
                            Step::Delete => e[at(i, j)] == cur,
                            Step::Insert => f[at(i, j)] == cur,
                        };
                        if hit {
                            taken = step;
                            break;
                        }
                    }
                    match taken {
                        Step::Diagonal => {
                            ops.push(if same { b'=' } else { b'X' });
                            i -= 1;
                            j -= 1;
                        }
                        step => chain = Some(step),
                    }
                }
            }
        }
        ops.extend(std::iter::repeat_n(b'D', j));
        ops.extend(std::iter::repeat_n(b'I', i));
        ops.reverse();
        encode(&ops).0
    }

    #[test]
    fn global_agrees_with_parasail_semiglobal() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            return;
        };
        let data = std::fs::read_to_string(&path).expect("cases");
        let (mut n, mut same, mut differ, mut shown) = (0usize, 0usize, 0usize, 0usize);
        for line in data
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let sc = Scoring {
                match_score: f[4].parse().unwrap(),
                mismatch: f[5].parse().unwrap(),
                open: f[6].parse().unwrap(),
                ext: f[7].parse().unwrap(),
            };
            let got = global(f[0].as_bytes(), f[1].as_bytes(), sc, TieBreak::PARASAIL);
            if got == f[2] {
                same += 1;
            } else {
                if shown < 3 {
                    shown += 1;
                    eprintln!(
                        "differs on {}x{}:\n  sg    : {}\n  global: {}",
                        f[0].len(),
                        f[1].len(),
                        f[2],
                        got
                    );
                }
                differ += 1;
            }
            n += 1;
        }
        eprintln!("global vs parasail sg: {same} identical, {differ} different, of {n}");
    }
}
// Removed on the way across from isONcorrect: the whole `blockalign` test module.
//
// It measured how close block-aligner's SIMD banded alignment came to this exact
// path, and whether the difference survived `guard::fix_correction`. Both halves
// are isONcorrect-specific: block-aligner was on its default path for the
// structural-overcorrection guard, and isONform has neither the guard nor the
// dependency. The exact semi-global path below is what isONform needs, and its
// own oracle test against recorded `parasail.sg_trace_scan_16` output is what
// verifies it.

// Removed on the way across from isONcorrect: the `one_case` probe module, a
// one-off that read two sequences from files and compared semi-global against
// `global_affine`. Both the comparison target and the guard question behind it
// are gone.
