//! The SIMD aligner: block-aligner behind a parasail-shaped call.
//!
//! # Why
//!
//! Profiling says alignment is the dominant cost of the whole tool, in *two*
//! stages, and the graph is not the problem:
//!
//! | stage | 1 000-read batch | 250-read batch, 7 343-node graph |
//! | --- | --- | --- |
//! | intervals | 2.33s | 0.25s |
//! | graph build | **0.35s** | **0.03s** |
//! | simplify (`pop_bubbles`) | 2.29s | **46.04s** |
//! | merge (`merge_consensuses`) | **58.05s** | 20.56s |
//!
//! Which of simplify and merge dominates depends on the graph, but sampling both
//! lands in `parasail::semiglobal_with` and `parasail::traceback`. So one aligner
//! swap addresses the whole pipeline.
//!
//! isONcorrect already made this swap and measured **8-99x**, on 1 400 verified
//! alignments. This port has 54 884 recorded real-parasail cases to check against,
//! which is 39x the evidence.
//!
//! # Free end gaps: what parasail does, and what block-aligner gives
//!
//! [`crate::parasail::semiglobal`] zeroes row 0 **and** column 0, and its end-cell
//! scan looks at the last row **and** the last column. Both sequences may dangle
//! at both ends --- an overlap alignment. That is what makes `align_to_merge`
//! work on isoforms whose lengths differ.
//!
//! block-aligner 0.5.1 parameterises `Block` as
//! `Block<TRACE, X_DROP, LOCAL_START, FREE_QUERY_START_GAPS, FREE_QUERY_END_GAPS>`
//! --- free gaps on the **query** only, not symmetric.
//!
//! With the query as the *shorter* sequence, freeing the query's start and end
//! gaps lets the shorter sequence sit anywhere inside the longer one, which is the
//! case `align_to_merge` cares about: it always compares a long consensus against
//! a short one and asks whether the short is contained. The asymmetric form is
//! therefore the right shape, *provided* the caller keeps the orientation. This
//! module enforces it rather than trusting callers --- see [`semiglobal_ops`].
//!
//! What it is **not** is bit-identical to parasail. Two reasons, and both are
//! measured rather than assumed by [`oracle`]:
//!
//! * the band is adaptive, so a score can in principle be sub-optimal;
//! * even at an equal score, a different equally-optimal path may be reported,
//!   which changes the CIGAR and therefore the diversity computation that decides
//!   a merge.
//!
//! Reference identity is no longer the goal, so a different-but-equally-optimal
//! path is acceptable. A *worse* score is not, and that is what the oracle gates.

use crate::align::CigarOp;
use crate::parasail::Scoring;

/// Upper bound on the block size, and therefore on the query length this can
/// handle. Block sizes must be powers of two and below 2^16.
///
/// # Free end gaps forbid banding
///
/// block-aligner asserts `min_size > query.len()` whenever
/// `FREE_QUERY_END_GAPS` is set, so the block must span the whole query and there
/// is no adaptive band. That is a real cost: isONcorrect's 8-99x came from SIMD
/// *and* banding on short segments, and only the SIMD half is available here.
/// Measured rather than assumed --- see `oracle` and the benchmark.
const MAX_BLOCK: usize = 32768;
/// Starting capacity and growth granularity, so a run of slightly longer
/// sequences does not reallocate every call.
const CAP_STEP: usize = 512;

/// The block size for a given query length: the smallest power of two strictly
/// greater than it, which is what the `FREE_QUERY_END_GAPS` assert demands.
fn block_size_for(query_len: usize) -> usize {
    (query_len + 1).next_power_of_two().clamp(32, MAX_BLOCK)
}

/// Semi-global alignment ops for a pair of sequences, SIMD.
///
/// **Orientation is enforced, not requested.** The shorter sequence becomes the
/// query, because block-aligner's free-gap flags apply to the query only. When the
/// arguments arrive the other way round the ops are transposed on the way out, so
/// callers see ops describing `(s1, s2)` in the order they passed them --- an
/// insertion stays an insertion relative to *their* first argument.
///
/// Returns `false` when the pair cannot be aligned this way (empty input, a score
/// that will not fit `i8`, or a truncated result), and the caller should fall back
/// to the exact path. Silently returning a wrong alignment is the one thing that
/// must not happen, since a merge decision hangs on it.
pub fn semiglobal_ops(s1: &[u8], s2: &[u8], sc: Scoring, out: &mut Vec<CigarOp>) -> Option<i32> {
    out.clear();
    if s1.is_empty() || s2.is_empty() {
        return None;
    }
    let swapped = s1.len() < s2.len();
    let (query, reference) = if swapped { (s1, s2) } else { (s2, s1) };
    let score = SCRATCH.with_borrow_mut(|s| s.align_into(query, reference, sc, out))?;
    // `out` describes (query, reference) = (shorter, longer). The caller passed
    // (s1, s2). Transpose when those disagree: an insertion in one orientation is
    // a deletion in the other.
    let caller_query_is_shorter = s1.len() <= s2.len();
    let produced_query_is_s1 = swapped;
    if produced_query_is_s1 != caller_query_is_shorter || !produced_query_is_s1 {
        for op in out.iter_mut() {
            op.1 = match op.1 {
                b'I' => b'D',
                b'D' => b'I',
                other => other,
            };
        }
    }
    Some(score)
}

thread_local! {
    /// Reusable allocations. **Not a micro-optimisation:** the exact path
    /// allocates and zeroes an `n * m` byte table per call --- ~1.9 MB at the
    /// measured mean consensus length of 1 364 bp --- and there are thousands of
    /// calls per batch. isONcorrect measured per-call setup at 18.84s against
    /// 11.54s reused on the same work.
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::new(Scratch::new());
}

struct Scratch {
    /// `<TRACE, X_DROP, LOCAL_START, FREE_QUERY_START_GAPS, FREE_QUERY_END_GAPS>`
    /// --- trace on (the CIGAR is needed), and the query's end gaps free, which is
    /// the semi-global shape parasail's `sg` provides. See the module docs.
    block: block_aligner::scan_block::Block<true, false, false, true, true>,
    q: block_aligner::scan_block::PaddedBytes,
    r: block_aligner::scan_block::PaddedBytes,
    cigar: block_aligner::cigar::Cigar,
    cap: usize,
    /// The block size the current allocations were built for.
    blk: usize,
}

impl Scratch {
    fn new() -> Self {
        use block_aligner::{scan_block::*, scores::*};
        let blk = block_size_for(CAP_STEP);
        Self {
            block: Block::new(CAP_STEP, CAP_STEP, blk),
            q: PaddedBytes::new::<NucMatrix>(CAP_STEP, blk),
            r: PaddedBytes::new::<NucMatrix>(CAP_STEP, blk),
            cigar: block_aligner::cigar::Cigar::new(CAP_STEP, CAP_STEP),
            cap: CAP_STEP,
            blk,
        }
    }

    /// `Block::new` takes *maximum* lengths and `align` the actual ones, so one
    /// block sized to the largest pair seen serves every smaller one. Grows only.
    fn reserve(&mut self, len: usize, query_len: usize) -> usize {
        use block_aligner::{scan_block::*, scores::*};
        let want_blk = block_size_for(query_len);
        if len <= self.cap && want_blk <= self.blk {
            return self.blk;
        }
        self.cap = len.max(self.cap).next_multiple_of(CAP_STEP);
        self.blk = want_blk.max(self.blk);
        self.block = Block::new(self.cap, self.cap, self.blk);
        self.q = PaddedBytes::new::<NucMatrix>(self.cap, self.blk);
        self.r = PaddedBytes::new::<NucMatrix>(self.cap, self.blk);
        self.cigar = block_aligner::cigar::Cigar::new(self.cap, self.cap);
        self.blk
    }

    fn align_into(
        &mut self,
        query: &[u8],
        reference: &[u8],
        sc: Scoring,
        out: &mut Vec<CigarOp>,
    ) -> Option<i32> {
        use block_aligner::{cigar::Operation, scores::*};

        // block-aligner takes i8 scores, and its gap penalties are negative.
        let (Ok(m), Ok(x), Ok(open), Ok(extend)) = (
            i8::try_from(sc.match_score),
            i8::try_from(sc.mismatch),
            i8::try_from(-sc.open),
            i8::try_from(-sc.ext),
        ) else {
            return None;
        };
        let matrix = NucMatrix::new_simple(m, x);
        let gaps = Gaps { open, extend };

        if query.len() >= MAX_BLOCK {
            return None; // beyond what a full-width block can cover
        }
        let blk = self.reserve(query.len().max(reference.len()), query.len());
        self.q.set_bytes::<NucMatrix>(query, blk);
        self.r.set_bytes::<NucMatrix>(reference, blk);
        // `min_size > query.len()` is required by FREE_QUERY_END_GAPS, so the
        // range is a single size: no band.
        let lo = block_size_for(query.len());
        self.block
            .align(&self.q, &self.r, &matrix, gaps, lo..=blk, 0);

        let res = self.block.res();
        // A result that did not consume the query is not the alignment asked for.
        if res.query_idx != query.len() {
            return None;
        }
        self.block.trace().cigar_eq(
            &self.q,
            &self.r,
            res.query_idx,
            res.reference_idx,
            &mut self.cigar,
        );

        out.reserve(self.cigar.len() + 4);
        for i in 0..self.cigar.len() {
            let ol = self.cigar.get(i);
            let op = match ol.op {
                Operation::Eq => b'=',
                Operation::X => b'X',
                Operation::I => b'I',
                Operation::D => b'D',
                Operation::M => b'M',
                _ => return None,
            };
            out.push((ol.len, op));
        }
        if out.is_empty() {
            return None;
        }

        // block-aligner's CIGAR covers only the aligned region: with free end
        // gaps the dangling prefix and suffix of the *reference* are simply not
        // in it. parasail's `sg` traceback does include them as I/D ops, and
        // `parse_cigar_diversity_isoform_level` needs them --- `mergeable_start`
        // and `mergeable_end` are computed from exactly those columns, which is
        // what `--delta_iso_len_5`/`_3` bound.
        //
        // So they are reconstructed here from the result indices rather than
        // left out. `res.query_idx`/`res.reference_idx` are the *end* positions,
        // and the CIGAR's own spans give the aligned extent, so the difference is
        // the leading dangle and the tail is whatever is left over.
        let (mut q_span, mut r_span) = (0usize, 0usize);
        for &(l, op) in out.iter() {
            match op {
                b'=' | b'X' | b'M' => {
                    q_span += l;
                    r_span += l;
                }
                b'I' => q_span += l,
                b'D' => r_span += l,
                _ => {}
            }
        }
        let q_lead = res.query_idx.saturating_sub(q_span);
        let r_lead = res.reference_idx.saturating_sub(r_span);
        let q_tail = query.len().saturating_sub(res.query_idx);
        let r_tail = reference.len().saturating_sub(res.reference_idx);
        if r_lead > 0 {
            out.insert(0, (r_lead, b'D'));
        }
        if q_lead > 0 {
            out.insert(0, (q_lead, b'I'));
        }
        if q_tail > 0 {
            out.push((q_tail, b'I'));
        }
        if r_tail > 0 {
            out.push((r_tail, b'D'));
        }
        Some(res.score)
    }
}

/// Does the SIMD aligner ever score *worse* than exact parasail?
///
/// The one question that matters. block-aligner's band is adaptive, so a
/// sub-optimal score is possible in principle; a different equally-optimal path
/// is expected and acceptable, since reference identity is no longer the goal.
/// So this gates the score and merely *reports* CIGAR differences.
///
/// ```text
/// PARASAIL_CASES=<dump>/parasail_cases.tsv \
///   cargo test --release --lib simd::oracle -- --nocapture
/// ```
///
/// isONcorrect validated the same swap on 1 400 alignments. This runs on 54 884.
#[cfg(test)]
pub mod oracle {
    use super::*;

    /// Does the SIMD aligner ever score *worse* than exact parasail?
    ///
    /// The one question that matters. block-aligner's band is adaptive, so a
    /// sub-optimal score is possible in principle; a different equally-optimal
    /// path is expected and acceptable now that reference identity is not the
    /// goal. So this gates the score and merely *reports* CIGAR differences.
    ///
    /// The recorded cases carry parasail's own score and scoring parameters, so
    /// the comparison is against what the real library returned rather than
    /// against a rescore of my own devising.
    ///
    /// isONcorrect validated the same swap on 1 400 alignments; this runs on
    /// 54 884.
    #[test]
    fn never_scores_worse_than_exact_parasail() {
        let Ok(path) = std::env::var("PARASAIL_CASES") else {
            eprintln!(
                "SKIPPED: set PARASAIL_CASES to a parasail_cases.tsv to check the \
                 SIMD aligner against the exact one. Without it this asserts nothing."
            );
            return;
        };
        let text = std::fs::read_to_string(&path).expect("cases readable");
        let (mut n, mut worse, mut better, mut refused) = (0u64, 0u64, 0u64, 0u64);
        let mut worst = 0i32;
        // Bucket by the longer sequence. If divergence concentrates in long
        // pairs it is a block-size / accumulator effect; if it is uniform the
        // objective differs. Guessing between those cost me a wrong conclusion
        // once already.
        let mut buckets = [[0u64; 3]; 5]; // [len bucket][worse, equal, better]
        let mut ops = Vec::new();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            let (s1, s2) = (f[0].as_bytes(), f[1].as_bytes());
            let want: i32 = match f[3].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            let sc = Scoring {
                match_score: f[4].parse().unwrap_or(2),
                mismatch: f[5].parse().unwrap_or(-2),
                open: f[6].parse().unwrap_or(12),
                ext: f[7].parse().unwrap_or(1),
            };
            if s1.is_empty() || s2.is_empty() {
                continue;
            }
            n += 1;
            match semiglobal_ops(s1, s2, sc, &mut ops) {
                None => refused += 1,
                Some(got) => {
                    let ln = s1.len().max(s2.len());
                    let bi = match ln {
                        0..=199 => 0,
                        200..=499 => 1,
                        500..=999 => 2,
                        1000..=1999 => 3,
                        _ => 4,
                    };
                    if got < want {
                        worse += 1;
                        worst = worst.min(got - want);
                        buckets[bi][0] += 1;
                    } else if got > want {
                        better += 1;
                        buckets[bi][2] += 1;
                    } else {
                        buckets[bi][1] += 1;
                    }
                }
            }
        }
        eprintln!(
            "simd oracle: {n} cases | {worse} scored worse (worst delta {worst}) | \
             {better} scored higher | {refused} refused (exact path handles those)"
        );
        for (bi, label) in ["<200", "200-499", "500-999", "1000-1999", "2000+"]
            .iter()
            .enumerate()
        {
            let [w, e, b] = buckets[bi];
            let tot = w + e + b;
            if tot > 0 {
                eprintln!(
                    "  len {label:>10}: {tot:6} cases | {w:6} worse | {e:6} equal | {b:6} better"
                );
            }
        }
        assert_eq!(
            worse, 0,
            "scored worse than exact parasail on {worse} of {n} cases (worst delta \
             {worst}); a different equally-optimal path is fine, a worse score is not"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ops_span(ops: &[CigarOp]) -> (usize, usize) {
        let mut q = 0;
        let mut r = 0;
        for &(l, op) in ops {
            match op {
                b'=' | b'X' | b'M' => {
                    q += l;
                    r += l;
                }
                b'I' => q += l,
                b'D' => r += l,
                _ => {}
            }
        }
        (q, r)
    }

    #[test]
    fn identical_sequences_align_end_to_end() {
        let s = b"ACGTACGTACGTTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(4);
        let mut ops = Vec::new();
        assert!(semiglobal_ops(&s, &s, Scoring::MERGE, &mut ops).is_some());
        let (q, r) = ops_span(&ops);
        assert_eq!((q, r), (s.len(), s.len()));
        assert!(
            ops.iter().all(|&(_, o)| o == b'=' || o == b'M'),
            "identical input should produce only matches, got {ops:?}"
        );
    }

    #[test]
    fn a_contained_sequence_aligns_inside_the_longer_one() {
        // The case `align_to_merge` exists for: free end gaps let the shorter
        // sequence sit inside the longer without paying for the overhang.
        let long = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGTACGTACGTACGT".repeat(8);
        assert!(long.len() > 260, "fixture must be long enough to slice");
        let short = long[40..200].to_vec();
        let mut ops = Vec::new();
        assert!(semiglobal_ops(&long, &short, Scoring::MERGE, &mut ops).is_some());
        let (q, r) = ops_span(&ops);
        // Ops are reported relative to the caller's (s1, s2) = (long, short).
        assert_eq!(q, long.len(), "s1 fully spanned");
        assert_eq!(r, short.len(), "s2 fully spanned");
    }

    #[test]
    fn the_orientation_is_symmetric() {
        // Swapping the arguments must transpose I and D, not change the alignment.
        let long = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(8);
        let short = long[30..150].to_vec();
        let (mut a, mut b) = (Vec::new(), Vec::new());
        assert!(semiglobal_ops(&long, &short, Scoring::MERGE, &mut a).is_some());
        assert!(semiglobal_ops(&short, &long, Scoring::MERGE, &mut b).is_some());
        let (qa, ra) = ops_span(&a);
        let (qb, rb) = ops_span(&b);
        assert_eq!((qa, ra), (long.len(), short.len()));
        assert_eq!(
            (qb, rb),
            (short.len(), long.len()),
            "reversing the arguments must reverse which side each op consumes"
        );
    }

    #[test]
    fn which_ends_are_free() {
        // Establish the objective empirically instead of reasoning from the flag
        // names. Query = shorter. If the *reference's* overhangs are free, adding
        // junk to the long side must not change the score.
        let core = b"ACGTTGCAAGGCTTAACCGGATTCAGGTACGATCGATCGGCTAACGT".repeat(6);
        let mut ops = Vec::new();
        let base = semiglobal_ops(&core, &core, Scoring::MERGE, &mut ops).unwrap();

        let mut ref_pad = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        ref_pad.extend_from_slice(&core);
        ref_pad.extend_from_slice(b"GGGGGGGGGGGGGGGGGGGG");
        let padded_ref = semiglobal_ops(&ref_pad, &core, Scoring::MERGE, &mut ops).unwrap();

        let mut q_pad = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        q_pad.extend_from_slice(&core);
        let padded_query = semiglobal_ops(&core, &q_pad, Scoring::MERGE, &mut ops).unwrap();

        eprintln!("  same/same          score {base}");
        eprintln!("  long side padded    score {padded_ref}  (free if == {base})");
        eprintln!("  short side padded   score {padded_query}  (free if == {base})");
        assert_eq!(padded_ref, base, "the long side's overhangs should be free");
    }

    #[test]
    fn empty_input_is_refused_rather_than_guessed() {
        let mut ops = Vec::new();
        assert!(semiglobal_ops(b"", b"ACGT", Scoring::MERGE, &mut ops).is_none());
        assert!(semiglobal_ops(b"ACGT", b"", Scoring::MERGE, &mut ops).is_none());
        assert!(ops.is_empty());
    }
}
