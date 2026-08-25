//! Bottom-*k* minhash sketches, to stop paying for alignments that cannot merge.
//!
//! # Why this exists
//!
//! Repairing cross-batch merging (`PORTING.md` finding 31) traded a correctness
//! bug for a performance cliff. The merge compares every isoform of every batch
//! against every isoform of every later batch, and on `droso_deep`'s deepest
//! cluster that is 26 batches × ~133 isoforms each: **325 batch pairs, ~5.7
//! million alignments, about 8 hours for one cluster.** Measured, not estimated
//! --- the run sat at 99% CPU for 93 minutes inside
//! `join_back_via_batch_merging` without finishing that cluster.
//!
//! A length-based prefilter was tried first and rejected on evidence: of 6 150
//! real pairs, the conditions a length test could decide (`min(len) < 100`,
//! shared region `< 100`) fired **zero** times, and the length differences of
//! merges and of rejects overlap almost completely --- largest difference among
//! merges 851, smallest among structural rejects 3. No sound length window
//! separates them.
//!
//! What does separate them is *sequence content*. 82% of pairs are rejected as
//! structural, and in a cluster of 25 000 reads many of those pairs are simply
//! different transcripts that happen to share a cluster. Two unrelated
//! consensuses share few k-mers; two isoforms of one transcript share most of
//! them even when they differ by a whole exon. A sketch comparison costs O(k)
//! against an alignment's O(n·m).
//!
//! # Not used. Rejected for lack of a false-negative bound.
//!
//! This module is kept for its reasoning and its measurements, not because
//! anything calls it. A minhash sketch gives no *guarantee*: its false-negative
//! rate can be measured on a corpus but not bounded, and a filter that decides
//! which pairs are **never** examined cannot rest on a rate observed once. The
//! bar is a window guarantee --- a scheme whose density and coverage are
//! provable, so the pairs it skips are ones that provably cannot merge.
//!
//! What survives here is the shape of the problem and the tests: the case a
//! filter must not reject is `an_isoform_missing_a_whole_exon_still_looks_similar`,
//! because alignment refuses that pair as structural and the filter runs *before*
//! alignment --- so skipping it would silently change which pairs are considered.
//! Any replacement should keep that test.
//!
//! # The measurement it was built with
//!
//! A sketch can be wrong in the direction that matters --- skipping a pair that
//! would have merged. So the filter ships with a verification mode
//! (`ISONFORM_SKETCH_VERIFY=1`) that aligns the pairs it would have skipped and
//! counts how many would have merged. The false-negative rate is a measurement,
//! not a hope, and `ISONFORM_SKETCH=0` turns the filter off entirely.

/// A bottom-*k* minhash sketch: the `K` smallest k-mer hashes, sorted.
///
/// Bottom-*k* rather than *k* independent hash functions: one pass, one hash per
/// k-mer, and the Jaccard estimate is the standard one.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Sketch {
    /// Sorted ascending, at most [`Sketch::SIZE`] entries.
    hashes: Vec<u64>,
    /// Distinct k-mers seen, so a short sequence with fewer than `SIZE` k-mers
    /// is compared exactly rather than by estimate.
    n_kmers: usize,
}

/// A 64-bit mix, so neighbouring k-mers do not land in neighbouring buckets.
#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    x
}

impl Sketch {
    /// How many hashes are kept. 128 keeps the Jaccard estimate's standard error
    /// near 1/sqrt(128) ≈ 9%, which is ample for a filter whose threshold is far
    /// from the interesting region.
    pub const SIZE: usize = 128;
    /// K-mer length. Short enough that a 15-mer survives the error rate of an
    /// ONT consensus, long enough to be specific in a transcriptome.
    pub const K: usize = 15;

    /// Sketch one sequence. O(len).
    pub fn of(seq: &[u8]) -> Self {
        if seq.len() < Self::K {
            return Self::default();
        }
        // Rolling 2-bit packing. A non-ACGT base breaks the window, which is the
        // conservative choice: it drops k-mers rather than inventing them.
        let mut hashes: Vec<u64> = Vec::with_capacity(seq.len());
        let mut code: u64 = 0;
        let mut run = 0usize;
        for &b in seq {
            let two = match b | 0x20 {
                b'a' => 0u64,
                b'c' => 1,
                b'g' => 2,
                b't' => 3,
                _ => {
                    run = 0;
                    continue;
                }
            };
            code = ((code << 2) | two) & ((1u64 << (2 * Self::K)) - 1);
            run += 1;
            if run >= Self::K {
                hashes.push(mix(code));
            }
        }
        hashes.sort_unstable();
        hashes.dedup();
        let n_kmers = hashes.len();
        hashes.truncate(Self::SIZE);
        Self { hashes, n_kmers }
    }

    /// Estimated Jaccard similarity of the two k-mer sets, in `0.0..=1.0`.
    ///
    /// The bottom-*k* estimator: merge the two sketches, take the smallest `SIZE`
    /// distinct values overall, and report the fraction of those present in both.
    /// When either side holds fewer than `SIZE` k-mers the sketch *is* the set and
    /// this is exact.
    pub fn jaccard(&self, other: &Self) -> f64 {
        if self.hashes.is_empty() || other.hashes.is_empty() {
            return 0.0;
        }
        let (mut i, mut j) = (0usize, 0usize);
        let (mut seen, mut shared) = (0usize, 0usize);
        while seen < Self::SIZE && (i < self.hashes.len() || j < other.hashes.len()) {
            let a = self.hashes.get(i).copied();
            let b = other.hashes.get(j).copied();
            match (a, b) {
                (Some(x), Some(y)) if x == y => {
                    i += 1;
                    j += 1;
                    shared += 1;
                    seen += 1;
                }
                (Some(x), Some(y)) if x < y => {
                    i += 1;
                    seen += 1;
                }
                (Some(_), Some(_)) => {
                    j += 1;
                    seen += 1;
                }
                (Some(_), None) => {
                    i += 1;
                    seen += 1;
                }
                (None, Some(_)) => {
                    j += 1;
                    seen += 1;
                }
                (None, None) => break,
            }
        }
        if seen == 0 {
            return 0.0;
        }
        shared as f64 / seen as f64
    }

    /// Distinct k-mers in the sequence, before truncation.
    pub fn n_kmers(&self) -> usize {
        self.n_kmers
    }
}

/// The filter's settings, from the environment. Read once.
#[derive(Debug, Clone, Copy)]
pub struct SketchFilter {
    /// Pairs below this estimated Jaccard skip alignment. `None` disables the
    /// filter, which is what `ISONFORM_SKETCH=0` does.
    pub min_jaccard: Option<f64>,
    /// Align the skipped pairs anyway and count how many would have merged.
    /// Turns the filter into a measurement.
    pub verify: bool,
}

impl SketchFilter {
    pub fn from_env() -> Self {
        static F: std::sync::OnceLock<SketchFilter> = std::sync::OnceLock::new();
        *F.get_or_init(|| {
            let raw = std::env::var("ISONFORM_SKETCH").ok();
            let min_jaccard = match raw.as_deref() {
                Some("0") | Some("off") => None,
                Some(v) => v.parse::<f64>().ok().or(Some(DEFAULT_MIN_JACCARD)),
                None => Some(DEFAULT_MIN_JACCARD),
            };
            SketchFilter {
                min_jaccard,
                verify: std::env::var_os("ISONFORM_SKETCH_VERIFY").is_some(),
            }
        })
    }

    /// Should this pair be aligned at all?
    pub fn admits(&self, a: &Sketch, b: &Sketch) -> bool {
        match self.min_jaccard {
            None => true,
            Some(t) => a.jaccard(b) >= t,
        }
    }
}

/// Deliberately far below the region where merges live, so the filter removes
/// unrelated pairs and nothing near the decision boundary. Chosen from the
/// measured distribution, not by taste --- see `PORTING.md`.
pub const DEFAULT_MIN_JACCARD: f64 = 0.02;

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(seed: u64, len: usize) -> Vec<u8> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = mix(s);
                b"ACGT"[(s >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn a_sequence_is_identical_to_itself() {
        let a = Sketch::of(&seq(1, 1000));
        assert_eq!(a.jaccard(&a), 1.0);
    }

    #[test]
    fn unrelated_sequences_share_almost_nothing() {
        let a = Sketch::of(&seq(1, 1400));
        let b = Sketch::of(&seq(999, 1400));
        let j = a.jaccard(&b);
        assert!(j < 0.01, "unrelated sequences estimated at {j}");
        assert!(
            !SketchFilter {
                min_jaccard: Some(DEFAULT_MIN_JACCARD),
                verify: false
            }
            .admits(&a, &b),
            "the filter should skip this pair"
        );
    }

    #[test]
    fn an_isoform_missing_a_whole_exon_still_looks_similar() {
        // The case the filter must NOT reject: two isoforms of one transcript
        // differing by a 200 bp internal exon. Alignment calls that structural
        // and refuses to merge, but the filter runs *before* alignment and must
        // let it through --- otherwise it changes which pairs get considered.
        let full = seq(7, 1400);
        let mut skipped = full[..600].to_vec();
        skipped.extend_from_slice(&full[800..]);
        let (a, b) = (Sketch::of(&full), Sketch::of(&skipped));
        let j = a.jaccard(&b);
        assert!(j > 0.5, "exon-skipping pair estimated at only {j}");
        assert!(SketchFilter {
            min_jaccard: Some(DEFAULT_MIN_JACCARD),
            verify: false
        }
        .admits(&a, &b));
    }

    #[test]
    fn a_noisy_copy_still_looks_similar() {
        // ONT consensus error: 2% of bases substituted.
        let a = seq(11, 1400);
        let mut b = a.clone();
        let mut s = 3u64;
        for _ in 0..28 {
            s = mix(s);
            let p = (s >> 33) as usize % b.len();
            b[p] = b"ACGT"[(s >> 11) as usize % 4];
        }
        let j = Sketch::of(&a).jaccard(&Sketch::of(&b));
        assert!(j > 0.3, "2% substituted estimated at only {j}");
    }

    #[test]
    fn a_short_sequence_sketches_exactly_or_not_at_all() {
        assert_eq!(Sketch::of(b"ACGT").n_kmers(), 0, "shorter than K");
        let s = Sketch::of(&seq(5, 40));
        assert!(s.n_kmers() > 0 && s.n_kmers() < Sketch::SIZE);
    }

    #[test]
    fn non_acgt_breaks_the_window_rather_than_inventing_kmers() {
        let clean = seq(13, 200);
        let mut dirty = clean.clone();
        dirty[100] = b'N';
        let a = Sketch::of(&clean);
        let b = Sketch::of(&dirty);
        assert!(b.n_kmers() < a.n_kmers(), "the N must drop its windows");
        assert!(a.jaccard(&b) > 0.5, "but the rest still matches");
    }

    #[test]
    fn the_filter_is_off_when_asked() {
        let a = Sketch::of(&seq(1, 500));
        let b = Sketch::of(&seq(2, 500));
        assert!(SketchFilter {
            min_jaccard: None,
            verify: false
        }
        .admits(&a, &b));
    }
}
