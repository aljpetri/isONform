//! Minimizer-pair anchors: `get_minimizers_and_positions`,
//! `minimizers_comb_iterator` and `get_minimizer_combinations_database`
//! from `main:118`.
//!
//! An *anchor* is a pair of minimizers `(m1, m2)` whose positions in a read are
//! between `x_low` and `x_high` apart. The database maps each pair to every
//! `(r_id, p1, p2)` that exhibits it, and the interval builder uses it to ask
//! "which other reads span the same pair of anchors, at about the same
//! distance?".
//!
//! Details of the reference that reach the output:
//!
//! * the inner scan **breaks** as soon as the span exceeds `x_high`, relying on
//!   minimizer positions being increasing;
//! * `minimizers_comb_iterator` yields each partner list **reversed**
//!   (`m1_curr_spans[::-1]`), and that order reaches the payload, which decides
//!   the order supporting reads are collected downstream. Not free to change;
//! * the outer loop is over `minimizers[:-1]`, so the last minimizer never
//!   starts a pair;
//! * a pair where both k-mers are poly-A (`'A' * k`) is skipped entirely;
//! * entries with three or fewer occurrences are dropped, as are entries more
//!   than three times as abundant as the read count.
//!
//! # Two places this differs from isONcorrect, which shares the rest
//!
//! 1. **The abundance filter.** isONcorrect drops a pair when it is more frequent
//!    than the read count *and* either the partner is poly-A or it exceeds ten
//!    times the read count, and prints the survivors it considers suspicious.
//!    isONform drops unconditionally at `> 3 * len(reads)` and prints nothing.
//! 2. **`del` is undone immediately.** The reference deletes a singleton entry
//!    and then, on the very next line, indexes `M2[m1][m2]` again — on a
//!    `defaultdict`, which *recreates* the key as an empty array. So the deleted
//!    entries come back empty rather than absent. This port stores only
//!    survivors: an empty entry and a missing entry both mean "no support", so
//!    the data is the same, and `get` returns an empty slice for either.

use rustc_hash::FxBuildHasher;
use std::collections::HashMap;

/// `(r_id, p1, p2)` --- one read exhibiting an anchor pair, and where.
pub type Occurrence = (u32, u32, u32);

/// A minimizer and its position in a read.
pub type Minimizer = (Vec<u8>, usize);

/// Minimizers per read, in `r_id` order.
pub type MinimizersByRead = [(u32, Vec<Minimizer>)];

/// The inner `m2 -> occurrences` map of `M2[m1]`.
pub type PartnerMap = HashMap<Vec<u8>, Vec<Occurrence>, FxBuildHasher>;

/// `M2`: `m1 -> m2 -> occurrences`.
///
/// Keyed by the k-mer bytes, as the reference keys by the k-mer string.
#[derive(Debug, Default)]
pub struct AnchorDb {
    map: HashMap<Vec<u8>, PartnerMap, FxBuildHasher>,
    /// How many entries the singleton rule removed, for reporting.
    pub singletons: usize,
    /// How many the abundance rule removed.
    pub too_abundant: usize,
}

impl AnchorDb {
    /// Occurrences for a pair, or an empty slice.
    ///
    /// Empty and absent are the same answer, which is what lets this store only
    /// survivors where the reference keeps resurrected empty entries.
    pub fn get(&self, m1: &[u8], m2: &[u8]) -> &[Occurrence] {
        self.map
            .get(m1)
            .and_then(|p| p.get(m2))
            .map_or(&[], |v| v.as_slice())
    }

    pub fn len(&self) -> usize {
        self.map.values().map(|p| p.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// `get_minimizers_and_positions`, for the `hash_fcn == "lex"` path `main` uses.
pub fn minimizers_by_read(
    reads: &[(u32, Vec<u8>)],
    w: usize,
    k: usize,
) -> Vec<(u32, Vec<Minimizer>)> {
    reads
        .iter()
        .map(|(r_id, seq)| (*r_id, crate::minimizers::minimizers_lex(seq, k, w)))
        .collect()
}

/// `minimizers_comb_iterator`: for each minimizer, the partners within
/// `(x_low, x_high]`, **reversed**.
///
/// Returns rather than yields; the reference's consumers materialise it anyway
/// (`read_min_comb = [... for ... in minimizers_comb_iterator(...)]`).
pub fn comb_iterator(
    minimizers: &[Minimizer],
    x_low: usize,
    x_high: usize,
) -> Vec<(usize, Vec<usize>)> {
    let mut out = Vec::new();
    if minimizers.is_empty() {
        return out;
    }
    // `minimizers[:-1]` --- the last minimizer never starts a pair.
    for i in 0..minimizers.len() - 1 {
        let p1 = minimizers[i].1;
        let mut spans = Vec::new();
        for (j, (_, p2)) in minimizers.iter().enumerate().skip(i + 1) {
            let d = p2.saturating_sub(p1);
            if x_low < d && d <= x_high {
                spans.push(j);
            } else if d > x_high {
                // Positions increase, so nothing later can qualify.
                break;
            }
        }
        // `m1_curr_spans[::-1]`. The order reaches the payload.
        spans.reverse();
        out.push((i, spans));
    }
    out
}

/// `get_minimizer_combinations_database`.
pub fn build(
    minimizers_by_read: &MinimizersByRead,
    k: usize,
    x_low: usize,
    x_high: usize,
    n_reads: usize,
) -> AnchorDb {
    let mut db = AnchorDb::default();
    let forbidden = vec![b'A'; k];

    for (r_id, minimizers) in minimizers_by_read {
        for (i, partners) in comb_iterator(minimizers, x_low, x_high) {
            let (m1, p1) = (&minimizers[i].0, minimizers[i].1);
            for j in partners {
                let (m2, p2) = (&minimizers[j].0, minimizers[j].1);
                // Both poly-A: skipped entirely.
                if *m2 == forbidden && *m1 == forbidden {
                    continue;
                }
                db.map
                    .entry(m1.clone())
                    .or_default()
                    .entry(m2.clone())
                    .or_default()
                    .push((*r_id, p1 as u32, p2 as u32));
            }
        }
    }

    // The reference walks `list(M2.keys())` and `list(M2[m1].keys())`, i.e.
    // snapshots, so deleting while iterating is safe. Same here by construction.
    let mut singletons = 0usize;
    let mut too_abundant = 0usize;
    for partners in db.map.values_mut() {
        partners.retain(|_, occ| {
            // Units matter here. The reference stores a flat `array("I")` of
            // three integers per occurrence, so its `len(...)` is
            // `3 * occurrences` and its `len(...) // 3` is `occurrences`. This
            // stores tuples, so `occ.len()` is already `occurrences`:
            //
            //   reference `len(...) > 3`            -> occurrences > 1
            //   reference `len(...) // 3 > 3 * n`   -> occurrences > 3 * n
            if occ.len() <= 1 {
                singletons += 1;
                return false;
            }
            if occ.len() > 3 * n_reads {
                too_abundant += 1;
                return false;
            }
            true
        });
    }
    db.map.retain(|_, partners| !partners.is_empty());
    db.singletons = singletons;
    db.too_abundant = too_abundant;
    db
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mins(spec: &[(&str, usize)]) -> Vec<Minimizer> {
        spec.iter()
            .map(|(m, p)| (m.as_bytes().to_vec(), *p))
            .collect()
    }

    #[test]
    fn partners_come_back_reversed() {
        // Positions 0, 5, 10, 15 with x_low=1, x_high=20: minimizer 0 pairs with
        // all three later ones, and the reference yields them furthest-first.
        let m = mins(&[("AA", 0), ("CC", 5), ("GG", 10), ("TT", 15)]);
        let got = comb_iterator(&m, 1, 20);
        assert_eq!(got[0].0, 0);
        assert_eq!(
            got[0].1,
            vec![3, 2, 1],
            "reversed --- `m1_curr_spans[::-1]`, and the order reaches the payload"
        );
    }

    #[test]
    fn the_last_minimizer_never_starts_a_pair() {
        let m = mins(&[("AA", 0), ("CC", 5)]);
        let got = comb_iterator(&m, 1, 20);
        assert_eq!(got.len(), 1, "`minimizers[:-1]` --- only index 0 starts");
    }

    #[test]
    fn the_span_window_is_half_open_low_and_closed_high() {
        // `x_low < p2 - p1 and p2 - p1 <= x_high`: a span exactly x_low is out,
        // exactly x_high is in.
        let m = mins(&[("AA", 0), ("CC", 5), ("GG", 10)]);
        let got = comb_iterator(&m, 5, 10);
        assert_eq!(got[0].1, vec![2], "span 5 excluded, span 10 included");
    }

    #[test]
    fn the_scan_breaks_rather_than_skipping_past_x_high() {
        // The break relies on positions increasing. A partner beyond x_high stops
        // the scan even though a *later* one could not qualify anyway --- the
        // observable difference is only cost, but the structure is the
        // reference's.
        let m = mins(&[("AA", 0), ("CC", 3), ("GG", 100), ("TT", 4)]);
        let got = comb_iterator(&m, 1, 10);
        assert_eq!(
            got[0].1,
            vec![1],
            "stopped at position 100, never reached the out-of-order 4"
        );
    }

    #[test]
    fn a_pair_of_polya_kmers_is_skipped_but_a_half_polya_pair_is_not() {
        let k = 3;
        let poly = "AAA";
        let by_read = vec![
            (1u32, mins(&[(poly, 0), (poly, 5), ("CGT", 10)])),
            (2u32, mins(&[(poly, 0), (poly, 5), ("CGT", 10)])),
            (3u32, mins(&[(poly, 0), (poly, 5), ("CGT", 10)])),
            (4u32, mins(&[(poly, 0), (poly, 5), ("CGT", 10)])),
        ];
        let db = build(&by_read, k, 1, 20, by_read.len());
        assert!(
            db.get(poly.as_bytes(), poly.as_bytes()).is_empty(),
            "poly-A against poly-A is skipped"
        );
        assert!(
            !db.get(poly.as_bytes(), b"CGT").is_empty(),
            "poly-A against a real k-mer is kept"
        );
    }

    #[test]
    fn singletons_are_dropped_and_an_absent_entry_reads_as_no_support() {
        // One read exhibiting a pair gives three array entries, which is `<= 3`.
        let by_read = vec![(1u32, mins(&[("AC", 0), ("GT", 5)]))];
        let db = build(&by_read, 2, 1, 20, 1);
        assert!(db.is_empty());
        assert_eq!(db.singletons, 1);
        // The reference's `del` is undone by a defaultdict read, leaving the key
        // present but empty. Absent and empty are the same answer here.
        assert_eq!(db.get(b"AC", b"GT"), &[] as &[Occurrence]);
    }

    #[test]
    fn two_or_more_occurrences_survive_and_keep_read_order() {
        let by_read = vec![
            (1u32, mins(&[("AC", 0), ("GT", 5)])),
            (2u32, mins(&[("AC", 1), ("GT", 7)])),
        ];
        let db = build(&by_read, 2, 1, 20, 2);
        assert_eq!(db.get(b"AC", b"GT"), &[(1, 0, 5), (2, 1, 7)]);
    }

    #[test]
    fn an_over_abundant_pair_is_dropped() {
        // Threshold is `occurrences > 3 * n_reads`. Five occurrences against a
        // declared read count of 1 is over it; against 5 it is not.
        let by_read: Vec<(u32, Vec<Minimizer>)> = (1..=5u32)
            .map(|r| (r, mins(&[("AC", 0), ("GT", 5)])))
            .collect();
        let dropped = build(&by_read, 2, 1, 20, 1);
        assert!(dropped.is_empty(), "5 > 3 * 1");
        assert_eq!(dropped.too_abundant, 1);
        let kept = build(&by_read, 2, 1, 20, 5);
        assert!(!kept.is_empty(), "5 <= 3 * 5");
    }
}
