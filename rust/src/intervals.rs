//! Interval selection: `find_most_supported_span` and the per-read loop of
//! `main` that turns it into `all_intervals_for_graph`.
//!
//! This is the last stage before graph construction. For each read it finds
//! every anchor pair other reads also span at about the same distance, turns
//! those into weighted candidate intervals, and hands them to
//! [`crate::wis`] to pick a non-overlapping set. Reads that yield none are
//! written to the skip file and never reach the graph.
//!
//! # `find_most_supported_span` does *not* transfer from isONcorrect
//!
//! Despite `PORTING.md`'s reconnaissance note, this one is isONform's own and
//! much simpler. isONcorrect's compares the actual subsequences, tracks quality
//! values and edit distances, and memoises across reads. isONform's keeps a
//! supporting read whenever the two spans differ in *length* by less than
//! `delta_len`, and never looks at a base. Fifteen lines against sixty, and
//! nothing in common but the name and the `grouper(relevant_reads, 3)` walk.
//!
//! # Roughly a third of the reference's loop is dead
//!
//! `previously_corrected_regions` is created as an empty `defaultdict(list)`,
//! read in four places, and `del`-ed --- **nothing ever appends to it**
//! (`main:423-495`). `prev_visited_intervals` is created empty and never
//! appended to either. So on every path through the reference:
//!
//! * `read_previously_considered_positions` is always the empty set, which makes
//!   `not_prev_corrected_spans` equal to `m1_curr_spans` --- the filter is a
//!   no-op;
//! * `pos_group` is always empty, so `not_prev_corrected_spans2` is always empty;
//! * `all_intervals.extend(prev_visited_intervals)` extends by nothing;
//! * both `if previously_corrected_regions[r_id]:` branches never run.
//!
//! That is inherited scaffolding: in isONcorrect the same machinery is live,
//! because correction happens in rounds and later rounds skip what earlier ones
//! fixed. isONform runs one pass and never populates it. This port therefore
//! implements the reachable path only, and `PORTING.md` finding 27 records the
//! dead code rather than the port silently dropping it.

use crate::anchors::AnchorDb;
use crate::wis::{self, Interval, WisOpts};
use rustc_hash::FxHashMap;

/// One read's chosen intervals, as `generateGraphfromIntervals` wants them.
///
/// `support` is the number of reads spanning it, including this one.
/// `instance` is the flat `(r_id, p1, p2)` triple list the reference carries in
/// an `array("I")`, in the order it collected them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChosenInterval {
    pub start: usize,
    pub stop: usize,
    pub support: usize,
    pub instance: Vec<(u32, u32, u32)>,
}

/// `find_most_supported_span`: candidate intervals for one `(m1, p1)`.
///
/// Appends to `out`, as the reference appends to `all_intervals`.
#[allow(clippy::too_many_arguments)]
/// Optional sequence check on occurrences, the piece isONcorrect has and
/// isONform's Python stripped.
///
/// isONform admits an occurrence on span *length* alone
/// (`abs((p2-p1)-(pos2-pos1)) < delta_len`) and never compares the sequence.
/// isONcorrect required `read_seq == ref_seq`, or an edit distance under a
/// quality-derived threshold. `ISONFORM_SEQ_CHECK` restores a simplified form:
/// the occurrence's sequence must be within the given divergence of the current
/// read's, measured as Hamming over the overlap plus the length difference.
///
/// `ISONFORM_SEQ_CHECK=0.0` means exact match only; unset means no check, which
/// is isONform's behaviour. Measured on real data: 83% of admitted occurrences
/// already match exactly, so this is expected to prune the tail rather than most
/// of the list --- which is why it is worth measuring rather than assuming.
pub fn seq_check_threshold() -> Option<f64> {
    static T: std::sync::OnceLock<Option<f64>> = std::sync::OnceLock::new();
    *T.get_or_init(|| {
        std::env::var("ISONFORM_SEQ_CHECK")
            .ok()
            .and_then(|v| v.trim().parse::<f64>().ok())
    })
}

#[allow(clippy::too_many_arguments)]
pub fn find_most_supported_span(
    r_id: u32,
    m1: &[u8],
    p1: usize,
    partners: &[(Vec<u8>, usize)],
    db: &AnchorDb,
    k: usize,
    delta_len: i64,
    seqs: Option<&[(u32, Vec<u8>)]>,
    out: &mut Vec<ChosenInterval>,
) {
    for (m2, p2) in partners {
        let relevant = db.get(m1, m2);
        // `len(relevant_reads) // 3 >= 3`: at least three occurrences.
        if relevant.len() < 3 {
            continue;
        }
        // The read itself goes in first, before any check.
        let mut instance = vec![(r_id, p1 as u32, *p2 as u32)];
        let this_span = *p2 as i64 - p1 as i64;
        for &(other, pos1, pos2) in relevant {
            if other == r_id {
                continue;
            }
            let other_span = pos2 as i64 - pos1 as i64;
            if (this_span - other_span).abs() < delta_len {
                // isONcorrect's check, restored behind `ISONFORM_SEQ_CHECK`.
                if let (Some(tol), Some(all)) = (seq_check_threshold(), seqs) {
                    let mine = all.iter().find(|(r, _)| *r == r_id).map(|(_, s)| s);
                    let theirs = all.iter().find(|(r, _)| *r == other).map(|(_, s)| s);
                    if let (Some(a), Some(b)) = (mine, theirs) {
                        let (x0, x1) = (p1.min(a.len()), (*p2 + k).min(a.len()));
                        let (y0, y1) = (
                            (pos1 as usize).min(b.len()),
                            (pos2 as usize + k).min(b.len()),
                        );
                        if x0 < x1 && y0 < y1 {
                            let (x, y) = (&a[x0..x1], &b[y0..y1]);
                            let n = x.len().min(y.len());
                            let ham = (0..n).filter(|&i| x[i] != y[i]).count();
                            let div = (ham + x.len().abs_diff(y.len())) as f64 / n.max(1) as f64;
                            if div > tol {
                                SEQ_REJECTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                continue;
                            }
                        }
                    }
                }
                instance.push((other, pos1, pos2));
            }
        }
        let support = instance.len();
        out.push(ChosenInterval {
            // `(p1 + k_size, p2, len(seqs) // 3, seqs)`. The interval starts
            // *after* the first anchor and ends at the start of the second.
            start: p1 + k,
            stop: *p2,
            support,
            instance,
        });
    }
}

/// What one batch's front half produces.
/// Instrumentation for the missing sequence check --- see `ISONFORM_SEQ_CHECK_STATS`.
/// Occurrences rejected by `ISONFORM_SEQ_CHECK`.
/// Candidate intervals seen before WIS, and how many have a near-duplicate of
/// equal support --- see `ISONFORM_TIE_STATS`.
pub static TIE_CANDIDATES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static TIE_AMBIGUOUS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static SEQ_REJECTED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static SEQ_TOTAL: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static SEQ_EXACT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static SEQ_DIV: [std::sync::atomic::AtomicU64; 5] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];

#[derive(Debug, Default)]
pub struct BatchIntervals {
    /// `all_intervals_for_graph`, keyed by the 1-based `graph_id` the reference
    /// assigns in read order.
    pub by_graph_id: FxHashMap<u32, Vec<ChosenInterval>>,
    /// `new_all_reads`: `graph_id -> r_id`, so the caller can map back.
    pub read_of_graph_id: FxHashMap<u32, u32>,
    /// Reads that produced no intervals, in the order the reference writes them
    /// to the skip file.
    pub skipped: Vec<u32>,
}

/// The per-read loop of `main`, for one batch.
///
/// `reads` must be in ascending `r_id` order: the reference iterates
/// `sorted(reads)` and `graph_id` is assigned in that order.
pub fn build_batch(
    reads: &[(u32, Vec<u8>)],
    w: usize,
    k: usize,
    x_low: usize,
    x_high: usize,
    delta_len: i64,
    opts: WisOpts,
) -> BatchIntervals {
    let mins = crate::anchors::minimizers_by_read(reads, w, k);
    let db = crate::anchors::build(&mins, k, x_low, x_high, reads.len());

    let mut out = BatchIntervals::default();
    let mut graph_id = 1u32;

    for (r_id, minimizers) in &mins {
        let mut all_intervals: Vec<ChosenInterval> = Vec::new();
        for (i, partners) in crate::anchors::comb_iterator(minimizers, x_low, x_high) {
            // The reference guards on `if not_prev_corrected_spans:`, which is
            // `m1_curr_spans` on every reachable path --- see the module docs.
            if partners.is_empty() {
                continue;
            }
            let (m1, p1) = (&minimizers[i].0, minimizers[i].1);
            // Note it passes `m1_curr_spans`, not the filtered list.
            let spans: Vec<(Vec<u8>, usize)> = partners
                .iter()
                .map(|&j| (minimizers[j].0.clone(), minimizers[j].1))
                .collect();
            find_most_supported_span(
                *r_id,
                m1,
                p1,
                &spans,
                &db,
                k,
                delta_len,
                seq_check_threshold().map(|_| reads),
                &mut all_intervals,
            );
        }

        if all_intervals.is_empty() {
            out.skipped.push(*r_id);
            continue;
        }

        // Tie structure among the *candidates*, before WIS selects. A candidate
        // is "ambiguous" when another candidate sits within 20 bases on both
        // ends with the same support: WIS then has no principled reason to
        // prefer either, and two reads can legitimately select different members
        // --- which would break node sharing without either read being wrong.
        if std::env::var_os("ISONFORM_TIE_STATS").is_some() {
            use std::sync::atomic::Ordering::Relaxed;
            TIE_CANDIDATES.fetch_add(all_intervals.len() as u64, Relaxed);
            for (i, a) in all_intervals.iter().enumerate() {
                let mut amb = false;
                for (j, b) in all_intervals.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    if a.start.abs_diff(b.start) <= 20
                        && a.stop.abs_diff(b.stop) <= 20
                        && a.support == b.support
                    {
                        amb = true;
                        break;
                    }
                }
                if amb {
                    TIE_AMBIGUOUS.fetch_add(1, Relaxed);
                }
            }
        }

        // `all_intervals.sort(key=lambda x: x[1])` --- by stop. Python's sort is
        // stable, so ties keep discovery order, and that order reaches the
        // selection.
        all_intervals.sort_by_key(|iv| iv.stop);

        let flat: Vec<Interval> = all_intervals
            .iter()
            .map(|iv| Interval {
                start: iv.start,
                stop: iv.stop,
                support: iv.support,
            })
            .collect();
        let mut chosen = wis::solve(&flat, opts);
        // `ISONFORM_TRACE_GRAPH_ID=<id>`: every candidate this read offered WIS,
        // with the weight `(support - 1) * (span + eps)` that decided it, and
        // whether WIS took it. Diffing one read's trace across two batch sizes
        // shows exactly which candidate displaced which when the graph spirals
        // (finding 45: the onset is between 240 and 250 reads).
        if let Ok(want) = std::env::var("ISONFORM_TRACE_GRAPH_ID") {
            if want.parse::<u32>() == Ok(graph_id) {
                let taken: std::collections::HashSet<usize> = chosen.iter().copied().collect();
                for (j, iv) in flat.iter().enumerate() {
                    let span = iv.stop as f64 - iv.start as f64 + 0.0001;
                    eprintln!(
                        "trace gid={} cand={j} start={} stop={} span={} support={} \
                         weight={:.1} taken={}",
                        graph_id,
                        iv.start,
                        iv.stop,
                        iv.stop - iv.start,
                        iv.support,
                        (iv.support as f64 - 1.0) * span,
                        taken.contains(&j)
                    );
                }
            }
        }
        // `opt_indicies[::-1]` --- the caller reverses before looking them up.
        chosen.reverse();

        if chosen.is_empty() {
            out.skipped.push(*r_id);
            continue;
        }
        let picked: Vec<ChosenInterval> =
            chosen.iter().map(|&j| all_intervals[j].clone()).collect();
        // Would isONcorrect's sequence check have admitted these occurrences?
        // isONform keeps only `abs((p2-p1)-(pos2-pos1)) < delta_len` --- a
        // *length* test with no sequence comparison. isONcorrect additionally
        // required `read_seq == ref_seq` or an edit distance under a
        // quality-derived threshold. This measures, on the intervals actually
        // chosen, how similar the co-occurring sequences really are.
        if std::env::var_os("ISONFORM_SEQ_CHECK_STATS").is_some() {
            use std::sync::atomic::Ordering::Relaxed;
            let me = reads.iter().find(|(r, _)| r == r_id).map(|(_, s)| s);
            for iv in picked.iter() {
                for &(r, pos1, pos2) in iv.instance.iter() {
                    if r == *r_id {
                        continue;
                    }
                    let (Some(mine), Some(other)) =
                        (me, reads.iter().find(|(rr, _)| *rr == r).map(|(_, s)| s))
                    else {
                        continue;
                    };
                    // `iv.start` is already `p1 + k` (see `find_most_supported_span`),
                    // while the occurrence's `pos1` is a raw minimizer position.
                    // isONcorrect slices `seq[p1 : p2 + k]`, so undo the +k here or
                    // the two slices are offset by k and every comparison fails.
                    let a = iv.start.saturating_sub(k);
                    let b = (iv.stop + k).min(mine.len());
                    let c = pos1 as usize;
                    let d = (pos2 as usize + k).min(other.len());
                    if a >= b || c >= d {
                        continue;
                    }
                    SEQ_TOTAL.fetch_add(1, Relaxed);
                    let (x, y) = (&mine[a..b], &other[c..d]);
                    if x == y {
                        SEQ_EXACT.fetch_add(1, Relaxed);
                    }
                    // A cheap divergence proxy: length delta plus Hamming over
                    // the overlap. Not edit distance, but enough to say whether
                    // these are "the same sequence" or not.
                    let n = x.len().min(y.len());
                    let ham = (0..n).filter(|&i| x[i] != y[i]).count();
                    let dv = ham + x.len().abs_diff(y.len());
                    let frac = dv as f64 / n.max(1) as f64;
                    let bucket = if x == y {
                        0
                    } else if frac <= 0.05 {
                        1
                    } else if frac <= 0.15 {
                        2
                    } else if frac <= 0.40 {
                        3
                    } else {
                        4
                    };
                    SEQ_DIV[bucket].fetch_add(1, Relaxed);
                }
            }
        }
        out.by_graph_id.insert(graph_id, picked);
        out.read_of_graph_id.insert(graph_id, *r_id);
        graph_id += 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::Minimizer;

    fn mins(spec: &[(&str, usize)]) -> Vec<Minimizer> {
        spec.iter()
            .map(|(m, p)| (m.as_bytes().to_vec(), *p))
            .collect()
    }

    fn db_from(by_read: &[(u32, Vec<Minimizer>)], k: usize, n: usize) -> AnchorDb {
        crate::anchors::build(by_read, k, 1, 100, n)
    }

    #[test]
    fn the_read_itself_is_always_the_first_instance_entry() {
        let by_read: Vec<(u32, Vec<Minimizer>)> = (1..=3u32)
            .map(|r| (r, mins(&[("AC", 0), ("GT", 10)])))
            .collect();
        let db = db_from(&by_read, 2, 3);
        let mut out = Vec::new();
        find_most_supported_span(
            2,
            b"AC",
            0,
            &[(b"GT".to_vec(), 10)],
            &db,
            2,
            5,
            None,
            &mut out,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].instance[0].0, 2, "the read being built comes first");
        assert!(!out[0].instance[1..].iter().any(|(r, _, _)| *r == 2));
    }

    #[test]
    fn the_interval_starts_after_the_first_anchor() {
        // `(p1 + k_size, p2, ...)`, not `(p1, p2)`.
        let by_read: Vec<(u32, Vec<Minimizer>)> = (1..=3u32)
            .map(|r| (r, mins(&[("ACGT", 5), ("TTTT", 40)])))
            .collect();
        let db = db_from(&by_read, 4, 3);
        let mut out = Vec::new();
        find_most_supported_span(
            1,
            b"ACGT",
            5,
            &[(b"TTTT".to_vec(), 40)],
            &db,
            4,
            5,
            None,
            &mut out,
        );
        assert_eq!((out[0].start, out[0].stop), (9, 40));
    }

    #[test]
    fn a_supporting_read_is_kept_only_when_the_span_lengths_are_close() {
        // Spans: read 1 has 10, read 2 has 10, read 3 has 40. With delta_len 5
        // only reads 1 and 2 agree; the length is all that is compared, never
        // the sequence.
        let by_read = vec![
            (1u32, mins(&[("AC", 0), ("GT", 10)])),
            (2u32, mins(&[("AC", 3), ("GT", 13)])),
            (3u32, mins(&[("AC", 0), ("GT", 40)])),
        ];
        let db = db_from(&by_read, 2, 3);
        let mut out = Vec::new();
        find_most_supported_span(
            1,
            b"AC",
            0,
            &[(b"GT".to_vec(), 10)],
            &db,
            2,
            5,
            None,
            &mut out,
        );
        let ids: Vec<u32> = out[0].instance.iter().map(|(r, _, _)| *r).collect();
        assert_eq!(ids, vec![1, 2], "read 3's span is 40, too far from 10");
        assert_eq!(out[0].support, 2);
    }

    #[test]
    fn fewer_than_three_occurrences_produces_no_interval() {
        // `len(relevant_reads) // 3 >= 3`. Two reads exhibiting the pair is not
        // enough, even though the anchor database keeps the entry.
        let by_read: Vec<(u32, Vec<Minimizer>)> = (1..=2u32)
            .map(|r| (r, mins(&[("AC", 0), ("GT", 10)])))
            .collect();
        let db = db_from(&by_read, 2, 2);
        assert_eq!(db.get(b"AC", b"GT").len(), 2, "the entry survives");
        let mut out = Vec::new();
        find_most_supported_span(
            1,
            b"AC",
            0,
            &[(b"GT".to_vec(), 10)],
            &db,
            2,
            5,
            None,
            &mut out,
        );
        assert!(out.is_empty(), "but two occurrences yield no interval");
    }

    #[test]
    fn a_read_with_no_intervals_is_skipped_and_gets_no_graph_id() {
        // Three identical reads share anchors; a fourth unrelated one does not.
        let shared = "ACGTACGTAAACCCGGGTTTACGTACGTTTTTGGGGCCCCAAAA";
        let lonely = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let reads: Vec<(u32, Vec<u8>)> = vec![
            (1, shared.into()),
            (2, shared.into()),
            (3, shared.into()),
            (4, lonely.into()),
        ];
        let got = build_batch(&reads, 10, 5, 1, 40, 5, WisOpts::reference());
        assert!(
            got.skipped.contains(&4),
            "the unrelated read is skipped, got {:?}",
            got.skipped
        );
        assert!(
            !got.read_of_graph_id.values().any(|r| *r == 4),
            "and never receives a graph id"
        );
    }

    #[test]
    fn graph_ids_are_assigned_in_read_order_starting_at_one() {
        let shared = "ACGTACGTAAACCCGGGTTTACGTACGTTTTTGGGGCCCCAAAA";
        let reads: Vec<(u32, Vec<u8>)> = (1..=3u32)
            .map(|r| (r, shared.as_bytes().to_vec()))
            .collect();
        let got = build_batch(&reads, 10, 5, 1, 40, 5, WisOpts::reference());
        assert!(!got.by_graph_id.is_empty(), "these reads share anchors");
        let mut ids: Vec<u32> = got.read_of_graph_id.keys().copied().collect();
        ids.sort_unstable();
        assert_eq!(ids[0], 1, "graph ids are 1-based");
        for (gid, r) in &got.read_of_graph_id {
            // Assigned in ascending read order, so the two rank the same.
            let rank = ids.iter().position(|x| x == gid).unwrap();
            let mut rs: Vec<u32> = got.read_of_graph_id.values().copied().collect();
            rs.sort_unstable();
            assert_eq!(rs[rank], *r);
        }
    }

    #[test]
    fn chosen_intervals_come_back_in_increasing_order() {
        // `solve_WIS` returns decreasing indices and the caller reverses, so what
        // reaches the graph is ascending. The graph builder depends on it.
        let shared = "ACGTACGTAAACCCGGGTTTACGTACGTTTTTGGGGCCCCAAAAGGTTCCAATTGGCCAA";
        let reads: Vec<(u32, Vec<u8>)> = (1..=4u32)
            .map(|r| (r, shared.as_bytes().to_vec()))
            .collect();
        let got = build_batch(&reads, 10, 5, 1, 40, 5, WisOpts::reference());
        for ivs in got.by_graph_id.values() {
            for pair in ivs.windows(2) {
                assert!(
                    pair[0].stop <= pair[1].stop,
                    "expected ascending by stop, got {:?} then {:?}",
                    pair[0],
                    pair[1]
                );
            }
        }
    }
}
