//! Weighted interval scheduling: `solve_WIS`, `fill_p2` and
//! `get_intervals_to_correct` from `main:250`.
//!
//! Given every candidate interval a read could be built from, pick the
//! highest-weight set that does not overlap. Weight is
//! `(support - 1) * (stop - start + epsilon)`: `support` counts the read itself,
//! so a span nothing else supports is worth zero.
//!
//! Two details of the reference reach the output:
//!
//! * the intervals are sorted by `stop` before this runs (`main:487`), and the
//!   whole method depends on it;
//! * the caller reverses the result (`opt_indicies[::-1]`) before turning it into
//!   intervals, so the order here is decreasing and that is not incidental.
//!
//! # `fill_p2` is off by one, and it is isONform's bug alone now
//!
//! `fill_p2` builds `p[j]`, the largest index whose interval finishes at or
//! before interval `j` starts. It stores a **0-based** interval index into a
//! table that `solve_WIS` then indexes with the **1-based** `OPT` array
//! (`OPT[p[j]]`). Two consequences:
//!
//! 1. every predecessor is shifted down by one, so compatible earlier intervals
//!    are treated as incompatible;
//! 2. `j_max` starts at 0 and means both "interval 0" and "nothing precedes
//!    this", so an interval with no predecessor is credited with interval 0's
//!    optimum.
//!
//! Both make the selection *conservative* rather than invalid --- fewer intervals
//! than the optimum, never an overlapping set --- which is why it went unnoticed.
//!
//! **This was found, proven and fixed in the isONcorrect port**, on both the Rust
//! and Python sides there. isONform still carries the original. So the situation
//! is unusual and worth stating plainly: the fix is already validated by
//! exhaustive search, but it changes this tool's output, so [`WisOpts::fix_p2`]
//! is **off** by default and the reference's behaviour is what runs. See
//! `PORTING.md` finding 26.
//!
//! Re-established here rather than inherited, on isONform's own code path: over
//! 3 000 random instances checked against exhaustive maximum-weight independent
//! set, the fix is optimal on **3 000 of 3 000** and the reference is worse on
//! **2 085** and better on **none**. (isONcorrect measured 2 040 on its own
//! generator, which is the same story from a different draw.)

/// One candidate interval. `support` counts the read being built *and* the
/// reads that agree with it, which is why the weight subtracts one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: usize,
    pub stop: usize,
    pub support: usize,
}

/// `epsilon` in `solve_WIS`, so a zero-length span still carries some weight.
const EPSILON: f64 = 0.0001;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct WisOpts {
    /// Store 1-based predecessor indices, as the maximum-weight independent set
    /// actually requires. **Off by default**: the reference is off by one and
    /// that is what the recorded goldens contain. See the module docs.
    pub fix_p2: bool,
}

/// `p[j]`: the index whose interval finishes at or before interval `j` starts,
/// or `0` for none. 1-based, with a leading placeholder, matching `p = [None]`.
fn fill_p2(intervals: &[Interval], opts: WisOpts) -> Vec<usize> {
    let mut p = vec![0usize; intervals.len() + 1];
    if intervals.is_empty() {
        return p;
    }
    // `all_intervals_sorted_by_finish[-1][1]` --- the last interval's stop, which
    // is the largest only because the caller sorted by stop.
    let max_stop = intervals[intervals.len() - 1].stop;

    // `{stop: j for j, ... if start < stop}`. A later interval with the same stop
    // overwrites an earlier one, as the dict comprehension does.
    let mut stop_to_max_j: Vec<Option<usize>> = vec![None; max_stop + 1];
    for (j, iv) in intervals.iter().enumerate() {
        if iv.start < iv.stop && iv.stop <= max_stop {
            stop_to_max_j[iv.stop] = Some(j);
        }
    }

    let mut coord: Vec<usize> = vec![0; max_stop + 1];
    let mut j_max = 0usize;
    for (i, slot) in coord.iter_mut().enumerate() {
        if let Some(j) = stop_to_max_j[i] {
            // The one-character difference. The reference stores `j`; the
            // independent set needs `j + 1`, because `OPT` is 1-based and `0`
            // has to keep meaning "no predecessor".
            j_max = if opts.fix_p2 { j + 1 } else { j };
        }
        *slot = j_max;
    }
    for (j, iv) in intervals.iter().enumerate() {
        // The reference indexes `all_choord_to_max_j[start]` directly. `start`
        // cannot exceed `max_stop` for a well-formed interval, but clamping
        // rather than panicking keeps a malformed one out of the panic path.
        p[j + 1] = coord[iv.start.min(max_stop)];
    }
    p
}

/// `(support - 1) * (stop - start + epsilon)`, exactly as the reference has it.
fn weight(iv: &Interval) -> f64 {
    let w = iv.support as f64 - 1.0;
    let span = (iv.stop as f64 - iv.start as f64) + EPSILON;
    w * span
}

/// `solve_WIS`: chosen interval indices, in the reference's **decreasing** order.
///
/// `intervals` must already be sorted by `stop`; the caller does that.
pub fn solve(intervals: &[Interval], opts: WisOpts) -> Vec<usize> {
    let mut opt_indices = Vec::new();
    if intervals.is_empty() {
        return opt_indices;
    }
    let p = fill_p2(intervals, opts);

    // v and OPT are 1-based with a leading placeholder, as in the reference.
    let n = intervals.len();
    let mut v = vec![0.0f64; n + 1];
    for (j, iv) in intervals.iter().enumerate() {
        v[j + 1] = weight(iv);
    }
    let mut opt = vec![0.0f64; n + 1];
    for j in 1..=n {
        opt[j] = (v[j] + opt[p[j]]).max(opt[j - 1]);
    }

    let mut j = n;
    while j > 0 {
        if v[j] + opt[p[j]] > opt[j - 1] {
            // "we have shifted all indices forward by one so we need to reduce
            // to j - 1 because of indexing in python works"
            opt_indices.push(j - 1);
            j = p[j];
        } else {
            j -= 1;
        }
    }
    opt_indices
}

/// `get_intervals_to_correct`: look the chosen indices back up.
///
/// The reference passes `opt_indicies[::-1]`, so callers reverse first. Kept as
/// its own function because it is where the reversal becomes visible.
pub fn intervals_to_correct(opt_indices: &[usize], intervals: &[Interval]) -> Vec<Interval> {
    opt_indices.iter().map(|&j| intervals[j]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: usize, stop: usize, support: usize) -> Interval {
        Interval {
            start,
            stop,
            support,
        }
    }

    /// Every subset, checked for overlap, scored. The ground truth the
    /// off-by-one is measured against.
    fn brute_force(intervals: &[Interval]) -> f64 {
        let n = intervals.len();
        let mut best = 0.0f64;
        for mask in 0u32..(1u32 << n) {
            let chosen: Vec<&Interval> = (0..n)
                .filter(|b| mask & (1 << b) != 0)
                .map(|b| &intervals[b])
                .collect();
            let mut ok = true;
            for a in 0..chosen.len() {
                for b in (a + 1)..chosen.len() {
                    // Sorted by stop, so `a` finishes first.
                    if chosen[b].start < chosen[a].stop {
                        ok = false;
                    }
                }
            }
            if ok {
                let total: f64 = chosen.iter().map(|x| weight(x)).sum();
                if total > best {
                    best = total;
                }
            }
        }
        best
    }

    fn score(intervals: &[Interval], picked: &[usize]) -> f64 {
        picked.iter().map(|&j| weight(&intervals[j])).sum()
    }

    #[test]
    fn the_reference_is_suboptimal_and_the_fix_is_not() {
        // Three intervals that do not overlap at all, so the optimum takes all
        // three. The off-by-one predecessor table cannot see that.
        let ivs = [iv(0, 10, 5), iv(10, 20, 5), iv(20, 30, 5)];
        let best = brute_force(&ivs);

        let fixed = solve(&ivs, WisOpts { fix_p2: true });
        assert_eq!(fixed.len(), 3, "all three are compatible");
        assert!(
            (score(&ivs, &fixed) - best).abs() < 1e-9,
            "the fix reaches the optimum"
        );

        let reference = solve(&ivs, WisOpts::default());
        assert!(
            score(&ivs, &reference) < best,
            "the reference does not: got {}, optimum {best}",
            score(&ivs, &reference)
        );
    }

    #[test]
    fn the_reference_never_returns_an_overlapping_set() {
        // The defect makes the selection conservative, not invalid --- which is
        // why nobody noticed. Worth pinning: whatever else changes, the
        // reference's output stays a legal independent set.
        let mut seed = 12345u64;
        let mut rand = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as usize
        };
        for _ in 0..500 {
            let n = 2 + rand() % 6;
            let mut ivs: Vec<Interval> = (0..n)
                .map(|_| {
                    let start = rand() % 50;
                    let len = 1 + rand() % 20;
                    iv(start, start + len, 1 + rand() % 6)
                })
                .collect();
            ivs.sort_by_key(|x| x.stop);
            for opts in [WisOpts::default(), WisOpts { fix_p2: true }] {
                let picked = solve(&ivs, opts);
                let mut chosen: Vec<Interval> = picked.iter().map(|&j| ivs[j]).collect();
                chosen.sort_by_key(|x| x.stop);
                for pair in chosen.windows(2) {
                    assert!(
                        pair[1].start >= pair[0].stop,
                        "overlapping selection {pair:?} with {opts:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_fix_is_optimal_and_the_reference_is_never_better() {
        // The claim the isONcorrect port established by exhaustive search,
        // re-established here on isONform's own code path rather than inherited.
        let mut seed = 99u64;
        let mut rand = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) as usize
        };
        let (mut fixed_optimal, mut ref_worse, mut ref_better) = (0, 0, 0);
        for _ in 0..3000 {
            let n = 2 + rand() % 7;
            let mut ivs: Vec<Interval> = (0..n)
                .map(|_| {
                    let start = rand() % 40;
                    let len = 1 + rand() % 15;
                    iv(start, start + len, 1 + rand() % 5)
                })
                .collect();
            ivs.sort_by_key(|x| x.stop);
            let best = brute_force(&ivs);
            let f = score(&ivs, &solve(&ivs, WisOpts { fix_p2: true }));
            let r = score(&ivs, &solve(&ivs, WisOpts::default()));
            if (f - best).abs() < 1e-9 {
                fixed_optimal += 1;
            }
            if r < f - 1e-9 {
                ref_worse += 1;
            }
            if r > f + 1e-9 {
                ref_better += 1;
            }
        }
        assert_eq!(fixed_optimal, 3000, "the fix is optimal on every instance");
        assert_eq!(ref_better, 0, "the reference is never better");
        eprintln!("MEASURED: reference suboptimal on {ref_worse} of 3000");
        assert!(
            ref_worse > 0,
            "and it is worse often enough to matter: {ref_worse} of 3000"
        );
    }

    #[test]
    fn indices_come_back_in_decreasing_order() {
        // The caller reverses (`opt_indicies[::-1]`), so this order is part of
        // the contract rather than an implementation detail.
        let ivs = [iv(0, 10, 5), iv(12, 20, 5), iv(22, 30, 5)];
        let picked = solve(&ivs, WisOpts { fix_p2: true });
        assert!(
            picked.windows(2).all(|w| w[0] > w[1]),
            "expected decreasing, got {picked:?}"
        );
    }

    #[test]
    fn a_span_nothing_else_supports_is_worth_nothing() {
        // support == 1 means only the read itself, so `w - 1` is zero.
        assert_eq!(weight(&iv(0, 100, 1)), 0.0);
        assert!(solve(&[iv(0, 100, 1)], WisOpts::default()).is_empty());
    }

    #[test]
    fn no_intervals_is_empty_rather_than_a_panic() {
        assert!(solve(&[], WisOpts::default()).is_empty());
        assert!(intervals_to_correct(&[], &[]).is_empty());
    }
}
