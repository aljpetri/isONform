//! Minimizer extraction, matching `main:81`'s `get_kmer_minimizers`.
//!
//! `hash_fcn` is hardcoded to `"lex"` in `main`, so minimizers are selected by
//! **lexicographic order on the k-mer string**, not by a hash. That is the
//! determinism fix of `PORTING.md` finding 1; before it the reference hashed the
//! k-mer with CPython's seed-randomised `str.__hash__` and every output file
//! varied run to run. The `rev_lex` maximizer variant is reachable only by
//! editing the source, and is ported here for completeness.
//!
//! Four details in the reference are easy to get wrong and all of them change
//! output:
//!
//! 1. **Ties resolve to the last position in the window** (`rindex`), not the
//!    first. The reference's own comment says so twice: "we now want the last
//!    occurrence of the smallest kmer not the first anymore".
//! 2. The rescan only happens when the outgoing k-mer *was* the minimizer *and*
//!    the minimizer has actually left the window (`minimizer_pos < i - w`).
//!    Comparing by value alone would rescan too eagerly.
//! 3. **`minimizer_pos` is not updated when a new minimum arrives mid-window.**
//!    The `elif new_kmer < curr_min` branch sets `curr_min` and emits, but leaves
//!    `minimizer_pos` pointing at the *previous* minimizer's position, which then
//!    feeds condition 2. isONcorrect's port of the same function adds
//!    `minimizer_pos = i` there, with no comment and no counterpart in either
//!    Python. Measured before copying it: the two agree on **3 000 of 3 000**
//!    random reads at `k=9, w=20`, so the extra line is unobservable rather than
//!    wrong. This port still omits it, because matching the reference needs no
//!    justification and diverging does.
//! 4. The loop is `range(w+1, len(seq) - k)`, whose upper bound is **exclusive**
//!    --- so the final k-mer at `len(seq) - k` is never examined. That is an
//!    off-by-one in the reference, and it is part of the contract.

use std::collections::VecDeque;

/// `seq[i:i+k]` with Python's clamping slice semantics.
///
/// Python never panics on an out-of-range slice, it truncates. Short reads
/// therefore compare truncated k-mers, and reproducing that matters more than it
/// looks: the comparison is lexicographic, and a truncated k-mer sorts before any
/// string it is a prefix of.
fn kmer(seq: &[u8], i: usize, k: usize) -> &[u8] {
    let start = i.min(seq.len());
    let end = i.saturating_add(k).min(seq.len());
    &seq[start..end]
}

/// `rindex(lst, value)`: the **last** index of `value`, not the first.
fn rindex(window: &VecDeque<&[u8]>, value: &[u8]) -> usize {
    window
        .iter()
        .rposition(|x| *x == value)
        .expect("value came from this window")
}

/// Lexicographic minimizers as `(kmer, position)`, in the reference's order.
pub fn minimizers_lex(seq: &[u8], k: usize, w_size: usize) -> Vec<(Vec<u8>, usize)> {
    extrema(seq, k, w_size, true)
}

/// `get_kmer_maximizers`: the same walk with the comparison reversed.
///
/// Reachable only via `hash_fcn == "rev_lex"`, which `main` never sets. Ported
/// because it is fifteen lines and because leaving one of a pair out is how the
/// two drift apart.
pub fn maximizers_lex(seq: &[u8], k: usize, w_size: usize) -> Vec<(Vec<u8>, usize)> {
    extrema(seq, k, w_size, false)
}

/// Both variants: `min` when `smallest`, `max` otherwise.
fn extrema(seq: &[u8], k: usize, w_size: usize, smallest: bool) -> Vec<(Vec<u8>, usize)> {
    let mut out = Vec::new();
    if k == 0 || w_size < k {
        return out;
    }
    let w = w_size - k;

    // `deque([seq[i:i+k] for i in range(w + 1)])` --- k-mer starts 0..=w.
    let mut window: VecDeque<&[u8]> = (0..=w).map(|i| kmer(seq, i, k)).collect();
    let pick = |win: &VecDeque<&[u8]>| -> Option<Vec<u8>> {
        if smallest {
            win.iter().min().map(|x| x.to_vec())
        } else {
            win.iter().max().map(|x| x.to_vec())
        }
    };
    let Some(mut curr_min) = pick(&window) else {
        return out;
    };
    let mut minimizer_pos = rindex(&window, &curr_min);
    out.push((kmer(seq, minimizer_pos, k).to_vec(), minimizer_pos));

    // Python: `range(w + 1, len(seq) - k)`. Upper bound exclusive --- detail 4.
    let upper = seq.len().saturating_sub(k);
    for i in (w + 1)..upper {
        let new_kmer = kmer(seq, i, k);
        let discarded = window.pop_front().expect("window is never empty");
        window.push_back(new_kmer);

        // `i - w` is the leftmost start position now in the window.
        if discarded == curr_min.as_slice() && minimizer_pos < i - w {
            curr_min = pick(&window).expect("window is never empty");
            minimizer_pos = rindex(&window, &curr_min) + i - w;
            out.push((kmer(seq, minimizer_pos, k).to_vec(), minimizer_pos));
        } else {
            let better = if smallest {
                new_kmer < curr_min.as_slice()
            } else {
                new_kmer > curr_min.as_slice()
            };
            if better {
                curr_min = new_kmer.to_vec();
                // No `minimizer_pos = i` here --- detail 3.
                out.push((kmer(seq, i, k).to_vec(), i));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lex(seq: &str, k: usize, w: usize) -> Vec<(String, usize)> {
        minimizers_lex(seq.as_bytes(), k, w)
            .into_iter()
            .map(|(m, p)| (String::from_utf8(m).unwrap(), p))
            .collect()
    }

    #[test]
    fn ties_go_to_the_last_position_in_the_window() {
        // Two identical minimal k-mers in the first window: `rindex` picks the
        // later one. With `index` this would be 0.
        let got = lex("ACACACGGGG", 2, 6);
        assert_eq!(got[0].0, "AC");
        assert_eq!(
            got[0].1, 4,
            "the last AC start in window 0..=4, not the first"
        );
    }

    #[test]
    fn the_final_kmer_is_never_examined() {
        // `range(w + 1, len(seq) - k)` excludes `len(seq) - k`, so a minimal
        // k-mer sitting exactly there is never emitted. Reference off-by-one,
        // reproduced.
        let seq = "TTTTTTTTAA"; // len 10, k 2 -> last k-mer start is 8
        let got = lex(seq, 2, 4);
        assert!(
            !got.iter().any(|(_, p)| *p == 8),
            "position 8 is out of the loop's range, got {got:?}"
        );
    }

    #[test]
    fn matches_the_reference_on_a_worked_example() {
        // Read off `get_kmer_minimizers` directly, not derived. Exercises both
        // emission branches and the stale-`minimizer_pos` path of detail 3.
        let got = lex("GGGGACCCCAAAAGGGGTTTT", 3, 6);
        let want: Vec<(String, usize)> = vec![
            ("GAC".into(), 3),
            ("ACC".into(), 4),
            ("CAA".into(), 8),
            ("AAA".into(), 9),
            ("AAA".into(), 10),
            ("AAG".into(), 11),
            ("AGG".into(), 12),
            ("GGG".into(), 14),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn a_sequence_shorter_than_k_yields_one_empty_kmer_rather_than_nothing() {
        // Not the obvious answer, and not a guess --- read off the reference:
        //
        //   get_kmer_minimizers("AC", 5, 8) -> [('', 3)]
        //   get_kmer_minimizers("",   3, 5) -> [('', 2)]
        //
        // Python's slice truncates rather than raising, so the initial window is
        // all-empty-or-truncated strings, `min` picks the empty one, and the
        // function emits it. An implementation that returned nothing here would
        // silently drop a read's only minimizer instead of a useless one, which
        // is a different graph.
        assert_eq!(lex("AC", 5, 8), vec![(String::new(), 3)]);
        assert_eq!(lex("", 3, 5), vec![(String::new(), 2)]);
        // And a read shorter than the window but longer than k still works.
        assert_eq!(lex("ACGT", 3, 5), vec![("ACG".to_string(), 0)]);
    }

    #[test]
    fn maximizers_are_the_same_walk_with_the_comparison_reversed() {
        let seq = b"ACGTACGTACGT";
        let mins = minimizers_lex(seq, 3, 6);
        let maxs = maximizers_lex(seq, 3, 6);
        assert!(!mins.is_empty() && !maxs.is_empty());
        assert_ne!(mins, maxs, "min and max should not agree on this input");
        // The first emission is the window extremum either way.
        assert!(mins[0].0 <= maxs[0].0);
    }
}
