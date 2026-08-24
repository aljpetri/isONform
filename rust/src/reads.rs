//! Read preparation: `remove_read_polyA_ends` and `batch` from `main:20`.
//!
//! Everything here runs before the front half proper. Reads are loaded, their
//! poly-A/poly-T tails are shortened, and they are cut into batches of
//! `--max_seqs`; each batch then gets its own minimizers, anchors, intervals and
//! graph.

/// `remove_read_polyA_ends`: shorten a long homopolymer run of A or T near the
/// end of the read to `to_len` bases.
///
/// Only the last `min(len / 2, 100)` bases are examined, and only `A` and `T`
/// runs longer than `threshold_len` are shortened; every other run is copied
/// through at its original length.
///
/// # The window is not what it looks like when the read is shorter than two bases
///
/// `end_length_window` is `min(len(seq) // 2, 100)`, which is **0** for a read of
/// length 0 or 1. The reference then does:
///
/// ```python
/// seq_list = [ seq[:-end_length_window] ]      # seq[:-0]  == seq[:0] == ""
/// for ch, g in itertools.groupby(seq[-end_length_window:]):  # seq[-0:] == seq
/// ```
///
/// Python's `-0` is `0`, so the "everything before the window" half becomes
/// **empty** and the window itself becomes the **whole read** — both halves flip
/// at once. The two mistakes cancel and the result is the read unchanged, which
/// is why it has never mattered. Reproduced rather than tidied, because "the
/// result happens to be right" is not the same as "the code does what it says",
/// and a future edit to either half would break the cancellation.
pub fn remove_polya_ends(seq: &[u8], threshold_len: usize, to_len: usize) -> Vec<u8> {
    let end_window = std::cmp::min(seq.len() / 2, 100);
    let mut out: Vec<u8> = Vec::with_capacity(seq.len());

    // `seq[:-end_window]` and `seq[-end_window:]`, including the `-0` behaviour.
    let (head, tail): (&[u8], &[u8]) = if end_window == 0 {
        (&[], seq)
    } else {
        seq.split_at(seq.len() - end_window)
    };
    out.extend_from_slice(head);

    // `itertools.groupby` over the window: consecutive runs of one character.
    let mut i = 0usize;
    while i < tail.len() {
        let ch = tail[i];
        let mut run = 1usize;
        while i + run < tail.len() && tail[i + run] == ch {
            run += 1;
        }
        let n = if run > threshold_len && (ch == b'A' || ch == b'T') {
            to_len
        } else {
            run
        };
        out.extend(std::iter::repeat_n(ch, n));
        i += run;
    }
    out
}

/// `batch`: cut the reads into groups of `size`, preserving order.
///
/// # This drops a read, and the last batch is short by one
///
/// The reference walks `enumerate(dictionary.items())` and starts a new
/// sub-dictionary whenever `i % size == 0`, putting the current read into the
/// *new* batch. That part is fine. What is not:
///
/// ```python
/// if i / size != 0:
///     sub_dict[acc] = seq        # re-inserts the LAST read, already present
///     batches.append(sub_dict)
/// elif len(dictionary) == 1:
///     batches.append(sub_dict)
/// ```
///
/// `i` is the loop variable surviving the loop, so this is "if the final index
/// over `size` is non-zero" — true whenever there is more than one read, and
/// using true division so it is false only when `i == 0`. The re-insert is a
/// no-op because the read is already in `sub_dict`. So the behaviour is: append
/// the final partial batch, unless there was exactly one read, in which case the
/// `elif` catches it.
///
/// The one real consequence: with **exactly `size` reads plus one**, the trailing
/// batch holds a single read, and a batch of one read is skipped entirely by the
/// driver (`len(reads) == 1` writes it to the skip file). That is the reference's
/// behaviour and it is reproduced.
pub fn batch<T: Clone>(items: &[T], size: usize) -> Vec<Vec<T>> {
    let mut batches: Vec<Vec<T>> = Vec::new();
    if items.is_empty() {
        // `i` is unbound in the reference, so this raises `UnboundLocalError`.
        // Returning nothing is the only sane reading; the driver never calls it
        // with an empty read set because an empty fastq stops earlier.
        return batches;
    }
    if size == 0 {
        // `i % 0` raises ZeroDivisionError. Not reachable: `--max_seqs` defaults
        // to 1000 and is an int the CLI does not validate, so this is the honest
        // boundary rather than a panic.
        return vec![items.to_vec()];
    }
    let mut sub: Vec<T> = Vec::new();
    for (i, item) in items.iter().enumerate() {
        if i > 0 && i % size == 0 {
            batches.push(std::mem::take(&mut sub));
        }
        sub.push(item.clone());
    }
    // `if i / size != 0` with `i` the final index: true unless the final index is
    // 0, i.e. unless there was exactly one read. The `elif len(...) == 1` then
    // catches that case, so in both branches the trailing batch is appended.
    batches.push(sub);
    batches
}

#[cfg(test)]
mod tests {
    use super::*;

    fn polya(s: &str) -> String {
        String::from_utf8(remove_polya_ends(s.as_bytes(), 12, 1)).unwrap()
    }

    #[test]
    fn a_long_polya_tail_is_shortened_to_one_base() {
        // Read off the reference: 'ACGT' + 'A'*30 -> 'ACGT' + 'A'*14.
        // The window is the last 17 bases (len 34 / 2), all A, so that run
        // collapses to 1 --- but the 13 A's *before* the window are untouched.
        assert_eq!(
            polya(&format!("ACGT{}", "A".repeat(30))),
            "ACGTAAAAAAAAAAAAAA"
        );
    }

    #[test]
    fn a_long_polyt_tail_is_shortened_too_but_other_bases_are_not() {
        assert_eq!(
            polya(&format!("ACGT{}", "T".repeat(30))),
            "ACGTTTTTTTTTTTTTTT"
        );
        let g = format!("ACGT{}", "G".repeat(30));
        assert_eq!(polya(&g), g, "only A and T runs are shortened");
    }

    #[test]
    fn a_run_at_or_below_the_threshold_is_left_alone() {
        // 11 A's in the window, threshold 12: `h_len > threshold_len` is false.
        let s = "ACGTACGTACAAAAAAAAAAAAA";
        assert_eq!(polya(s), s);
    }

    #[test]
    fn a_read_shorter_than_two_bases_survives_the_minus_zero_window() {
        // `end_length_window` is 0, so `seq[:-0]` is empty and `seq[-0:]` is the
        // whole read. Both halves flip and the mistakes cancel. Values read off
        // the reference, not derived.
        assert_eq!(polya(""), "");
        assert_eq!(polya("A"), "A");
        // Length 2 is the first case where the window is a real suffix.
        assert_eq!(polya("AC"), "AC");
        assert_eq!(polya("ACGT"), "ACGT");
        assert_eq!(polya("AAAA"), "AAAA");
    }

    #[test]
    fn batches_are_size_long_and_keep_order() {
        let items: Vec<u32> = (1..=10).collect();
        let got = batch(&items, 4);
        assert_eq!(got, vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8], vec![9, 10]]);
    }

    #[test]
    fn size_plus_one_reads_leave_a_batch_of_one() {
        // The consequence worth pinning: the driver skips a batch of one read
        // entirely, writing it to the skip file. So 5 reads at size 4 loses one.
        let items: Vec<u32> = (1..=5).collect();
        let got = batch(&items, 4);
        assert_eq!(got.len(), 2);
        assert_eq!(
            got[1],
            vec![5],
            "a lone trailing read, which the driver skips"
        );
    }

    #[test]
    fn one_read_still_produces_one_batch() {
        // The reference needs its `elif len(dictionary) == 1` for this, because
        // `i / size != 0` is false when the only index is 0.
        assert_eq!(batch(&[7u32], 4), vec![vec![7]]);
    }

    #[test]
    fn an_exact_multiple_does_not_produce_a_trailing_empty_batch() {
        let items: Vec<u32> = (1..=8).collect();
        let got = batch(&items, 4);
        assert_eq!(got, vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]]);
    }
}
