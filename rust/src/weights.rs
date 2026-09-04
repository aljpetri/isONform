//! Per-read weights, for feeding *consensuses* back through the pipeline.
//!
//! # Why
//!
//! The cross-batch merge is `O(batch pairs x isoforms per batch^2)` and proves
//! non-mergeability one alignment at a time. Running isONform's own front end
//! over the surviving consensuses instead lets the graph do the merging, which is
//! not quadratic in the batch count. Measured on `droso_deep` cluster 3: one such
//! pass folds 386 consensuses to 233 in 20.3s, against 47 merges for 55 643
//! aligned pairs in 24.1s, and it converges (233 -> 215 -> 211) on the 215 that
//! processing all 9 261 reads as a single batch produces.
//!
//! # What breaks without weights
//!
//! Support is a *count of read accessions* everywhere it matters, so on a
//! recursive pass it becomes a count of *consensuses* --- typically one to three.
//! `--iso_abundance 5` then discards nearly everything, and interval selection
//! sees a fraction of the true support. A consensus standing for 300 reads has to
//! say so.
//!
//! # How
//!
//! The weight rides on the read accession as a `:w=<n>` suffix, so it survives
//! every stage that carries accessions without threading a parameter through
//! them. Absent suffix means weight 1, which is what an ordinary read is. The
//! whole thing is off unless `ISONFORM_WEIGHTS=1`, so default behaviour is
//! untouched --- with it off, `of` and `sum_accs` are exactly the counts the
//! unweighted code computed.
//!
//! Composing passes needs no new output format: `mapping{batch}` already lists
//! each isoform's supporting accessions, and those carry their own `:w=`, so the
//! next pass's weight is their sum.

use std::cell::{Cell, RefCell};

thread_local! {
    /// On for this thread. Initialised from `ISONFORM_WEIGHTS`, and forced on by
    /// the recursive merge stage, which needs weights regardless of the flag.
    static ON: Cell<bool> = const { Cell::new(false) };
    /// Weights by `r_id` for the run *this thread* is executing.
    ///
    /// Thread-local rather than global because the parallel driver runs one
    /// cluster per thread and each has its own read set --- a process-wide table
    /// would be whichever cluster installed it first.
    static TABLE: RefCell<Vec<u32>> = const { RefCell::new(Vec::new()) };
    static INIT: Cell<bool> = const { Cell::new(false) };
}

/// Is weighting on for this thread? `ISONFORM_WEIGHTS=1`, or forced by
/// [`enable_for_thread`].
pub fn enabled() -> bool {
    INIT.with(|i| {
        if !i.get() {
            i.set(true);
            ON.with(|on| {
                on.set(matches!(
                    std::env::var("ISONFORM_WEIGHTS").ok().as_deref(),
                    Some("1") | Some("on")
                ))
            });
        }
    });
    ON.with(|on| on.get())
}

/// Turn weighting on for this thread whatever the environment says.
pub fn enable_for_thread() {
    enabled();
    ON.with(|on| on.set(true));
}

/// The weight encoded in an accession, or 1.
///
/// Parsed from the last `:w=` in the accession, so an accession that already
/// carries one from an earlier pass is replaced rather than accumulated --- the
/// suffix is the weight, not a history.
pub fn parse_acc(acc: &str) -> u32 {
    if !enabled() {
        return 1;
    }
    parse_suffix(acc)
}

/// The parse itself, independent of whether weighting is on.
fn parse_suffix(acc: &str) -> u32 {
    acc.rfind(":w=")
        .and_then(|i| acc[i + 3..].split_whitespace().next())
        .and_then(|v| v.parse::<u32>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(1)
}

/// Record the weights for the run this thread is about to execute, replacing any
/// previous one.
pub fn install(reads: &[(u32, String, Vec<u8>)]) {
    if !enabled() {
        return;
    }
    let max = reads.iter().map(|(r, _, _)| *r).max().unwrap_or(0) as usize;
    let mut t = vec![1u32; max + 1];
    for (r, acc, _) in reads {
        t[*r as usize] = parse_suffix(acc);
    }
    TABLE.with(|c| *c.borrow_mut() = t);
}

/// The weight of one read. 1 when weighting is off or the id is unknown.
pub fn of(r_id: u32) -> u32 {
    if !enabled() {
        return 1;
    }
    TABLE.with(|c| c.borrow().get(r_id as usize).copied().unwrap_or(1))
}

/// Total weight behind a set of accessions --- the support an isoform really has.
///
/// With weighting off this is `accs.len()`, which is what the unweighted code
/// compared against `--iso_abundance`.
pub fn sum_accs<S: AsRef<str>>(accs: &[S]) -> usize {
    if !enabled() {
        return accs.len();
    }
    accs.iter().map(|a| parse_suffix(a.as_ref()) as usize).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_suffix_parses_only_when_it_is_a_positive_integer() {
        assert_eq!(parse_suffix("r1"), 1, "no suffix is weight 1");
        assert_eq!(parse_suffix("r1:w=42"), 42);
        assert_eq!(parse_suffix("r1:w=0"), 1, "zero is not a weight");
        assert_eq!(parse_suffix("r1:w=x"), 1, "junk is not a weight");
        assert_eq!(parse_suffix("r1:w=2:w=9"), 9, "the last suffix wins");
    }

    #[test]
    fn forcing_it_on_makes_support_a_sum_and_leaves_other_threads_alone() {
        // A fresh thread, so the flag cannot leak from another test.
        std::thread::spawn(|| {
            assert_eq!(sum_accs(&["a:w=10", "b:w=20"]), 2, "off: a count");
            enable_for_thread();
            assert_eq!(sum_accs(&["a:w=10", "b:w=20"]), 30, "on: a sum");
            install(&[(1, "a:w=10".into(), vec![])]);
            assert_eq!(of(1), 10);
            assert_eq!(of(99), 1, "unknown ids weigh 1");
        })
        .join()
        .unwrap();
        // This thread is untouched unless the environment turned it on.
        if !enabled() {
            assert_eq!(sum_accs(&["a:w=10"]), 1);
        }
    }
}
