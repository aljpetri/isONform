//! Command-line parity with the Python implementation.
//!
//! Flag names, defaults, validation order, stderr text and exit codes are
//! copied from `main::__main__` and `isONform_parallel::__main__`. Where this
//! file and the reference disagree, the reference is right — see `PORTING.md`.
//!
//! # What "parity" covers, and what it does not
//!
//! Same choice as the isONcorrect port: parity means **flag names, defaults,
//! validation order, stderr text and exit codes**, asserted by the tests below
//! and by `bench/equivalence.sh`. It does **not** mean reproducing argparse's
//! `--help` layout byte for byte — clap's formatter wraps and orders
//! differently, and nothing downstream parses help text. The reference's help
//! output is recorded in `bench/golden/cli/` so the *content* can be diffed by
//! eye when a flag changes.
//!
//! # Traps in the reference this file reproduces on purpose
//!
//! * `--xmin` is silently raised to `2 * k` **before** the no-arguments check,
//!   so `main` with no arguments writes `xmin set to 40` to stderr and *then*
//!   prints help. Order matters and is tested.
//! * `main`'s `--xmin`/`--xmax` help strings are swapped ("Upper interval
//!   length" documents `--xmin`). `isONform_parallel` has them the right way
//!   round. Reproduced; fixing it is its own commit.
//! * `--parallel` is `type=bool` in argparse, which applies `bool()` to the
//!   *string*, so `--parallel False` is **true**. See [`ParallelFlag`].
//! * The window-size message says "smaller than 100" but the check admits
//!   `--w 100`.
//! * `isONform_parallel` prints `len(sys.argv)` to stdout before doing anything
//!   else. That is a debug leftover and it is in the observable contract; see
//!   `PORTING.md` finding 4.
//!
//! Four of `main`'s flags do nothing in the reference (`--T`,
//! `--disable_numpy`, `--exact`, `--exact_instance_limit` — `PORTING.md`
//! finding 3). They are accepted here with the same names and defaults, because
//! existing pipeline scripts pass them — `isONform_parallel` itself passes
//! `--exact` on every child invocation — and rejecting them would break those
//! scripts for no gain. They are marked in the struct so nothing downstream
//! quietly starts depending on them.

use std::ffi::OsString;
use std::path::PathBuf;

use clap::{ArgAction, Parser};

/// The version string both entry points advertise. argparse formats it as
/// `%(prog)s 0.3.9`, so `--version` prints `main 0.3.9` /
/// `isONform_parallel 0.3.9`.
pub const VERSION: &str = "0.3.9";

/// Exit code the reference uses for its own validation failures (`sys.exit(1)`).
pub const EXIT_VALIDATION: i32 = 1;

/// Exit code argparse uses for an unrecognised or malformed argument.
pub const EXIT_USAGE: i32 = 2;

/// The reference's window-size diagnostic, verbatim. Note "smaller than 100"
/// where the check is `100 < w` — the message is off by one, not the check.
pub const MSG_BAD_WINDOW: &str =
    "Please specify a window of size larger or equal to k, and smaller than 100.";

/// `--parallel`'s value, reproducing argparse `type=bool`.
///
/// argparse applies `type` to the raw string, and Python's `bool` on a `str` is
/// "is it non-empty". So every one of `--parallel True`, `--parallel False`,
/// `--parallel 0` and `--parallel no` enables the flag, and only
/// `--parallel ""` disables it. `isONform_parallel` passes `"True"`, which
/// works by accident.
///
/// This is `PORTING.md` finding 2. Reproduced rather than fixed, because a
/// behaviour change does not land in the same commit as the port.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParallelFlag(pub bool);

impl ParallelFlag {
    pub fn is_set(&self) -> bool {
        self.0
    }
}

impl From<&str> for ParallelFlag {
    fn from(s: &str) -> Self {
        ParallelFlag(!s.is_empty())
    }
}

impl std::str::FromStr for ParallelFlag {
    // argparse's `type=bool` cannot fail: `bool()` accepts every string.
    type Err = std::convert::Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(ParallelFlag::from(s))
    }
}

/// `main` — reconstruct isoforms for a single cluster.
#[derive(Debug, Clone, PartialEq, Parser)]
#[command(
    name = "main",
    version = VERSION,
    about = "De novo error correction of long-read transcriptome reads",
    disable_help_subcommand = true,
    // argparse accepts any unambiguous prefix of a long option (`allow_abbrev`
    // defaults to true), and it is measurable: `main --comp` sets
    // `--compression`, `main --disable_num` sets `--disable_numpy`,
    // `isONform_parallel --fastq d` sets `--fastq_folder`. Existing pipeline
    // scripts may well rely on it, so the port matches. clap's inference
    // prefers an exact match, as argparse's does, which is what keeps `--k`
    // from colliding with `--keep_old`.
    infer_long_args = true
)]
pub struct MainArgs {
    /// Path to input fastq file with reads
    #[arg(long)]
    pub fastq: Option<PathBuf>,

    /// Kmer size
    #[arg(long, default_value_t = 20)]
    pub k: usize,

    /// Window size
    #[arg(long, default_value_t = 31)]
    pub w: usize,

    // The reference's help text for these two is swapped. Kept as-is.
    /// Upper interval length
    #[arg(long, default_value_t = 18)]
    pub xmin: usize,

    /// Lower interval length
    #[arg(long, default_value_t = 80)]
    pub xmax: usize,

    /// Minimum fraction keeping substitution
    ///
    /// INERT in the reference: `args.T` is read nowhere. `PORTING.md` finding 3.
    #[arg(long = "T", default_value_t = 0.1)]
    pub t_threshold: f64,

    /// Get exact solution for WIS for every read
    ///
    /// INERT in the reference: its only effect is to reset
    /// `previously_corrected_regions`, which is never populated.
    #[arg(long, action = ArgAction::SetTrue)]
    pub exact: bool,

    /// Do not require numpy to be installed
    ///
    /// INERT in the reference: numpy is never imported. The help text's "about
    /// 1.5x slower" describes nothing.
    #[arg(long = "disable_numpy", action = ArgAction::SetTrue)]
    pub disable_numpy: bool,

    /// Maximum length difference between two reads intervals for which they
    /// would still be merged
    #[arg(long = "delta_len", default_value_t = 5)]
    pub delta_len: usize,

    /// Maximum number of seqs to spoa
    #[arg(long = "max_seqs_to_spoa", default_value_t = 200)]
    pub max_seqs_to_spoa: usize,

    /// Maximum number of seqs to correct at a time (in case of large clusters)
    #[arg(long = "max_seqs", default_value_t = 1000)]
    pub max_seqs: usize,

    /// Activates slower exact mode for instance smaller than this limit
    ///
    /// INERT in the reference: its only effect is `args.exact = True`, and
    /// `--exact` is itself inert.
    #[arg(long = "exact_instance_limit", default_value_t = 0)]
    pub exact_instance_limit: usize,

    /// Set w = k + max(2*k, floor(cluster_size/1000))
    #[arg(long = "set_w_dynamically", action = ArgAction::SetTrue)]
    pub set_w_dynamically: bool,

    /// Print various developer stats
    #[arg(long, action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Use homopolymer compressed reads (deprecated)
    ///
    /// Live and load-bearing, unlike in isONcorrect where it was unreachable:
    /// output differs with it set.
    #[arg(long, action = ArgAction::SetTrue)]
    pub compression: bool,

    /// The outfolder of isONform, into which the algorithm will write the results
    #[arg(long)]
    pub outfolder: Option<PathBuf>,

    /// Cutoff parameter: maximum length difference at 3prime end, for which
    /// subisoforms are still merged into longer isoforms
    #[arg(long = "delta_iso_len_3", default_value_t = 30)]
    pub delta_iso_len_3: usize,

    /// Cutoff parameter: maximum length difference at 5prime end, for which
    /// subisoforms are still merged into longer isoforms
    #[arg(long = "delta_iso_len_5", default_value_t = 50)]
    pub delta_iso_len_5: usize,

    /// Indicates whether we run the parallelization wrapper script
    ///
    /// Takes a value, and every non-empty value is true. See [`ParallelFlag`].
    #[arg(long, default_value = "")]
    pub parallel: ParallelFlag,

    /// EXPERIMENTAL: use the slow mode for graph simplification (every bubble
    /// gets popped)
    #[arg(long, action = ArgAction::SetTrue)]
    pub slow: bool,
}

impl MainArgs {
    /// The reference's `print("ARGS", args)` line, verbatim (`main:343`).
    ///
    /// This is a debug leftover — `PORTING.md` finding 4 — but it lands on
    /// stdout on every run, so it is observable contract and the port emits it.
    /// It also pays for itself: one line carrying every resolved argument is a
    /// complete differential oracle for argument parsing, so
    /// `bench/equivalence.sh` catches a wrong default or a mis-parsed value
    /// without a case per flag.
    ///
    /// Field order is argparse's `add_argument` order, which is what
    /// `Namespace.__repr__` iterates. `--version` and `-h` do not appear.
    ///
    /// Reproduced faithfully includes the odd parts:
    ///
    /// * `--fastq` defaults to `False`, not `None`, because the reference
    ///   declares `default=False` on a `type=str` argument. `--outfolder`
    ///   really does default to `None`.
    /// * `--T` is a float and prints as Python `repr(float)`, so an integral
    ///   value is `1.0` rather than `1`.
    /// * `--parallel` prints as a genuine `True`/`False`, since argparse stores
    ///   the result of `bool(str)`.
    ///
    /// When the print is deleted upstream, this method and the call in
    /// `bin/main.rs` go with it.
    pub fn args_line(&self) -> String {
        format!(
            "ARGS Namespace(fastq={}, k={}, w={}, xmin={}, xmax={}, T={}, exact={}, \
             disable_numpy={}, delta_len={}, max_seqs_to_spoa={}, max_seqs={}, \
             exact_instance_limit={}, set_w_dynamically={}, verbose={}, compression={}, \
             outfolder={}, delta_iso_len_3={}, delta_iso_len_5={}, parallel={}, slow={})",
            // argparse's declared default here is `False`, not `None`.
            opt_path_or(self.fastq.as_deref(), "False"),
            self.k,
            self.w,
            self.xmin,
            self.xmax,
            py_float(self.t_threshold),
            py_bool(self.exact),
            py_bool(self.disable_numpy),
            self.delta_len,
            self.max_seqs_to_spoa,
            self.max_seqs,
            self.exact_instance_limit,
            py_bool(self.set_w_dynamically),
            py_bool(self.verbose),
            py_bool(self.compression),
            opt_path_or(self.outfolder.as_deref(), "None"),
            self.delta_iso_len_3,
            self.delta_iso_len_5,
            py_bool(self.parallel.is_set()),
            py_bool(self.slow),
        )
    }
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

/// `repr(float)`. Rust prints `1` where Python prints `1.0`; everything else
/// agrees for the values a CLI produces.
fn py_float(v: f64) -> String {
    if v.fract() == 0.0 && v.is_finite() {
        format!("{v:.1}")
    } else {
        format!("{v}")
    }
}

/// `repr(str)` for a path, or the argparse default when absent.
///
/// Python prefers single quotes and switches to double quotes for a string
/// containing a single quote but no double quote. Paths like that are pathological
/// and the reference would be the odd one out either way, so only the single-quote
/// form is reproduced; a path containing a quote is the one case where this line
/// may differ.
fn opt_path_or(p: Option<&std::path::Path>, absent: &str) -> String {
    match p {
        Some(p) => format!("'{}'", p.display()),
        None => absent.to_string(),
    }
}

/// `isONform_parallel` — fan out over a folder of per-cluster fastqs.
#[derive(Debug, Clone, PartialEq, Parser)]
#[command(
    name = "isONform_parallel",
    version = VERSION,
    about = "De novo reconstruction of long-read transcriptome reads",
    disable_help_subcommand = true,
    // See the note on `MainArgs`.
    infer_long_args = true
)]
pub struct ParallelArgs {
    /// Path to input fastq folder with reads in clusters
    #[arg(long = "fastq_folder")]
    pub fastq_folder: Option<PathBuf>,

    /// Number of cores allocated for clustering
    #[arg(long = "t", default_value_t = 8)]
    pub nr_cores: usize,

    /// Kmer size
    #[arg(long, default_value_t = 20)]
    pub k: usize,

    /// Window size
    #[arg(long, default_value_t = 31)]
    pub w: usize,

    // Correct here, unlike in `main`.
    /// Lower interval length
    #[arg(long, default_value_t = 18)]
    pub xmin: usize,

    /// Upper interval length
    #[arg(long, default_value_t = 80)]
    pub xmax: usize,

    /// Do exact correction for clusters under this threshold
    ///
    /// Forwarded to the child, where it is inert (`PORTING.md` finding 3). Note
    /// the default differs from `main`'s, which is 0.
    #[arg(long = "exact_instance_limit", default_value_t = 50)]
    pub exact_instance_limit: usize,

    /// Do not recompute previous results if corrected_reads.fq is found and has
    /// the same number of reads as the input file
    #[arg(long = "keep_old", action = ArgAction::SetTrue)]
    pub keep_old: bool,

    /// Set w = k + max(2*k, floor(cluster_size/1000))
    #[arg(long = "set_w_dynamically", action = ArgAction::SetTrue)]
    pub set_w_dynamically: bool,

    /// Maximum number of seqs to correct at a time (in case of large clusters)
    ///
    /// Reaches file splitting only: it is commented out of the child command
    /// line at `isONform_parallel:79`.
    #[arg(long = "max_seqs", default_value_t = 1000)]
    pub max_seqs: usize,

    /// Process reads per batch (of max_seqs sequences) instead of per cluster
    #[arg(long = "split_wrt_batches", action = ArgAction::SetTrue)]
    pub split_wrt_batches: bool,

    /// Outfolder with all corrected reads
    #[arg(long)]
    pub outfolder: Option<PathBuf>,

    /// Maximum length difference between two reads intervals for which they
    /// would still be merged
    #[arg(long = "delta_len", default_value_t = 5)]
    pub delta_len: usize,

    /// Diversity rate used to compare sequences
    ///
    /// Reaches batch merging only — never forwarded to the child, which
    /// hard-codes `delta = 0.15` (`main:560`). `PORTING.md` deferred
    /// improvements.
    #[arg(long, default_value_t = 0.1)]
    pub delta: f64,

    /// Maximum number of seqs to spoa
    ///
    /// Reaches batch merging only; not forwarded to the child.
    #[arg(long = "max_seqs_to_spoa", default_value_t = 200)]
    pub max_seqs_to_spoa: usize,

    /// Print various developer stats
    #[arg(long, action = ArgAction::SetTrue)]
    pub verbose: bool,

    /// Cutoff parameter: abundance of reads that have to support an isoform to
    /// show in results
    #[arg(long = "iso_abundance", default_value_t = 5)]
    pub iso_abundance: usize,

    /// Cutoff parameter: maximum length difference at 3prime end, for which
    /// subisoforms are still merged into longer isoforms
    #[arg(long = "delta_iso_len_3", default_value_t = 30)]
    pub delta_iso_len_3: usize,

    /// Cutoff parameter: maximum length difference at 5prime end, for which
    /// subisoforms are still merged into longer isoforms
    #[arg(long = "delta_iso_len_5", default_value_t = 50)]
    pub delta_iso_len_5: usize,

    /// OPTIONAL: absolute path to a custom folder in which to store temporary
    /// files
    #[arg(long)]
    pub tmpdir: Option<PathBuf>,

    /// Output the final transcriptome as fastq rather than fasta
    #[arg(long = "write_fastq", action = ArgAction::SetTrue)]
    pub write_fastq: bool,
}

/// What the caller should do after argument handling.
#[derive(Debug, PartialEq, Eq)]
pub enum Action<T> {
    /// Run with these arguments.
    Run(Box<T>),
    /// Print help to stdout and exit 0. The reference reaches this when
    /// `len(sys.argv) == 1`.
    HelpAndExit,
    /// Write `message` to stderr and exit with `code`.
    Fail { message: String, code: i32 },
}

/// Anything the caller must emit before acting, in order. Keeping these out of
/// the `Action` means the validation order is expressed once, here, and the
/// tests can assert on it without capturing process output.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Emissions {
    pub stdout: Vec<String>,
    pub stderr: Vec<String>,
}

/// Reproduce `main`'s post-parse handling, in the reference's order.
///
/// From `main:628-639`:
///
/// 1. raise `xmin` to `2 * k` if below it, announcing it on stderr;
/// 2. **then** if there were no arguments at all, print help and exit 0;
/// 3. create `--outfolder` if it does not exist (the caller does this — it is
///    the only step with a side effect on the filesystem);
/// 4. reject a window outside `[k, 100]` with exit 1.
///
/// Steps 1 and 2 really are in that order, which is why `main` with no
/// arguments writes `xmin set to 40` before its own help text.
pub fn resolve_main(mut args: MainArgs, argc: usize) -> (Emissions, Action<MainArgs>) {
    let mut out = Emissions::default();

    if args.xmin < 2 * args.k {
        args.xmin = 2 * args.k;
        out.stderr.push(format!("xmin set to {}", args.xmin));
    }

    if argc == 1 {
        return (out, Action::HelpAndExit);
    }

    if 100 < args.w || args.w < args.k {
        return (
            out,
            Action::Fail {
                message: MSG_BAD_WINDOW.to_string(),
                code: EXIT_VALIDATION,
            },
        );
    }

    (out, Action::Run(Box::new(args)))
}

/// Reproduce `isONform_parallel`'s post-parse handling, in the reference's
/// order.
///
/// From `isONform_parallel:369-378`. The first thing it does is
/// `print(len(sys.argv))` — a debug leftover on stdout, reproduced because it
/// is observable. There is no window-size check here.
pub fn resolve_parallel(args: ParallelArgs, argc: usize) -> (Emissions, Action<ParallelArgs>) {
    let mut out = Emissions::default();
    out.stdout.push(argc.to_string());

    if argc == 1 {
        return (out, Action::HelpAndExit);
    }

    (out, Action::Run(Box::new(args)))
}

/// `argc` as the reference sees it: `len(sys.argv)`, so the program name counts.
pub fn argc_from<I, T>(argv: I) -> usize
where
    I: IntoIterator<Item = T>,
    T: Into<OsString>,
{
    argv.into_iter().count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_main(extra: &[&str]) -> (Emissions, Action<MainArgs>) {
        let mut argv = vec!["main"];
        argv.extend_from_slice(extra);
        let argc = argv.len();
        let args = MainArgs::try_parse_from(&argv).expect("parses");
        resolve_main(args, argc)
    }

    fn parse_parallel(extra: &[&str]) -> (Emissions, Action<ParallelArgs>) {
        let mut argv = vec!["isONform_parallel"];
        argv.extend_from_slice(extra);
        let argc = argv.len();
        let args = ParallelArgs::try_parse_from(&argv).expect("parses");
        resolve_parallel(args, argc)
    }

    fn run_main(extra: &[&str]) -> MainArgs {
        match parse_main(extra).1 {
            Action::Run(a) => *a,
            other => panic!("expected Run, got {other:?}"),
        }
    }

    // --- defaults, checked against `python main --help` -------------------

    #[test]
    fn main_defaults_match_the_reference() {
        let a = run_main(&["--fastq", "x.fq"]);
        assert_eq!(a.k, 20);
        assert_eq!(a.w, 31);
        // 18 in argparse, but resolve_main raises it to 2*k before anything
        // else can see it. The *declared* default is asserted separately below.
        assert_eq!(a.xmin, 40);
        assert_eq!(a.xmax, 80);
        assert_eq!(a.t_threshold, 0.1);
        assert_eq!(a.delta_len, 5);
        assert_eq!(a.max_seqs_to_spoa, 200);
        assert_eq!(a.max_seqs, 1000);
        assert_eq!(a.exact_instance_limit, 0);
        assert_eq!(a.delta_iso_len_3, 30);
        assert_eq!(a.delta_iso_len_5, 50);
        assert!(!a.exact);
        assert!(!a.disable_numpy);
        assert!(!a.set_w_dynamically);
        assert!(!a.verbose);
        assert!(!a.compression);
        assert!(!a.slow);
        assert!(!a.parallel.is_set());
        assert_eq!(a.outfolder, None);
    }

    #[test]
    fn main_declared_xmin_default_is_18() {
        // Before validation. argparse's own default, as printed by --help.
        let a = MainArgs::try_parse_from(["main", "--fastq", "x.fq"]).unwrap();
        assert_eq!(a.xmin, 18);
    }

    #[test]
    fn parallel_defaults_match_the_reference() {
        let a = match parse_parallel(&["--fastq_folder", "d"]).1 {
            Action::Run(a) => *a,
            other => panic!("expected Run, got {other:?}"),
        };
        assert_eq!(a.nr_cores, 8);
        assert_eq!(a.k, 20);
        assert_eq!(a.w, 31);
        // No 2*k adjustment on this entry point — 18 survives.
        assert_eq!(a.xmin, 18);
        assert_eq!(a.xmax, 80);
        // Deliberately different from `main`'s 0.
        assert_eq!(a.exact_instance_limit, 50);
        assert_eq!(a.max_seqs, 1000);
        assert_eq!(a.delta_len, 5);
        assert_eq!(a.delta, 0.1);
        assert_eq!(a.max_seqs_to_spoa, 200);
        assert_eq!(a.iso_abundance, 5);
        assert_eq!(a.delta_iso_len_3, 30);
        assert_eq!(a.delta_iso_len_5, 50);
        assert!(!a.keep_old);
        assert!(!a.set_w_dynamically);
        assert!(!a.split_wrt_batches);
        assert!(!a.verbose);
        assert!(!a.write_fastq);
        assert_eq!(a.tmpdir, None);
    }

    #[test]
    fn the_two_entry_points_disagree_about_exact_instance_limit() {
        // Not a typo in either place — the reference really does default this
        // to 0 in `main` and 50 in the driver.
        let m = run_main(&["--fastq", "x.fq"]);
        let p = match parse_parallel(&["--fastq_folder", "d"]).1 {
            Action::Run(a) => *a,
            other => panic!("{other:?}"),
        };
        assert_eq!(m.exact_instance_limit, 0);
        assert_eq!(p.exact_instance_limit, 50);
    }

    // --- validation order -------------------------------------------------

    #[test]
    fn xmin_is_raised_to_twice_k_and_announced() {
        let (em, act) = parse_main(&["--k", "25", "--fastq", "x.fq"]);
        assert_eq!(em.stderr, vec!["xmin set to 50"]);
        assert!(matches!(act, Action::Run(_)));
    }

    #[test]
    fn xmin_above_twice_k_is_left_alone_and_silent() {
        let (em, act) = parse_main(&["--k", "10", "--xmin", "60", "--fastq", "x.fq"]);
        assert!(em.stderr.is_empty());
        match act {
            Action::Run(a) => assert_eq!(a.xmin, 60),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn xmin_exactly_twice_k_is_left_alone() {
        // The check is `<`, not `<=`.
        let (em, _) = parse_main(&["--k", "20", "--xmin", "40", "--fastq", "x.fq"]);
        assert!(em.stderr.is_empty());
    }

    #[test]
    fn no_arguments_announces_xmin_before_printing_help() {
        // The order is the point: the reference adjusts xmin, *then* notices
        // there were no arguments. Verified against
        // `bench/golden/cli/main.noargs.stderr`, which contains exactly this.
        let (em, act) = parse_main(&[]);
        assert_eq!(em.stderr, vec!["xmin set to 40"]);
        assert!(em.stdout.is_empty());
        assert_eq!(act, Action::HelpAndExit);
    }

    #[test]
    fn window_smaller_than_k_is_rejected_after_the_no_args_check() {
        let (em, act) = parse_main(&["--k", "20", "--w", "10", "--fastq", "x.fq"]);
        // The xmin line still comes first.
        assert_eq!(em.stderr, vec!["xmin set to 40"]);
        assert_eq!(
            act,
            Action::Fail {
                message: MSG_BAD_WINDOW.to_string(),
                code: EXIT_VALIDATION,
            }
        );
    }

    #[test]
    fn window_larger_than_100_is_rejected() {
        let (_, act) = parse_main(&["--w", "101", "--fastq", "x.fq"]);
        assert!(matches!(act, Action::Fail { code: 1, .. }));
    }

    #[test]
    fn window_of_exactly_100_is_accepted_despite_the_message() {
        // The check is `100 < w`, so 100 passes, even though the diagnostic
        // says "smaller than 100". Off by one in the message, not the check.
        let (_, act) = parse_main(&["--w", "100", "--fastq", "x.fq"]);
        assert!(matches!(act, Action::Run(_)));
    }

    #[test]
    fn window_equal_to_k_is_accepted() {
        let (_, act) = parse_main(&["--k", "31", "--w", "31", "--fastq", "x.fq"]);
        assert!(matches!(act, Action::Run(_)));
    }

    #[test]
    fn parallel_driver_has_no_window_check() {
        // `isONform_parallel` never validates w, so a window the child will
        // reject is accepted here and fails one level down.
        let (_, act) = parse_parallel(&["--k", "20", "--w", "5", "--fastq_folder", "d"]);
        assert!(matches!(act, Action::Run(_)));
    }

    // --- the driver's stray stdout ---------------------------------------

    #[test]
    fn parallel_prints_argc_to_stdout_first() {
        // `print(len(sys.argv))` at isONform_parallel:370. Observable, so it
        // is contract. See PORTING.md finding 4.
        let (em, act) = parse_parallel(&["--fastq_folder", "d"]);
        assert_eq!(em.stdout, vec!["3"]);
        assert!(matches!(act, Action::Run(_)));
    }

    #[test]
    fn parallel_with_no_arguments_prints_one_then_help() {
        // Matches bench/golden/cli/isONform_parallel.noargs.stdout, whose
        // first line is `1`.
        let (em, act) = parse_parallel(&[]);
        assert_eq!(em.stdout, vec!["1"]);
        assert_eq!(act, Action::HelpAndExit);
    }

    // --- the type=bool bug ------------------------------------------------

    #[test]
    fn parallel_flag_is_true_for_every_non_empty_value() {
        // argparse `type=bool` means `bool(str)`, i.e. "non-empty". Measured
        // against the reference: `--parallel True` and `--parallel False`
        // produce byte-identical output, and both differ from the default.
        for v in ["True", "False", "0", "no", "off", "None"] {
            let a = run_main(&["--fastq", "x.fq", "--parallel", v]);
            assert!(a.parallel.is_set(), "--parallel {v} should be true");
        }
    }

    #[test]
    fn parallel_flag_is_false_only_for_the_empty_string() {
        let a = run_main(&["--fastq", "x.fq", "--parallel", ""]);
        assert!(!a.parallel.is_set());
    }

    #[test]
    fn parallel_flag_requires_a_value() {
        // It is not a store_true: bare `--parallel` is a usage error in
        // argparse too.
        assert!(MainArgs::try_parse_from(["main", "--parallel"]).is_err());
    }

    // --- flag names -------------------------------------------------------

    #[test]
    fn multi_word_flags_keep_their_underscores() {
        // clap rewrites `field_name` to `--field-name` unless told otherwise,
        // which would silently break every pipeline script. Method point 1.
        for f in [
            "--disable_numpy",
            "--delta_len",
            "--max_seqs_to_spoa",
            "--max_seqs",
            "--exact_instance_limit",
            "--set_w_dynamically",
            "--delta_iso_len_3",
            "--delta_iso_len_5",
        ] {
            let mut argv = vec!["main", "--fastq", "x.fq", f];
            // The valued ones need a value; the switches must reject one.
            if f.starts_with("--delta") || f.contains("max_seqs") || f.contains("limit") {
                argv.push("7");
            }
            assert!(
                MainArgs::try_parse_from(&argv).is_ok(),
                "{f} not accepted as spelled"
            );
        }
    }

    #[test]
    fn parallel_multi_word_flags_keep_their_underscores() {
        for (f, val) in [
            ("--fastq_folder", Some("d")),
            ("--exact_instance_limit", Some("7")),
            ("--keep_old", None),
            ("--set_w_dynamically", None),
            ("--max_seqs", Some("7")),
            ("--split_wrt_batches", None),
            ("--delta_len", Some("7")),
            ("--max_seqs_to_spoa", Some("7")),
            ("--iso_abundance", Some("7")),
            ("--delta_iso_len_3", Some("7")),
            ("--delta_iso_len_5", Some("7")),
            ("--write_fastq", None),
        ] {
            let mut argv = vec!["isONform_parallel", f];
            if let Some(v) = val {
                argv.push(v);
            }
            assert!(
                ParallelArgs::try_parse_from(&argv).is_ok(),
                "{f} not accepted as spelled"
            );
        }
    }

    #[test]
    fn hyphenated_spellings_are_rejected() {
        // Guards against a future `long = "..."` being dropped: if clap starts
        // accepting --delta-len, the underscore spelling may have moved.
        assert!(MainArgs::try_parse_from(["main", "--delta-len", "7"]).is_err());
        assert!(
            ParallelArgs::try_parse_from(["isONform_parallel", "--fastq-folder", "d"]).is_err()
        );
    }

    #[test]
    fn t_is_capital_on_main_and_lowercase_on_the_driver() {
        // `--T` is a float threshold on `main`; `--t` is the core count on the
        // driver. Two different flags, distinguished only by case.
        assert!(MainArgs::try_parse_from(["main", "--T", "0.5"]).is_ok());
        assert!(ParallelArgs::try_parse_from(["isONform_parallel", "--t", "4"]).is_ok());
        // And neither accepts the other's spelling.
        assert!(MainArgs::try_parse_from(["main", "--t", "4"]).is_err());
        assert!(ParallelArgs::try_parse_from(["isONform_parallel", "--T", "0.5"]).is_err());
    }

    #[test]
    fn driver_only_flags_are_absent_from_main() {
        for f in [
            "--keep_old",
            "--split_wrt_batches",
            "--tmpdir",
            "--iso_abundance",
            "--write_fastq",
        ] {
            assert!(
                MainArgs::try_parse_from(["main", f, "x"]).is_err(),
                "{f} should not exist on main"
            );
        }
    }

    #[test]
    fn main_only_flags_are_absent_from_the_driver() {
        // Note `--fastq` is deliberately *not* in this list: with prefix
        // inference it resolves to `--fastq_folder` on the driver, exactly as
        // argparse resolves it. See `argparse_prefix_abbreviations_are_accepted`.
        for f in [
            "--exact",
            "--disable_numpy",
            "--compression",
            "--slow",
            "--parallel",
        ] {
            assert!(
                ParallelArgs::try_parse_from(["isONform_parallel", f, "x"]).is_err(),
                "{f} should not exist on isONform_parallel"
            );
        }
    }

    // --- argparse's prefix abbreviation -----------------------------------

    #[test]
    fn argparse_prefix_abbreviations_are_accepted() {
        // Measured against the reference: each of these parses there.
        let a = MainArgs::try_parse_from(["main", "--comp"]).expect("--comp");
        assert!(a.compression);
        let a = MainArgs::try_parse_from(["main", "--disable_num"]).expect("--disable_num");
        assert!(a.disable_numpy);
        let a = MainArgs::try_parse_from(["main", "--set_w"]).expect("--set_w");
        assert!(a.set_w_dynamically);
        let a = MainArgs::try_parse_from(["main", "--fast", "x.fq"]).expect("--fast");
        assert_eq!(a.fastq, Some(PathBuf::from("x.fq")));

        // The driver's `--fastq` is the one that bites: it is an unambiguous
        // prefix of `--fastq_folder`, so `isONform_parallel --fastq d` silently
        // means `--fastq_folder d`. The reference does this too.
        let a =
            ParallelArgs::try_parse_from(["isONform_parallel", "--fastq", "d"]).expect("--fastq");
        assert_eq!(a.fastq_folder, Some(PathBuf::from("d")));
    }

    #[test]
    fn an_exact_match_beats_a_longer_flag_it_prefixes() {
        // `--k` prefixes `--keep_old` on the driver, but it is also a flag in
        // its own right, and an exact match wins in argparse and in clap alike.
        let a = ParallelArgs::try_parse_from(["isONform_parallel", "--k", "15"]).expect("--k");
        assert_eq!(a.k, 15);
        assert!(!a.keep_old);
    }

    #[test]
    fn ambiguous_prefixes_are_still_rejected() {
        // `--delta` is exact on the driver, but on `main` there is no `--delta`
        // and it prefixes --delta_len, --delta_iso_len_3 and --delta_iso_len_5.
        assert!(MainArgs::try_parse_from(["main", "--delta", "7"]).is_err());
        // `--m` prefixes --max_seqs and --max_seqs_to_spoa.
        assert!(MainArgs::try_parse_from(["main", "--m", "7"]).is_err());
    }

    #[test]
    fn unknown_flags_are_a_usage_error() {
        let err = MainArgs::try_parse_from(["main", "--bogus"]).unwrap_err();
        assert_eq!(err.exit_code(), EXIT_USAGE);
    }

    #[test]
    fn version_matches_setup_py() {
        assert_eq!(VERSION, "0.3.9");
    }

    // --- the ARGS Namespace line -----------------------------------------
    //
    // Expected strings are captured from the reference verbatim:
    //   PYTHONHASHSEED=0 python main <args> 2>/dev/null | head -1

    #[test]
    fn args_line_matches_the_reference_for_defaults() {
        let a = run_main(&["--k", "25"]);
        assert_eq!(
            a.args_line(),
            "ARGS Namespace(fastq=False, k=25, w=31, xmin=50, xmax=80, T=0.1, exact=False, \
             disable_numpy=False, delta_len=5, max_seqs_to_spoa=200, max_seqs=1000, \
             exact_instance_limit=0, set_w_dynamically=False, verbose=False, compression=False, \
             outfolder=None, delta_iso_len_3=30, delta_iso_len_5=50, parallel=False, slow=False)"
        );
    }

    #[test]
    fn args_line_shows_parallel_as_a_real_bool() {
        let a = run_main(&["--parallel", "True"]);
        assert_eq!(
            a.args_line(),
            "ARGS Namespace(fastq=False, k=20, w=31, xmin=40, xmax=80, T=0.1, exact=False, \
             disable_numpy=False, delta_len=5, max_seqs_to_spoa=200, max_seqs=1000, \
             exact_instance_limit=0, set_w_dynamically=False, verbose=False, compression=False, \
             outfolder=None, delta_iso_len_3=30, delta_iso_len_5=50, parallel=True, slow=False)"
        );
    }

    #[test]
    fn args_line_matches_the_reference_with_everything_set() {
        let a = run_main(&[
            "--fastq",
            "/x.fq",
            "--outfolder",
            "/tmp/of",
            "--T",
            "0.5",
            "--slow",
            "--compression",
            "--verbose",
            "--exact",
            "--disable_numpy",
            "--set_w_dynamically",
            "--max_seqs",
            "7",
            "--max_seqs_to_spoa",
            "8",
            "--exact_instance_limit",
            "9",
            "--delta_len",
            "3",
            "--delta_iso_len_3",
            "11",
            "--delta_iso_len_5",
            "12",
            "--xmax",
            "99",
        ]);
        assert_eq!(
            a.args_line(),
            "ARGS Namespace(fastq='/x.fq', k=20, w=31, xmin=40, xmax=99, T=0.5, exact=True, \
             disable_numpy=True, delta_len=3, max_seqs_to_spoa=8, max_seqs=7, \
             exact_instance_limit=9, set_w_dynamically=True, verbose=True, compression=True, \
             outfolder='/tmp/of', delta_iso_len_3=11, delta_iso_len_5=12, parallel=False, slow=True)"
        );
    }

    #[test]
    fn integral_floats_print_python_style() {
        // repr(1.0) is "1.0" in Python; Rust's Display gives "1".
        assert_eq!(py_float(1.0), "1.0");
        assert_eq!(py_float(0.1), "0.1");
        assert_eq!(py_float(0.5), "0.5");
        let a = run_main(&["--T", "1"]);
        assert!(a.args_line().contains("T=1.0,"), "{}", a.args_line());
    }
}
