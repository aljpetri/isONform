//! `main` — the single-instance isoform reconstructor.
//!
//! Named `main` because that is what the Python ships (`setup.py`
//! `scripts=['isONform_parallel','main']`) and what `isONform_parallel` execs.
//!
//! The whole pipeline now runs: reads are trimmed and batched, and each batch
//! goes through intervals, graph construction, simplification and isoform
//! generation. Each stage is verified against the reference by its own oracle;
//! this binary is the wiring plus the file IO.
//!
//! # The intermediates are not Python pickles
//!
//! The reference writes `{id}batch`, `spoa{id}` and `mapping{id}` with
//! `pickle.dump`, for `isONform_parallel` to read back. This port writes the same
//! *content* in a plain text format instead --- see [`isonform::driver`] for why,
//! and for the consequence: **the Rust and Python halves cannot be mixed.** The
//! user-visible outputs are byte-compatible; these internal handoff files are not.

use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isonform::cli::{self, Action, MainArgs};

fn main() -> ExitCode {
    let argv: Vec<std::ffi::OsString> = std::env::args_os().collect();
    let argc = argv.len();

    // clap handles --help/--version itself, exiting 0, and unknown flags,
    // exiting 2 — the same codes argparse uses.
    let args = match MainArgs::try_parse_from(&argv) {
        Ok(a) => a,
        Err(e) => e.exit(),
    };

    let (emissions, action) = cli::resolve_main(args, argc);

    for line in &emissions.stdout {
        println!("{line}");
    }
    for line in &emissions.stderr {
        eprintln!("{line}");
    }

    match action {
        Action::HelpAndExit => {
            // argparse's `parser.print_help(); sys.exit()` — stdout, exit 0.
            MainArgs::command()
                .print_help()
                .expect("writing help to stdout");
            println!();
            ExitCode::SUCCESS
        }
        Action::Fail { message, code } => {
            eprintln!("{message}");
            ExitCode::from(code as u8)
        }
        Action::Run(args) => {
            if let Some(outfolder) = &args.outfolder {
                // The reference does this before the window check; we do it
                // after, because it is the only step here with a side effect
                // and the window check is pure. Observable difference: a bad
                // `--w` with a nonexistent `--outfolder` leaves the folder
                // uncreated. Noted rather than reproduced; revisit if a
                // pipeline turns out to depend on the empty folder.
                if let Err(e) = std::fs::create_dir_all(outfolder) {
                    eprintln!("could not create outfolder {}: {e}", outfolder.display());
                    return ExitCode::from(cli::EXIT_VALIDATION as u8);
                }
            }

            // `print("ARGS", args)` at the top of the reference's `main()`,
            // after the outfolder is created. A debug leftover, but it is on
            // stdout, so it is contract — and it doubles as the differential
            // oracle for argument parsing. See `MainArgs::args_line`.
            println!("{}", args.args_line());

            match run(args.as_ref()) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("isONform: {e}");
                    ExitCode::from(1)
                }
            }
        }
    }
}

/// `main.main(args)`: read, trim, batch, and run each batch through the stages.
fn run(args: &MainArgs) -> std::io::Result<()> {
    let Some(fastq) = &args.fastq else {
        return Ok(());
    };
    let outfolder = args
        .outfolder
        .clone()
        .unwrap_or_else(|| std::path::PathBuf::from("."));

    let text = std::fs::read_to_string(fastq)?;
    let records = isonform::fastq::read_fastq(&text);

    // `if len(all_reads) <= args.exact_instance_limit: args.exact = True` --- the
    // reference mutates its own args here. `exact` only clears
    // `previously_corrected_regions`, which is never populated (finding 27), so
    // it has no observable effect and is not threaded through.
    let params = isonform::driver::RunParams {
        k: args.k,
        w: args.w,
        xmin: args.xmin,
        xmax: args.xmax,
        delta_len: args.delta_len as i64,
        max_seqs: args.max_seqs,
        max_seqs_to_spoa: args.max_seqs_to_spoa,
        set_w_dynamically: args.set_w_dynamically,
        slow: args.slow,
        // `main:568` hard-codes `delta = 0.15`; there is no `--delta` on this
        // entry point, unlike `isONform_parallel` which defaults it to 0.1.
        delta: 0.15,
        delta_iso_len_3: args.delta_iso_len_3,
        delta_iso_len_5: args.delta_iso_len_5,
        wis: Default::default(),
        build: Default::default(),
    };

    let batches = isonform::driver::run_cluster(&records, &params);
    for b in &batches {
        write_batch_outputs(
            &outfolder,
            b,
            args.parallel.is_set(),
            parallel_batch_id(fastq).as_deref(),
        )?;
    }
    Ok(())
}

/// `p_batch_id`: the batch id `--parallel` takes from the *filename*.
///
/// ```python
/// filename = args.fastq.split("/")[-1]
/// tmp_filename = filename.split("_")
/// tmp_lastpart = tmp_filename[-1].split(".")
/// p_batch_id = tmp_lastpart[0]
/// ```
///
/// So `.../3.fastq` and `.../cluster_3.fastq` both give `"3"`. It is a *string*,
/// never parsed, and it reaches the isoform ids in the final output, so a file
/// whose stem is not a number produces ids containing that stem.
fn parallel_batch_id(fastq: &std::path::Path) -> Option<String> {
    let filename = fastq.to_str()?.rsplit('/').next()?;
    let last = filename.rsplit('_').next()?;
    Some(last.split('.').next()?.to_string())
}

/// The four per-batch files, in the reference's names.
///
/// # Under `--parallel` every batch writes to the *same* four filenames
///
/// The reference derives all of them from `p_batch_id`, which comes from the
/// input filename and does not change as the batch loop runs (`main:373-389`,
/// and `batch_id = p_batch_id` again at `main:572` before the isoform step). A
/// cluster large enough to be split into two batches therefore has its first
/// batch's four files overwritten by its second, and only the last batch
/// survives on disk --- silently, since every open is `'w'`. Without
/// `--parallel` the loop index is used and each batch keeps its own files.
///
/// Reproduced rather than repaired: writing the batches in order and letting
/// `File::create` truncate gives the same surviving files. `PORTING.md`
/// finding 33.
fn write_batch_outputs(
    outfolder: &std::path::Path,
    b: &isonform::driver::BatchOutput,
    parallel: bool,
    p_batch_id: Option<&str>,
) -> std::io::Result<()> {
    use std::io::Write as _;

    // The id in the filenames, and the suffix on the batch file: `"{id}_batch"`
    // under `--parallel`, `"{id}batch"` without it.
    let (id, batch_name) = match (parallel, p_batch_id) {
        (true, Some(p)) => (p.to_string(), format!("{p}_batch")),
        _ => (b.batch_id.to_string(), format!("{}batch", b.batch_id)),
    };

    // `skip{id}.fa` --- the one the reference writes as text, so this is
    // byte-compatible with it.
    let mut skip = std::fs::File::create(outfolder.join(format!("skip{id}.fa")))?;
    for (acc, seq) in &b.skipped {
        writeln!(skip, ">{acc}")?;
        skip.write_all(seq)?;
        writeln!(skip)?;
    }

    // The three the reference pickles. Same content, plain text --- see the
    // module docs.
    let mut batch = std::fs::File::create(outfolder.join(&batch_name))?;
    for (acc, seq) in &b.reads {
        write!(batch, "{acc}\t")?;
        batch.write_all(seq)?;
        writeln!(batch)?;
    }

    let mut spoa = std::fs::File::create(outfolder.join(format!("spoa{id}")))?;
    for (id, seq) in &b.isoforms {
        write!(spoa, "{id}\t")?;
        spoa.write_all(seq)?;
        writeln!(spoa)?;
    }

    let mut mapping = std::fs::File::create(outfolder.join(format!("mapping{id}")))?;
    for (id, accs) in &b.mapping {
        writeln!(mapping, "{id}\t{}", accs.join(","))?;
    }
    Ok(())
}
