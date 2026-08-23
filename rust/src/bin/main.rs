//! `main` — the single-instance isoform reconstructor.
//!
//! Named `main` because that is what the Python ships (`setup.py`
//! `scripts=['isONform_parallel','main']`) and what `isONform_parallel` execs.
//!
//! Only argument handling is ported so far. Reaching the algorithm exits 3 with
//! a message saying so, rather than exiting 0 and emitting nothing — a port that
//! silently produces no isoforms is worse than one that says it cannot yet.

use std::io::Write;
use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isonform::cli::{self, Action, MainArgs};

/// Prefix on every line of the not-yet-implemented notice, so a harness can
/// strip it and compare the stderr the reference genuinely produced.
const PORT_NOTICE: &str = "[port-unimplemented] ";

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

            // Every line carries the sentinel so bench/equivalence.sh can strip
            // the notice and compare the stderr the reference actually produced.
            // Nothing downstream parses it; it exists so a half-ported binary
            // cannot be mistaken for one that ran and found no isoforms.
            let mut err = std::io::stderr();
            for line in [
                "isONform (Rust): argument handling is ported; the algorithm is not.",
                "Nothing was written. Use the Python entry point from the",
                "isonform-ref environment for real runs \u{2014} see PORTING.md.",
            ] {
                let _ = writeln!(err, "{PORT_NOTICE}{line}");
            }
            ExitCode::from(3)
        }
    }
}
