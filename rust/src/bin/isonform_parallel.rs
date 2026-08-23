//! `isONform_parallel` — fans out over a folder of per-cluster fastqs.
//!
//! Only argument handling is ported so far; see `main.rs` for why reaching the
//! algorithm exits 3.

use std::io::Write;
use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use isonform::cli::{self, Action, ParallelArgs};

/// Prefix on every line of the not-yet-implemented notice, so a harness can
/// strip it and compare the stderr the reference genuinely produced.
const PORT_NOTICE: &str = "[port-unimplemented] ";

fn main() -> ExitCode {
    let argv: Vec<std::ffi::OsString> = std::env::args_os().collect();
    let argc = argv.len();

    let args = match ParallelArgs::try_parse_from(&argv) {
        Ok(a) => a,
        Err(e) => e.exit(),
    };

    let (emissions, action) = cli::resolve_parallel(args, argc);

    // The first of these is `print(len(sys.argv))` — a debug leftover in the
    // reference that lands on stdout, so it is part of the contract. See
    // PORTING.md finding 4; removing it is its own commit.
    for line in &emissions.stdout {
        println!("{line}");
    }
    for line in &emissions.stderr {
        eprintln!("{line}");
    }

    match action {
        Action::HelpAndExit => {
            ParallelArgs::command()
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
                if let Err(e) = std::fs::create_dir_all(outfolder) {
                    eprintln!("could not create outfolder {}: {e}", outfolder.display());
                    return ExitCode::from(cli::EXIT_VALIDATION as u8);
                }
            }

            // Every line carries the sentinel so bench/equivalence.sh can strip
            // the notice and compare the stderr the reference actually produced.
            // Nothing downstream parses it; it exists so a half-ported binary
            // cannot be mistaken for one that ran and found no isoforms.
            let mut err = std::io::stderr();
            for line in [
                "isONform_parallel (Rust): argument handling is ported; the algorithm is not.",
                "Nothing was written. Use the Python entry point from the",
                "isonform-ref environment for real runs \u{2014} see PORTING.md.",
            ] {
                let _ = writeln!(err, "{PORT_NOTICE}{line}");
            }
            ExitCode::from(3)
        }
    }
}
