//! `isONform_parallel` — fans out over a folder of per-cluster fastqs.
//!
//! One `main` process per fastq, `--t` of them at a time, then the outputs are
//! read back, filtered by `--iso_abundance`, and concatenated into the three
//! `transcriptome*` files. The stages themselves live in [`isonform::parallel`].
//!
//! Like the reference, this **spawns the `main` executable** rather than calling
//! the library in-process: that is what produces the per-cluster `stderr.txt` and
//! `stdout{cl}_{batch}.txt`, and it keeps `--t` meaning the same thing. `main` is
//! looked up next to this executable, as the reference looks for it next to its
//! own `__file__`.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use clap::{CommandFactory, Parser};
use isonform::cli::{self, Action, ParallelArgs};
use isonform::parallel;

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
            match run(&args) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("isONform_parallel: {e}");
                    ExitCode::from(1)
                }
            }
        }
    }
}

fn run(args: &ParallelArgs) -> std::io::Result<()> {
    let globstart = Instant::now();
    let Some(directory) = &args.fastq_folder else {
        // `args.fastq_folder` defaults to `False` in the reference, and
        // `os.listdir(False)` is a TypeError. Reported instead.
        return Err(std::io::Error::other("--fastq_folder is required"));
    };
    let outfolder = args.outfolder.clone().unwrap_or_else(|| PathBuf::from("."));

    // Rewrites `directory` in place when it holds only subdirectories. See
    // `restructure_isoncorrect_output` — this is destructive, and deliberately so.
    parallel::restructure_isoncorrect_output(directory)?;

    // `write_low_abundance = False`, a local that nothing ever assigns again:
    // no flag reaches it. PORTING.md finding 35.
    let write_low_abundance = false;

    // Known bugs are fixed by default; `ISONFORM_BUG_COMPAT` puts named ones
    // back. Read here as well as in `main` because cross-batch merging (finding
    // 31) lives on this side of the fork.
    let compat = isonform::driver::bug_compat_from_env().map_err(std::io::Error::other)?;
    let reproduced = compat.describe();
    if !reproduced.is_empty() {
        eprintln!(
            "isONform_parallel: reproducing reference bugs: {}",
            reproduced.join(",")
        );
    }

    let mut split_tmp: Option<PathBuf> = None;
    let split_directory = if args.split_wrt_batches {
        let tmp_work_dir = match &args.tmpdir {
            Some(d) => {
                std::fs::create_dir_all(d)?;
                d.clone()
            }
            None => {
                let d = std::env::temp_dir().join(format!("isonform_{}", std::process::id()));
                std::fs::create_dir_all(&d)?;
                d
            }
        };
        println!("Temporary workdirectory: {}", tmp_work_dir.display());
        let d = parallel::split_cluster_in_batches_clust(directory, &tmp_work_dir, args.max_seqs)?;
        split_tmp = Some(d.clone());
        d
    } else {
        println!("SplitDIR {}", directory.display());
        directory.clone()
    };

    let instances = match parallel::discover_instances(&split_directory, &outfolder)? {
        Ok(v) => v,
        Err(msg) => return Err(std::io::Error::other(msg)),
    };
    println!("Printing instances");
    for i in &instances {
        println!(
            "('{}', '{}', '{}', '{}')",
            i.fastq.display(),
            i.outfolder.display(),
            i.batch_id,
            i.cl_id
        );
    }
    println!("Using {0} cores.", args.nr_cores);

    let exe = main_executable()?;
    let start_multi = Instant::now();
    run_instances(&exe, &instances, args)?;
    println!(
        "Time elapsed multiprocessing: {:?}",
        start_multi.elapsed().as_secs_f64()
    );

    println!("Merging...");
    let file_handling = Instant::now();
    // `delta` and `max_seqs_to_spoa` reach cross-batch merging and nothing else,
    // which is why finding 37 lists them as inert --- they were inert only
    // because the merge never ran. Now they matter.
    let merge_opts = isonform::isoforms::MergeOpts {
        delta: args.delta,
        delta_len: args.delta_len,
        delta_iso_len_3: args.delta_iso_len_3,
        delta_iso_len_5: args.delta_iso_len_5,
        max_seqs_to_spoa: args.max_seqs_to_spoa,
        cigar_diversity_counts_runs: compat.cigar_diversity_counts_runs,
    };
    parallel::join_back_via_batch_merging(
        &outfolder,
        args.iso_abundance,
        args.write_fastq,
        write_low_abundance,
        merge_opts,
        compat.batch_merge_no_op,
    )?;
    parallel::generate_full_output(&outfolder, args.write_fastq, write_low_abundance)?;
    parallel::remove_folders(&outfolder)?;
    if let Some(d) = split_tmp {
        println!("Removed the split directory");
        std::fs::remove_dir_all(d)?;
    }
    println!(
        "Joined back batched files in: {}",
        file_handling.elapsed().as_secs_f64()
    );
    println!(
        "Finished full algo after : {}",
        globstart.elapsed().as_secs_f64()
    );
    Ok(())
}

/// `os.path.join(isONform_location, "main")`, where `isONform_location` is the
/// directory holding the entry point. Here that is the directory holding *this*
/// executable, so a built `target/release` tree works without installation.
fn main_executable() -> std::io::Result<PathBuf> {
    let exe = std::env::current_exe()?;
    let dir = exe.parent().unwrap_or(Path::new("."));
    let cand = dir.join("main");
    if cand.is_file() {
        return Ok(cand);
    }
    Err(std::io::Error::other(format!(
        "cannot find the `main` executable next to {} — build both binaries",
        exe.display()
    )))
}

/// `wccount`: the reference shells out to `wc -l`.
fn line_count(path: &Path) -> usize {
    std::fs::read_to_string(path)
        .map(|t| t.lines().count())
        .unwrap_or(0)
}

/// Has this instance already produced complete output? (`--keep_old`.)
///
/// True when the batch file holds one record per read in the input fastq. A
/// partial file --- an interrupted run --- has fewer and is recomputed.
fn already_complete(inst: &parallel::Instance) -> bool {
    let reads = line_count(&inst.fastq) / 4;
    if reads == 0 {
        return false;
    }
    // The batch file's name depends on how many batches `main` wrote; any of
    // them being complete for this instance is what matters.
    let Ok(entries) = std::fs::read_dir(&inst.outfolder) else {
        return false;
    };
    for e in entries.flatten() {
        let name = e.file_name().to_string_lossy().to_string();
        if name.ends_with("_batch") && line_count(&e.path()) == reads {
            return true;
        }
    }
    false
}

/// `Pool(processes=nr_cores).imap_unordered(isONform, instances)`.
///
/// Work is taken off a shared queue by `--t` threads, each spawning one `main`
/// process and waiting for it — the same concurrency the reference gets from a
/// process pool, and the same unordered completion reporting.
fn run_instances(
    exe: &Path,
    instances: &[parallel::Instance],
    args: &ParallelArgs,
) -> std::io::Result<()> {
    let next = AtomicUsize::new(0);
    let start = Instant::now();
    let failure: Mutex<Option<String>> = Mutex::new(None);
    let n_threads = args.nr_cores.max(1).min(instances.len().max(1));

    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                let Some(inst) = instances.get(i) else { break };
                if failure.lock().unwrap().is_some() {
                    break;
                }
                match run_one(exe, inst, args) {
                    Ok(true) => {
                        println!(
                            "{} (Time elapsed: {}s)",
                            inst.batch_id,
                            start.elapsed().as_secs()
                        );
                    }
                    Ok(false) => {}
                    Err(e) => {
                        *failure.lock().unwrap() = Some(e.to_string());
                        break;
                    }
                }
            });
        }
    });

    match failure.into_inner().unwrap() {
        Some(e) => Err(std::io::Error::other(e)),
        None => Ok(()),
    }
}

/// `isONform(data)`: one `main` invocation, with its stdout and stderr captured
/// into the cluster's folder.
///
/// Returns whether it ran — `--keep_old` can skip it.
fn run_one(exe: &Path, inst: &parallel::Instance, args: &ParallelArgs) -> std::io::Result<bool> {
    std::fs::create_dir_all(&inst.outfolder)?;

    if args.keep_old {
        // Finding 37, fixed. The reference looks for `isoforms.fastq`, a name
        // that appears exactly once in the whole codebase --- on the line that
        // looks for it. Nothing writes it, so the test is always false and
        // `--keep_old` has never skipped anything.
        //
        // What it *meant* to ask is "did this instance already finish?", and the
        // answer this port uses is: the batch file exists and holds one record
        // per read in the input. That is the closest analogue of the reference's
        // "same number of reads as the input, i.e. complete" --- note it cannot
        // be a line-for-line comparison as the reference attempts, because a
        // fastq has four lines per read and the batch file has one.
        if already_complete(inst) {
            println!(
                "Skipping batch_id:{}.{} --- already complete",
                inst.cl_id, inst.batch_id
            );
            return Ok(false);
        }
    }

    // `print(..., end=' ')` in the reference, so "Done with batch_id" lands on
    // the same line.
    print!(
        "Running isONform batch_id:{}.{}... ",
        inst.cl_id, inst.batch_id
    );
    std::io::stdout().flush()?;

    let err_file = std::fs::File::create(inst.outfolder.join("stderr.txt"))?;
    let out_file = std::fs::File::create(
        inst.outfolder
            .join(format!("stdout{}_{}.txt", inst.cl_id, inst.batch_id)),
    )?;

    // Exactly the argument list the reference builds. `--max_seqs` is commented
    // out there, and `--set_w_dynamically` and `--slow` are collected into the
    // params dict but never passed — so they do not reach `main`. This argument
    // list is why seven of this entry point's flags do nothing: finding 37.
    let status = std::process::Command::new(exe)
        .arg("--fastq")
        .arg(&inst.fastq)
        .arg("--outfolder")
        .arg(&inst.outfolder)
        .arg("--exact_instance_limit")
        .arg(args.exact_instance_limit.to_string())
        .arg("--k")
        .arg(args.k.to_string())
        .arg("--w")
        .arg(args.w.to_string())
        .arg("--xmin")
        .arg(args.xmin.to_string())
        .arg("--xmax")
        .arg(args.xmax.to_string())
        .arg("--delta_len")
        .arg(args.delta_len.to_string())
        .arg("--exact")
        .arg("--parallel")
        .arg("True")
        .arg("--delta_iso_len_3")
        .arg(args.delta_iso_len_3.to_string())
        .arg("--delta_iso_len_5")
        .arg(args.delta_iso_len_5.to_string())
        .stdout(out_file)
        .stderr(err_file)
        .status()?;

    if !status.success() {
        // `subprocess.check_call` raises, and the pool propagates it.
        return Err(std::io::Error::other(format!(
            "main exited {} on {}",
            status.code().unwrap_or(-1),
            inst.fastq.display()
        )));
    }
    println!("Done with batch_id:{}.{}", inst.cl_id, inst.batch_id);
    Ok(true)
}
