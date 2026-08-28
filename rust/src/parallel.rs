//! `isONform_parallel`: fan out over a folder of clusters, then join back.
//!
//! Covers `isONform_parallel`'s own body plus
//! `modules/Parallelization_side_functions.py` and the file-reading half of
//! `batch_merging_parallel.join_back_via_batch_merging`. The merging half is in
//! [`crate::batch_merge`], where it is a no-op (`PORTING.md` finding 31).
//!
//! The shape is: discover one instance per fastq, run `main` on each, then read
//! every cluster's intermediates back, decide what to write, and concatenate the
//! per-cluster files into the three `transcriptome*` files.
//!
//! # `--iso_abundance` silently discards, because `write_low_abundance` is a local
//!
//! `main(args)` opens with `write_low_abundance = False` and never assigns it
//! again; there is no flag that reaches it. So every isoform below
//! `--iso_abundance` (default 5) is dropped with no file recording it, and the
//! whole low-abundance half of `write_final_output` is unreachable from this entry
//! point --- including finding 32's always-1 support count. `PORTING.md`
//! finding 35.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use crate::batch_merge::{format_consensus, select_output, Destination, Isoform};

/// One unit of work: one fastq, run by one `main` invocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Instance {
    pub fastq: PathBuf,
    /// `outfolder/{cl_id}` --- shared by every instance of the same cluster.
    pub outfolder: PathBuf,
    pub cl_id: String,
    pub batch_id: String,
}

/// What `restructure_isoncorrect_output` decided the input folder was.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputStructure {
    /// Top-level fastq files, one per cluster. Nothing is moved.
    IsonClust,
    /// Subdirectories holding `corrected_reads.fastq`. **The input folder is
    /// rewritten in place** --- see [`restructure_isoncorrect_output`].
    IsonCorrect,
}

/// `restructure_isoncorrect_output`: flatten isONcorrect's output in place.
///
/// If the folder contains no *files* at the top level, the reference treats it as
/// isONcorrect output: every subdirectory's `corrected_reads.fastq` is **moved**
/// up to `{subdirectory}.fastq`, and then every subdirectory is **deleted**.
///
/// # This destroys the input folder, and that is the reference's behaviour
///
/// Not a copy --- `shutil.move` followed by `remove_folders`. A folder of
/// subdirectories that happens to contain no top-level files is flattened and its
/// subdirectories removed, whatever else was in them. It is reproduced because
/// pipelines depend on it (this is how isONcorrect's output is fed in), but it is
/// worth knowing before pointing either implementation at a directory you care
/// about. `PORTING.md` finding 36.
pub fn restructure_isoncorrect_output(directory: &Path) -> std::io::Result<InputStructure> {
    let mut any_file = false;
    for e in std::fs::read_dir(directory)? {
        if e?.file_type()?.is_file() {
            any_file = true;
            break;
        }
    }
    if any_file {
        println!("isONclust structure");
        return Ok(InputStructure::IsonClust);
    }
    println!("isONcorrect structure");
    for e in std::fs::read_dir(directory)? {
        let e = e?;
        if !e.file_type()?.is_dir() {
            continue;
        }
        let source = e.path().join("corrected_reads.fastq");
        let name = e.file_name();
        let target = directory.join(format!("{}.fastq", name.to_string_lossy()));
        if source.exists() {
            std::fs::rename(&source, &target)?;
            println!(
                "Moved and renamed {} to {}",
                source.display(),
                target.display()
            );
        } else {
            println!("File {} does not exist.", source.display());
        }
    }
    remove_folders(directory)?;
    Ok(InputStructure::IsonCorrect)
}

/// `splitfile`: cut one fastq into `{cl_id}_{i}.{ext}` chunks of `chunksize`
/// **lines**.
///
/// The reference counts lines, not records, and the caller always passes
/// `4 * max_seqs` --- so a fastq whose records are not exactly four lines each
/// (a wrapped sequence) is cut mid-record. Reproduced; it does not arise on
/// isONclust output, which is always four lines per record.
///
/// The empty trailing chunk that appears when the file divides evenly is removed,
/// as the reference does.
pub fn splitfile(
    infile: &Path,
    tmp_outdir: &Path,
    chunksize: usize,
    cl_id: &str,
    ext: &str,
) -> std::io::Result<()> {
    let text = std::fs::read_to_string(infile)?;
    // `infile.readline()` yields lines *with* their terminator, and the last line
    // of a file with no trailing newline still counts. Splitting inclusively
    // reproduces that.
    let mut lines: Vec<String> = Vec::new();
    let mut cur = String::new();
    for ch in text.chars() {
        cur.push(ch);
        if ch == '\n' {
            lines.push(std::mem::take(&mut cur));
        }
    }
    if !cur.is_empty() {
        lines.push(cur);
    }

    let mut i = 0usize;
    let mut pos = 0usize;
    loop {
        let out = tmp_outdir.join(format!("{cl_id}_{i}.{ext}"));
        let end = (pos + chunksize).min(lines.len());
        let chunk: String = lines[pos..end].concat();
        std::fs::write(&out, chunk.as_bytes())?;
        let written = end > pos;
        // "Corner case: Original file is even multiple of max_seqs, hence the
        // last file becomes empty. Remove this."
        if std::fs::metadata(&out)?.len() == 0 {
            std::fs::remove_file(&out)?;
        }
        pos = end;
        if !written {
            break;
        }
        i += 1;
    }
    Ok(())
}

/// `split_cluster_in_batches_clust`: one file per cluster becomes one file per
/// batch, under `{tmp_work_dir}/split_in_batches`.
///
/// A cluster with more than `4 * max_seqs` lines is cut up; a smaller one is
/// symlinked to `{cl_id}_0.fastq`. Returns the directory the instances are then
/// discovered from.
///
/// Reached only with `--split_wrt_batches`. Without it the input folder is used
/// directly and the cluster id is the whole filename stem.
pub fn split_cluster_in_batches_clust(
    indir: &Path,
    tmp_work_dir: &Path,
    max_seqs: usize,
) -> std::io::Result<PathBuf> {
    let out = tmp_work_dir.join("split_in_batches");
    std::fs::create_dir_all(&out)?;

    // `Path(indir).rglob('*.fastq')` --- recursive.
    let mut files = Vec::new();
    collect_fastqs(indir, &mut files)?;
    for path in files {
        let fname = path.file_name().unwrap().to_string_lossy().to_string();
        // `fastq_file.rsplit('.', 1)` --- split on the *last* dot.
        let (cl_id, ext) = match fname.rsplit_once('.') {
            Some((a, b)) => (a.to_string(), b.to_string()),
            None => (fname.clone(), String::new()),
        };
        let n_lines = BufReader::new(std::fs::File::open(&path)?).lines().count();
        if n_lines > 4 * max_seqs {
            splitfile(&path, &out, 4 * max_seqs, &cl_id, &ext)?;
        } else {
            let link = out.join(format!("{cl_id}_0.{ext}"));
            println!("SYMLINK {}", link.display());
            symlink_force(&path, &link)?;
        }
    }
    Ok(out)
}

fn collect_fastqs(dir: &Path, out: &mut Vec<PathBuf>) -> std::io::Result<()> {
    for e in std::fs::read_dir(dir)? {
        let e = e?;
        let p = e.path();
        if e.file_type()?.is_dir() {
            collect_fastqs(&p, out)?;
        } else if p.extension().is_some_and(|x| x == "fastq") {
            out.push(p);
        }
    }
    Ok(())
}

/// `symlink_force`: create a symlink, replacing a broken one.
///
/// The reference replaces the link only when it is *broken*; an existing link
/// that resolves is left alone, so a stale-but-valid link to a different file
/// survives. Reproduced.
pub fn symlink_force(target: &Path, link: &Path) -> std::io::Result<()> {
    match std::os::unix::fs::symlink(target, link) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
            if std::fs::metadata(link).is_err() {
                println!("path {} is a broken symlink", link.display());
                std::fs::remove_file(link)?;
                symlink_force(target, link)
            } else {
                Ok(())
            }
        }
        Err(e) => Err(e),
    }
}

/// Discover one [`Instance`] per `.fastq` in `split_directory`.
///
/// The cluster id is everything before the first `_` of the stem, and the batch
/// id is what follows it (`0` when there is none). Both come from the *filename*,
/// which is why `--parallel` names its outputs after the file rather than the
/// loop index (finding 34).
///
/// Sorted by `(int(cl_id), int(batch_id))`, so a cluster whose filename stem is
/// not a number is a hard error in the reference (`ValueError` out of the sort
/// key). Kept: such a file would also break the output ids downstream.
pub fn discover_instances(
    split_directory: &Path,
    outfolder: &Path,
) -> std::io::Result<Result<Vec<Instance>, String>> {
    let mut instances = Vec::new();
    for e in std::fs::read_dir(split_directory)? {
        let e = e?;
        let name = e.file_name().to_string_lossy().to_string();
        if !name.ends_with(".fastq") {
            continue;
        }
        // `read_fastq_file.split(".")[0]` --- the *first* dot, not the last.
        let stem = name.split('.').next().unwrap_or(&name);
        let mut parts = stem.splitn(2, '_');
        let cl_id = parts.next().unwrap_or("").to_string();
        let batch_id = parts.next().unwrap_or("0").to_string();
        instances.push(Instance {
            fastq: split_directory.join(&name),
            outfolder: outfolder.join(&cl_id),
            cl_id,
            batch_id,
        });
    }
    for i in &instances {
        if i.cl_id.parse::<i64>().is_err() || i.batch_id.parse::<i64>().is_err() {
            return Ok(Err(format!(
                "cluster id {:?} or batch id {:?} is not a number, from {}",
                i.cl_id,
                i.batch_id,
                i.fastq.display()
            )));
        }
    }
    instances.sort_by_key(|i| {
        (
            i.cl_id.parse::<i64>().unwrap(),
            i.batch_id.parse::<i64>().unwrap(),
        )
    });
    Ok(Ok(instances))
}

/// One cluster's intermediates, read back from its output folder.
#[derive(Debug, Default, Clone)]
pub struct ClusterBatches {
    /// `(batch_id, [(isoform_id, isoform)])`, in the order the batch files were
    /// listed --- which is the order records reach the output files.
    ///
    /// The batch id is a **string** because it is not always a single number:
    /// one `main` invocation that writes several batches names them
    /// `{p_batch_id}_{batch_id}` so they cannot overwrite each other
    /// (finding 34). Ordering uses [`crate::batch_merge::batch_sort_key`].
    pub batches: Vec<(String, Vec<(u32, Isoform)>)>,
}

/// Read one cluster's `{id}_batch`, `spoa{id}` and `mapping{id}` files.
///
/// # These are this port's own intermediates, not the reference's
///
/// `main` in the reference pickles all three; this port writes the same content
/// as TSV, because a Rust binary emitting Python pickles would be absurd and
/// nothing outside a single run reads them. The consequence is that the two
/// implementations cannot be mixed *within* one pipeline run --- the Python
/// `isONform_parallel` cannot join back this port's batches, and vice versa. The
/// user-visible outputs are byte-compatible; only the intermediates differ.
pub fn read_cluster(cl_dir: &Path) -> std::io::Result<ClusterBatches> {
    let mut out = ClusterBatches::default();
    // `for batchfile in os.listdir(cl_dir)` --- directory order, and the batch
    // ids are inserted into `all_infos_dict` in exactly that order.
    for e in std::fs::read_dir(cl_dir)? {
        let name = e?.file_name().to_string_lossy().to_string();
        let Some(prefix) = name.strip_suffix("_batch") else {
            continue;
        };
        // The reference does `int(batchfile.split('_')[0])`, which cannot tell
        // `3_0_batch` from `3_1_batch`. The whole prefix is the id here; it must
        // be non-empty and all-numeric per component.
        if prefix.is_empty()
            || !prefix
                .split('_')
                .all(|c| !c.is_empty() && c.bytes().all(|b| b.is_ascii_digit()))
        {
            continue;
        }
        let batch_id = prefix.to_string();
        let spoa = read_tsv(&cl_dir.join(format!("spoa{batch_id}")))?;
        let mapping = read_tsv(&cl_dir.join(format!("mapping{batch_id}")))?;
        let map: BTreeMap<&str, &str> = mapping
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();

        let mut isoforms = Vec::new();
        for (cons_id, seq) in &spoa {
            let Ok(id) = cons_id.parse::<u32>() else {
                continue;
            };
            let reads: Vec<String> = match map.get(cons_id.as_str()) {
                Some(v) if !v.is_empty() => v.split(',').map(|s| s.to_string()).collect(),
                _ => Vec::new(),
            };
            isoforms.push((
                id,
                Isoform {
                    sequence: seq.as_bytes().to_vec(),
                    reads,
                    // `Read(seq, reads, False)` --- nothing arrives merged.
                    merged: false,
                },
            ));
        }
        out.batches.push((batch_id, isoforms));
    }
    Ok(out)
}

fn read_tsv(path: &Path) -> std::io::Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let Ok(text) = std::fs::read_to_string(path) else {
        return Ok(out);
    };
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let (k, v) = line.split_once('\t').unwrap_or((line, ""));
        out.push((k.to_string(), v.to_string()));
    }
    Ok(out)
}

/// `repr()` of a Python list of `str`, which is what lands in the mapping file.
///
/// `mapping_file.write(">{0}\n{1}\n".format(new_id, infos.reads))` formats a
/// *list*, so the file contains `['a', 'b']` --- brackets, single quotes, and a
/// comma-space between entries. Not a format anything reads back, but it is the
/// bytes the reference produces, so it is the bytes this produces.
pub fn python_list_repr(items: &[String]) -> String {
    let mut s = String::from("[");
    for (i, x) in items.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        s.push('\'');
        s.push_str(x);
        s.push('\'');
    }
    s.push(']');
    s
}

/// `write_final_output`, with its file handles.
///
/// Called **once per batch** by the reference, each call reopening every file
/// with `'w'` and rewriting the whole thing (finding 32). Since the content does
/// not depend on which batch triggered the call, this writes once; the surviving
/// bytes are identical.
pub fn write_final_output(
    batches: &[(String, Vec<(u32, Isoform)>)],
    outfolder: &Path,
    cl_id: &str,
    iso_abundance: usize,
    write_fastq: bool,
    write_low_abundance: bool,
) -> std::io::Result<()> {
    let records = select_output(batches, cl_id, iso_abundance, write_low_abundance);

    let ext = if write_fastq { "fq" } else { "fa" };
    let mut consensus =
        std::fs::File::create(outfolder.join(format!("cluster{cl_id}_merged.{ext}")))?;
    let mut mapping = std::fs::File::create(outfolder.join(format!("cluster{cl_id}_mapping.txt")))?;
    let mut support = std::fs::File::create(outfolder.join(format!("support_{cl_id}.txt")))?;

    // The low-abundance files are only *created* when they are written to ---
    // the reference opens them inside `if write_low_abundance:`. Since
    // `isONform_parallel` hardcodes that to false, they never exist there.
    let mut low: Option<(std::fs::File, std::fs::File, std::fs::File)> = if write_low_abundance {
        Some((
            std::fs::File::create(
                outfolder.join(format!("cluster{cl_id}_merged_low_abundance.{ext}")),
            )?,
            std::fs::File::create(
                outfolder.join(format!("cluster{cl_id}_mapping_low_abundance.txt")),
            )?,
            std::fs::File::create(outfolder.join(format!("support_{cl_id}low_abundance.txt")))?,
        ))
    } else {
        None
    };

    for rec in &records {
        match rec.destination {
            Destination::Main => {
                writeln!(mapping, ">{}", rec.id)?;
                writeln!(mapping, "{}", python_list_repr(&rec.reads))?;
                consensus.write_all(&format_consensus(rec, write_fastq))?;
                writeln!(support, "{}: {}", rec.id, rec.support)?;
            }
            Destination::LowAbundance => {
                if let Some((c, m, s)) = low.as_mut() {
                    c.write_all(&format_consensus(rec, write_fastq))?;
                    writeln!(m, ">{}", rec.id)?;
                    writeln!(m, "{}", python_list_repr(&rec.reads))?;
                    writeln!(s, "{}: {}", rec.id, rec.support)?;
                }
            }
            Destination::Dropped => {}
        }
    }
    Ok(())
}

/// `join_back_via_batch_merging`: read every cluster back, merge across its
/// batches, and write its output.
///
/// The merging step is repaired here rather than reproduced --- see
/// [`crate::batch_merge::actual_merging_process`] and `PORTING.md` finding 31.
/// `batch_merge_no_op` puts the reference's do-nothing behaviour back.
pub fn join_back_via_batch_merging(
    outdir: &Path,
    iso_abundance: usize,
    write_fastq: bool,
    write_low_abundance: bool,
    opts: crate::isoforms::MergeOpts,
    batch_merge_no_op: bool,
    threads: usize,
) -> std::io::Result<()> {
    println!("Batch Merging");
    let dirs = subdirectories(outdir)?;

    // One cluster per thread. Clusters are independent --- `actual_merging_process`
    // only ever sees one cluster's batches, and `write_final_output` writes that
    // cluster's own `cluster{id}_merged.*` --- so this is embarrassingly parallel.
    //
    // The reference closes its process pool *before* calling this, so it merges
    // one cluster at a time on one core while the rest idle. That is faithful and
    // slow: the within-batch half of the deep-12 run uses 8 threads and takes
    // 179s, and this half had one core for hours. Nothing about the merge's
    // *result* depends on cluster order, so the port parallelises it.
    let n = threads.max(1).min(dirs.len().max(1));
    let next = std::sync::atomic::AtomicUsize::new(0);
    let failed: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);
    std::thread::scope(|scope| {
        for _ in 0..n {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(cl_dir) = dirs.get(i) else { break };
                    if let Err(e) = merge_one_cluster(
                        cl_dir,
                        outdir,
                        iso_abundance,
                        write_fastq,
                        write_low_abundance,
                        opts,
                        batch_merge_no_op,
                    ) {
                        let mut slot = failed.lock().unwrap();
                        if slot.is_none() {
                            *slot = Some(e);
                        }
                        break;
                    }
                }
            });
        }
    });
    if let Some(e) = failed.lock().unwrap().take() {
        return Err(e);
    }
    if std::env::var_os("ISONFORM_MERGE_STATS").is_some() {
        use crate::batch_merge::{MERGES, PAIRS_ALIGNED, PAIRS_SEEN, PAIRS_SKIPPED, SKIPPED_WOULD_MERGE};
        use std::sync::atomic::Ordering::Relaxed;
        let (seen, skipped, aligned) = (
            PAIRS_SEEN.load(Relaxed),
            PAIRS_SKIPPED.load(Relaxed),
            PAIRS_ALIGNED.load(Relaxed),
        );
        eprintln!(
            "merge-stats pairs_seen={seen} skipped={skipped} ({:.1}%) aligned={aligned} \
             merges={} skipped_would_merge={}",
            if seen > 0 { 100.0 * skipped as f64 / seen as f64 } else { 0.0 },
            MERGES.load(Relaxed),
            SKIPPED_WOULD_MERGE.load(Relaxed)
        );
    }
    Ok(())
}

/// One cluster's cross-batch merge and output. Split out of
/// [`join_back_via_batch_merging`] so the cluster loop can be threaded.
fn merge_one_cluster(
    cl_dir: &Path,
    outdir: &Path,
    iso_abundance: usize,
    write_fastq: bool,
    write_low_abundance: bool,
    opts: crate::isoforms::MergeOpts,
    batch_merge_no_op: bool,
) -> std::io::Result<()> {
    {
        let cl_id = cl_dir
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let mut cluster = read_cluster(cl_dir)?;
        let mut engine = crate::isoforms::SpoaParasailMerge;
        crate::batch_merge::actual_merging_process(
            &mut engine,
            &mut cluster.batches,
            opts,
            batch_merge_no_op,
        );
        // The reference calls `write_final_output` inside its per-batch loop, so
        // a cluster with no batch files at all is never written. Same here.
        if cluster.batches.is_empty() {
            return Ok(());
        }
        write_final_output(
            &cluster.batches,
            outdir,
            &cl_id,
            iso_abundance,
            write_fastq,
            write_low_abundance,
        )?;
    }
    Ok(())
}

/// `[f.path for f in os.scandir(outfolder) if f.is_dir()]` --- directory order.
///
/// Not sorted. The order reaches the output: it is the order clusters are
/// concatenated into `transcriptome.fasta`. `readdir(3)` order is what both
/// implementations get, so they agree on the same filesystem.
fn subdirectories(dir: &Path) -> std::io::Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    for e in std::fs::read_dir(dir)? {
        let e = e?;
        if e.file_type()?.is_dir() {
            out.push(e.path());
        }
    }
    Ok(out)
}

/// `remove_folders`: delete every subdirectory of `outfolder`.
pub fn remove_folders(outfolder: &Path) -> std::io::Result<()> {
    for d in subdirectories(outfolder)? {
        std::fs::remove_dir_all(d)?;
    }
    Ok(())
}

/// `generate_full_output`: concatenate the per-cluster files.
///
/// Which clusters are concatenated is decided by the *subdirectories* that still
/// exist, not by the files themselves — so this must run before
/// [`remove_folders`], and a cluster whose subdirectory name is not all digits is
/// skipped even if its output files are there.
pub fn generate_full_output(
    outfolder: &Path,
    write_fastq: bool,
    write_low_abundance: bool,
) -> std::io::Result<()> {
    let ext = if write_fastq { "fq" } else { "fa" };
    if write_fastq {
        println!("Generating transcriptome.fastq");
    } else {
        println!("Generating transcriptome.fasta");
    }
    let t_ext = if write_fastq { "fastq" } else { "fasta" };
    concat_per_cluster(
        outfolder,
        &format!("transcriptome.{t_ext}"),
        |c| format!("cluster{c}_merged.{ext}"),
        false,
    )?;
    concat_per_cluster(
        outfolder,
        "transcriptome_mapping.txt",
        |c| format!("cluster{c}_mapping.txt"),
        false,
    )?;
    concat_per_cluster(
        outfolder,
        "transcriptome_support.txt",
        |c| format!("support_{c}.txt"),
        false,
    )?;
    if write_low_abundance {
        concat_per_cluster(
            outfolder,
            &format!("transcriptome_low.{t_ext}"),
            |c| format!("cluster{c}_merged_low_abundance.{ext}"),
            true,
        )?;
        concat_per_cluster(
            outfolder,
            "transcriptome_mapping_low.txt",
            |c| format!("cluster{c}_mapping_low_abundance.txt"),
            false,
        )?;
    }
    Ok(())
}

/// The shared body of the five `generate_*` functions.
///
/// `stamp_headers` is `generate_low_abundance_output`'s quirk and only its own:
/// it does `line = line + str(actual_folder)` on every line starting with `@` or
/// `>`. The line still carries its newline, so the cluster id lands at the
/// *start of the following line* rather than the end of the header --- it
/// corrupts the sequence. Unreachable from `isONform_parallel`, which hardcodes
/// `write_low_abundance = False`; reproduced so it stays visible.
fn concat_per_cluster(
    outfolder: &Path,
    target: &str,
    per_cluster: impl Fn(&str) -> String,
    stamp_headers: bool,
) -> std::io::Result<()> {
    let mut out = std::fs::File::create(outfolder.join(target))?;
    for d in subdirectories(outfolder)? {
        let name = d
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        // `str.isdigit()` on an empty string is false.
        if name.is_empty() || !name.bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        let src = outfolder.join(per_cluster(&name));
        let Ok(text) = std::fs::read_to_string(&src) else {
            continue;
        };
        if !stamp_headers {
            out.write_all(text.as_bytes())?;
            continue;
        }
        for line in text.split_inclusive('\n') {
            if line.starts_with('@') || line.starts_with('>') {
                write!(out, "{line}{name}")?;
            } else {
                write!(out, "{line}")?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("isonform_par_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn iso(seq: &str, reads: &[&str]) -> Isoform {
        Isoform {
            sequence: seq.as_bytes().to_vec(),
            reads: reads.iter().map(|s| s.to_string()).collect(),
            merged: false,
        }
    }

    #[test]
    fn the_mapping_file_holds_a_python_list_repr() {
        // `"{1}".format(reads)` on a list, so the file really does contain
        // brackets and single quotes.
        assert_eq!(
            python_list_repr(&["87".into(), "37".into()]),
            "['87', '37']"
        );
        assert_eq!(python_list_repr(&[]), "[]");
        assert_eq!(python_list_repr(&["only".into()]), "['only']");
    }

    #[test]
    fn instances_sort_numerically_on_cluster_then_batch() {
        let d = tmp("discover");
        for f in [
            "10_2.fastq",
            "2_0.fastq",
            "10_0.fastq",
            "2_10.fastq",
            "notes.txt",
        ] {
            std::fs::write(d.join(f), "").unwrap();
        }
        let got = discover_instances(&d, Path::new("/out")).unwrap().unwrap();
        let keys: Vec<(String, String)> = got
            .iter()
            .map(|i| (i.cl_id.clone(), i.batch_id.clone()))
            .collect();
        // Numeric, not lexicographic --- "10" sorts after "2".
        assert_eq!(
            keys,
            vec![
                ("2".into(), "0".into()),
                ("2".into(), "10".into()),
                ("10".into(), "0".into()),
                ("10".into(), "2".into()),
            ]
        );
        assert_eq!(got[0].outfolder, Path::new("/out/2"));
    }

    #[test]
    fn a_bare_cluster_file_gets_batch_id_zero() {
        let d = tmp("bare");
        std::fs::write(d.join("7.fastq"), "").unwrap();
        let got = discover_instances(&d, Path::new("/out")).unwrap().unwrap();
        assert_eq!(got[0].cl_id, "7");
        assert_eq!(got[0].batch_id, "0", "`else 0` in the reference");
    }

    #[test]
    fn a_non_numeric_cluster_id_is_an_error_not_a_silent_skip() {
        // The reference's sort key is `int(x[5])`, so this is a ValueError
        // traceback there. Reported rather than reproduced as a panic.
        let d = tmp("nonnum");
        std::fs::write(d.join("sample_a.fastq"), "").unwrap();
        assert!(discover_instances(&d, Path::new("/out")).unwrap().is_err());
    }

    #[test]
    fn splitfile_cuts_on_lines_and_drops_the_empty_tail() {
        let d = tmp("split");
        let src = d.join("3.fastq");
        // Eight lines: exactly two chunks of four, so the reference's third
        // chunk is empty and gets removed.
        std::fs::write(&src, "a\nb\nc\nd\ne\nf\ng\nh\n").unwrap();
        let out = d.join("o");
        std::fs::create_dir_all(&out).unwrap();
        splitfile(&src, &out, 4, "3", "fastq").unwrap();
        let mut names: Vec<String> = std::fs::read_dir(&out)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
            .collect();
        names.sort();
        assert_eq!(
            names,
            vec!["3_0.fastq", "3_1.fastq"],
            "no empty third chunk"
        );
        assert_eq!(
            std::fs::read_to_string(out.join("3_0.fastq")).unwrap(),
            "a\nb\nc\nd\n"
        );
        assert_eq!(
            std::fs::read_to_string(out.join("3_1.fastq")).unwrap(),
            "e\nf\ng\nh\n"
        );
    }

    #[test]
    fn splitfile_keeps_a_short_final_chunk() {
        let d = tmp("split2");
        let src = d.join("3.fastq");
        std::fs::write(&src, "a\nb\nc\nd\ne\n").unwrap();
        let out = d.join("o");
        std::fs::create_dir_all(&out).unwrap();
        splitfile(&src, &out, 4, "3", "fastq").unwrap();
        assert_eq!(
            std::fs::read_to_string(out.join("3_1.fastq")).unwrap(),
            "e\n"
        );
    }

    #[test]
    fn the_final_output_names_and_formats_match_the_reference() {
        let d = tmp("final");
        let batches = vec![("0".to_string(), vec![(72u32, iso("ACGT", &["87", "37"]))])];
        write_final_output(&batches, &d, "0", 1, false, false).unwrap();
        assert_eq!(
            std::fs::read_to_string(d.join("cluster0_merged.fa")).unwrap(),
            ">0_0_72\nACGT\n"
        );
        assert_eq!(
            std::fs::read_to_string(d.join("cluster0_mapping.txt")).unwrap(),
            ">0_0_72\n['87', '37']\n"
        );
        assert_eq!(
            std::fs::read_to_string(d.join("support_0.txt")).unwrap(),
            "0_0_72: 2\n"
        );
        assert!(
            !d.join("cluster0_merged_low_abundance.fa").exists(),
            "the low-abundance files are not created when they are not written"
        );
    }

    #[test]
    fn below_the_cutoff_and_without_low_abundance_output_nothing_records_it() {
        // finding 35: `isONform_parallel` hardcodes write_low_abundance = False,
        // so an isoform under --iso_abundance vanishes with no trace.
        let d = tmp("cutoff");
        let batches = vec![("0".to_string(), vec![(1u32, iso("ACGT", &["a"]))])];
        write_final_output(&batches, &d, "0", 5, false, false).unwrap();
        assert_eq!(
            std::fs::read_to_string(d.join("cluster0_merged.fa")).unwrap(),
            "",
            "written nowhere at all"
        );
        assert_eq!(
            std::fs::read_to_string(d.join("support_0.txt")).unwrap(),
            ""
        );
    }

    #[test]
    fn only_all_digit_subfolders_are_concatenated() {
        let d = tmp("concat");
        for c in ["0", "12", "Analysis"] {
            std::fs::create_dir_all(d.join(c)).unwrap();
            std::fs::write(
                d.join(format!("cluster{c}_merged.fa")),
                format!(">{c}\nAC\n"),
            )
            .unwrap();
        }
        generate_full_output(&d, false, false).unwrap();
        let got = std::fs::read_to_string(d.join("transcriptome.fasta")).unwrap();
        assert!(got.contains(">0\n") && got.contains(">12\n"));
        assert!(
            !got.contains(">Analysis"),
            "`actual_folder.isdigit()` gates it"
        );
    }

    #[test]
    fn a_cluster_with_no_output_file_is_skipped_rather_than_failing() {
        let d = tmp("missing");
        std::fs::create_dir_all(d.join("4")).unwrap();
        generate_full_output(&d, false, false).unwrap();
        assert_eq!(
            std::fs::read_to_string(d.join("transcriptome.fasta")).unwrap(),
            ""
        );
    }

    #[test]
    fn reading_a_cluster_back_pairs_consensuses_with_their_reads() {
        let d = tmp("readback");
        std::fs::write(d.join("3_batch"), "acc1\tACGT\n").unwrap();
        std::fs::write(d.join("spoa3"), "1\tACGT\n7\tTTTT\n").unwrap();
        std::fs::write(d.join("mapping3"), "1\tacc1,acc2\n7\tacc3\n").unwrap();
        let got = read_cluster(&d).unwrap();
        assert_eq!(got.batches.len(), 1);
        let (bid, isoforms) = &got.batches[0];
        assert_eq!(bid, "3");
        assert_eq!(isoforms.len(), 2);
        assert_eq!(isoforms[0].0, 1);
        assert_eq!(isoforms[0].1.reads, vec!["acc1", "acc2"]);
        assert!(!isoforms[0].1.merged, "nothing arrives merged");
        assert_eq!(isoforms[1].1.sequence, b"TTTT");
    }

    #[test]
    fn a_folder_with_top_level_files_is_left_alone() {
        let d = tmp("structure");
        std::fs::write(d.join("0.fastq"), "@r\nA\n+\nI\n").unwrap();
        std::fs::create_dir_all(d.join("sub")).unwrap();
        assert_eq!(
            restructure_isoncorrect_output(&d).unwrap(),
            InputStructure::IsonClust
        );
        assert!(d.join("sub").exists(), "nothing is moved or deleted");
    }

    #[test]
    fn a_folder_of_only_subdirectories_is_flattened_in_place() {
        // Destructive, and the reference's behaviour --- see the doc comment.
        let d = tmp("restructure");
        std::fs::create_dir_all(d.join("5")).unwrap();
        std::fs::write(d.join("5/corrected_reads.fastq"), "@r\nA\n+\nI\n").unwrap();
        assert_eq!(
            restructure_isoncorrect_output(&d).unwrap(),
            InputStructure::IsonCorrect
        );
        assert_eq!(
            std::fs::read_to_string(d.join("5.fastq")).unwrap(),
            "@r\nA\n+\nI\n"
        );
        assert!(!d.join("5").exists(), "the subdirectory is removed");
    }
}
