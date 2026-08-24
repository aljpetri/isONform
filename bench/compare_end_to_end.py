#!/usr/bin/env python3
"""End-to-end differential check: run the reference and the port on the same
clusters and diff every file `main` writes.

This is the last oracle, and the only one that exercises the whole pipeline at
once. The six stage oracles replay recorded reference *inputs* through one Rust
function; this runs both programs from a fastq and compares what lands on disk.
A stage oracle can pass while the wiring between stages is wrong --- this is what
catches that.

    bench/compare_end_to_end.py --fastq-folder bench/corpus/sirv_small \
        --workdir /tmp/e2e [--parallel]

Three of `main`'s four outputs are Python pickles (`{id}_batch`, `spoa{id}`,
`mapping{id}`); the port writes the same content as tab-separated text, because
a Rust program has no business emitting pickles and nothing downstream of a
finished run reads them. So the comparison unpickles the reference's side rather
than comparing bytes. `skip{id}.fa` is text in both and *is* compared as bytes.
"""
import argparse
import os
import pickle
import subprocess
import sys


def load_ref(d, name):
    with open(os.path.join(d, name), "rb") as fh:
        return pickle.load(fh)


def load_port(d, name):
    out = []
    with open(os.path.join(d, name)) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            k, _, v = line.partition("\t")
            out.append((k, v))
    return out


def compare(ref_dir, port_dir, batch_id):
    """Returns a list of problem strings, empty when the two agree."""
    problems = []

    # skip file: text on both sides, so compare bytes.
    for d in (ref_dir, port_dir):
        p = os.path.join(d, f"skip{batch_id}.fa")
        if not os.path.exists(p):
            problems.append(f"skip{batch_id}.fa missing from {d}")
    if not problems:
        a = open(os.path.join(ref_dir, f"skip{batch_id}.fa"), "rb").read()
        b = open(os.path.join(port_dir, f"skip{batch_id}.fa"), "rb").read()
        if a != b:
            problems.append(
                f"skip{batch_id}.fa differs: {a.count(b'>')} vs {b.count(b'>')} record(s)"
            )

    # The batch file. Its name depends on --parallel, so find whichever exists.
    for cand in (f"{batch_id}_batch", f"{batch_id}batch"):
        if os.path.exists(os.path.join(ref_dir, cand)):
            ref_batch, port_batch = load_ref(ref_dir, cand), load_port(port_dir, cand)
            ref_pairs = [
                (str(k), v[0] if isinstance(v, (list, tuple)) else v)
                for k, v in ref_batch.items()
            ]
            if ref_pairs != port_batch:
                problems.append(
                    f"{cand}: {len(ref_pairs)} vs {len(port_batch)} read(s)"
                    + first_diff(ref_pairs, port_batch)
                )
            break

    # Isoform consensuses.
    ref_spoa = load_ref(ref_dir, f"spoa{batch_id}")
    port_spoa = dict(load_port(port_dir, f"spoa{batch_id}"))
    ref_spoa_s = {str(k): v for k, v in ref_spoa.items()}
    if ref_spoa_s != port_spoa:
        only_ref = sorted(set(ref_spoa_s) - set(port_spoa))
        only_port = sorted(set(port_spoa) - set(ref_spoa_s))
        differing = [k for k in ref_spoa_s if k in port_spoa and ref_spoa_s[k] != port_spoa[k]]
        problems.append(
            f"spoa{batch_id}: {len(ref_spoa_s)} vs {len(port_spoa)} isoform(s); "
            f"{len(differing)} sequence(s) differ, {len(only_ref)} only in reference, "
            f"{len(only_port)} only in port"
        )

    # Read-to-isoform mapping.
    ref_map = load_ref(ref_dir, f"mapping{batch_id}")
    port_map = {k: (v.split(",") if v else []) for k, v in load_port(port_dir, f"mapping{batch_id}")}
    ref_map_s = {str(k): list(v) for k, v in ref_map.items()}
    if ref_map_s != port_map:
        differing = [k for k in ref_map_s if port_map.get(k) != ref_map_s[k]]
        problems.append(
            f"mapping{batch_id}: {len(ref_map_s)} vs {len(port_map)} entry(s), "
            f"{len(differing)} differ (first: {differing[:3]})"
        )
    return problems


def first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return f"; first difference at index {i}: {x[0]!r} vs {y[0]!r}"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-folder", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--port", default="rust/target/release/main")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--parallel", action="store_true",
                    help="pass --parallel True, as isONform_parallel does")
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--w", type=int, default=31)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--max-disagreeing", type=int, default=0,
        help="exit 0 while at most this many clusters disagree. Use ONLY to pin a "
             "divergence that is already measured and written down --- today that "
             "means finding 28's set-order renumbering. It is a regression pin, not "
             "a mute button: raising it needs a matching entry in PORTING.md. A run "
             "that fails to run at all is never tolerated.",
    )
    ap.add_argument(
        "--extra", nargs=argparse.REMAINDER, default=[],
        help="passed through to both programs. `--extra --max_seqs 25` forces a "
             "cluster to split into several batches, which is the only way to "
             "exercise per-batch versus per-file state (PORTING.md finding 33).",
    )
    args = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    fastqs = sorted(
        f for f in os.listdir(args.fastq_folder) if f.endswith((".fastq", ".fq"))
    )
    if args.limit:
        fastqs = fastqs[: args.limit]
    if not fastqs:
        sys.exit(f"no fastq in {args.fastq_folder}")

    ok = bad = failed = 0
    report = []
    for name in fastqs:
        stem = name.rsplit(".", 1)[0]
        src = os.path.join(args.fastq_folder, name)
        ref_dir = os.path.join(args.workdir, stem, "ref")
        port_dir = os.path.join(args.workdir, stem, "port")
        for d in (ref_dir, port_dir):
            os.makedirs(d, exist_ok=True)

        common = ["--fastq", src, "--k", str(args.k), "--w", str(args.w), "--exact"]
        if args.parallel:
            common += ["--parallel", "True"]
        common += list(args.extra)

        r1 = subprocess.run(
            [args.python, os.path.join(root, "main"), *common, "--outfolder", ref_dir],
            capture_output=True, cwd=root,
            env={**os.environ, "PYTHONHASHSEED": "0"},
        )
        r2 = subprocess.run(
            [os.path.join(root, args.port), *common, "--outfolder", port_dir],
            capture_output=True, cwd=root,
        )
        if r1.returncode or r2.returncode:
            failed += 1
            report.append(
                f"\n=== {name} === run failed: reference exit {r1.returncode}, "
                f"port exit {r2.returncode}\n  {r2.stderr.decode()[-400:]}"
            )
            continue

        # Which batch ids exist. Under --parallel every batch overwrites one
        # set of files, so there is exactly one.
        ids = sorted(
            f[4:-3] for f in os.listdir(ref_dir)
            if f.startswith("skip") and f.endswith(".fa")
        )
        problems = []
        for bid in ids:
            problems += compare(ref_dir, port_dir, bid)
        if problems:
            bad += 1
            report.append(f"\n=== {name} === ({len(ids)} batch(es))")
            report += [f"  {p}" for p in problems[:8]]
        else:
            ok += 1

    print(
        f"end-to-end: {len(fastqs)} cluster(s), {ok} agree, {bad} disagree, "
        f"{failed} failed to run"
    )
    if report:
        print("\n".join(report))
    if failed:
        sys.exit(1)
    if bad > args.max_disagreeing:
        print(
            f"FAIL: {bad} disagreeing, at most {args.max_disagreeing} expected",
            file=sys.stderr,
        )
        sys.exit(1)
    if bad:
        print(f"({bad} disagreeing, within the {args.max_disagreeing} allowed)")
    sys.exit(0)


if __name__ == "__main__":
    main()
