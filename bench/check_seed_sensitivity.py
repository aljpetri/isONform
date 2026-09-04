#!/usr/bin/env python3
"""Run the reference under many `PYTHONHASHSEED` values and report what changes.

Why this exists, and why it runs before anything else is recorded
----------------------------------------------------------------

`main` selects minimizers by **`hash()` of the k-mer string** (`get_kmer_minimizers`,
`main:85`):

    window_kmers = deque([hash(seq[i:i+k_size]) for i in range(w + 1)])

CPython's `str.__hash__` is SipHash keyed by `PYTHONHASHSEED`, which is random
per process unless the variable is set in the environment *before* the
interpreter starts. So every minimizer, every minimizer combination, every
interval, the graph, the popped bubbles and the emitted isoforms depend on the
seed. isONcorrect's `get_kmer_minimizers_lex` compares the k-mer **strings**;
isONform is otherwise character-for-character the same function, so the `hash()`
is the whole difference.

`main` does try to pin it, at `main:634`:

    os.environ['PYTHONHASHSEED'] = '0'

That line cannot work. CPython reads `PYTHONHASHSEED` once at startup, so
assigning it from inside a running process affects only child processes — and
`main` spawns none that care. `isONform_parallel` has the same intent as a
commented-out `#PYTHONHASHSEED = 0` at module scope, which is a no-op twice
over: it never executed, and a module-level assignment would not have been read
by anything.

Consequences, all measured on `bench/corpus/sirv_small` (see PORTING.md):

  * 8 seeds gave 8 different `mapping0`, `skip0.fa` and `spoa0`.
  * The divergence starts at the first number the tool prints — minimizer
    combinations generated ranged 10 549 to 12 532 across those 8 seeds.
  * Reads skipped for low interval abundance ranged 14 to 17 of 100.
  * `transcriptome.fasta` from `isONform_parallel` differs between plain runs.
  * Cluster 0 and cluster 1 of `sirv_small` are byte-identical inputs, and in
    one run they produced *different* isoforms — each child process draws its
    own seed. That is a self-contained proof needing no second implementation.

So: goldens are only meaningful with the seed pinned, every bench invocation
pins it, and the port cannot reproduce seeded SipHash. Making minimizer
selection deterministic is an upstream fix, not a porting decision.

Usage
-----

    bench/check_seed_sensitivity.py --fastq bench/corpus/sirv_small/0.fastq
    bench/check_seed_sensitivity.py --fastq-folder bench/corpus/sirv_small --seeds 24

Exits non-zero when output varies across seeds, so it works as a regression
check once minimizer selection is deterministic.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def digest(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def outputs(outdir: str) -> dict[str, str]:
    """Every emitted file, relative path -> digest.

    Deliberately not a filter on known names: a stage that starts emitting
    something new should show up here rather than be silently ignored.
    """
    got = {}
    for root, _dirs, files in os.walk(outdir):
        for name in files:
            if name.startswith(".bench-"):
                continue
            full = os.path.join(root, name)
            got[os.path.relpath(full, outdir)] = digest(full)
    return got


def run_once(entry: str, target: str, seed: int | None, workdir: str, extra: list[str]) -> dict[str, str]:
    outdir = os.path.join(workdir, "out")
    os.makedirs(outdir, exist_ok=True)

    if entry == "main":
        cmd = [sys.executable, "main", "--fastq", target, "--outfolder", outdir]
    else:
        # The parallel driver mutates its input directory (restructure_isoncorrect_output
        # moves files and removes folders), so it always gets a private copy.
        indir = os.path.join(workdir, "in")
        shutil.copytree(target, indir)
        cmd = [sys.executable, "isONform_parallel", "--fastq_folder", indir,
               "--outfolder", outdir, "--t", "1"]
    cmd += extra

    env = dict(os.environ)
    if seed is None:
        env.pop("PYTHONHASHSEED", None)
    else:
        env["PYTHONHASHSEED"] = str(seed)

    proc = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"reference failed (seed={seed}, exit={proc.returncode}):\n{proc.stderr[-2000:]}\n")
        sys.exit(2)
    return outputs(outdir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--fastq", help="single cluster fastq; drives `main`")
    src.add_argument("--fastq-folder", help="folder of cluster fastqs; drives `isONform_parallel`")
    ap.add_argument("--seeds", type=int, default=8,
                    help="number of seeds to try, 0..N-1 (default 8). "
                         "isONcorrect needed 24 before its own defect stopped "
                         "looking like a porting bug; prefer more than feels necessary.")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="extra flags passed through to the entry point")
    args = ap.parse_args()

    entry = "main" if args.fastq else "isONform_parallel"
    target = os.path.abspath(args.fastq or args.fastq_folder)

    # `values[relpath][digest] = [seeds]`
    values: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    seen_files: set[str] = set()

    print(f"entry point: {entry}")
    print(f"input:       {target}")
    print(f"seeds:       0..{args.seeds - 1}")
    print()

    for seed in range(args.seeds):
        with tempfile.TemporaryDirectory(prefix="isonform-seed-") as workdir:
            got = run_once(entry, target, seed, workdir, args.extra)
        seen_files |= set(got)
        for rel, dig in got.items():
            values[rel][dig].append(seed)
        print(f"  seed {seed:<3} " + "  ".join(f"{r}={d}" for r, d in sorted(got.items())))

    print()
    unstable = {rel: v for rel, v in values.items() if len(v) > 1}
    missing = {rel for rel in seen_files if sum(len(s) for s in values[rel].values()) != args.seeds}

    for rel in sorted(seen_files):
        n = len(values[rel])
        note = "STABLE" if n == 1 else f"{n} DISTINCT VALUES over {args.seeds} seeds"
        extra = "  (not emitted by every seed)" if rel in missing else ""
        print(f"  {rel:<28} {note}{extra}")

    print()
    if unstable or missing:
        print(f"SEED-DEPENDENT: {len(unstable)} of {len(seen_files)} output files vary with PYTHONHASHSEED.")
        print("Goldens are only comparable with the seed pinned. See the module docstring.")
        return 1

    print(f"deterministic across {args.seeds} seeds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
