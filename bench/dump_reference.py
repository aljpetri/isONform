#!/usr/bin/env python3
"""Record the graph-construction stage's inputs and outputs, from the live driver.

This is the differential oracle for `GraphGeneration.generateGraphfromIntervals`.
The Rust port replays these recorded inputs and diffs its graph against the
recorded output, which localises a disagreement to one stage instead of leaving
"the transcriptome differs" as the only signal.

# Why it wraps rather than edits

`modules/GraphGeneration.py` is not modified. This script imports it, replaces
`generateGraphfromIntervals` with a recording wrapper, and then runs `main`
through `runpy`. Two consequences that matter:

* **The dump comes from the real driver**, not from a standalone harness that
  reconstructs plausible inputs. That distinction caught the one genuine bug in
  the isONcorrect port: it passed every oracle built from recorded inputs and
  was found only by diffing a dump taken from the running driver, because a port
  that follows a different trajectory still satisfies replayed inputs.
* Inputs are **deep-copied on entry**, before the reference touches them. That is
  not defensive habit: `convert_array_to_hash` *mutates* the interval arrays it
  is given, popping their first three elements
  (`GraphGeneration.py:129-131`). Recording after the call would record
  something the function never saw.

# Format

One text file per call, `graph_<batch>_<call>.txt`, with sectioned records so a
diff points at a line rather than a byte offset in a pickle. All sections are
sorted deterministically; nothing here depends on dict iteration order.

    # params k=<k> delta_len=<d>
    R <r_id> <sequence>                     one per read the graph sees
    L <r_id> <length>                       read_len_dict
    I <r_id> <pos> <start> <end> <weight> <a,b,c,...>    intervals, in call order
    N <node> <end_mini_seq> <r_id:start:end:orig,...>    nodes, sorted
    E <u> -> <v> <length> <r_id,...>                     edges, sorted
    F <r_id,...>                            reads_for_isoforms
    S <n_nodes> <n_edges>

Node names are the reference's own (`"<start>, <end>, <r_id>"`, plus `s` and
`t`), because the port has to reproduce the *graph*, not the naming: keeping the
reference's names in the dump is what lets a port using integer node ids diff
against it after mapping.

# Usage

    bench/dump_reference.py --fastq bench/corpus/sirv_small/0.fastq --outdir /tmp/d
    bench/dump_reference.py --fastq-folder bench/corpus/sirv_small --outdir /tmp/d

`PYTHONHASHSEED` is exported by this script before it re-executes itself.
Minimizer selection is lexicographic now so the reference no longer depends on
it, but `convert_array_to_hash` still calls `hash()` on a tuple of ints, and
pinning the seed costs nothing and keeps the dumps comparable to older ones.
"""

from __future__ import annotations

import argparse
import copy
import os
import runpy
import sys

# CPython reads PYTHONHASHSEED once, at startup, so it cannot be set from inside
# a running interpreter --- re-exec instead. This is the mistake that sat in
# isONcorrect's dump script for a while and made every dump silently
# randomly-seeded.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def fmt_reads_attr(reads):
    """A node's `reads` dict -> a stable string.

    `Read_infos` is a namedtuple `(start_mini_end, end_mini_start,
    original_support)`.
    """
    parts = []
    for r_id in sorted(reads):
        ri = reads[r_id]
        parts.append(f"{r_id}:{ri.start_mini_end}:{ri.end_mini_start}:{int(bool(ri.original_support))}")
    return ",".join(parts) if parts else "-"


def dump_call(path, k, delta_len, intervals, read_len_dict, all_reads, dg, reads_for_isoforms):
    with open(path, "w") as fh:
        fh.write(f"# params k={k} delta_len={delta_len}\n")

        for r_id in sorted(all_reads):
            fh.write(f"R {r_id} {all_reads[r_id][1]}\n")
        for r_id in sorted(read_len_dict):
            fh.write(f"L {r_id} {read_len_dict[r_id]}\n")

        # Interval order within a read is the order the reference iterates, and
        # it decides which node is created first, so it is preserved verbatim
        # rather than sorted.
        for r_id in sorted(intervals):
            for pos, inter in enumerate(intervals[r_id]):
                arr = ",".join(str(x) for x in inter[3])
                fh.write(f"I {r_id} {pos} {inter[0]} {inter[1]} {inter[2]} {arr}\n")

        for node in sorted(dg.nodes()):
            attrs = dg.nodes[node]
            reads = attrs.get("reads", {}) or {}
            ems = attrs.get("end_mini_seq", "")
            fh.write(f"N {node} {ems if ems else '-'} {fmt_reads_attr(reads)}\n")

        for u, v in sorted(dg.edges()):
            d = dg[u][v]
            length = d.get("length", "NA")
            supp = d.get("edge_supp", []) or []
            fh.write(f"E {u} -> {v} {length} {','.join(str(x) for x in sorted(supp))}\n")

        fh.write("F " + ",".join(str(x) for x in reads_for_isoforms) + "\n")
        fh.write(f"S {dg.number_of_nodes()} {dg.number_of_edges()}\n")


def install(outdir):
    """Replace `generateGraphfromIntervals` with a recording wrapper."""
    from modules import GraphGeneration as GG

    original = GG.generateGraphfromIntervals
    state = {"n": 0}

    def wrapper(all_intervals_for_graph, k, delta_len, read_len_dict, all_reads):
        # Deep-copy before the call: convert_array_to_hash mutates the interval
        # arrays in place.
        snap_intervals = copy.deepcopy(all_intervals_for_graph)
        snap_lens = dict(read_len_dict)
        snap_reads = {r: tuple(v) for r, v in all_reads.items()}

        result = original(all_intervals_for_graph, k, delta_len, read_len_dict, all_reads)
        dg, reads_for_isoforms = result

        state["n"] += 1
        path = os.path.join(outdir, f"graph_{state['n']:04d}.txt")
        dump_call(path, k, delta_len, snap_intervals, snap_lens, snap_reads, dg, reads_for_isoforms)
        print(f"[dump] {path}: {dg.number_of_nodes()} nodes, "
              f"{dg.number_of_edges()} edges, {len(snap_intervals)} reads", file=sys.stderr)
        return result

    GG.generateGraphfromIntervals = wrapper

    # `main` does `from modules import ... GraphGeneration ...` and then calls
    # `GraphGeneration.generateGraphfromIntervals(...)`, i.e. it looks the
    # attribute up on the module at call time. Patching the module attribute is
    # therefore enough; there is no separately-bound name to chase.
    return original


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--fastq", help="single cluster; drives `main` in-process")
    src.add_argument(
        "--fastq-folder",
        help="folder of clusters. Runs `main` once per cluster in-process, "
        "rather than going through isONform_parallel, because the driver "
        "spawns children and a monkeypatch does not survive a fork+exec.",
    )
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--w", type=int, default=31)
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sys.path.insert(0, REPO_ROOT)
    install(args.outdir)

    if args.fastq:
        targets = [args.fastq]
    else:
        import glob
        targets = sorted(glob.glob(os.path.join(args.fastq_folder, "*.fastq")))
        if not targets:
            sys.exit(f"no *.fastq under {args.fastq_folder}")

    for i, fq in enumerate(targets):
        workdir = os.path.join(args.outdir, f"run_{i}")
        os.makedirs(workdir, exist_ok=True)
        sys.argv = [
            "main", "--fastq", os.path.abspath(fq), "--outfolder", workdir,
            "--k", str(args.k), "--w", str(args.w),
        ] + list(args.extra)
        print(f"[dump] running main on {fq}", file=sys.stderr)
        cwd = os.getcwd()
        try:
            os.chdir(REPO_ROOT)
            runpy.run_path(os.path.join(REPO_ROOT, "main"), run_name="__main__")
        except SystemExit as e:
            if e.code not in (0, None):
                print(f"[dump] main exited {e.code} on {fq}", file=sys.stderr)
        finally:
            os.chdir(cwd)

    n = len([f for f in os.listdir(args.outdir) if f.startswith("graph_")])
    print(f"[dump] {n} graph-stage records in {args.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
