#!/usr/bin/env python3
"""Record a stage's inputs and outputs, from the live driver.

Two stages, both recorded from the same run:

* **graph construction** --- `GraphGeneration.generateGraphfromIntervals`, written
  to `graph_*.txt`;
* **simplification** --- `SimplifyGraph.simplifyGraph`, written to
  `simplify_*.txt` as the graph before and the graph after, since that stage
  mutates in place and returns nothing.
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

One text file per call, `graph_<call>.txt`, with sectioned records so a diff
points at a line rather than a byte offset in a pickle. Records are ordered
deterministically, and *within* a record, read lists and edge-support lists keep
their **insertion** order rather than being sorted. That is stricter than the
reference currently requires — the order is not observable downstream, see
`rust/src/graph.rs` — but a stricter dump costs nothing and catches a divergence
the moment one matters.

    # params k=<k> delta_len=<d>
    R <r_id> <sequence>                     one per read the graph sees
    L <r_id> <length>                       read_len_dict
    I <r_id> <pos> <start> <end> <weight> <a,b,c,...>    intervals, in call order
    N <node> <end_mini_seq> <r_id:start:end:orig,...>    nodes sorted by name;
                                                         reads in INSERTION order
    E <u> -> <v> <length> <r_id,...>                     edges sorted by (u,v);
                                                         support in INSERTION order
    T <node>|<node>|...                     nx.topological_sort order (see below)
    F <r_id,...>                            reads_for_isoforms
    S <n_nodes> <n_edges>

and for `simplify_*.txt`:

    # params k=<k> delta_len=<d> mode=<slow>
    R <r_id> <sequence>
    BN / BE / BT                            the graph on entry (nodes/edges/topo)
    AN / AE / AT                            the graph on exit
    S <before_n> <before_e> <after_n> <after_e>

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
    """A node's `reads` dict -> a string, **in insertion order**.

    Not sorted, deliberately — but for a weaker reason than first claimed here.
    Insertion order of this map turns out **not** to be observable in the current
    reference: every consumer either sorts it or aggregates commutatively (see
    `rust/src/graph.rs` for the audit of all five). Recording it unsorted is kept
    anyway, because it makes the oracle strictly stricter at no cost, and because
    a stricter dump catches a divergence early if a future consumer does start
    depending on the order.

    `Read_infos` is a namedtuple `(start_mini_end, end_mini_start,
    original_support)`.
    """
    parts = []
    for r_id, ri in reads.items():
        parts.append(f"{r_id}:{ri.start_mini_end}:{ri.end_mini_start}:{int(bool(ri.original_support))}")
    return ",".join(parts) if parts else "-"


def write_nodes(fh, dg, tag="N"):
    for node in sorted(dg.nodes()):
        attrs = dg.nodes[node]
        reads = attrs.get("reads", {}) or {}
        ems = attrs.get("end_mini_seq", "")
        fh.write(f"{tag} {node} {ems if ems else '-'} {fmt_reads_attr(reads)}\n")


def write_edges(fh, dg, tag="E"):
    for u, v in sorted(dg.edges()):
        d = dg[u][v]
        length = d.get("length", "NA")
        supp = d.get("edge_supp", []) or []
        # Insertion order again, not sorted, for the same reason as above:
        # stricter than currently necessary, and free.
        fh.write(f"{tag} {u} -> {v} {length} {','.join(str(x) for x in supp)}\n")


def write_topo(fh, dg, tag="T"):
    """The topological order the next stage will actually use.

    Recorded because a topological order is NOT unique and the reference compares
    the resulting indices to decide which node pairs are candidate bubbles
    (SimplifyGraph.py:115) --- so a port producing a different valid order finds
    different bubbles. networkx 2.8.4's order is generation-based Kahn seeded in
    node insertion order; this line is what proves the port reproduces it rather
    than merely producing *a* valid order.
    """
    try:
        import networkx as _nx
        fh.write(f"{tag} " + "|".join(_nx.topological_sort(dg)) + "\n")
    except Exception as exc:  # noqa: BLE001 --- a cycle is informative, not fatal
        fh.write(f"{tag} CYCLIC {exc}\n")


def dump_simplify_call(path, k_size, delta_len, mode, all_reads, before, after):
    """One `simplifyGraph` call: the graph in, the graph out.

    Simplification mutates `DG` in place, so `before` is a deep copy taken on
    entry. The B/A prefixes keep the two graphs apart in one file so a diff shows
    what the stage did rather than only what it produced.

    **This stage runs spoa.** `generate_consensus_path` shells out to it whenever a
    bubble path has more than two supporting reads, and the consensus decides how
    the bubble is linearised — so an exact oracle for simplification is an oracle
    for spoa too. Reproducing it is `poa.rs`/spoars' job (the invocation is the
    same `-l 0 -r 0 -g -2` the isONcorrect port already matches). Bubbles with
    exactly two supporting reads skip spoa entirely and just take the longer
    subsequence, so a corpus of those is the place to start.
    """
    with open(path, "w") as fh:
        fh.write(f"# params k={k_size} delta_len={delta_len} mode={mode}\n")
        for r_id in sorted(all_reads):
            fh.write(f"R {r_id} {all_reads[r_id][1]}\n")
        write_nodes(fh, before, "BN")
        write_edges(fh, before, "BE")
        write_topo(fh, before, "BT")
        write_nodes(fh, after, "AN")
        write_edges(fh, after, "AE")
        write_topo(fh, after, "AT")
        fh.write(
            f"S {before.number_of_nodes()} {before.number_of_edges()} "
            f"{after.number_of_nodes()} {after.number_of_edges()}\n"
        )


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

        write_nodes(fh, dg, "N")
        write_edges(fh, dg, "E")

        write_topo(fh, dg, "T")

        fh.write("F " + ",".join(str(x) for x in reads_for_isoforms) + "\n")
        fh.write(f"S {dg.number_of_nodes()} {dg.number_of_edges()}\n")


def install_simplify(outdir):
    """Replace `simplifyGraph` with a recording wrapper.

    `simplifyGraph` mutates `DG` in place and returns nothing, so the "before"
    graph has to be deep-copied on entry --- there is no other way to see what the
    stage changed.
    """
    from modules import SimplifyGraph as SG

    original = SG.simplifyGraph
    state = {"n": 0}

    def wrapper(DG, all_reads, work_dir, k_size, delta_len, mode):
        before = copy.deepcopy(DG)
        snap_reads = {r: tuple(v) for r, v in all_reads.items()}
        result = original(DG, all_reads, work_dir, k_size, delta_len, mode)
        state["n"] += 1
        path = os.path.join(outdir, f"simplify_{state['n']:04d}.txt")
        dump_simplify_call(path, k_size, delta_len, mode, snap_reads, before, DG)
        print(
            f"[dump] {path}: {before.number_of_nodes()}n/{before.number_of_edges()}e "
            f"-> {DG.number_of_nodes()}n/{DG.number_of_edges()}e",
            file=sys.stderr,
        )
        return result

    SG.simplifyGraph = wrapper
    return original


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
    ap.add_argument(
        "--stage",
        choices=("graph", "simplify", "both"),
        default="both",
        help="which stage(s) to record. `graph` writes graph_*.txt, `simplify` "
        "writes simplify_*.txt. Both by default: they come from the same run, "
        "so recording one alone wastes the other.",
    )
    ap.add_argument("--k", type=int, default=20)
    ap.add_argument("--w", type=int, default=31)
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[])
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sys.path.insert(0, REPO_ROOT)
    if args.stage in ("graph", "both"):
        install(args.outdir)
    if args.stage in ("simplify", "both"):
        install_simplify(args.outdir)

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

    n_graph = len([f for f in os.listdir(args.outdir) if f.startswith("graph_")])
    n_simpl = len([f for f in os.listdir(args.outdir) if f.startswith("simplify_")])
    print(
        f"[dump] {n_graph} graph-stage and {n_simpl} simplify-stage records "
        f"in {args.outdir}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
