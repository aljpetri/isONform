#!/usr/bin/env python3
"""Record a stage's inputs and outputs, from the live driver.

Two stages, both recorded from the same run:

* **graph construction** --- `GraphGeneration.generateGraphfromIntervals`, written
  to `graph_*.txt`;
* **simplification** --- `SimplifyGraph.simplifyGraph`, written to
  `simplify_*.txt` as the graph before and the graph after, since that stage
  mutates in place and returns nothing.

and, with `--record-spoa`, one thing that is not a stage:

* **every `spoa` call**, written to `spoa_cases.tsv` as `consensus<TAB>seq<TAB>…`.
  spoa is reached from both simplification and isoform generation, so it is
  recorded independently of `--stage`. This is what lets `crate::poa` be checked
  against the `spoa` binary on *isONform's* inputs; the equivalence it shipped
  with was measured on isONcorrect's correction intervals, which are different
  sequences from a different stage.
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

The B/A records are in networkx **insertion** order, not sorted: the oracle
rebuilds the graph from them, and node and adjacency insertion order decide the
topological order, which decides which bubbles get found.
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


def write_nodes(fh, dg, tag="N", order="sorted"):
    """Node records.

    `order="insertion"` emits them in networkx's own order --- i.e. the order they
    were added --- which the simplification oracle *needs*, because it rebuilds the
    graph from these records and node insertion order decides the topological
    order, which decides which bubbles are found. Sorted output would make the
    rebuild produce a different (still valid) topological order and every case
    would disagree for the wrong reason.
    """
    nodes = dg.nodes() if order == "insertion" else sorted(dg.nodes())
    for node in nodes:
        attrs = dg.nodes[node]
        reads = attrs.get("reads", {}) or {}
        ems = attrs.get("end_mini_seq", "")
        fh.write(f"{tag} {node} {ems if ems else '-'} {fmt_reads_attr(reads)}\n")


def write_edges(fh, dg, tag="E", order="sorted"):
    """Edge records. See `write_nodes` for why order matters.

    networkx yields edges grouped by source node in node-insertion order, each
    group in edge-insertion order --- so replaying these in sequence reproduces
    each node's adjacency order, which is what the successor-order tie-breaks
    depend on.
    """
    edges = dg.edges() if order == "insertion" else sorted(dg.edges())
    for u, v in edges:
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
        # Insertion order, not sorted --- the oracle rebuilds from these.
        write_nodes(fh, before, "BN", order="insertion")
        write_edges(fh, before, "BE", order="insertion")
        write_topo(fh, before, "BT")
        write_nodes(fh, after, "AN", order="insertion")
        write_edges(fh, after, "AE", order="insertion")
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


def read_fasta_in_order(path):
    """Sequences from a two-line-per-record fasta, **in file order**.

    Order is load-bearing and not a stylistic choice: partial-order alignment is
    order-sensitive, so replaying these through `crate::poa` in any other order
    would compare a different computation. Every writer in the reference emits
    `">{name}\\n{seq}\\n"`, so this does not need a general fasta parser --- but it
    does need to not sort, not deduplicate, and not skip blanks.
    """
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.startswith(">"):
                seqs.append(line)
    return seqs


def install_spoa(outdir):
    """Record every `spoa` call the run makes, as replayable cases.

    Writes `spoa_cases.tsv` in the format `rust/src/poa.rs`'s `oracle` test
    already reads (`SPOA_CASES`): one line per call, `consensus<TAB>seq<TAB>seq…`,
    `#` lines ignored. So this closes the open caveat on `crate::poa` --- verified
    against isONcorrect's correction intervals, never against isONform's inputs
    --- without a new test needing to be written.

    Patching `IsoformGeneration.run_spoa` catches **all** of it. That is worth
    stating rather than assuming, because there are four live call sites in three
    modules: `SimplifyGraph.py:570` (bubble-path consensus), `IsoformGeneration.py`
    at 190, 413 and 436, and `batch_merging_parallel.py:34`. Two reach it as
    `IsoformGeneration.run_spoa(...)` and two as a bare `run_spoa(...)` inside
    `IsoformGeneration` itself, and both forms resolve the name through the module
    dict at call time, so one attribute patch covers every one.
    (`modules/consensus.py:245` and `create_augmented_reference.py` also name a
    `run_spoa`, but those take a third `spoa_path` argument and belong to a module
    nothing imports --- see PORTING.md, reconnaissance correction 2.)

    Identical inputs are recorded once. spoa is deterministic, so a duplicate
    exercises nothing new; the counts printed at the end report how many calls
    collapsed, so the deduplication is visible rather than silent.
    """
    from modules import IsoformGeneration as IG

    original = IG.run_spoa
    path = os.path.join(outdir, "spoa_cases.tsv")
    fh = open(path, "w")
    fh.write(
        "# One `spoa -l 0 -r 0 -g -2` call per line: "
        "<consensus>\\t<seq>\\t<seq>...\n"
        "# Sequence order is the order spoa saw them --- POA is order-sensitive.\n"
    )
    state = {"path": path, "fh": fh, "calls": 0, "unique": 0, "empty_seq": 0,
             "sites": {}, "seen": set()}

    def wrapper(reads, spoa_out_file):
        consensus = original(reads, spoa_out_file)
        state["calls"] += 1
        try:
            seqs = read_fasta_in_order(reads)
        except OSError:
            return consensus
        if not seqs:
            return consensus

        key = (consensus, tuple(seqs))
        if key in state["seen"]:
            return consensus
        state["seen"].add(key)
        state["unique"] += 1
        if any(not s for s in seqs):
            state["empty_seq"] += 1

        # Provenance, as a comment the oracle skips: which call site produced
        # this case. Kept per case rather than aggregated so that a mismatch
        # localises to a stage without regenerating anything.
        frame = sys._getframe(1)
        site = f"{os.path.basename(frame.f_code.co_filename)}:{frame.f_lineno}"
        state["sites"][site] = state["sites"].get(site, 0) + 1

        fh.write(f"# {site}\n")
        fh.write(consensus + "\t" + "\t".join(seqs) + "\n")
        fh.flush()
        return consensus

    IG.run_spoa = wrapper
    return state


def install_parasail(outdir):
    """Record every `parasail_alignment` call as a replayable case.

    Writes `parasail_cases.tsv` in the format `rust/src/parasail.rs`'s `oracle`
    module already reads (`PARASAIL_CASES`): one line per call,
    `s1<TAB>s2<TAB>cigar<TAB>score<TAB>match<TAB>mismatch<TAB>open<TAB>ext`.

    This exists because a claim in `PORTING.md` turned out to be wrong. The
    parasail port was described as needing no re-verification on isONform's
    inputs --- "exact by construction and verified at the CIGAR level" --- where
    the spoa port did. The score is exact by construction; the **CIGAR is not**,
    because a semi-global alignment usually has several optimal paths and which
    one you get is a tie-break. `parse_cigar_diversity` reads the CIGAR, not the
    score, so a tie-break difference changes whether a bubble pops. One did.

    Patching `consensus.parasail_alignment` catches all three live call sites:
    `SimplifyGraph.py:657` (bubble poppability, match=2/mismatch=-8),
    `IsoformGeneration.py:381` (`align_to_merge`, 2/-2) and
    `consensus.py:147` (2/-2 by default). The first two name the module, the
    third calls the bare name inside `consensus` itself; both forms resolve
    through the module dict at call time.
    """
    from modules import consensus as CONS

    original = CONS.parasail_alignment
    path = os.path.join(outdir, "parasail_cases.tsv")
    fh = open(path, "w")
    fh.write(
        "# s1\ts2\tcigar\tscore\tmatch\tmismatch\topen\text --- one "
        "parasail_alignment call per line\n"
    )
    state = {"path": path, "fh": fh, "calls": 0, "unique": 0,
             "sites": {}, "seen": set()}

    def wrapper(s1, s2, match_score=2, mismatch_penalty=-2,
                opening_penalty=12, gap_ext=1):
        res = original(s1, s2, match_score, mismatch_penalty,
                       opening_penalty, gap_ext)
        s1_aln, s2_aln, cigar_string, cigar_tuples, score = res
        state["calls"] += 1

        key = (s1, s2, match_score, mismatch_penalty, opening_penalty, gap_ext)
        if key in state["seen"]:
            return res
        state["seen"].add(key)
        state["unique"] += 1

        frame = sys._getframe(1)
        site = f"{os.path.basename(frame.f_code.co_filename)}:{frame.f_lineno}"
        state["sites"][site] = state["sites"].get(site, 0) + 1

        fh.write(f"# {site}\n")
        fh.write(
            f"{s1}\t{s2}\t{cigar_string}\t{score}\t{match_score}\t"
            f"{mismatch_penalty}\t{opening_penalty}\t{gap_ext}\n"
        )
        fh.flush()
        return res

    CONS.parasail_alignment = wrapper
    return state


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
    ap.add_argument(
        "--record-parasail",
        action="store_true",
        help="also record every `parasail_alignment` call as a replayable case "
        "in parasail_cases.tsv, for rust/src/parasail.rs's PARASAIL_CASES "
        "oracle and its tie-break sweep.",
    )
    ap.add_argument(
        "--record-spoa",
        action="store_true",
        help="also record every `spoa` call as a replayable case in "
        "spoa_cases.tsv, for rust/src/poa.rs's SPOA_CASES oracle. Orthogonal to "
        "--stage: spoa is called from simplification and from isoform "
        "generation, so this is not a property of one stage.",
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
    spoa_state = install_spoa(args.outdir) if args.record_spoa else None
    para_state = install_parasail(args.outdir) if args.record_parasail else None

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
    if spoa_state is not None:
        spoa_state["fh"].close()
        sites = ", ".join(
            f"{s}={n}" for s, n in sorted(spoa_state["sites"].items())
        )
        print(
            f"[dump] {spoa_state['calls']} spoa call(s), "
            f"{spoa_state['unique']} unique -> {spoa_state['path']}"
            + (f" [{sites}]" if sites else ""),
            file=sys.stderr,
        )
        if spoa_state["empty_seq"]:
            print(
                f"[dump] NOTE {spoa_state['empty_seq']} case(s) contain an empty "
                "sequence; recorded as-is rather than skipped, so the oracle "
                "judges them rather than this script.",
                file=sys.stderr,
            )
    if para_state is not None:
        para_state["fh"].close()
        sites = ", ".join(
            f"{s}={n}" for s, n in sorted(para_state["sites"].items())
        )
        print(
            f"[dump] {para_state['calls']} parasail call(s), "
            f"{para_state['unique']} unique -> {para_state['path']}"
            + (f" [{sites}]" if sites else ""),
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
