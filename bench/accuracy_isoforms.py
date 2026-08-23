#!/usr/bin/env python3
"""Isoform-reconstruction accuracy against a known transcriptome.

The companion to isONcorrect's `bench/accuracy.py`, but the question is a
different one and so is the metric. isONcorrect emits *corrected reads*, one per
input read, so accuracy is per-read error against the read's source transcript.
isONform emits *isoforms* — a reconstructed transcriptome whose size is an output,
not an input. Per-read error does not apply; what matters is:

* did it find the transcripts that are really there,  (**recall**)
* is what it emitted real,                            (**precision**)
* and is each emitted sequence *accurate*.            (**per-isoform identity**)

All three are needed together. Any one alone is trivially gamed: emit one
high-confidence isoform and precision and identity look perfect while recall
collapses; emit hundreds and recall rises while precision collapses.

# How truth is assigned

Each reconstructed isoform is aligned against every reference transcript with
edlib infix mode (`HW`), which finds the isoform's best placement inside a
transcript, and is assigned to its best match. Ties break by transcript name so
the choice is deterministic.

**Assignment is per-isoform and independent per set, and unlike `accuracy.py`
that is safe here.** There, letting each implementation pick its own
best-matching transcript would flatter every one of them, because a correction
that dragged a read toward the wrong isoform would score *better*. Here the unit
being scored is not shared between sets — each set emits its own isoforms — so
there is no paired comparison to corrupt. What replaces that safeguard is
reporting recall and precision alongside identity, so a set cannot win by
assigning itself an easy subset.

# The two thresholds, and why both are reported

A reconstructed isoform "matches" a transcript when identity over the aligned
region is at least `--min-identity` (default 0.95) **and** its length is within
`--len-tol` of the transcript's (default 0.10, i.e. 10%). The length condition is
not optional: infix alignment happily places a 300 bp fragment inside a 2 000 bp
transcript at 99% identity, and calling that a recovered transcript would be
wrong. Results are reported at a strict and a lenient threshold both, because a
single cutoff hides whether a change moved sequences or moved them across a line.

# Metrics

    recall            fraction of *expressed* reference transcripts matched by
                      at least one reconstructed isoform
    precision         fraction of reconstructed isoforms that match some
                      reference transcript
    redundancy        reconstructed isoforms per matched transcript; 1.0 is
                      ideal, higher means the same transcript reported more
                      than once
    identity          median and mean identity of matching isoforms
    length ratio      reconstructed length / transcript length, to expose
                      systematic truncation

"Expressed" matters: a SIRV transcript with no reads in the corpus cannot be
recovered and counting it against recall measures the corpus, not the tool. Pass
`--expressed-from` with the *input* reads (simulated headers name their source
transcript) to restrict the denominator; without it, every reference transcript
counts and recall is a lower bound.

# Usage

    conda activate isonform-ref
    bench/accuracy_isoforms.py \\
        --transcriptome sirv_transcriptome.fasta \\
        --isoforms lex=/tmp/out_lex/transcriptome.fasta \\
        --isoforms hash=/tmp/out_hash/transcriptome.fasta \\
        --expressed-from /path/to/reads_10k_err7%.fastq

Each `--isoforms NAME=PATH` is a fasta (or a folder containing
`transcriptome.fasta`). Multiple sets are compared side by side.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import statistics
import sys

try:
    import edlib
except ImportError:  # pragma: no cover
    sys.exit(
        "error: edlib is required.\n"
        "  The isonform-ref environment does not include it, because isONform\n"
        "  itself never imports it. Install it for evaluation only:\n"
        "      pip install edlib\n"
        "  (isONcorrect's isoncorrect-ref environment already has it.)"
    )


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------


def read_fasta(path):
    """name -> sequence. Names are truncated at the first whitespace."""
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line.upper())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def read_fastq_headers(path):
    """Header lines of a fastq, without the leading '@'."""
    out = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                out.append(line[1:].rstrip())
    return out


def resolve_isoform_file(path):
    """Accept a fasta, or a folder that contains isONform's transcriptome.fasta."""
    if os.path.isfile(path):
        return path
    cand = os.path.join(path, "transcriptome.fasta")
    if os.path.isfile(cand):
        return cand
    hits = sorted(glob.glob(os.path.join(path, "*.fa"))) + sorted(
        glob.glob(os.path.join(path, "*.fasta"))
    )
    if len(hits) == 1:
        return hits[0]
    sys.exit(f"error: could not find a single fasta under {path} (found {len(hits)})")


# --------------------------------------------------------------------------
# Assignment
# --------------------------------------------------------------------------


def best_match(query, refs, ref_names):
    """Best (identity, ref_name, aligned_len) for `query` over `refs`.

    edlib `HW` mode is infix: it finds the placement of `query` inside a
    reference that minimises edit distance, leaving the reference's flanks free.
    That is the right mode here — a reconstructed isoform may legitimately be a
    little short at either end — and it is exactly why the length check in
    `classify` is not optional.

    Ties in edit distance break by reference name, so the result does not depend
    on dict order.
    """
    best = None
    for name in ref_names:
        r = refs[name]
        res = edlib.align(query, r, mode="HW", task="locations")
        ed = res["editDistance"]
        if ed < 0:
            continue
        locs = res.get("locations") or [(0, len(r) - 1)]
        start, end = locs[0]
        start = 0 if start is None else start
        aligned_len = end - start + 1
        # Identity over the alignment: the query is placed in full, so the
        # denominator is the longer of query and the reference window.
        denom = max(len(query), aligned_len)
        identity = 1.0 - (ed / denom) if denom else 0.0
        key = (-identity, name)
        if best is None or key < best[0]:
            best = (key, identity, name, aligned_len)
    if best is None:
        return 0.0, None, 0
    return best[1], best[2], best[3]


def classify(identity, iso_len, ref_len, min_identity, len_tol):
    """Does this isoform count as a reconstruction of that transcript?

    Both conditions, deliberately. Infix alignment will place a short fragment
    inside a long transcript at high identity, and counting that as a recovered
    transcript would make truncation look like success.
    """
    if identity < min_identity:
        return False
    if ref_len == 0:
        return False
    return abs(iso_len - ref_len) / ref_len <= len_tol


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def evaluate(isoforms, refs, expressed, min_identity, len_tol):
    ref_names = sorted(refs)
    matched_refs, rows = set(), []

    for iso_name in sorted(isoforms):
        iso = isoforms[iso_name]
        identity, ref, aligned_len = best_match(iso, refs, ref_names)
        ref_len = len(refs[ref]) if ref else 0
        is_match = classify(identity, len(iso), ref_len, min_identity, len_tol)
        if is_match:
            matched_refs.add(ref)
        rows.append(
            {
                "isoform": iso_name,
                "len": len(iso),
                "ref": ref,
                "ref_len": ref_len,
                "identity": identity,
                "match": is_match,
            }
        )

    n_iso = len(rows)
    n_match = sum(1 for r in rows if r["match"])
    denom_refs = expressed if expressed else set(ref_names)
    # A match against a transcript outside the expressed set still counts as a
    # true positive for precision — the transcript is real, our estimate of
    # which ones are expressed is what is approximate.
    recall = len(matched_refs & denom_refs) / len(denom_refs) if denom_refs else 0.0
    precision = n_match / n_iso if n_iso else 0.0
    redundancy = n_match / len(matched_refs) if matched_refs else 0.0

    ids = [r["identity"] for r in rows if r["match"]]
    ratios = [r["len"] / r["ref_len"] for r in rows if r["match"] and r["ref_len"]]

    return {
        "n_isoforms": n_iso,
        "n_matching": n_match,
        "n_refs_hit": len(matched_refs & denom_refs),
        "n_refs_denom": len(denom_refs),
        "recall": recall,
        "precision": precision,
        "redundancy": redundancy,
        "identity_median": statistics.median(ids) if ids else float("nan"),
        "identity_mean": statistics.fmean(ids) if ids else float("nan"),
        "len_ratio_median": statistics.median(ratios) if ratios else float("nan"),
        "rows": rows,
    }


def f1(p, r):
    return 2 * p * r / (p + r) if (p + r) else 0.0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--transcriptome", required=True, help="reference transcripts, fasta")
    ap.add_argument(
        "--isoforms",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="a set of reconstructed isoforms; repeat to compare sets",
    )
    ap.add_argument(
        "--expressed-from",
        help="fastq of the INPUT reads. Simulated headers name their source "
        "transcript (read_1_from_SIRV612), which restricts the recall "
        "denominator to transcripts that are actually present. Without it "
        "every reference transcript counts and recall is a lower bound.",
    )
    ap.add_argument("--min-identity", type=float, default=0.95)
    ap.add_argument("--len-tol", type=float, default=0.10)
    ap.add_argument(
        "--lenient",
        nargs=2,
        type=float,
        default=(0.90, 0.20),
        metavar=("IDENTITY", "LEN_TOL"),
        help="second, looser threshold reported alongside the strict one "
        "(default 0.90 0.20). Reporting both shows whether a change moved "
        "sequences or only moved them across a line.",
    )
    ap.add_argument("--per-isoform", action="store_true", help="dump every isoform's best match")
    args = ap.parse_args()

    refs = read_fasta(args.transcriptome)
    if not refs:
        sys.exit(f"error: no sequences in {args.transcriptome}")

    expressed = set()
    if args.expressed_from:
        pat = re.compile(r"_from_(\S+)")
        unmatched = 0
        for h in read_fastq_headers(args.expressed_from):
            m = pat.search(h)
            if m:
                expressed.add(m.group(1))
            else:
                unmatched += 1
        if not expressed:
            print(
                f"warning: no '_from_<transcript>' headers in {args.expressed_from}; "
                "recall is against the whole transcriptome and is a lower bound.",
                file=sys.stderr,
            )
        elif unmatched:
            print(
                f"note: {unmatched} headers carried no source transcript and were ignored.",
                file=sys.stderr,
            )
        # Keep only names we have a reference sequence for.
        missing = expressed - set(refs)
        if missing:
            print(
                f"note: {len(missing)} source transcripts are absent from the "
                f"transcriptome and were dropped from the denominator "
                f"(e.g. {sorted(missing)[:3]}).",
                file=sys.stderr,
            )
        expressed &= set(refs)

    sets = []
    for spec in args.isoforms:
        if "=" not in spec:
            sys.exit(f"error: --isoforms wants NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        sets.append((name, resolve_isoform_file(path)))

    print(f"reference:  {args.transcriptome}  ({len(refs)} transcripts)")
    if expressed:
        print(f"expressed:  {len(expressed)} transcripts, from {args.expressed_from}")
    else:
        print("expressed:  unknown — recall is against the full transcriptome (lower bound)")
    print()

    strict = (args.min_identity, args.len_tol)
    lenient = tuple(args.lenient)

    for label, (mi, lt) in (("strict", strict), ("lenient", lenient)):
        print(f"== {label}: identity >= {mi:.2f}, length within {lt:.0%} ==")
        print(
            f"  {'set':<14} {'isoforms':>8} {'matching':>8} {'recall':>8} "
            f"{'precis':>8} {'F1':>7} {'redund':>7} {'identity':>9} {'len.ratio':>9}"
        )
        for name, path in sets:
            isoforms = read_fasta(path)
            r = evaluate(isoforms, refs, expressed, mi, lt)
            print(
                f"  {name:<14} {r['n_isoforms']:>8} {r['n_matching']:>8} "
                f"{r['recall']:>7.1%} {r['precision']:>7.1%} "
                f"{f1(r['precision'], r['recall']):>7.3f} "
                f"{r['redundancy']:>7.2f} {r['identity_median']:>9.4f} "
                f"{r['len_ratio_median']:>9.3f}"
            )
            if label == "strict":
                print(
                    f"                 ({r['n_refs_hit']} of {r['n_refs_denom']} "
                    f"reference transcripts recovered)"
                )
                if args.per_isoform:
                    for row in r["rows"]:
                        flag = "match" if row["match"] else "     "
                        print(
                            f"      {flag} {row['isoform']:<16} len={row['len']:>5} "
                            f"-> {str(row['ref']):<12} ref_len={row['ref_len']:>5} "
                            f"id={row['identity']:.4f}"
                        )
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
