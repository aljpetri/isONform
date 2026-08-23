#!/usr/bin/env python3
"""Isoform-reconstruction accuracy against a genome, via spliced alignment.

The companion to `accuracy_isoforms.py`. That one needs a trustworthy
transcriptome and scores each reconstructed isoform against one transcript with
edlib; this one needs only a **genome** and lets `minimap2 -ax splice` place each
isoform across exons. Use it for real ONT cDNA from a species whose annotation
you do not want to depend on — the Drosophila set, here.

    bench/accuracy_isoforms_genome.py --genome fruitfly.fa \\
        --isoforms lex=/tmp/droso_lex/transcriptome.fasta \\
        --isoforms hash=/tmp/droso_hash/transcriptome.fasta

With `--annotation` it additionally classifies every isoform into the SQANTI
structural categories (FSM / ISM / NIC / NNC / genic / antisense / fusion /
intergenic), which is a much sharper instrument than canonical-splice fraction:

    bench/accuracy_isoforms_genome.py --genome fruitfly.fa \\
        --annotation dmel-all-r6.68.gff \\
        --isoforms lex=/tmp/droso_lex/transcriptome.fasta

See bench/annotation.py for what each category means and why FSM counts are
reported as absolute numbers plus a per-detected-gene rate rather than as a
single "recall" percentage.

# What is measured, and why each number is here

Per isoform, from the primary alignment's CIGAR and `NM` tag:

* **error rate** = `NM / aligned_isoform_bases`. Sequence accuracy of the
  consensus.
* **aligned fraction** = `aligned_isoform_bases / isoform_length`. Soft clips do
  not count as aligned.

**Both, never one.** Error rate alone is trivially gamed by aligning less: an
isoform whose ends were mangled would soft-clip them away and *improve* its
error rate. Together they cannot move that way.

Then, because this is isoform reconstruction and not read correction, three more
that per-read metrics have no analogue for:

* **canonical splice fraction** — of the introns each isoform implies (CIGAR `N`
  operations), how many have canonical `GT..AG` dinucleotides on the genome (or
  `CT..AC`, the reverse-strand form). This is the strongest annotation-free
  quality signal available: a correctly reconstructed junction lands on a real
  splice site, a junction invented by mis-assembly generally does not. **It is
  also the one metric here that cannot be gamed by emitting less**, because it is
  a rate over junctions actually claimed.
* **distinct junction chains** — the number of unique intron chains across all
  emitted isoforms, against the number of isoforms. Two isoforms with the same
  chain are the same splice structure reported twice; the ratio is redundancy.
* **multi-exon fraction** — how many isoforms span at least one intron at all.
  A tool that collapsed everything to single-exon fragments would score well on
  error rate and aligned fraction and be useless.

# Caveats to state whenever these numbers are quoted

`NM` over a genome includes real biological difference — SNPs against the
reference assembly, unannotated sites minimap2 mis-places — so the absolute
error rate **overestimates** sequencing error. The inflation is identical across
sets, so comparisons hold where the absolute figure does not.

Sets are aligned independently. That is fair here, unlike the per-read
transcriptome case where each implementation picking its own best-matching
transcript would flatter all of them: a genome locus is essentially unambiguous,
so there is no menu of near-identical targets. What is still a risk is isoforms
dropping out of the comparison, so **unaligned isoforms are counted and
reported** rather than silently excluded — a set that emits garbage minimap2
cannot place must not thereby look clean.

Requires `minimap2` on PATH. The genome is indexed once per invocation by
minimap2 itself; a `.fai` is not needed, but the genome is read into memory to
check splice dinucleotides.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import statistics
from collections import Counter
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import annotation as annot  # noqa: E402  (path set above)

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# Canonical intron boundaries. A forward-strand intron begins GT and ends AG; the
# same intron read from the other strand is CT..AC. minimap2's `ts`/`ax splice`
# output does not always let us know the transcript strand, so both count.
CANONICAL = {("GT", "AG"), ("CT", "AC")}


def read_fasta(path):
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


def resolve_isoform_file(path):
    if os.path.isfile(path):
        return path
    cand = os.path.join(path, "transcriptome.fasta")
    if os.path.isfile(cand):
        return cand
    sys.exit(f"error: no transcriptome.fasta under {path}")


def run_minimap2(genome, fasta, threads, extra):
    """Primary alignments only, as (qname, ref, pos, cigar, nm) tuples."""
    cmd = [
        "minimap2",
        "-ax",
        "splice",
        "-uf",  # transcripts are forward-strand relative to the mRNA
        "-t",
        str(threads),
        "--secondary=no",
        genome,
        fasta,
    ] + list(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"minimap2 failed:\n{proc.stderr[-3000:]}")

    rows = []
    for line in proc.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x100 or flag & 0x800:  # secondary / supplementary
            continue
        qname, ref, pos, cigar = f[0], f[2], int(f[3]) - 1, f[5]
        if ref == "*" or cigar == "*":
            continue
        nm = None
        for tag in f[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag[5:])
                break
        rows.append((qname, ref, pos, cigar, nm, flag))
    return rows


def parse_cigar(cigar, ref_start, genome_seq):
    """Aligned query bases, query length, and the intron list.

    Introns are `N` operations, returned as (start, end) half-open genome
    coordinates so the flanking dinucleotides can be looked up.
    """
    aligned_q = 0
    q_len = 0
    ref_pos = ref_start
    introns = []
    for num, op in CIGAR_RE.findall(cigar):
        n = int(num)
        if op in "M=X":
            aligned_q += n
            q_len += n
            ref_pos += n
        elif op == "I":
            aligned_q += n
            q_len += n
        elif op == "D":
            ref_pos += n
        elif op == "N":
            introns.append((ref_pos, ref_pos + n))
            ref_pos += n
        elif op == "S":
            q_len += n
        elif op == "H":
            # Hard clip: those bases are not in this record's SEQ, but they are
            # part of the isoform, so they count toward its length.
            q_len += n
    return aligned_q, q_len, introns


def intron_is_canonical(genome_seq, start, end):
    if start < 2 or end > len(genome_seq):
        return False
    donor = genome_seq[start : start + 2]
    acceptor = genome_seq[end - 2 : end]
    return (donor, acceptor) in CANONICAL


def evaluate(genome, fasta, threads, extra, genome_seqs, ann=None):
    isoforms = read_fasta(fasta)
    rows = run_minimap2(genome, fasta, threads, extra)
    by_q = {}
    for qname, ref, pos, cigar, nm, flag in rows:
        by_q.setdefault(qname, (ref, pos, cigar, nm, flag))

    err_rates, aln_fracs = [], []
    chains = set()
    n_multi = 0
    n_introns = 0
    n_canonical = 0
    n_aligned = 0

    # SQANTI-style accounting, only when an annotation was supplied.
    cats = Counter()
    subs = Counter()
    fsm_tx = set()
    genes_touched = set()
    # For each NNC isoform: how many of its junctions are unannotated. An NNC
    # with exactly one bad junction and the rest right is a different kind of
    # failure from one that is wrong throughout, and the distinction points at
    # different code.
    nnc_novel = []

    for name in sorted(isoforms):
        if name not in by_q:
            continue
        ref, pos, cigar, nm, flag = by_q[name]
        gseq = genome_seqs.get(ref, "")
        aligned_q, q_len, introns = parse_cigar(cigar, pos, gseq)
        if aligned_q == 0:
            continue
        n_aligned += 1
        if nm is not None:
            err_rates.append(nm / aligned_q)
        # Prefer the true isoform length from the fasta over the CIGAR's view.
        true_len = len(isoforms[name]) or q_len
        aln_fracs.append(min(aligned_q / true_len, 1.0) if true_len else 0.0)
        if introns:
            n_multi += 1
            chains.add((ref, tuple(introns)))
            for s_, e_ in introns:
                n_introns += 1
                if intron_is_canonical(gseq, s_, e_):
                    n_canonical += 1

        if ann is not None:
            # `minimap2 -ax splice -uf` orients the alignment to the transcript,
            # so the SAM strand flag is the transcript's strand. That is what the
            # classifier needs: an isoform matching a gene on the other strand is
            # an antisense call, not a match.
            strand = "-" if (flag & 0x10) else "+"
            span_start = pos
            span_end = pos + sum(
                int(n) for n, op in CIGAR_RE.findall(cigar) if op in "M=XDN"
            )
            cat, sub, gene, tx = annot.classify(
                tuple(introns), ref, span_start, span_end, strand, ann
            )
            cats[cat] += 1
            if sub:
                subs[sub] += 1
            if cat == "FSM" and tx:
                fsm_tx.add(tx)
            if gene:
                genes_touched.update(gene.split(","))
            if cat == "NNC" and introns:
                known = set()
                for g_, _st, _gs, _ge in ann.genes_overlapping(ref, span_start, span_end):
                    known |= ann.gene_junctions[g_]
                nnc_novel.append((sum(1 for j in introns if j not in known), len(introns)))

    n_iso = len(isoforms)
    result = {
        "n_isoforms": n_iso,
        "n_aligned": n_aligned,
        "n_unaligned": n_iso - n_aligned,
        "err_median": statistics.median(err_rates) if err_rates else float("nan"),
        "err_mean": statistics.fmean(err_rates) if err_rates else float("nan"),
        "aln_median": statistics.median(aln_fracs) if aln_fracs else float("nan"),
        "n_multi_exon": n_multi,
        "n_introns": n_introns,
        "canonical_frac": (n_canonical / n_introns) if n_introns else float("nan"),
        "n_chains": len(chains),
    }
    if ann is not None:
        # Of the annotated transcripts belonging to genes this set produced
        # anything for, how many did it reconstruct exactly? A fair denominator,
        # unlike "all annotated transcripts" --- see bench/annotation.py.
        denom = set()
        for g in genes_touched:
            denom |= ann.gene_tx.get(g, set())
        result.update(
            {
                "cats": cats,
                "subs": subs,
                "n_fsm_tx": len(fsm_tx),
                "n_genes": len(genes_touched),
                "n_tx_in_genes": len(denom),
                "nnc_novel": nnc_novel,
            }
        )
    return result


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--genome", required=True)
    ap.add_argument(
        "--isoforms",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help="a set of reconstructed isoforms (fasta, or a folder holding "
        "transcriptome.fasta); repeat to compare sets",
    )
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument(
        "--annotation",
        help="GFF3 of reference gene models. Enables SQANTI-style structural "
        "classification (FSM/ISM/NIC/NNC/...). The FlyBase whole-genome GFF "
        "works but is mostly alignment evidence rather than annotation --- "
        "6.7 GB for r6.68 --- so pre-filtering it to transcript and exon "
        "records is worth doing once. See bench/corpus/README.md.",
    )
    ap.add_argument(
        "--minimap2-extra",
        nargs=argparse.REMAINDER,
        default=[],
        help="extra minimap2 flags, passed through verbatim",
    )
    args = ap.parse_args()

    if shutil.which("minimap2") is None:
        sys.exit("error: minimap2 not found on PATH")

    sets = []
    for spec in args.isoforms:
        if "=" not in spec:
            sys.exit(f"error: --isoforms wants NAME=PATH, got {spec!r}")
        name, path = spec.split("=", 1)
        sets.append((name, resolve_isoform_file(path)))

    print(f"genome: {args.genome}", file=sys.stderr)
    print("reading genome for splice-site checks...", file=sys.stderr)
    genome_seqs = read_fasta(args.genome)
    print(f"  {len(genome_seqs)} sequences", file=sys.stderr)

    ann = None
    if args.annotation:
        print(f"reading annotation {args.annotation}...", file=sys.stderr)
        ann = annot.Annotation.from_gff(args.annotation)
    print(file=sys.stderr)

    print(f"genome:  {args.genome}  ({len(genome_seqs)} sequences)")
    if ann is not None:
        print(f"annot:   {args.annotation}  "
              f"({len(ann.tx_chain)} transcripts, {len(ann.gene_tx)} genes)")
    print()

    results = []
    print(
        f"  {'set':<12} {'isof':>5} {'algn':>5} {'unalgn':>6} {'err.med':>8} "
        f"{'aln.frac':>8} {'multiexon':>9} {'introns':>7} {'canonical':>9} {'chains':>6}"
    )
    for name, path in sets:
        r = evaluate(args.genome, path, args.threads, args.minimap2_extra, genome_seqs, ann)
        results.append((name, r))
        print(
            f"  {name:<12} {r['n_isoforms']:>5} {r['n_aligned']:>5} "
            f"{r['n_unaligned']:>6} {r['err_median']:>8.4f} {r['aln_median']:>8.3f} "
            f"{r['n_multi_exon']:>9} {r['n_introns']:>7} "
            f"{r['canonical_frac']:>9.3f} {r['n_chains']:>6}"
        )

    print()
    print("err.med   = median NM / aligned isoform bases. Inflated by real")
    print("            biological difference from the assembly; comparable")
    print("            across sets, not quotable as an absolute error rate.")
    print("aln.frac  = median aligned isoform bases / isoform length. Reported")
    print("            because error rate alone is gamed by aligning less.")
    print("canonical = fraction of implied introns with GT..AG (or CT..AC).")
    print("            The annotation-free junction-quality signal, and the one")
    print("            number here that emitting fewer isoforms cannot improve.")
    print("chains    = distinct intron chains; chains < multiexon means the same")
    print("            splice structure was reported more than once.")

    if ann is None:
        return 0

    # --- SQANTI structural categories ---------------------------------------
    order = [
        "FSM", "ISM", "NIC", "NNC",
        "genic", "genic_intron", "antisense", "fusion", "intergenic",
    ]
    print()
    print("== SQANTI-style structural categories ==")
    print()
    hdr = "  {:<12}".format("set") + "".join(f"{c:>13}" for c in order)
    print(hdr)
    for name, r in results:
        cats = r["cats"]
        total = sum(cats.values()) or 1
        row = "  {:<12}".format(name)
        for c in order:
            n = cats.get(c, 0)
            row += f"{n:>6} {n/total:>5.0%}"
        print(row)

    print()
    print("  {:<12} {:>10} {:>12} {:>10} {:>18}".format(
        "set", "genes", "FSM tx", "tx in those", "FSM / tx in those"))
    for name, r in results:
        frac = r["n_fsm_tx"] / r["n_tx_in_genes"] if r["n_tx_in_genes"] else float("nan")
        print("  {:<12} {:>10} {:>12} {:>10} {:>17.1%}".format(
            name, r["n_genes"], r["n_fsm_tx"], r["n_tx_in_genes"], frac))

    print()
    for name, r in results:
        ir = r["subs"].get("intron_retention", 0)
        mono = r["subs"].get("mono_exon", 0)
        kj = r["subs"].get("known_junctions", 0)
        ks = r["subs"].get("known_sites", 0)
        print(f"  {name:<12} intron_retention={ir}  mono_exon={mono}  "
              f"NIC(known junctions)={kj}  NIC(known sites)={ks}")

    print()
    print("  How wrong is an NNC isoform? novel junctions / total junctions:")
    for name, r in results:
        nn = r["nnc_novel"]
        if not nn:
            print(f"  {name:<12} no NNC isoforms")
            continue
        dist = Counter(novel for novel, _tot in nn)
        only_one = dist.get(1, 0)
        print(f"  {name:<12} {len(nn)} NNC; " +
              ", ".join(f"{k} novel: {v}" for k, v in sorted(dist.items())) +
              f"   ({only_one}/{len(nn)} are wrong in exactly one junction)")

    print()
    print("FSM       = junction chain identical to an annotated transcript's.")
    print("            End positions are deliberately excluded: ONT reads are")
    print("            truncated and end position is not what isoform")
    print("            reconstruction is being judged on.")
    print("ISM       = chain is a contiguous sub-chain of an annotated one, i.e.")
    print("            a fragment. Skipping a middle junction is NIC, not ISM.")
    print("NIC       = every junction, or at least every splice site, is")
    print("            annotated but the combination is not. Real novel isoform")
    print("            or a chimera of known parts --- this metric cannot tell.")
    print("NNC       = at least one unannotated splice site. The strongest")
    print("            error signal here.")
    print()
    print("'FSM / tx in those' has a fair denominator on purpose: annotated")
    print("transcripts of genes this set produced anything for. A percentage")
    print("against all", len(ann.tx_chain), "annotated transcripts would measure")
    print("sequencing depth, not the tool.")
    print()
    print("An NNC wrong in exactly one junction is the signature of a single")
    print("mis-resolved bubble rather than a wholesale mis-assembly, so that")
    print("split is worth reading before deciding where to look in the graph code.")
    print()
    print("intron_retention cross-cuts ISM and NIC rather than being a category:")
    print("it fires when an exon swallows an annotated intron whole. For a")
    print("graph-based assembler that means a bubble that did not get popped, so")
    print("its rate is diagnostic rather than merely descriptive.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
