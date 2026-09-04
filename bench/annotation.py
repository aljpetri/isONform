#!/usr/bin/env python3
"""SQANTI-style structural classification of reconstructed isoforms.

Used by `accuracy_isoforms_genome.py --annotation`. Standalone so the
classification logic can be unit-tested without minimap2 or a genome; run this
file directly to execute those tests.

# Why this and not just canonical-splice fraction

`accuracy_isoforms_genome.py` without an annotation reports how many implied
introns carry `GT..AG`. That is a real signal and it needs no annotation, but it
cannot distinguish a *correctly reconstructed known transcript* from a plausible
novel one, and it cannot measure recall at all. With an annotation, every
reconstructed isoform can be placed in one of the SQANTI structural categories,
which is the vocabulary the field already uses for exactly this question.

# The categories, and what each one means for isONform specifically

Comparison is on the **splice junction chain** — the ordered list of introns the
isoform implies — in genomic coordinates. 5'/3' end positions are deliberately
not part of FSM, because ONT reads are truncated and end position is not what
isoform reconstruction is being judged on.

| category | meaning | what it says about the tool |
| --- | --- | --- |
| `FSM` | Full Splice Match: junction chain identical to a reference transcript's | got it exactly right |
| `ISM` | Incomplete Splice Match: chain is a *contiguous* sub-chain of a reference transcript's | right as far as it goes; a fragment |
| `NIC` | Novel In Catalog: every junction, or at least every splice site, is annotated, but the combination is not | either a real novel isoform or a mis-assembled chimera of known parts |
| `NNC` | Novel Not in Catalog: at least one splice donor or acceptor is unannotated | most likely wrong; the strongest error signal here |
| `genic` | overlaps a gene but fits none of the above (typically mono-exonic) | uninformative output |
| `genic_intron` | lies entirely inside an annotated intron | usually an artifact |
| `antisense` | overlaps a gene only on the opposite strand | strand assignment or assembly error |
| `fusion` | junctions span two or more genes | mis-assembly across loci |
| `intergenic` | overlaps no annotated gene | novel locus, or garbage |

`intron_retention` is reported separately rather than as a category, because it
cross-cuts ISM and NIC: it fires when one of the isoform's exons fully spans an
annotated intron of the gene it matched. It is called out because for a
graph-based assembler it has a specific meaning — a bubble that did not get
popped — so its rate is diagnostic rather than merely descriptive.

# On recall, and the denominator problem

`FSM` counts give a genuine recall measure that the annotation-free version
cannot: how many annotated transcripts were reconstructed exactly. But the
denominator cannot be "all 35 000 annotated transcripts", because a few hundred
isoforms from 13 000 reads could not possibly cover them and the resulting
percentage would measure sequencing depth, not the tool.

So two numbers are reported and neither is called "recall" on its own:

* **absolute** count of distinct reference transcripts with at least one FSM;
* the same restricted to **genes the isoform set actually touched**, which asks
  "of the loci we produced anything for, how much of their annotated structure
  did we recover" — a fair question with a fair denominator.

# Coordinates

Everything internal is **0-based half-open**, matching what a CIGAR walk
produces. GFF is 1-based inclusive and is converted on read. An intron between
exons `[s1,e1)` and `[s2,e2)` is `[e1, s2)`.
"""

from __future__ import annotations

import bisect
import gzip
import re
import sys
from collections import defaultdict

# FlyBase and Ensembl both use these; the ones that carry a Parent pointing at a
# gene and have exon children. `mRNA` covers protein-coding; the rest are the
# non-coding classes FlyBase emits as first-class transcript types.
TRANSCRIPT_TYPES = {
    "mRNA",
    "ncRNA",
    "pre_miRNA",
    "miRNA",
    "tRNA",
    "snRNA",
    "snoRNA",
    "rRNA",
    "pseudogene",
    "transcript",
    "lnc_RNA",
}

ATTR_RE = re.compile(r"([^=;]+)=([^;]*)")


def _attrs(field):
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(field)}


class Annotation:
    """Transcript structures and per-gene splice-site catalogues."""

    def __init__(self):
        self.tx_exons = defaultdict(list)          # tx -> [(start, end)] 0-based half-open
        self.tx_gene = {}                          # tx -> gene
        self.tx_chrom = {}
        self.tx_strand = {}
        self.gene_span = {}                        # gene -> (chrom, start, end, strand)
        self.tx_chain = {}                         # tx -> ((istart, iend), ...)
        self.chain_index = defaultdict(list)       # chain -> [tx]
        self.gene_junctions = defaultdict(set)
        self.gene_sites = defaultdict(set)         # every donor and acceptor position
        self.gene_introns = defaultdict(set)
        self.gene_tx = defaultdict(set)
        self._chrom_genes = defaultdict(list)      # chrom -> [(start, end, gene)] sorted
        self._chrom_starts = {}

    # -- loading -----------------------------------------------------------

    @classmethod
    def from_gff(cls, path, verbose=True):
        """Parse a GFF3.

        Tolerates the FlyBase whole-genome GFF, which is mostly *not*
        annotation: r6.68 is 6.7 GB of which the gene models are a rounding
        error, the bulk being `match`/`match_part` alignment evidence. Only
        transcript and exon records are kept, so a pre-filtered file (same nine
        columns, or the 7-column projection this module also accepts) loads far
        faster and is worth making once.
        """
        self = cls()
        opener = gzip.open if str(path).endswith(".gz") else open
        n_lines = 0
        with opener(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                # Accept the full 9-column GFF or a 7-column projection
                # (chrom, source, type, start, end, strand, attrs).
                if len(f) == 9:
                    chrom, _src, ftype, start, end, _sc, strand, _ph, attr = f
                elif len(f) == 7:
                    chrom, _src, ftype, start, end, strand, attr = f
                else:
                    continue
                if ftype not in TRANSCRIPT_TYPES and ftype != "exon" and ftype != "gene":
                    continue
                n_lines += 1
                if verbose and n_lines % 200_000 == 0:
                    print(f"  ...{n_lines} annotation records", file=sys.stderr)

                a = _attrs(attr)
                s, e = int(start) - 1, int(end)   # -> 0-based half-open

                if ftype == "gene":
                    gid = a.get("ID")
                    if gid:
                        self.gene_span[gid] = (chrom, s, e, strand)
                elif ftype in TRANSCRIPT_TYPES:
                    tid = a.get("ID")
                    if not tid:
                        continue
                    parent = a.get("Parent", "").split(",")[0]
                    self.tx_gene[tid] = parent or tid
                    self.tx_chrom[tid] = chrom
                    self.tx_strand[tid] = strand
                else:  # exon
                    # An exon may list several transcript parents.
                    for parent in a.get("Parent", "").split(","):
                        if parent:
                            self.tx_exons[parent].append((s, e))

        self._finalise(verbose)
        return self

    def _finalise(self, verbose=True):
        dropped = 0
        for tid, exons in self.tx_exons.items():
            if tid not in self.tx_gene:
                dropped += 1
                continue
            exons.sort()
            chain = tuple(
                (exons[i][1], exons[i + 1][0])
                for i in range(len(exons) - 1)
                if exons[i + 1][0] > exons[i][1]      # skip zero-length/overlapping
            )
            self.tx_chain[tid] = chain
            gene = self.tx_gene[tid]
            self.gene_tx[gene].add(tid)
            if chain:
                self.chain_index[chain].append(tid)
                for js, je in chain:
                    self.gene_junctions[gene].add((js, je))
                    self.gene_introns[gene].add((js, je))
                    self.gene_sites[gene].add(js)
                    self.gene_sites[gene].add(je)
            # A gene record is not guaranteed to precede its transcripts, and
            # some annotations omit it; derive the span from exons if needed.
            if gene not in self.gene_span:
                self.gene_span[gene] = (
                    self.tx_chrom[tid],
                    exons[0][0],
                    exons[-1][1],
                    self.tx_strand[tid],
                )
            else:
                chrom, gs, ge, strand = self.gene_span[gene]
                self.gene_span[gene] = (chrom, min(gs, exons[0][0]), max(ge, exons[-1][1]), strand)

        for gene, (chrom, s, e, strand) in self.gene_span.items():
            if gene in self.gene_tx:          # genes with no transcripts are not useful here
                self._chrom_genes[chrom].append((s, e, gene, strand))
        for chrom in self._chrom_genes:
            self._chrom_genes[chrom].sort()
            self._chrom_starts[chrom] = [g[0] for g in self._chrom_genes[chrom]]

        if verbose:
            n_multi = sum(1 for c in self.tx_chain.values() if c)
            print(
                f"  {len(self.tx_chain)} transcripts ({n_multi} multi-exon) "
                f"in {len(self.gene_tx)} genes"
                + (f"; {dropped} exon groups had no transcript record" if dropped else ""),
                file=sys.stderr,
            )

    # -- queries -----------------------------------------------------------

    def genes_overlapping(self, chrom, start, end):
        """Genes whose span intersects [start, end). Any strand."""
        starts = self._chrom_starts.get(chrom)
        if not starts:
            return []
        # Genes are shorter than a megabase in practice; scanning back a bounded
        # window is simpler than an interval tree and fast enough at this scale.
        i = bisect.bisect_right(starts, end)
        out = []
        genes = self._chrom_genes[chrom]
        j = i - 1
        while j >= 0:
            gs, ge, gene, strand = genes[j]
            if ge <= start:
                # Spans vary, so do not break immediately; bound the scan.
                if start - gs > 1_000_000:
                    break
            else:
                out.append((gene, strand, gs, ge))
            j -= 1
        return out


def _is_contiguous_subchain(query, ref):
    """Is `query` a contiguous run of junctions inside `ref`, and shorter?

    ISM means the isoform is a *fragment*: it must agree with the reference over
    a stretch and simply stop, not skip a junction in the middle. A
    non-contiguous match is a novel combination, i.e. NIC — which is why this is
    a sublist test and not a subset test.
    """
    n, m = len(query), len(ref)
    if n == 0 or n >= m:
        return False
    for i in range(m - n + 1):
        if ref[i : i + n] == query:
            return True
    return False


def classify(chain, chrom, start, end, strand, ann):
    """SQANTI-style structural category for one aligned isoform.

    `chain` is the ordered tuple of introns from the CIGAR, all coordinates
    0-based half-open. Returns `(category, subcategory, gene, matched_tx)` where
    subcategory is `"intron_retention"` or `""`.
    """
    overlaps = ann.genes_overlapping(chrom, start, end)
    if not overlaps:
        return "intergenic", "", None, None

    same = [o for o in overlaps if o[1] == strand]
    if not same:
        return "antisense", "", overlaps[0][0], None

    # Fusion: the isoform's junctions come from two or more distinct genes.
    #
    # Two conditions, and the first version of this had neither, which made it
    # fire on 13% of a real Drosophila set --- FlyBase has 1 683 overlapping
    # same-strand gene pairs and 3 305 fully nested ones, so *any* test based on
    # gene-span overlap reports a fusion whenever an intron happens to span a
    # nested gene. Measured; that is where the bogus 13% came from, and it was
    # swallowing the entire NIC category.
    #
    #  1. a gene only counts if at least one of the isoform's junctions is
    #     *annotated in that gene*, not merely inside its span;
    #  2. the contributing genes must occupy disjoint spans, so a nested or
    #     overlapping gene pair cannot look like a fusion.
    if chain:
        contributors = {}
        for gene, _gstrand, gs, ge in same:
            known = ann.gene_junctions[gene]
            if any(j in known for j in chain):
                contributors[gene] = (gs, ge)
        if len(contributors) > 1:
            spans = sorted(contributors.values())
            disjoint = all(spans[i][1] <= spans[i + 1][0] for i in range(len(spans) - 1))
            if disjoint:
                return "fusion", "", ",".join(sorted(contributors)), None

    # Score against each candidate gene and keep the best category. Genes are
    # ranked by how much of the isoform they explain, so a gene nested inside
    # another's span cannot win by accident.
    def gene_key(o):
        gs, ge = o[2], o[3]
        return -(min(end, ge) - max(start, gs))

    best = None
    for gene, _gstrand, gs, ge in sorted(same, key=gene_key):
        if not chain:
            # Mono-exonic. FSM only against a mono-exonic reference transcript
            # whose single exon this substantially covers.
            for tid in ann.gene_tx[gene]:
                if ann.tx_chain.get(tid):
                    continue
                exons = ann.tx_exons.get(tid) or []
                if not exons:
                    continue
                ts, te = exons[0][0], exons[-1][1]
                ov = min(end, te) - max(start, ts)
                if ov > 0 and ov >= 0.8 * (te - ts) and ov >= 0.8 * (end - start):
                    return "FSM", "mono_exon", gene, tid
            # Inside an intron of every transcript that spans it?
            if any(js <= start and end <= je for js, je in ann.gene_introns[gene]):
                cand = ("genic_intron", "", gene, None)
            else:
                cand = ("genic", "mono_exon", gene, None)
            best = best or cand
            continue

        # Multi-exonic.
        retention = "intron_retention" if _has_retained_intron(chain, start, end, gene, ann) else ""

        hits = ann.chain_index.get(chain)
        if hits:
            same_gene = [t for t in hits if ann.tx_gene.get(t) == gene]
            if same_gene:
                return "FSM", retention, gene, sorted(same_gene)[0]

        for tid in sorted(ann.gene_tx[gene]):
            if _is_contiguous_subchain(chain, ann.tx_chain.get(tid, ())):
                return "ISM", retention, gene, tid

        known_j = ann.gene_junctions[gene]
        known_s = ann.gene_sites[gene]
        if all(j in known_j for j in chain):
            cand = ("NIC", retention or "known_junctions", gene, None)
        elif all(js in known_s and je in known_s for js, je in chain):
            cand = ("NIC", retention or "known_sites", gene, None)
        else:
            cand = ("NNC", retention, gene, None)
        # NIC beats NNC if another candidate gene explains it better.
        if best is None or (best[0] == "NNC" and cand[0] == "NIC"):
            best = cand

    return best or ("genic", "", same[0][0], None)


def _has_retained_intron(chain, start, end, gene, ann):
    """Does one of this isoform's exons swallow an annotated intron whole?

    For a graph-based assembler this is a specific failure — a bubble that was
    not popped — so it is worth detecting rather than folding into ISM.
    """
    # Exons are the gaps between consecutive introns, plus the two flanks.
    bounds = [start]
    for js, je in chain:
        bounds.append(js)
        bounds.append(je)
    bounds.append(end)
    exons = [(bounds[i], bounds[i + 1]) for i in range(0, len(bounds) - 1, 2)]
    for es, ee in exons:
        for js, je in ann.gene_introns[gene]:
            if es < js and je < ee:
                return True
    return False


# --------------------------------------------------------------------------
# Self-tests. `python bench/annotation.py` runs them; no genome or minimap2
# needed, which is the point of keeping this module standalone.
# --------------------------------------------------------------------------


def _toy():
    """A two-gene toy annotation.

    geneA on '+' at 1000-2000 with two transcripts:
        txA1 exons 1000-1100, 1300-1400, 1700-2000  -> introns (1100,1300) (1400,1700)
        txA2 exons 1000-1100, 1700-2000             -> introns (1100,1700)
    geneB on '-' at 5000-6000, one transcript:
        txB1 exons 5000-5100, 5500-6000             -> introns (5100,5500)
    geneC on '+' at 8000-8500, mono-exonic txC1.
    geneD on '+' at 6500-7500, one transcript:
        txD1 exons 6500-6600, 7000-7500             -> introns (6600,7000)

    geneD exists so fusion can be tested on the *same* strand, which is the only
    way this implementation reports it: a query overlapping a neighbour on the
    opposite strand is antisense or NNC, not a fusion, and conflating the two
    would make the fusion count meaningless.
    """
    ann = Annotation()
    def add(tid, gene, chrom, strand, exons):
        ann.tx_gene[tid] = gene
        ann.tx_chrom[tid] = chrom
        ann.tx_strand[tid] = strand
        ann.tx_exons[tid] = list(exons)
    add("txA1", "geneA", "2L", "+", [(1000, 1100), (1300, 1400), (1700, 2000)])
    add("txA2", "geneA", "2L", "+", [(1000, 1100), (1700, 2000)])
    add("txB1", "geneB", "2L", "-", [(5000, 5100), (5500, 6000)])
    add("txC1", "geneC", "2L", "+", [(8000, 8500)])
    add("txD1", "geneD", "2L", "+", [(6500, 6600), (7000, 7500)])
    ann.gene_span["geneA"] = ("2L", 1000, 2000, "+")
    ann.gene_span["geneB"] = ("2L", 5000, 6000, "-")
    ann.gene_span["geneC"] = ("2L", 8000, 8500, "+")
    ann.gene_span["geneD"] = ("2L", 6500, 7500, "+")
    ann._finalise(verbose=False)
    return ann


def _toy_with_nested_gene():
    """`_toy` plus geneE nested inside geneA's first intron, same strand.

    FlyBase is full of these --- 3 305 fully nested gene pairs in r6.68 --- so a
    classifier that cannot handle them is not usable on real annotation.
    """
    ann = _toy()
    ann.tx_gene["txE1"] = "geneE"
    ann.tx_chrom["txE1"] = "2L"
    ann.tx_strand["txE1"] = "+"
    ann.tx_exons["txE1"] = [(1150, 1180), (1220, 1260)]
    ann.gene_span["geneE"] = ("2L", 1150, 1260, "+")
    # Rebuild the derived indices from scratch; _finalise is idempotent over
    # tx_exons but the per-gene sets are additive, so start clean.
    fresh = Annotation()
    fresh.tx_exons = ann.tx_exons
    fresh.tx_gene = ann.tx_gene
    fresh.tx_chrom = ann.tx_chrom
    fresh.tx_strand = ann.tx_strand
    fresh.gene_span = {
        g: v for g, v in ann.gene_span.items()
    }
    fresh._finalise(verbose=False)
    return fresh


def _test():
    ann = _toy()
    fails = []

    def check(name, got, want):
        if got != want:
            fails.append(f"{name}: got {got!r}, want {want!r}")

    # FSM: exact chain of txA1.
    cat, sub, gene, tx = classify(((1100, 1300), (1400, 1700)), "2L", 1000, 2000, "+", ann)
    check("FSM txA1", (cat, gene, tx), ("FSM", "geneA", "txA1"))

    # FSM even with different ends — end position is deliberately not part of it.
    cat, _, _, tx = classify(((1100, 1300), (1400, 1700)), "2L", 1050, 1900, "+", ann)
    check("FSM truncated ends", (cat, tx), ("FSM", "txA1"))

    # FSM against the two-exon transcript.
    cat, _, _, tx = classify(((1100, 1700),), "2L", 1000, 2000, "+", ann)
    check("FSM txA2", (cat, tx), ("FSM", "txA2"))

    # ISM: a contiguous prefix of txA1's chain.
    cat, _, _, tx = classify(((1100, 1300),), "2L", 1000, 1400, "+", ann)
    check("ISM prefix", (cat, tx), ("ISM", "txA1"))

    # ISM: contiguous suffix.
    cat, _, _, tx = classify(((1400, 1700),), "2L", 1300, 2000, "+", ann)
    check("ISM suffix", (cat, tx), ("ISM", "txA1"))

    # NIC: both junctions known, combination is not. Pair txA2's long intron
    # with txA1's second — no reference transcript has that chain.
    cat, sub, _, _ = classify(((1100, 1700), (1400, 1700)), "2L", 1000, 2000, "+", ann)
    check("NIC known junctions", cat, "NIC")

    # NIC from known sites, novel junction: donor 1100 with acceptor 1700 exists,
    # but donor 1400 paired to acceptor 1300 does not, and both sites are known.
    cat, sub, _, _ = classify(((1100, 1400),), "2L", 1000, 2000, "+", ann)
    check("NIC known sites", (cat, sub), ("NIC", "known_sites"))

    # NNC: an unannotated donor.
    cat, _, _, _ = classify(((1150, 1300),), "2L", 1000, 1400, "+", ann)
    check("NNC novel site", cat, "NNC")

    # Antisense: geneB is on '-', query on '+'.
    cat, _, _, _ = classify(((5100, 5500),), "2L", 5000, 6000, "+", ann)
    check("antisense", cat, "antisense")

    # Fusion: junctions reaching into geneA and geneD, both on '+'.
    cat, _, gene, _ = classify(((1100, 1300), (6600, 7000)), "2L", 1000, 7500, "+", ann)
    check("fusion", (cat, gene), ("fusion", "geneA,geneD"))

    # A query spanning a gene on the opposite strand is NOT a fusion. geneB is
    # on '-', so only geneA is a same-strand candidate and its second junction is
    # unannotated there.
    cat, _, _, _ = classify(((1100, 1300), (5100, 5500)), "2L", 1000, 6000, "+", ann)
    check("cross-strand span is not fusion", cat, "NNC")

    # Regression: a NESTED same-strand gene must not produce a fusion call.
    # geneE sits inside geneA's intron on the same strand. An isoform that is a
    # perfect FSM for txA1 spans it, and the first version of the fusion test
    # called that a fusion --- which on real Drosophila data fired on 13% of
    # isoforms and emptied the NIC category.
    nested = _toy_with_nested_gene()
    cat, _, _, tx = classify(((1100, 1300), (1400, 1700)), "2L", 1000, 2000, "+", nested)
    check("nested gene does not create a fusion", (cat, tx), ("FSM", "txA1"))

    # And a genuine fusion still is one: junctions annotated in two genes whose
    # spans are disjoint.
    cat, _, gene, _ = classify(((1100, 1300), (6600, 7000)), "2L", 1000, 7500, "+", nested)
    check("real fusion survives the fix", (cat, gene), ("fusion", "geneA,geneD"))

    # Intergenic.
    cat, _, _, _ = classify(((3100, 3300),), "2L", 3000, 3400, "+", ann)
    check("intergenic", cat, "intergenic")

    # Mono-exonic FSM against the mono-exonic geneC transcript.
    cat, sub, _, tx = classify((), "2L", 8000, 8500, "+", ann)
    check("mono FSM", (cat, tx), ("FSM", "txC1"))

    # Mono-exonic inside an intron of geneA.
    cat, _, _, _ = classify((), "2L", 1450, 1650, "+", ann)
    check("genic_intron", cat, "genic_intron")

    # Intron retention: chain of txA2, but the single exon spans txA1's
    # (1100,1300) intron... constructed as an isoform that keeps intron
    # (1400,1700) unspliced while splicing (1100,1300).
    cat, sub, _, _ = classify(((1100, 1300),), "2L", 1000, 2000, "+", ann)
    check("intron retention flagged", sub, "intron_retention")

    # A contiguous-subchain test that must NOT be ISM: skipping a middle
    # junction is a novel combination, not a fragment.
    assert not _is_contiguous_subchain(((1100, 1300),), ((1100, 1300),)), "equal chains are FSM, not ISM"
    assert _is_contiguous_subchain(((1400, 1700),), ((1100, 1300), (1400, 1700)))
    assert not _is_contiguous_subchain(
        ((1100, 1300), (9000, 9100)), ((1100, 1300), (1400, 1700), (9000, 9100))
    ), "non-contiguous match must not be ISM"

    if fails:
        print("FAILED:", file=sys.stderr)
        for f in fails:
            print("  -", f, file=sys.stderr)
        return 1
    print("annotation.py: all self-tests pass")
    return 0


if __name__ == "__main__":
    sys.exit(_test())
