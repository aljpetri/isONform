# bench

Differential harness for the isONform Rust port. The method is carried over
from the isONcorrect port; `PORTING.md` at the repository root has the reasoning
and the measurements, this file is the operator's guide.

## Start here

```bash
bench/setup_reference_env.sh      # pinned conda env `isonform-ref`
bench/equivalence.sh cli          # port vs reference, argument handling
```

`bench/equivalence.sh` builds the Rust port itself (always, so a stale release
binary cannot pass), locates the reference env, and runs both against each other.

## The scripts

| script | what it does |
| --- | --- |
| `setup_reference_env.sh` | creates/verifies the pinned `isonform-ref` env and writes `env/resolved-<platform>.txt`, which is the audit trail for every golden |
| `check_seed_sensitivity.py` | runs the reference under N hash seeds and reports which output files vary. The regression check for the minimizer determinism fix |
| `equivalence.sh` | the differential harness: port against reference. Three layers, `cli`, `invariants` and `stages` |
| `evaluate.sh` | builds the three evaluation corpora and scores isONform on them. `corpora` / `run` / `score` / `compare` |
| `accuracy_isoforms.py` | isoform accuracy against a known transcriptome (SIRV): recall, precision, redundancy, per-isoform identity |
| `accuracy_isoforms_genome.py` | isoform accuracy against a genome (Drosophila) via `minimap2 -ax splice`: error rate, aligned fraction, canonical splice fraction, intron-chain redundancy. With `--annotation`, also SQANTI structural categories |
| `annotation.py` | GFF3 parsing and SQANTI-style classification (FSM/ISM/NIC/NNC/...). Run it directly to execute its self-tests — no genome or minimap2 needed |
| `dump_reference.py` | records a stage’s inputs and outputs from the live driver: `--stage intervals` writes `intervals_*.txt` (the front half of `main`), `--stage isoforms` writes `isoforms_*.txt`, `--stage graph` writes `graph_*.txt`, `--stage simplify` writes the graph before and after each `simplifyGraph` call, `--record-spoa` and `--record-parasail` write every `spoa` / `parasail_alignment` call to `spoa_cases.tsv` / `parasail_cases.tsv` |

## Two different questions

Keep them apart, because they need different tools and they fail differently.

**"Do the two implementations agree?"** — `equivalence.sh`'s `cli` and `stages`
layers. Byte-level, no tolerance, and how the port is proven faithful stage by
stage. Its `invariants` layer answers a third question that neither of the other
two can — *which side is wrong* — by checking properties the reference must
satisfy regardless of implementation.

**"Is the output any good?"** — `evaluate.sh` and the two `accuracy_*` scripts.
This is the question that matters whenever the *reference itself* changes: a
determinism fix, a deliberate divergence in the port, a performance tradeoff.
Equivalence cannot answer it, by construction — two implementations can agree
perfectly on a wrong answer.

isONform's accuracy question is also not isONcorrect's. isONcorrect emits one
corrected read per input read, so accuracy is per-read error against the read's
source transcript. isONform emits *isoforms*, and how many is an output rather
than an input, so the metric has to be recall, precision and per-isoform identity
**together**. Any one alone is trivially gamed: emit one confident isoform and
precision and identity look perfect while recall collapses; emit hundreds and the
reverse.

## `check_seed_sensitivity.py` — what it used to say, and what it says now

`main` used to select minimizers by `hash()` of the k-mer string, and CPython's
string hash is SipHash seeded randomly per process. So the shipped tool produced
a different transcriptome on every run. Measured on `corpus/sirv_small` over 8
seeds: all 9 of `isONform_parallel`'s output files varied and
`transcriptome.fasta` took 5 distinct values.

That is fixed — minimizer selection compares the k-mers lexicographically, as
isONcorrect's `get_kmer_minimizers_lex` does — so this script is now a
**regression check** and exits 0:

```bash
bench/check_seed_sensitivity.py --fastq-folder bench/corpus/sirv_small --seeds 24
```

The harness still exports `PYTHONHASHSEED=0` everywhere, and that is **not** just
belt and braces any more. `PORTING.md` finding 14 documents a second, independent
seed dependency — `find_connecting_edges` returns a `set` of string tuples and
`prepare_adding_edges` takes `conn_list[0]` — which survives the minimizer fix and
produces two distinct graphs across eight seeds on real data. It is rare (one
occurrence in 19 831 calls) and `corpus/sirv_small` never reaches it, which is
exactly why this script still exits 0. So a pass here means "minimizer selection
is deterministic", not "the tool is".

Note the sample size. In isONcorrect, six seeds made six records look like
porting bugs; at 24 they were all reference behaviour. Prefer more seeds than
feels necessary.

## The oracles

Six, all replaying recorded reference behaviour through the Rust port:

```bash
bench/dump_reference.py --fastq-folder bench/corpus/sirv_small --outdir /tmp/d \
    --record-spoa --record-parasail

ISONFORM_GRAPH_DUMPS=/tmp/d    cargo test --manifest-path rust/Cargo.toml \
    --release --test graph_oracle -- --nocapture
ISONFORM_SIMPLIFY_DUMPS=/tmp/d cargo test --manifest-path rust/Cargo.toml \
    --release --test simplify_oracle -- --nocapture
SPOA_CASES=/tmp/d/spoa_cases.tsv cargo test --manifest-path rust/Cargo.toml \
    --release --lib poa::oracle -- --nocapture
PARASAIL_CASES=/tmp/d/parasail_cases.tsv cargo test --manifest-path rust/Cargo.toml \
    --release --lib parasail::oracle -- --nocapture
```

All six **skip loudly** without their variable rather than passing silently,
and CI sets all six.

The **interval** oracle is the odd one out and the most load-bearing:

```bash
bench/dump_reference.py --fastq-folder bench/corpus/sirv_small --outdir /tmp/i --stage intervals
ISONFORM_INTERVAL_DUMPS=/tmp/i cargo test --manifest-path rust/Cargo.toml \
    --release --test intervals_oracle -- --nocapture
```

It covers the front half of `main` — minimizer selection, the anchor database, span support and
weighted interval scheduling — and it is the only oracle whose input is a **read** rather than a
recorded intermediate. That matters: the graph oracle replays the *reference's* intervals, so a wrong
minimizer would have produced a graph it happily agreed with. It compares the chosen intervals, the
`graph_id` assignment and each interval's full instance list in order.

Its dump is also produced differently. The other stages are recorded by replacing an importable
function with a wrapper; the front half's values (`w`, the `x` bounds, the batch's reads *including*
the ones it skips) exist only as locals inside `main.main`, which `runpy` executes. So
`--stage intervals` compiles `main`'s own source with a single recorder call injected just before
`generateGraphfromIntervals`. Same code, same driver, one extra line.

The next two replay a *stage*: recorded inputs in, diff the outputs. The graph
oracle passes on 72 recorded calls covering 47 963 nodes. The simplification
oracle rebuilds the entry graph from the dump — in **insertion order**, because
node and adjacency order decide `nx.topological_sort` and therefore which node
pairs are candidate bubbles — and asserts that order matches the reference's
before it will judge anything downstream.

The third replays a *dependency* rather than a stage. `crate::poa` wraps `spoars`
in place of the `spoa` binary, and simplification's poppability decisions run
through it, so it needs its own evidence. `--record-spoa` wraps
`IsoformGeneration.run_spoa` — the single function all four live call sites
resolve through — and writes `consensus<TAB>seq<TAB>seq…` lines, one per call,
deduplicated (spoa is deterministic, so a repeat exercises nothing) with the
call site kept as a `#` comment above each case so a mismatch localises without
regenerating anything.

The fourth does the same for `crate::parasail`, and needed doing for the same
reason even though this file previously said otherwise: the *score* is exact by
construction, but the **CIGAR** is a tie-break among equally-optimal paths, and
`parse_cigar_diversity` reads the CIGAR, not the score. Recording isONform's own
calls found 12 outright score errors and 136 CIGAR errors in 54 884
(`PORTING.md` finding 25). Both are now zero across **56 549** recorded calls and
both are gated exactly.

`PARASAIL_SWEEP` re-runs the tie-break sweep and `PARASAIL_HARD_OUT` writes just
the mismatching cases out. Two warnings on those, both learned the hard way. Do
**not** sweep on the failing subset and take the winner — it picked a setting that
fixed all 112 of them and made the full corpus worse (54 772 → 54 522). And a
sweep that fits nothing bounds the *parameter space*, not the problem: the last
CIGAR errors were not a `TieBreak` value at all but a rule outside the
parameterisation, and "no setting fits, therefore structural" sat in `PORTING.md`
for a commit before that turned out to be wrong.

That third oracle is what makes the second one unconditional. The simplification
oracle used to report-but-not-fail any disagreement that had called spoa, since
attributing one would have been claiming evidence that did not exist. It now
exists — 5 368 of 5 368 recorded isONform calls match the binary, see
`PORTING.md` — so every disagreement fails.

When a simplification case does disagree, the report gives three things before you reach for a
debugger: the port's **iteration and pop counts** (the reference prints its own, so the useful
comparison is the sequence rather than the total), node `reads` differences **naming** the reads and
marking synthetic ones (`*` = `original_support == false`, i.e. invented by
`additional_node_support`), and edge support compared as a **multiset** (`<read>x<surplus>`, because
`edge_supp` is a Python list and a read can legitimately appear twice — a set difference reports "no
difference" on exactly that case).

`ISONFORM_TRACE_POPS=1` goes one level deeper: one line per pop, with iteration, branch, bubble
endpoints and both path supports. `ISONFORM_TRACE_DECIDE=<comma-separated read ids>` goes deeper
still, dumping the inputs, both consensus sequences and the verdict for the one bubble whose path
carries exactly that read set --- the companion question, *why* did this one pop. That pair is what
found findings 24 and 25: in both, the two sides computed byte-identical consensus sequences and
disagreed anyway, which is what moved the search out of the stage and into the aligner. Diffing that against the reference's equivalent and finding the
first surplus or missing pop has localised two of the three bugs found here. It is worth knowing that
the *aggregate* counts have twice beaten a more precise-looking local signal — in finding 24 the
"only synthetic reads differ" observation was true and pointed at the wrong function, and what
corrected it was 94 pops against 92.

The **isoform** oracle has the most nuanced gate, and the nuance is the point:

```bash
bench/dump_reference.py --fastq-folder bench/corpus/sirv_small --outdir /tmp/o --stage isoforms
ISONFORM_ISOFORM_DUMPS=/tmp/o cargo test --manifest-path rust/Cargo.toml \
    --release --test isoforms_oracle -- --nocapture
```

It reports three outcomes, not one. A wrong **partition** — reads in the wrong
groups — fails. A **merging** disagreement on a case whose grouping and order both
matched fails. A difference that is *only* CPython **set-iteration order** — same
reads, same groups, different representative id or different order within a group —
is reported and counted but does **not** fail, because reproducing it means
modelling CPython's set internals and that is an open decision (`PORTING.md`
finding 28) rather than a settled one. On 114 real cases: 0 wrong partitions,
0 merging failures, 28 set-order differences.

That split is worth copying rather than collapsing. "28 of 114 disagree" and
"0 of 114 disagree on anything this port decides" are both true, and only the
second one tells you where to look.

One boundary worth holding onto when reading a failure: the spoa result licenses
"same inputs ⇒ same consensus", not "the port computes the same inputs". A
divergence in span extraction hands spoa different sequences and gets a
faithfully different consensus back. Still a port bug — hence the unconditional
gate — but it means `spoa_calls` in the report is now a *diagnostic* (did this
failure go through the consensus path at all?) rather than a verdict.

## Scoring Drosophila against the annotation

```bash
# once: the FlyBase whole-genome GFF is 6.7 GB and mostly alignment evidence
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=9 && ($3=="gene" || $3=="mRNA" || $3=="exon" \
    || $3=="ncRNA" || $3=="pseudogene" || $3=="pre_miRNA" || $3=="miRNA" || $3=="tRNA" \
    || $3=="snRNA" || $3=="snoRNA" || $3=="rRNA") {print $1,$2,$3,$4,$5,$7,$9}' \
    ~/data/annotations/dmel-all-r6.68.gff > ~/work/dmel-r6.68-transcripts.tsv

bench/accuracy_isoforms_genome.py --genome ~/data/genomes/fruitfly.fa \
    --annotation ~/work/dmel-r6.68-transcripts.tsv \
    --isoforms lex=$ISONFORM_WORK/eval/droso__a/out/transcriptome.fasta
```

`annotation.py` reads either the full 9-column GFF or that 7-column projection.
The filtered file is 18 MB and loads in seconds; the raw one works but is not
worth re-reading.

**Read `bench/annotation.py`'s docstring before interpreting the output.** Two
things in there are easy to get wrong and both were, at first: the recall
denominator (a percentage against all 34 989 annotated transcripts measures
sequencing depth, not the tool) and the fusion test (a span-based one reports 13%
fusions on real data, because FlyBase has 3 305 fully nested gene pairs).

## Evaluating a change to the reference

```bash
bench/evaluate.sh corpora                                    # once, ~10 min
git archive <commit> | tar -x -C /tmp/before                 # the "before" tree
bench/evaluate.sh compare /path/to/repo /tmp/before          # after vs before
```

`compare` runs the "before" implementation at **three** hash seeds rather than
one, deliberately: with a seeded hash, a single run of it says nothing about what
the tool does, and the spread across seeds *is* the finding.

Extracting the old tree with `git archive` rather than checking it out is the
point — the working tree is never disturbed, and both implementations are live at
once.

## The corpus

`corpus/sirv_small` is not checked in — `corpus/README.md` says how to rebuild
it, and why the reference needs *corrected* reads to do anything at all.

## Reading a mismatch

A mismatch is evidence, not a verdict: the first question is which side is
wrong. Answering it needs a check independent of both implementations —
brute force over small inputs, replay under many hash seeds, or an invariant the
reference should satisfy regardless. Graph code is good hunting ground for the
third: acyclicity, connectivity, and node and edge counts are all checkable
without trusting either implementation. `corpus/sirv_small` carries one such
invariant for free (its two clusters are the same reads, so they must yield the
same isoforms), and that invariant is what proved finding 1 without reference to
the port at all.
