# Corpora

Four inputs. Everything here is derived, and this file is how each one is
rebuilt. `bench/evaluate.sh corpora` builds the three large ones.

`sirv_small` is the one exception that **is** committed — 88 KB, and CI's
determinism gate needs it. A gate that silently skips because its input is
absent is worse than no gate, so it is tracked rather than rebuilt in CI.

The other three land in `$ISONFORM_WORK` (default `~/work/isonform-corpus`),
outside the repository: they are gigabytes of derived fastq and committing any of
that would be a mistake.

## Why there is more than one

PORTING.md method point 4 exists because simulated and spike-in data gave the
wrong answer *twice* during the isONcorrect port — endorsing an alignment band
that real data broke, and reversing the sign of an accuracy effect. So the
corpora are chosen to fail in different directions, and a change is only
believed when it moves all of them the same way.

| corpus | reads | clusters | truth | what it is for |
| --- | --- | --- | --- | --- |
| `sirv_small` | 100 ×2 | 2 (identical) | none | smoke test, oracles, CI |
| `sirv_sim_gene` | 10 000 simulated | 7, 713–2242 reads | exact, per read | recall with a real denominator |
| `sirv_real` | 9 968 real ONT | 26, 12–4595 reads | transcriptome | real error profile, known answer |
| `droso` | 13 164 real ONT | 561, 5–488 reads | genome only | transcriptome scale, real biology |

Every one of them needs isONcorrect run over it first. That is not optional
preprocessing: on raw isONclust output isONform's interval-abundance filter skips
essentially every read and the graph comes out with 2 nodes and 0 edges. The
corrected reads are the input contract.

---

## `sirv_small` — smoke test

Two clusters, 100 reads each, ~190 bp: isONcorrect's own
`test_data/isoncorrect/{0,1}.fastq` pushed through isONcorrect.

```bash
for c in 0 1; do
  isONcorrect --fastq $ISONCORRECT/test_data/isoncorrect/$c.fastq --outfolder /tmp/c$c
  cp /tmp/c$c/corrected_reads.fastq bench/corpus/sirv_small/$c.fastq
done
```

### It is one cluster twice, and that is useful

`test_data/isoncorrect/0.fastq`, `1.fastq` and `test_data/100reads.fq` are
**byte-identical** (sha1 `ef0f78f8…`). Kept deliberately, because it buys an
invariant that needs neither implementation to be trusted: **clusters 0 and 1
must yield the same isoforms modulo the cluster id.**

That invariant is what proved the seed-dependence defect (PORTING.md finding 1)
without reference to the port at all — before the fix, one run emitted `0_0_71`
for cluster 0 and `1_1_2` for cluster 1 from identical input, because each child
process drew its own hash seed. After the fix the two differ only in the id.

It says nothing about cross-cluster behaviour, which is what the other three are
for. **Small and structurally thin: a smoke test, not evidence.**

---

## `sirv_sim_gene` — simulated, grouped by gene

10 000 simulated SIRV reads at 7% error, grouped into 7 clusters by **gene**, so
isoforms of one gene share a cluster and reads differ by whole exons. That
grouping is the point: it is the shape isONform exists to handle.

```bash
bench/evaluate.sh corpora        # or, just this corpus:
make_sirv_corpus.py --fastq reads_10k_err7%.fastq --outdir $W/sirv_sim_gene/clusters --group gene
run_isoncorrect --fastq_folder $W/sirv_sim_gene/clusters --outfolder $W/sirv_sim_gene/corrected --t 8 --split_wrt_batches
```

Cluster sizes 713 / 878 / 953 / 1551 / 1717 / 1946 / 2242 — the largest cross
`--max_seqs 1000`, so batching runs for real, and all 7 are far above
`--exact_instance_limit`.

**Truth is exact and free.** Every read header names its source transcript
(`@read_1_from_SIRV612`), so `accuracy_isoforms.py --expressed-from` can restrict
the recall denominator to transcripts actually present. Without that, recall is
measuring the corpus rather than the tool.

Reproducible from the fastq alone — no clustering run in the loop — which is what
a verification corpus wants.

**Weakness, and it is the important one:** errors are uniform and i.i.d., there
are no chimeras, and every read is full length. That is exactly the profile that
misled the isONcorrect port twice. Never conclude from this corpus alone.

---

## `sirv_real` — real ONT reads, known answer

9 968 real ONT SIRV spike-in reads through isONclust: 26 clusters, largest 4 595
reads, 20 of them at 50 reads or more.

```bash
bench/evaluate.sh corpora
```

This is the primary accuracy corpus, because SIRV is the one setting where a real
error profile and a trustworthy ground-truth transcriptome coexist. Real reads
bring truncation, chimeras and the ONT indel bias that the simulator does not
model.

Scored with `accuracy_isoforms.py` and no `--expressed-from`: real reads carry no
source transcript in their headers, so every one of the 68 reference transcripts
counts and **recall is a lower bound**. The script says so on every run.

**Weakness:** 68 transcripts from one synthetic spike-in set. It is a correctness
test, not a transcriptome-scale one.

---

## `droso` — real ONT Drosophila, genome truth

13 164 real ONT Drosophila cDNA reads (pychopper'd full-length output) through
isONclust: 561 clusters, 71 of them at 50 reads or more.

```bash
bench/evaluate.sh corpora
```

Transcriptome scale and real biology, and the only corpus here whose clusters number in the
hundreds. There is no trustworthy per-isoform truth, so it is scored against the **genome** with
`accuracy_isoforms_genome.py`, in two modes.

**Without an annotation** it reports error rate (`NM` / aligned bases) **and** aligned fraction —
never one alone, since error rate by itself is gamed by aligning less — plus **canonical splice
fraction**, the share of implied introns with `GT..AG`. That last is the annotation-free
junction-quality signal and the only metric there that emitting fewer isoforms cannot improve.

**With `--annotation`** it additionally classifies every isoform into the SQANTI structural
categories — FSM / ISM / NIC / NNC / genic / antisense / fusion / intergenic — which is a much
sharper instrument: it distinguishes a correctly reconstructed known transcript from a merely
plausible novel one, and it gives a real recall number. See `bench/annotation.py`.

```bash
# once: FlyBase's whole-genome GFF is 6.7 GB and mostly alignment evidence
# (19 M match_part records against 85 601 exons), so filter it first.
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=9 && ($3=="gene" || $3=="mRNA" || $3=="exon" \
    || $3=="ncRNA" || $3=="pseudogene" || $3=="pre_miRNA" || $3=="miRNA" || $3=="tRNA" \
    || $3=="snRNA" || $3=="snoRNA" || $3=="rRNA") {print $1,$2,$3,$4,$5,$7,$9}' \
    ~/data/annotations/dmel-all-r6.68.gff > ~/work/dmel-r6.68-transcripts.tsv
```

139 209 records, 18 MB. `annotation.py` reads either form.

`NM` against a genome includes real biological difference — SNPs against the assembly, sites
minimap2 mis-places — so the absolute error rate is an overestimate. The inflation is identical
across compared sets, so comparisons hold where the absolute number does not. Do not quote it as
"the error rate".

**Two traps, both of which were sprung before being fixed.** A recall percentage against all 34 989
annotated transcripts measures sequencing depth rather than the tool, so the reported denominator is
the transcripts of genes the isoform set actually touched. And a fusion test based on gene-span
overlap reports 13% fusions here, because FlyBase r6.68 has 1 683 overlapping same-strand gene pairs
and 3 305 fully nested ones. Both are covered in `bench/annotation.py`'s docstring and its
self-tests.

---

## Inputs these are derived from

Not in the repository. `bench/evaluate.sh` expects them at these paths and every
one is overridable in the environment.

| variable | default |
| --- | --- |
| `SIRV_SIM_FASTQ` | `~/data/lrRNA-seq/sirv/reads_10k_err7%.fastq` |
| `SIRV_REAL_FASTQ` | `~/data/lrRNA-seq/sirv/SIRV_real_10k.fastq` |
| `DROSO_FASTQ` | `~/data/lrRNA-seq/droso/full_length_output_first_20k.fq` |
| `SIRV_TRANSCRIPTOME` | `~/source/isONcorrect/test_data/sirv_transcriptome.fasta` (68 transcripts) |
| `DROSO_GENOME` | `~/data/genomes/fruitfly.fa` |
| `DROSO_ANNOTATION` | `~/work/dmel-r6.68-transcripts.tsv`, filtered from `~/data/annotations/dmel-all-r6.68.gff` |

Larger versions of the read sets exist beside them (`SIRV_real_100k.fastq`,
`SIRV_real_full.fastq`, `full_length_output_first_100k.fq`,
`full_length_output.fq`) for when a run is worth the wall clock. The 10k/20k
subsets are the routine ones.
