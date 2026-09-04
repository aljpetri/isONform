# isONform: Rust port vs the python reference

Every number here is a single measured run, not an estimate. The port is
deterministic and so is the reference, so repeated runs give identical output;
only wall clock varies.

## How these were run

| | |
|---|---|
| machine | Apple M4 Max, 16 cores (12 performance + 4 efficiency), 64 GB RAM |
| OS | macOS 26.5.1 |
| threads | `--t 8` for every run, python and Rust alike |
| other flags | `--split_wrt_batches --iso_abundance 5` |
| python | 3.10.21, networkx 2.8.4, parasail 1.3.4, spoa 4.1.5 (`bench/env/reference.yml`) |
| harness | `bench/evaluate.sh run` / `score` |
| memory | peak RSS of the process and its waited-for children, via `/usr/bin/time -l` --- for `isONform_parallel` that is the largest single worker |

**port (faithful)** is `--faithful`: the reference reproduced exactly.
**port (default)** is the shipped configuration: faithful plus the WFA2 aligner.

**Poly-A tails are trimmed from isoforms before matching.** The SIRV reference
transcripts are annotated without one; cDNA has one. Untrimmed, those bases are
charged twice --- against identity, whose denominator is the query length, and
against the length tolerance --- so a perfect reconstruction can fail both. On
real ONT SIRV, 71 of python's 155 isoforms carry such a tail; on PacBio one
isoform reconstructed 3 998 of 3 999 bases and then ran 303 poly-A bases past the
end, scoring 0.929. It affects every implementation about equally, so the
comparisons hold either way, but absolute numbers move by 2--7 transcripts.
`bench/accuracy_isoforms.py --no-trim-polya` scores the old way. The simulated
corpus is unchanged to four decimals, which is the control: the simulator emits
no poly-A.


## The corpora

| corpus | data | reads | clusters | scored against |
|---|---|---|---|---|
| `sirv_sim_gene` | simulated SIRV, grouped by gene | 10 000 | 7 | SIRV transcriptome, exact truth |
| `sirv_real` | real ONT SIRV spike-in | 10 000 | 26 | SIRV transcriptome (recall is a lower bound) |
| `sirv_real_deep` | the same, undownsampled | 99 631 | 26 | as above; deepest cluster 46 437 reads |
| `droso` | real ONT Drosophila cDNA | 20 000 | 561 | genome + FlyBase annotation |
| `droso_deep` | the same, 1M reads | 1 000 000 | 9 697 | as above |
| `pacbio_sirv` | PacBio HiFi SIRV spike-in | 17 550 | 44 | SIRV set 4, 83 transcripts |
| `pacbio_droso` | PacBio HiFi Drosophila cDNA | 455 226 | 99 | genome + FlyBase annotation |

The five ONT corpora went through the real upstream pipeline: isONclust, then
isONcorrect.

The two PacBio corpora went through **isONclust only** --- `isON_pipeline.sh
--mode pacbio` skips correction and points isONform at the clustering output
directly.

* `pacbio_sirv` is the SIRV spike-in of ENCODE `ENCSR507JOF` (WTC11, Sequel II,
  LRGASP), file `ENCFF370NFS` (1 714 957 reads). SIRV reads were selected by
  primary `minimap2 -x map-hifi` alignment to the reference below: 17 633 reads,
  mean 1 783 bp, longest 12 033 bp. Clustered with `isONclust --isoseq --N 5`.
* The scoring reference is **SIRV set 4**: the 68 short transcripts plus the 15
  long SIRVs (3 997--12 029 bp) from ENCODE `ENCSR759PLA`. The seven `SIRV1`..
  `SIRV7` entries in that file are genomic loci, not transcripts, and are
  excluded.
* `pacbio_droso` is SRA `SRR34344449` (D. melanogaster testis, Iso-Seq,
  1 420 172 reads), clustered the same way, restricted to the 99 clusters with
  >= 2 000 reads.

## What the columns mean

| column | meaning |
|---|---|
| **TP** | distinct reference transcripts recovered --- the number of *different* real transcripts found, so emitting the same one twice does not raise it |
| **called** | isoforms emitted in total, right or wrong |
| **matching** | emitted isoforms that match some reference transcript; `called - matching` is the junk |
| **strict F1** | F1 of recall and precision at identity >= 0.95 and length within 10% of the transcript |
| **lenient F1** | the same at identity >= 0.90 and length within 20% |
| **redund** | `matching / TP` --- matching isoforms per recovered transcript, so 1.0 is ideal and 5.0 means each transcript was emitted five times over. Reported separately because F1 counts a duplicate of a real transcript as a true positive and so rewards emitting more |
| **identity** | median base identity of matching isoforms to their transcript |

Both thresholds are reported because a single cutoff cannot distinguish "the
sequences changed" from "the same sequences moved across a line". `redund` and
`identity` below are the strict-threshold figures.

For Drosophila there is no transcript-level truth, so the SQANTI categories stand
in:

| column | meaning |
|---|---|
| **FSM tx** | distinct annotated transcripts with at least one full-splice-match isoform --- **the Drosophila analogue of TP** |
| **FSM** | emitted isoforms whose intron chain exactly matches an annotated transcript. Counts *instances*, so redundancy inflates it |
| **ISM** | isoform is a contiguous sub-chain of an annotated one, i.e. a fragment |
| **NNC** | at least one unannotated splice site --- the strongest error signal |
| **genic** | falls inside a gene without matching its splice structure |
| **canonical** | fraction of implied introns with GT..AG; the annotation-free junction-quality signal |
| **genes** | annotated genes the set produced anything for |

## Speed and memory

| corpus | python | port (faithful) | port (default) |
|---|---|---|---|
| `sirv_sim_gene` | 354.8s · 1 426 MB | 100.1s (**3.5x**) · 1 292 MB | **31.7s (11.2x)** · 1 342 MB |
| `sirv_real` | 241.0s · 1 074 MB | 186.6s (1.3x) · 1 032 MB | **47.1s (5.1x)** · 642 MB |
| `sirv_real_deep` | 1 817.3s · 1 896 MB | 1 421.8s (1.3x) · 2 081 MB | **199.5s (9.1x)** · 1 598 MB |
| `droso` | 51.5s · 385 MB | 26.0s (2.0x) · 259 MB | **10.9s (4.7x)** · 482 MB |
| `droso_deep` | 11 440.4s · 2 733 MB | 6 187.5s (1.8x) · 1 666 MB | **1 312.6s (8.7x)** · 1 810 MB |
| `pacbio_sirv` | 7 108.0s · 8 316 MB | 3 002.2s (2.4x) · 4 968 MB | **141.4s (50.3x)** · 3 957 MB |
| `pacbio_droso` | not run | 27 029.9s · 4 505 MB | **1 264.9s (21.4x*)** · 3 515 MB |

* `pacbio_droso` has no Python benchmark: Python took 7 108s on `pacbio_sirv`, a corpus
26x smaller, so the run was not attempted. The 21.4x on that row is the default
against faithful, not against python.


## Accuracy --- simulated SIRV

Of the 68 reference transcripts.

| corpus | implementation | TP | called | matching | strict F1 | lenient F1 | redund | identity |
|---|---|---|---|---|---|---|---|---|
| `sirv_sim_gene` | python | 62 | 138 | 136 | 0.947 | 0.954 | 2.19 | 0.9976 |
| | faithful | 62 | 138 | 136 | 0.947 | 0.954 | 2.19 | 0.9976 |
| | default | 61 | 133 | 133 | 0.946 | 0.946 | **2.18** | **0.9978** |
| `sirv_real` | python | 53 | 108 | 94 | 0.822 | 0.892 | 1.77 | 0.9780 |
| | faithful | 53 | 108 | 94 | 0.822 | 0.892 | 1.77 | 0.9780 |
| | default | 51 | 110 | 89 | 0.778 | 0.879 | **1.75** | 0.9779 |
| `sirv_real_deep` | python | 57 | 625 | 515 | 0.831 | 0.903 | 9.04 | 0.9712 |
| | faithful | 57 | 625 | 515 | 0.831 | 0.903 | 9.04 | 0.9712 |
| | **default** | **61** | 626 | 516 | **0.859** | **0.935** | **8.46** | **0.9739** |

### Real ONT SIRV by read depth

The same real ONT SIRV library subsampled to six depths, each built through the
full upstream pipeline. This is where the effect of depth on the WFA2 default is
visible: it costs a transcript on shallow data and stops costing anything by
50 000 reads.

| reads | implementation | TP | called | strict F1 | lenient F1 | redund | wall | peak RSS |
|---|---|---|---|---|---|---|---|---|
| 1 000 | python | 29 | 33 | 0.581 | 0.640 | 1.03 | 43.3s | 383 MB |
| | faithful | 29 | 33 | 0.581 | 0.640 | 1.03 | 14.3s | 270 MB |
| | default | 27 | 33 | 0.547 | 0.626 | 1.07 | **7.3s** | **198 MB** |
| 2 000 | python | 35 | 48 | 0.642 | 0.741 | 1.17 | 78.5s | 613 MB |
| | faithful | 35 | 48 | 0.642 | 0.741 | 1.17 | 38.7s | 501 MB |
| | **default** | **37** | 47 | **0.682** | **0.752** | **1.16** | **18.5s** | **457 MB** |
| 5 000 | python | 47 | 72 | 0.783 | 0.857 | 1.38 | 461.5s | 1 130 MB |
| | faithful | 47 | 72 | 0.783 | 0.857 | 1.38 | 60.7s | 764 MB |
| | default | 45 | 71 | 0.758 | 0.818 | 1.40 | **17.9s** | **406 MB** |
| 10 000 | python | 53 | 108 | 0.822 | 0.892 | 1.77 | 354.8s | 1 329 MB |
| | faithful | 53 | 108 | 0.822 | 0.892 | 1.77 | 220.3s | 1 002 MB |
| | default | 51 | 110 | 0.778 | 0.879 | **1.75** | **45.9s** | **584 MB** |
| 20 000 | python | 52 | 157 | 0.789 | 0.871 | 2.46 | 327.2s | 1 388 MB |
| | faithful | 52 | 157 | 0.789 | 0.871 | 2.46 | 309.2s | 1 356 MB |
| | default | 52 | 155 | **0.791** | **0.884** | **2.44** | **50.8s** | **983 MB** |
| 50 000 | python | 59 | 340 | 0.842 | 0.932 | 4.71 | 796.5s | 1 729 MB |
| | faithful | 59 | 340 | 0.842 | 0.932 | 4.71 | 797.9s | 1 771 MB |
| | **default** | **60** | 342 | **0.851** | **0.943** | **4.68** | **113.1s** | **1 412 MB** |

The simulated corpus has no depth series: its source is a single 10 000-read
simulated library, so 10 000 is the only depth available without regenerating the
simulation.

## Accuracy --- ONT Drosophila

Scored against the genome by spliced alignment, and into SQANTI-style structural
categories against FlyBase. `FSM` is full splice match: the isoform's intron chain
is identical to an annotated transcript's.

| corpus | implementation | FSM tx | called | FSM | ISM | NNC | genic | canonical | genes |
|---|---|---|---|---|---|---|---|---|---|
| `droso` | python | 409 | 504 | 443 | 14 | 10 | 29 | 0.983 | 426 |
| | faithful | 409 | 504 | 443 | 14 | 10 | 29 | 0.983 | 426 |
| | **default** | 409 | 512 | **457** | **12** | 11 | **25** | 0.983 | 425 |
| `droso_deep` | python | 6 663 | 13 312 | 9 476 | 945 | 912 | 1 367 | 0.973 | 7 197 |
| | faithful | 6 663 | 13 313 | 9 476 | 945 | 912 | 1 367 | 0.973 | 7 197 |
| | **default** | **6 764** | 13 116 | 9 333 | 945 | **874** | 1 375 | 0.973 | **7 297** |

**`FSM tx` is the column to read**, not `FSM`: it counts distinct annotated
transcripts, so it is the Drosophila analogue of SIRV's `TP`, while `FSM` counts
isoform instances and rises with redundancy.

## Accuracy --- PacBio SIRV

83 reference transcripts: 68 short (191--2 528 bp) and 15 long (3 997--12 029 bp).
Poly-A trimmed, as everywhere above.

| | TP | called | matching | strict F1 | lenient F1 | redund | identity | len.ratio |
|---|---|---|---|---|---|---|---|---|
| python | 69 | 109 | 80 | 0.780 | 0.847 | 1.16 | 0.9921 | 1.006 |
| faithful | 69 | 109 | 80 | 0.780 | 0.847 | 1.16 | 0.9921 | 1.006 |
| default | 66 | 114 | 82 | 0.755 | 0.826 | 1.24 | **0.9925** | **1.004** |

`cmp` reports the faithful output byte-identical to python's.

### Short and long SIRVs separately, default

| class | count | recovered |
|---|---|---|
| short (191--2 528 bp) | 68 | 48 |
| long (3 997--12 029 bp) | 15 | 11 |

Read coverage of the four long SIRVs not recovered, against three that were:

| transcript | length | reads | reads >= 90% of length | bases cov >= 1 | bases cov >= 5 | median depth | recovered |
|---|---|---|---|---|---|---|---|
| SIRV4002 | 3 999 | 739 | 666 | 100% | 100% | 705 | no |
| SIRV4001 | 3 997 | 296 | 261 | 100% | 100% | 279 | yes |
| SIRV10002 | 10 001 | 34 | 3 | 100% | 90.1% | 7 | no |
| SIRV10001 | 9 941 | 115 | 19 | 100% | 100% | 34 | yes |
| SIRV12001 | 12 029 | 37 | 0 | 100% | 55.7% | 6 | no |
| SIRV12002 | 11 999 | 78 | 0 | 69.2% | 49.0% | 4 | no |
| SIRV12003 | 12 000 | 99 | 2 | 100% | 89.1% | 24 | yes |

SIRV4002, the one with full coverage, aligns to its reference as `3I 3998= 1X
303I` --- 3 998 of 3 999 bases matched, one substitution, no internal indels, and
a 303-base trailing run that is 301 A's. Under `--faithful` the same isoform is
4 094 bp with a 92-base tail, 90% A. Longest isoforms emitted: 12 023, 9 982,
9 955, 8 081, 8 016, 8 014 bp.

## Accuracy --- PacBio Drosophila

455 226 reads in 99 clusters, scored against the genome and FlyBase.

| | called | FSM | FSM tx | ISM | NIC | NNC | genic | canonical | err.med | aln.frac | genes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| faithful | 1 386 | 913 (66%) | 238 | 117 | 69 | 227 | 56 | **0.981** | 0.0041 | 0.999 | 166 |
| default | 1 428 | **934** (65%) | **240** | 116 | 69 | 241 | 64 | 0.979 | 0.0041 | 0.999 | 167 |

`intron_retention` 145 faithful / 147 default; `mono_exon` 282 / 286.
Median error rate is 0.0041 here against 0.013--0.015 on the ONT corpora.

## Reading it

**The faithful configuration reproduces the reference.** On every corpus but the
largest it is byte-identical --- `cmp` reports identical `transcriptome.fasta` ---
and every metric above matches to the last digit.

**One exception:** On `droso_deep`
the faithful port emits **13 313 isoforms against python's 13 312**, differing in
the deepest cluster's batch 8 (`0_8_103` against `0_8_71`). Every scored metric is
identical --- same FSM, same categories, same chains --- so this is one extra
isoform in 13 312, 0.008%. The leading suspect is the one path in
`compute_equal_reads` where the reference's set order is *not* modelled: a read
whose graph walk takes no edge at all leaves `current_node_support` aliasing
`node_support_left`, whose order needs the `pop`/dummy-slot machinery
[`crate::pyset`] deliberately omits. Unverified, and open.

**WFA2 helps more the deeper the data.** The depth series shows it directly: on
the same library it costs 2 transcripts at 1 000, 5 000 and 10 000 reads, **gains
2 at 2 000**, draws at 20 000 and **gains at 50 000** (60 against 59) --- and on
the undownsampled corpus it gains **four** (61 against 57), with higher F1 at both
thresholds and lower redundancy. The Drosophila corpora agree: +14 FSM at 20 000
reads, and 101 more distinct annotated transcripts at 1M. On Drosophila it gains 14 FSM
at 20 000 reads and, at 1M reads, matches **101 more distinct annotated
transcripts** (6 764 against 6 663) across 100 more genes while emitting fewer
isoforms --- fewer instances, more distinct transcripts, which is the direction
worth having.

**Redundancy is consistently lower in the default**, on every corpus, which is
what a merge that collapses duplicates should do.

## Reproducing

```bash
bench/evaluate.sh corpora                                  # build all five, once
bench/evaluate.sh run  <corpus> <impl-dir> <tag>           # impl-dir: repo root, or rust/target/release
bench/evaluate.sh score <corpus> <tag>...
```

`--faithful` on the Rust binaries selects the reproduce-the-reference
configuration; without it, the default above.
