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


## The corpora

| corpus | data | reads | clusters | scored against |
|---|---|---|---|---|
| `sirv_sim_gene` | simulated SIRV, grouped by gene | 10 000 | 7 | SIRV transcriptome, exact truth |
| `sirv_real` | real ONT SIRV spike-in | 10 000 | 26 | SIRV transcriptome (recall is a lower bound) |
| `sirv_real_deep` | the same, undownsampled | 99 631 | 26 | as above; deepest cluster 46 437 reads |
| `droso` | real ONT Drosophila cDNA | 20 000 | 561 | genome + FlyBase annotation |
| `droso_deep` | the same, 1M reads | 1 000 000 | 9 697 | as above |

Both SIRV corpora and both Drosophila corpora went through the real upstream
pipeline (isONclust, then isONcorrect).

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

The default is **4.7--11.2x faster than python**; the faithful configuration,
which gives up the faster aligner in exchange for reproducing the reference, is
1.3--3.5x. Peak memory is comparable to python throughout and lower on the two
largest corpora.

## Accuracy --- SIRV

Of the 68 reference transcripts.

| corpus | implementation | TP | called | matching | strict F1 | lenient F1 | redund | identity |
|---|---|---|---|---|---|---|---|---|
| `sirv_sim_gene` | python | 62 | 138 | 136 | 0.947 | 0.954 | 2.19 | 0.9976 |
| | faithful | 62 | 138 | 136 | 0.947 | 0.954 | 2.19 | 0.9976 |
| | default | 61 | 133 | 133 | 0.946 | 0.946 | **2.18** | **0.9978** |
| `sirv_real` | python | 49 | 108 | 90 | 0.773 | 0.892 | 1.84 | 0.9697 |
| | faithful | 49 | 108 | 90 | 0.773 | 0.892 | 1.84 | 0.9697 |
| | default | 48 | 110 | 84 | 0.734 | 0.879 | **1.75** | **0.9710** |
| `sirv_real_deep` | python | 55 | 625 | 490 | 0.796 | 0.902 | 8.91 | 0.9663 |
| | faithful | 55 | 625 | 490 | 0.796 | 0.902 | 8.91 | 0.9663 |
| | **default** | **60** | 626 | 474 | **0.815** | **0.934** | **7.90** | **0.9667** |

### Real SIRV by read depth

The same real ONT SIRV library subsampled to six depths, each built through the
full upstream pipeline. This is where the effect of depth on the WFA2 default is
visible: it costs a transcript on shallow data and stops costing anything by
50 000 reads.

| reads | implementation | TP | called | strict F1 | lenient F1 | redund | wall | peak RSS |
|---|---|---|---|---|---|---|---|---|
| 1 000 | python | 28 | 33 | 0.561 | 0.640 | 1.04 | 43.3s | 383 MB |
| | faithful | 28 | 33 | 0.561 | 0.640 | 1.04 | 14.3s | 270 MB |
| | default | 26 | 33 | 0.527 | 0.626 | 1.08 | **7.3s** | **198 MB** |
| 2 000 | python | 35 | 48 | 0.642 | 0.741 | 1.17 | 78.5s | 613 MB |
| | faithful | 35 | 48 | 0.642 | 0.741 | 1.17 | 38.7s | 501 MB |
| | default | 35 | 47 | **0.647** | **0.752** | **1.15** | **18.5s** | **457 MB** |
| 5 000 | python | 46 | 72 | 0.763 | 0.857 | 1.37 | 461.5s | 1 130 MB |
| | faithful | 46 | 72 | 0.763 | 0.857 | 1.37 | 60.7s | 764 MB |
| | default | 44 | 71 | 0.738 | 0.818 | 1.39 | **17.9s** | **406 MB** |
| 10 000 | python | 49 | 108 | 0.773 | 0.892 | 1.84 | 354.8s | 1 329 MB |
| | faithful | 49 | 108 | 0.773 | 0.892 | 1.84 | 220.3s | 1 002 MB |
| | default | 48 | 110 | 0.734 | 0.879 | **1.75** | **45.9s** | **584 MB** |
| 20 000 | python | 50 | 157 | 0.756 | 0.871 | 2.44 | 327.2s | 1 388 MB |
| | faithful | 50 | 157 | 0.756 | 0.871 | 2.44 | 309.2s | 1 356 MB |
| | default | 49 | 155 | 0.746 | **0.873** | 2.45 | **50.8s** | **983 MB** |
| 50 000 | python | 56 | 340 | 0.791 | 0.929 | 4.62 | 796.5s | 1 729 MB |
| | faithful | 56 | 340 | 0.791 | 0.929 | 4.62 | 797.9s | 1 771 MB |
| | default | 56 | 342 | 0.783 | 0.924 | 4.55 | **113.1s** | **1 412 MB** |



## Accuracy --- Drosophila

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

## Reading it

**The faithful configuration reproduces the reference.** On every corpus but the
largest it is byte-identical --- `cmp` reports identical `transcriptome.fasta` ---
and every metric above matches to the last digit.

**One exception, and it is recorded rather than smoothed over.** On `droso_deep`
the faithful port emits **13 313 isoforms against python's 13 312**, differing in
the deepest cluster's batch 8 (`0_8_103` against `0_8_71`). Every scored metric is
identical --- same FSM, same categories, same chains --- so this is one extra
isoform in 13 312, 0.008%. The leading suspect is the one path in
`compute_equal_reads` where the reference's set order is *not* modelled: a read
whose graph walk takes no edge at all leaves `current_node_support` aliasing
`node_support_left`, whose order needs the `pop`/dummy-slot machinery
[`crate::pyset`] deliberately omits. Unverified, and open.

**WFA2 helps more the deeper the data.** The depth series shows it directly: on
the same library it costs 2 transcripts at 1 000 reads, 2 at 5 000, 1 at 10 000
and 20 000, and **nothing at 50 000** --- and on the undownsampled corpus it gains
**five** (60 against 55), with higher F1 at both thresholds and lower redundancy. On Drosophila it gains 14 FSM
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
