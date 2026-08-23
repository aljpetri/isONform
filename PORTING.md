# isONform — Rust rewrite

## Status: CLI ported and verified, reference determinism fixed, evaluation corpora built

This document is the plan, the reconnaissance, and the method carried over from the isONcorrect
port, which is finished and released as v0.2.0.

Done so far:

- **Pinned reference environment** (`bench/env/reference.yml`, `bench/setup_reference_env.sh`,
  resolved versions in `bench/env/resolved-macosx-11.0-arm64.txt`). `conda env` named
  `isonform-ref`: python 3.10.21, networkx 2.8.4, parasail 1.3.4, spoa 4.1.5.
- **The reference runs**, on a corpus built through the real upstream pipeline
  (`bench/corpus/sirv_small`, see `bench/corpus/README.md`). Both entry points, exit 0, isoforms
  emitted.
- **`--disable_numpy`, `--compression` and `PYTHONHASHSEED` checked before recording anything**,
  as the plan demanded. One of the three is a stop-the-line finding; see *Findings in the
  reference* below.
- **Reconnaissance corrected** — several load-bearing assumptions in the section below were wrong
  in the port's favour. See *Reconnaissance corrections*.
- **CLI contract captured** from the reference into `bench/golden/cli/` (`--help`, `--version`,
  no-args, and the validation paths).
- **CLI ported and differentially verified.** `rust/` builds two binaries under the contract names
  `main` and `isONform_parallel`; `cargo test` runs 33 unit tests over `cli.rs`, and
  `bench/equivalence.sh cli` runs **37 cases against the live reference, all passing** — exit codes,
  the reference's own stderr diagnostics verbatim, the stray stdout prints, argparse's prefix
  abbreviation, and the full resolved argument set. Clippy is `-D warnings` and `cargo fmt --check`
  is clean.
- **CI on day one**, `.github/workflows/rust.yml`: build, test, clippy, fmt and the entry-point
  name check across linux x86_64, linux aarch64, macos arm64 and macos x86_64, plus the
  differential harness in the pinned conda env on two of them.
- **The determinism defect is fixed in the Python** — minimizer selection is lexicographic, matching
  isONcorrect. `bench/check_seed_sensitivity.py` now exits 0 at 24 seeds for both entry points, and
  it is a hard gate in CI. Measured, not assumed: accuracy-neutral on simulated reads and a small
  consistent gain on real ONT reads (48 of 68 SIRV transcripts recovered against 46 for *both* hash
  seeds). This is what unblocks stage goldens. See *Finding 1* and *The determinism fix, measured*.
- **Three evaluation corpora**, built through the real upstream pipeline and scored by two new
  harnesses (`bench/accuracy_isoforms.py` against a transcriptome,
  `bench/accuracy_isoforms_genome.py` against a genome). `bench/evaluate.sh` builds and runs them.
  All three agree that the determinism fix costs nothing. They also disagree sharply with each other
  about how good isONform is — 89.7% recall on simulated reads, 70.6% on real ones — which is the
  entire reason there are three. See *Evaluating isONform*.

Not done: no isoform-generating Rust code yet. Stage goldens are now recordable and are the next
step.

### The argument-parity oracle, for free

`main:343` does `print("ARGS", args)` — a debug leftover that dumps the whole argparse `Namespace`
to stdout on every run. Since it is observable, the port reproduces it
(`MainArgs::args_line`), and it turns out to pay for itself: one line carrying every resolved
argument is a complete differential oracle for argument parsing, so a wrong default or a
mis-parsed value is caught without a harness case per flag. Three of `cli.rs`'s unit tests assert
it against strings captured verbatim from the reference.

That is method point 2 applied to the CLI, and it is worth remembering when the print is finally
deleted upstream: the oracle should move behind an env var rather than disappear.

> **Permission: granted.** This is a co-author's repository
> (<https://github.com/aljpetri/isONform>) and the owner approved the work on 2026-08-23. The
> *Working agreements* below still stand, and the never-commit rule matters more here than in our own
> repos, not less: leave work in the tree and propose the message.

## Goal

Port isONform from Python to Rust. **Identical CLI**, faster, lower memory. Same objectives as the
isONcorrect port, including the one that mattered most there: **the specification is accuracy, not
byte-identity**. Byte-identity is the tool used to prove the port faithful stage by stage; it is not
the thing being optimised, and where the reference is wrong or slow-by-compromise, the port may
diverge — measured, documented, and behind an env var that restores the reference path.

Two entry points must keep their names, flags and defaults:

- `main` — the single-instance corrector/assembler
- `isONform_parallel` — fans out over a folder of per-cluster fastqs

## Reconnaissance

Measured on `f0cd3c7` (2026-08-23), before any work:

| | |
| --- | --- |
| Python source | **~4 500 lines** (isONcorrect was ~3 300 before trimming) |
| dependencies | `networkx==2.8.4`, `parasail==1.2`, setuptools |
| external binaries | `spoa` (two different invocations), and optionally `minimap2`, `racon`, `medaka_consensus` |

Largest modules, which is where the work is:

| file | lines | what |
| --- | --- | --- |
| `modules/SimplifyGraph.py` | 1 166 | graph simplification |
| `modules/IsoformGeneration.py` | 631 | isoform calling |
| `modules/GraphGeneration.py` | 581 | graph construction |
| `main` | 639 | single-instance driver |
| `isONform_parallel` | 379 | batch driver |
| `modules/batch_merging_parallel.py` | 306 | merging across batches |
| `modules/consensus.py` | 247 | spoa/racon/medaka wrappers |

### What transfers directly from isONcorrect

The CLI overlap is large enough to be worth stating: `--fastq`, `--k`, `--w`, `--xmin`, `--xmax`,
`--exact`, `--max_seqs`, `--max_seqs_to_spoa`, `--exact_instance_limit`, `--set_w_dynamically`,
`--verbose`, `--outfolder`, `--compression`, `--disable_numpy` are all shared, and
`isONform_parallel` mirrors `run_isoncorrect` with `--fastq_folder`, `--t`, `--keep_old`,
`--split_wrt_batches`. So `cli.rs`, `params.rs` and `validate.rs` are a starting point rather than a
blank page.

| isONcorrect module | reusable here? |
| --- | --- |
| `parasail.rs` (exact affine semi-global + global, tie-break measured) | **yes** — isONform depends on `parasail==1.2` |
| `poa.rs` (spoars, `-l 0 -r 0 -g -2`) | **yes for one of the two spoa calls** — the invocation is character-for-character the same |
| `editdist.rs`, `align.rs`, `fastq.rs`, `profile.rs`, `simd.rs` | likely, pending a check of whether isONform uses edlib at all (it does not import it) |
| `bench/` harness: `dump_reference.py`, `equivalence.sh`, `accuracy.py`, `accuracy_genome.py`, `profile_rust.py` | **yes, as the method** — the scripts need isONform's stages, but the shape carries over |

### What is genuinely new

- **The graph — and it probably should not be a graph *library*.** isONform is graph-based where
  isONcorrect was not, but the networkx surface is small and entirely basic: `DiGraph`, add/remove
  node and edge, adjacency, `successors`, `predecessors`, `out_edges`, degrees, `has_edge`,
  `topological_sort`, `find_cycle`, and node/edge attributes. No PageRank, max-flow, isomorphism or
  community detection.

  **So do not reach for `petgraph` by default.** Nothing here needs it, and object-heavy graphs are
  where the time goes — networkx allocates a Python object per node plus a dict of attributes per
  node and per edge, which is a large part of why the reference is slow. Prefer the shape that won
  repeatedly in the isONcorrect port: **nodes as integer indices, CSR-style adjacency in flat `Vec`s,
  attributes as parallel arrays rather than per-node maps.** A flat buffer plus offsets beat
  `Vec<Vec<u8>>` there by 2x; one pass beat a `BTreeSet` per slot by 1.5x. `topological_sort` is
  Kahn's algorithm and `find_cycle` is a DFS three-colouring — a few lines each, no dependency.

  This is also worth telling the Python side: the same representation change would speed up the
  reference, and it is the kind of thing that can be measured before the port exists.

  **The work is isONform's own ~2 400 lines of algorithms over that graph, not the container.**
- **A second spoa configuration**: `-q <reads> -l 2 -r 2 -x -8 -m 10 -o ...`, which is not the
  `-l 0 -r 0 -g -2` that `poa.rs` reproduces. Resolve `-l 2`/`-r 2` against spoa's CLI the same way
  the isONcorrect port resolved the first one, and check whether `spoars` covers it before assuming.
- **`minimap2`, `racon`, `medaka_consensus`** subprocesses. Decide per call site whether each is on
  the default path or an optional extra, exactly as `--use_racon` was triaged in isONcorrect.
- **Isoform output**, which has no analogue in isONcorrect and needs its own equivalence contract.

### Expect the reference to be wrong

The owner's own expectation, stated up front: **the Python is likely bug-prone.** That is not a
criticism, it is a design input, and isONcorrect makes it credible — nine defects found there, and the
*two largest findings of the entire project* were reference bugs rather than porting bugs (a
`PYTHONHASHSEED`-dependent `set`, and an off-by-one predecessor table that made the tool correct fewer
regions than intended, worth more accuracy than both deliberate divergences combined).

The practical consequence changes how oracles are read. **A mismatch is evidence, not a verdict.**
When one fires, the first question is *which side is wrong*, and answering it needs a check
independent of both implementations:

- **brute force over small inputs** — how `solve_WIS` was proven suboptimal on 2040 of 3000 random
  instances, which is what justified fixing it rather than reproducing it;
- **replay under many `PYTHONHASHSEED` values** — how the `set` defect was found. Note the sample size
  mattered: at 6 seeds six records looked like porting bugs; at 24 they were all reference behaviour;
- **invariants the reference should satisfy** — a graph that should be acyclic, a partition that should
  be disjoint, counts that should agree between two stages.

Graph code is a good hunting ground for all three, since acyclicity, connectivity and node/edge
counts are all checkable without trusting either implementation.

### Two things to check early, because they were traps in isONcorrect

- **`--disable_numpy` and `--compression` appear here too.** In isONcorrect `--disable_numpy` had
  *never worked* — it raised `ValueError` on any input that reached it — and `--compression` was
  deprecated and unreachable from the batch driver. Test both here before deciding to port them.
- **Determinism.** isONcorrect had a defect where a Python `set` made corrected output depend on
  `PYTHONHASHSEED`. Graph code iterates containers constantly, so check for seed-dependence *before*
  recording any goldens: run the same input under several `PYTHONHASHSEED` values and diff. A harness
  that pins the seed in-process does not work — CPython reads it once at startup, so the harness must
  re-exec itself.

## Reconnaissance corrections

Everything in the section above was measured before any code was read closely. Five of its
load-bearing claims turned out to be wrong, all of them in the port's favour. Recorded here rather
than edited into the section above, because *how much the first estimate was off by* is itself
useful — the recon overstated the work.

### 1. `minimap2`, `racon` and `medaka_consensus` are unreachable. Not optional — unreachable.

`modules/consensus.py` does wrap them (`run_racon`, `run_medaka`, and `minimap2` inside
`run_racon`), and `polish_sequences` dispatches on `args.medaka` / `args.racon`. But **neither
entry point defines those flags**, and nothing outside `modules/consensus.py` calls any of those
functions:

| symbol in `modules/consensus.py` | call sites outside the module |
| --- | --- |
| `parasail_alignment` | 2 (`IsoformGeneration.py:381`, `SimplifyGraph.py:644`) |
| `run_spoa` | 0 |
| `run_medaka` | 0 |
| `run_racon` | 0 |
| `polish_sequences` | 0 |
| `form_draft_consensus` | 0 |
| `detect_reverse_complements` | 0 |

So the whole module reduces to one function, `parasail_alignment`, and there is no
"decide per call site whether each is on the default path" triage to do. **Three external binaries
drop out of the port entirely.**

### 2. There is no second spoa configuration.

The `-q <reads> -l 2 -r 2 -x -8 -m 10 -o …` invocation the recon flagged lives only in
`modules/create_augmented_reference.py` — five variants of it, plus three more commented out.
**That module is imported by nothing**; the single mention of its name anywhere in the repository
is a commented-out line in `modules/consensus.py:243`.

The only spoa invocation reachable from either entry point is `IsoformGeneration.run_spoa`
(`modules/IsoformGeneration.py:148`):

```python
subprocess.check_call(["spoa", reads, "-l", "0", "-r", "0", "-g", "-2"], ...)
```

which is character-for-character the invocation `poa.rs` already reproduces. **`poa.rs` transfers
whole, and the "resolve `-l 2`/`-r 2` against spoa's CLI" task does not exist.**

### 3. The dependency list is shorter than either of the two the repository publishes.

`requirements.txt` and `setup.py` disagree with each other, and both disagree with the code. What
is actually imported, across `main`, `isONform_parallel` and every reachable module, is
**`networkx` and `parasail`** and nothing else third-party.

| claimed | by | reality |
| --- | --- | --- |
| `networkx==2.8.4` | requirements.txt | used |
| `parasail==1.2` | requirements.txt | used — but `setup.py` demands `>=1.3.3`, so the repo contradicts itself, and 1.2 has no osx-arm64 build |
| `setuptools==78.1.1` | requirements.txt | build-time only |
| `edlib>=1.1.2` | setup.py | **never imported** |
| `recordclass>=0.17.2` | setup.py | **never imported** |
| `numpy` | `--disable_numpy`'s help text | **never imported** (see finding 3 below) |
| `ordered-set`, `matplotlib`, `pyinstrument`, `namedtuple` | README | **never imported** |

The only parasail surface used is `matrix_create` plus `sg_trace_scan_16` with a `_32` fallback on
saturation — the same call `isONcorrect` makes, at different scores (`match=2, mismatch=-2,
open=12, ext=1` here versus `4, -8, 12, 1` there). `parasail.rs` is documented as reproducing
exactly `sg_trace_scan_16`, so **it transfers as-is, parameterised by `Scoring`.**

`modules/consensus.py:13` and `modules/help_functions.py:41` are two copies of `cigar_to_seq`
differing only in whether they also return `cigar_tuples` (and a stray `print("error")`) — the same
duplication isONcorrect has, and `align.rs::cigar_to_seq` covers both.

### 4. Far more than the CLI transfers: the entire front half of `main` is isONcorrect's.

The recon framed the reuse as "`cli.rs`, `params.rs` and `validate.rs` are a starting point". It is
much larger than that. `main` lines 81–340 are isONcorrect's interval machinery, function for
function: `get_kmer_minimizers`, `get_kmer_maximizers`, `get_minimizers_and_positions`,
`get_minimizers_and_positions_compressed`, `get_minimizer_combinations_database`,
`minimizers_comb_iterator`, `fill_p2`, `solve_WIS`, `get_intervals_to_correct`, `batch`, `grouper`,
`find_most_supported_span`, `rindex`, `remove_read_polyA_ends`.

That is `minimizers.rs`, `anchors.rs`, `contexts.rs` and `wis.rs` — plus the `solve_WIS`
suboptimality already proven and fixed there by brute force. What isONform does differently is what
happens *after* the intervals are chosen: instead of correcting them it feeds them to a graph. **The
~2 400 lines of genuinely new algorithm are `GraphGeneration` → `SimplifyGraph` →
`IsoformGeneration` → `batch_merging_parallel`, and nothing before them.**

One difference inside that shared front half is not cosmetic, and it is the subject of finding 1
below.

### 5. The two entry points' CLIs overlap much less than stated.

The recon listed 14 shared flags. Per entry point, actually:

| flag | `main` | `isONform_parallel` |
| --- | --- | --- |
| `--k --w --xmin --xmax --max_seqs --max_seqs_to_spoa --exact_instance_limit --set_w_dynamically --verbose --outfolder --delta_len --delta_iso_len_3 --delta_iso_len_5` | yes | yes |
| `--fastq` | yes | — (`--fastq_folder`) |
| `--T --exact --disable_numpy --compression --parallel --slow` | yes | **no** |
| `--t --keep_old --split_wrt_batches --iso_abundance --delta --tmpdir --write_fastq` | **no** | yes |

Defaults differ where the flags are shared: `--exact_instance_limit` is `0` in `main` and `50` in
`isONform_parallel`. And `--xmin`/`--xmax` help text is *swapped* in `main` ("Upper interval length"
for `--xmin`) and correct in `isONform_parallel`; `--help` output is part of the CLI contract, so
the port reproduces the swap and a fix is a separate commit.

Two flags exist on `isONform_parallel` but cannot reach the per-cluster work, because
`isONform_parallel:79` builds the child command line by hand and does not forward them:
`--max_seqs_to_spoa` and `--delta` affect only batch merging, and `--max_seqs` only file splitting.
`main` hard-codes `delta = 0.15` (`main:560`) while `isONform_parallel --delta` defaults to `0.1`,
so the two halves of one run use different diversity rates. `--T` is commented out of that command
line too — and inert anyway (finding 3).

## Findings in the reference

Nine defects came out of the isONcorrect port and the two largest were reference bugs rather than
porting bugs. Same pattern here: the first day of reconnaissance produced these, before a line of
algorithm was ported.

### Finding 1 — isONform was non-deterministic run to run. Every output file, every run. **Fixed.**

**This blocked recording goldens and blocked the port, so it came first.** The fix is described in
*The determinism fix, measured* below; what follows is the defect as found.

`main:85`, inside `get_kmer_minimizers`:

```python
window_kmers = deque([hash(seq[i:i+k_size]) for i in range(w + 1)])
```

Minimizers are selected by **`hash()` of the k-mer string**. CPython's `str.__hash__` is SipHash
keyed by `PYTHONHASHSEED`, which is *random per process* unless set in the environment before the
interpreter starts. isONcorrect's `get_kmer_minimizers_lex` is otherwise the same function
character for character — same `rindex`, same `minimizer_pos < i - w` guard — but compares the
k-mer **strings**. The `hash()` is the entire difference, and it makes every minimizer, and
therefore everything downstream, seed-dependent.

The code twice tries to pin the seed and twice fails:

- `main:631`, `os.environ['PYTHONHASHSEED'] = '0'`, executed *after* interpreter startup. CPython
  reads the variable once, at startup. This affects child processes only, and `main` spawns none
  that care. (PORTING.md's method section already warned about exactly this mistake; it is in the
  reference.)
- `isONform_parallel:212`, `#PYTHONHASHSEED = 0`, commented out — and a module-level Python
  assignment would not have been read by anything even if it had executed.

Measured on `bench/corpus/sirv_small` with `bench/check_seed_sensitivity.py`:

| | |
| --- | --- |
| `main`, 8 seeds | `mapping0`, `skip0.fa`, `spoa0` each took **8 distinct values**; only the input-echo `0batch` was stable |
| `isONform_parallel`, 8 seeds | **all 9 output files** varied; `transcriptome.fasta` took **5 distinct values** |
| minimizer combinations generated | ranged **10 549 – 12 532** across those 8 seeds |
| reads skipped for low interval abundance | ranged **14 – 17** of 100 |
| `PYTHONHASHSEED` unset, 3 runs of `main` | 3 different outputs |
| `PYTHONHASHSEED` unset, 3 runs of `isONform_parallel` | `transcriptome.fasta` differed between runs |
| `PYTHONHASHSEED` pinned, 3 runs of each | identical |

The divergence is not a reordering. Under seed 0 read 3 joins the large isoform; under seed 1 it
forms a group with reads 73 and 48. Different isoforms, different read-to-isoform assignments, a
different number of reads dropped.

**An invariant proves it without reference to any second implementation.** Clusters 0 and 1 of
`sirv_small` are byte-identical inputs (they come from the same file — see
`bench/corpus/README.md`), so they must produce the same isoform. `isONform_parallel` runs them as
two child processes, each of which draws *its own* seed, and in one observed run cluster 0 emitted
`0_0_71` where cluster 1 emitted `1_1_2` — different sequences from identical input, in a single
run of the shipped tool. Nothing about the port is involved in that.

Consequences, in order:

1. Every bench invocation pins `PYTHONHASHSEED`, and still does as belt and braces.
   `bench/check_seed_sensitivity.py` is the regression check.
2. **The port could not reproduce this and should not have tried.** Matching seeded SipHash means
   shipping CPython's hash and threading a seed through, to reproduce a coin toss.
3. So minimizer selection had to become deterministic in the *Python* first — its own commit, ahead
   of the port. Done: see below.
4. It needed measuring, not just fixing. "How much does the answer move between seeds" is a
   reproducibility number this tool never had, and the accuracy comparison between seed-roulette and
   a fixed rule is what says whether the fix costs anything.

### Invariants, and one that turned out to be a real detector

`bench/equivalence.sh invariants` holds checks that trust neither implementation — the third of the
three independent techniques this document asks for. Two so far:

- **The identical-cluster invariant.** `sirv_small`'s two clusters are byte-identical inputs, so they
  must yield the same isoform sequences with accessions differing only in the cluster and batch ids.
  `isONform_parallel` runs them as separate child processes, so it catches per-instance state
  leakage and fan-out bugs.
- **Skipped-read accounting.** The count `main` prints must equal the number of records it wrote to
  `skip<batch>.fa`, and neither may exceed the input. A drift means reads are being lost between the
  interval stage and the output, which a golden would not necessarily reveal.

The first was validated against the defect it was written for rather than assumed to work: run on the
pre-fix tree with the seed unset it **fires on 2 of 3 runs**, and on the fixed tree it holds 3 of 3.
That it holds on one pre-fix run in three is the useful detail — the two children sometimes draw
seeds that happen to agree — and it is why an invariant like this belongs in a harness that runs
repeatedly rather than being checked once by hand.

A caution worth recording, because it cost time: the obvious way to write the id normalisation
(`sed 's/^>0_/>X_/'`) strips only the cluster id and leaves the batch id, so it reports a violation
on output that is in fact correct (`0_0_72` vs `1_1_72`). A check that cries wolf is worse than no
check. The scripted version normalises both fields.

Graph code is good hunting ground for more of these — acyclicity, connectivity, node and edge counts
agreeing between stages — and they should be added as the graph stages get ported.

### The determinism fix, measured

Out of numbering order on purpose: this sits next to the defect it closes rather than at
the end of the list.

The change is two lines in `get_kmer_minimizers` (`main:81`): compare the k-mer **strings** instead
of `hash()` of them, which is what isONcorrect's `get_kmer_minimizers_lex` has always done. The two
dead seed-pinning lines go with it — `os.environ['PYTHONHASHSEED'] = '0'` in `main` and the
commented-out `#PYTHONHASHSEED = 0` in `isONform_parallel` — since both existed only to (fail to)
contain this defect.

Determinism, after:

| check | result |
| --- | --- |
| `main`, 24 seeds, `sirv_small` | **all 4 output files identical** |
| `isONform_parallel`, 8 seeds, `sirv_small` | **all 9 output files identical** |
| the identical-cluster invariant | **holds** — cluster 0 emits `0_0_72`, cluster 1 emits `1_1_72`, sequences byte-identical. Before the fix these were different isoforms |

`bench/check_seed_sensitivity.py` exits 0 and is now a regression check.
`bench/equivalence.sh --record` no longer refuses outright: it runs the reference at two seeds and
compares before recording, so the gate cannot go stale the way a hardcoded refusal would.

#### What the fix costs: nothing, and on real data it gains a little

The question the plan asked was not "does it still work" but "what does the fix cost in accuracy".
Measured on two corpora, comparing lexicographic selection against the old hash at **two different
seeds** — because a single run of seeded code is one sample and says nothing about what the tool
does.

`sirv_sim_gene` (10 000 simulated reads, 7 gene-level clusters), strict thresholds:

| | isoforms | recall | precision | F1 | redundancy | identity |
| --- | --- | --- | --- | --- | --- | --- |
| **lexicographic** | 138 | 89.7% (61/68) | 97.8% | 0.936 | **2.21** | 0.9978 |
| hash, seed 0 | 137 | 89.7% (61/68) | 99.3% | 0.942 | 2.23 | 0.9977 |
| hash, seed 1 | 137 | 89.7% (61/68) | 100.0% | 0.946 | 2.25 | **0.9979** |

Indistinguishable. Recall is identical to the transcript in all three; lex falls *between* the two
hash seeds on identity and is best on redundancy; the precision spread between the two hash seeds
(99.3% → 100.0%) is the same order as the gap to lex. At the lenient threshold all three sit at
100.0%, so the strict-threshold differences are borderline lengths, not wrong sequences.

`sirv_real` (9 968 real ONT reads, 26 clusters), strict thresholds:

| | isoforms | recall | precision | F1 | redundancy | identity |
| --- | --- | --- | --- | --- | --- | --- |
| **lexicographic** | 106 | **70.6% (48/68)** | **82.1%** | **0.759** | **1.81** | 0.9673 |
| hash, seed 0 | 108 | 67.6% (46/68) | 80.6% | 0.735 | 1.89 | **0.9701** |
| hash, seed 1 | 101 | 67.6% (46/68) | 79.2% | 0.730 | 1.74 | 0.9696 |

Here lex is ahead of **both** seeds, not merely inside their spread: 48 reference transcripts
recovered against 46, with better precision and F1 than either. The same ordering holds at the
lenient threshold (79.4% recall against 77.9% for both seeds). Identity marginally favours hash, by
0.3 of a percentage point. So the fix is accuracy-neutral on simulated data and a small consistent
gain on real data — which is the direction one would hope for, since lexicographic selection also
turned out to sit closer to the ideal minimizer density.

**Note how much lower every real-data number is** — 70.6% recall and 0.967 identity against 89.7%
and 0.998 on simulated reads, same tool and same flags. That gap is the whole argument for method
point 4: a report built on the simulated corpus alone would have described a considerably better
tool than this one.

`droso` (13 164 real ONT Drosophila reads, 561 clusters), scored against the genome — no
transcriptome truth exists here, so this is the splice-level check:

| | isoforms | aligned | err. median | aligned frac. | multi-exon | introns | **canonical** | chains |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **lexicographic** | 505 | 504 | 0.0134 | 0.978 | 385 | 779 | **98.2%** | 359 |
| hash, seed 0 | 503 | 503 | 0.0129 | 0.978 | 391 | 788 | **98.4%** | 355 |

Indistinguishable again, and the canonical splice fraction is the number to read: **98.2% against
98.4%** of implied introns carrying `GT..AG`. Both are close to what real annotated junctions look
like, which is a useful independent signal that isONform is reconstructing genuine splice structures
rather than plausible-looking assembly artifacts — and that replacing the minimizer rule did not
disturb them. Aligned fraction is identical at 0.978, so neither set is winning on error rate by
aligning less.

`chains` (359) below `multi-exon` (385) says roughly 26 isoforms duplicate another's splice
structure, consistent with the ~1.8 redundancy measured on `sirv_real`.

**All three corpora move the same way, or rather none of them moves.** That is the test method
point 4 asks for: a change believed on the simulated corpus alone is a change fitted to simulated
error.

#### The instability the fix removes

The aggregate metrics above barely move. The *output* did, enormously. Comparing runs by exact
isoform sequence:

| comparison | corpus | shared | Jaccard |
| --- | --- | --- | --- |
| hash seed 0 vs hash seed 1 | `sirv_sim_gene` | 76 of 127 each | **0.427** |
| hash seed 0 vs hash seed 1 | `sirv_real` | 38 | 0.222 |
| lex vs hash seed 0 | `sirv_real` | 33 | 0.182 |
| lex vs hash seed 1 | `sirv_real` | 30 | 0.169 |

On the simulated corpus, where isoforms are near-perfect (identity 0.998) so an exact-sequence match
is a meaningful test, **fewer than half the reconstructed isoforms survived a change of hash seed**.
The old code was not less accurate on average; it was **differently wrong every run**, and the
aggregate metrics hid that completely. Any user diffing two runs, or any reviewer trying to reproduce
a figure, was looking at a coin toss.

The real-data rows also settle whether replacing the rule is a bigger change than re-rolling the
seed. It is not: two hash seeds disagree with each other (0.222) about as much as either disagrees
with lex (0.182, 0.169). Changing the seed was already as large a change as changing the rule — which
is the clearest possible statement of why the seed could not be left in charge.

Read the Jaccard figures on real data with care, though. Isoform identity there is ~0.967, so two
near-identical isoforms differing by a single base count as entirely different sequences; the metric
is harsh on real data and only really diagnostic on the simulated corpus.

#### Runtime: no claim either way

Recorded so it is not guessed at later: `sirv_real` took 242 s under lex against 200 s and 231 s for
the two hash seeds, and `droso` 51 s against 54 s. The `sirv_sim_gene` hash arms took 480 s and
475 s, and the lex arm there was not timed comparably. Mixed, within the noise of an 8-thread laptop
run, and **no runtime claim is made in either direction**. Lexicographic selection does yield slightly more minimizers per base (below),
so if anything it should cost a little more work.

#### The other `hash()` in the codebase is fine, and should be left alone

`modules/GraphGeneration.py:133`, in `convert_array_to_hash`, does `hash(tuple(info_array))`. That
one is **not** a defect: `info_array` is an `array("I")` of integers, and CPython randomises only
`str`/`bytes` hashing, not `int` or tuples of them. Verified directly — the same tuple hashes
identically at `PYTHONHASHSEED` 0, 1, 42 and 999 — and confirmed empirically by the 24-seed check
passing with it untouched. Worth recording so nobody "fixes" it.

#### Was lexicographic selection the right choice? The poly-A worry, measured

The obvious objection to lexicographic minimizers is that poly-A is lexicographically smallest, so
`AAAA…` should be selected disproportionately and minimizers should clump. There is even
circumstantial evidence the reference once feared this: `get_minimizer_combinations_database` carries
a `forbidden = 'A' * k` guard that skips a combination when both ends are pure poly-A — a guard that
makes sense only for lexicographic selection, which suggests the `hash()` came later.

Measured on 400 reads of `sirv_real` cluster 0 at `k=20, w=31`, comparing the two implementations
directly:

| | lexicographic | `hash()` |
| --- | --- | --- |
| minimizers per 100 bases | **16.80** | 14.98 |
| theoretical ideal for `w=31, k=20` | 16.67 | 16.67 |
| distinct k-mers, as a share of selections | 1.7% | 1.8% |
| minimizers that are ≥80% one base | **0** | 0 |
| pure poly-A minimizers | **0** | 0 |

So the worry does not materialise on real data, and lexicographic selection is in fact *closer* to
the ideal minimizer density than hashing was — hashing under-selected. The reason poly-A does not
dominate is that `remove_read_polyA_ends` already collapses terminal poly-A runs longer than 12 down
to 1 before minimizers are taken.

This is worth re-measuring on a corpus with genuine internal low-complexity, which SIRV does not
have.

### Finding 2 — `--parallel` accepts any value as true, including `False`.

`main:619`: `parser.add_argument('--parallel', type=bool, default=False, ...)`. argparse applies
`type` to the *string*, and `bool("False")` is `True`. `isONform_parallel:79` invokes the child with
`"--parallel", "True"`, which works by accident; `--parallel False` also enables it. Measured: with
`--parallel True` and `--parallel False` the output is identical and differs from the default
(the batch pickle is named `0_batch` rather than `0batch`).

The flag is internal — set only by `isONform_parallel` — so the fix is `action="store_true"`, or
dropping it from the public surface. Either is a behaviour change and needs its own commit.

### Finding 3 — four of `main`'s twenty flags do nothing.

Counted as uses of `args.<name>` in `main`, then confirmed by running each flag and diffing every
output file at a pinned seed:

| flag | uses | effect |
| --- | --- | --- |
| `--T` | **0** | none. Also commented out of `isONform_parallel`'s child command line |
| `--disable_numpy` | **0** | none — numpy is never imported, so its help text ("about 1.5x slower") describes nothing |
| `--exact` | 2 | **none.** Both uses reset `previously_corrected_regions`, which is *never assigned to* anywhere in `main` — only read and `del`eted. It is always empty, so `--exact` is a no-op. `isONform_parallel` passes `--exact` on every child invocation |
| `--exact_instance_limit` | 1 | **none** — its only effect is `args.exact = True` |
| `--compression` | 1 | **works** (output differs), unlike in isONcorrect where it was unreachable. Must be ported |
| `--slow` | 1 | works (output differs) |
| `--verbose` | 1 | stderr only |

So the isONcorrect traps did not repeat, but they inverted: there `--disable_numpy` raised
`ValueError` and `--compression` was unreachable; here `--disable_numpy` is silently inert and
`--compression` is live and load-bearing.

The dead `previously_corrected_regions` also means `main:438–461` — the `pos_group` construction —
is unreachable, along with the `set_not_prev` / `not_prev_corrected_spans2` filtering at
`main:467–476`. Roughly 40 lines that always take the same branch. Worth *not* porting, and worth
telling the Python side, but note it changes nothing observable, so it is a representation choice
rather than a behaviour change.

### Finding 4 — two stray debug prints are part of the observable output.

- `main:343`, `print("ARGS", args)` — dumps the whole argparse `Namespace` to stdout on every run.
- `isONform_parallel:370`, `print(len(sys.argv))` — prints the argument count to stdout before
  anything else, so `isONform_parallel` with no arguments prints `1` and then the help text.

Both are on stdout, so both are in the CLI contract and `bench/golden/cli/` records them. They are
one-line deletions upstream, and the port should not enshrine them — but the fix is a commit of its
own and the goldens move with it.

### Finding 5 — `isONform_parallel` rewrites its own input directory.

`main()` calls `restructure_isoncorrect_output(args.fastq_folder)` before anything else
(`isONform_parallel:217`). On an isONcorrect-shaped input that function `shutil.move`s every
`<cl_id>/corrected_reads.fastq` up to `<cl_id>.fastq` and then calls
`Parallelization_side_functions.remove_folders(directory)` on the input folder. So pointing the tool
at an isONcorrect output directory *destroys it* — the reads are still there under new names, but
the isONcorrect run's structure and anything else in those folders is gone, and a second run over
the same path takes a different branch.

Not a determinism issue and not a porting question, but it is why `bench/check_seed_sensitivity.py`
copies its input before every invocation, and the port must not inherit the behaviour silently.

### Finding 6 — argparse accepts unambiguous flag prefixes, and pipelines may rely on it

Not a defect, but a parity trap that only surfaced because the harness ran the two implementations
against each other. `argparse.ArgumentParser` defaults to `allow_abbrev=True`, so **any unambiguous
prefix of a long option works**. Measured against the reference:

| invocation | resolves to |
| --- | --- |
| `main --comp` | `--compression` |
| `main --disable_num` | `--disable_numpy` |
| `main --set_w` | `--set_w_dynamically` |
| `main --fast x.fq` | `--fastq x.fq` |
| `isONform_parallel --fastq d` | **`--fastq_folder d`** |

The last one is the interesting one: `--fastq` is not a flag on the driver at all, yet
`isONform_parallel --fastq d` works, because it is an unambiguous prefix of `--fastq_folder`. Any
pipeline script written that way is silently relying on abbreviation.

clap does not infer prefixes by default, so the port would have rejected all five. It now sets
`infer_long_args = true` on both entry points, which matches argparse's rule including its
preference for an exact match — that is what keeps the driver's `--k` from colliding with
`--keep_old`. Ambiguous prefixes stay rejected on both sides: `main --delta` prefixes three flags
and `main --m` prefixes two.

### Finding 7 — smaller things, recorded and not acted on

- `main`'s window check is `if 100 < args.w or args.w < args.k` but its message reads "smaller than
  100"; `--w 100` is accepted. Off-by-one between check and message.
- `--xmin`/`--xmax` help text is swapped in `main` (see *Reconnaissance corrections* 5).
- `isONform_parallel:13` imports `from ast import Param`, which is unused and is an IDE
  autocomplete artifact.
- `isONform_parallel:79` spawns `["python", ...]` rather than `[sys.executable, ...]`, so the child
  is whatever `python` resolves to on `PATH` — not necessarily the interpreter running the parent,
  and absent entirely on systems that only ship `python3`.
- `--xmin 200 --xmax 80` (xmin > xmax) is not validated and crashes with a traceback rather than a
  message.
- A nonexistent `--fastq` crashes with a `FileNotFoundError` traceback (exit 1) rather than a
  diagnostic.

## Evaluating isONform

Equivalence testing answers *do the two implementations agree*. It cannot answer *is the output any
good* — two implementations can agree perfectly on a wrong answer — and that second question is the
one that matters whenever the reference itself changes: a determinism fix, a deliberate divergence
in the port, a performance tradeoff. This is the machinery for it.

### The metric is not isONcorrect's

isONcorrect emits one corrected read per input read, so accuracy is per-read error against the
read's source transcript, and `bench/accuracy.py` there measures exactly that. **isONform emits
isoforms, and how many is an output rather than an input.** Per-read error has no meaning; the
metric has to be

| | |
| --- | --- |
| **recall** | fraction of reference transcripts recovered by at least one reconstructed isoform |
| **precision** | fraction of reconstructed isoforms that match some reference transcript |
| **redundancy** | reconstructed isoforms per matched transcript; 1.0 is ideal |
| **identity** | per-isoform identity against the matched transcript |
| **length ratio** | reconstructed length over transcript length, to expose systematic truncation |

**All of them, together, always.** Any one alone is trivially gamed: emit a single confident isoform
and precision and identity look perfect while recall collapses; emit hundreds and the reverse. A
change that improves one and is silent about the others has not been measured.

A match requires identity above a threshold **and** a length within tolerance. The length condition
is not decoration: edlib infix alignment will happily place a 300 bp fragment inside a 2 000 bp
transcript at 99% identity, and counting that as a recovered transcript would make truncation look
like success. Both a strict (0.95, 10%) and a lenient (0.90, 20%) threshold are reported, because a
single cutoff cannot distinguish "sequences moved" from "sequences moved across a line".

### Three corpora, because one lies

Method point 4 is in this document because simulated and spike-in data gave the wrong answer twice
during the isONcorrect port. So the corpora are chosen to fail in *different* directions, and a
result is believed when all three move the same way.

| corpus | reads | clusters | truth | what it is for |
| --- | --- | --- | --- | --- |
| `sirv_sim_gene` | 10 000 simulated, 7% error | 7, 713–2242 reads | exact, in every read header | recall with a real denominator |
| `sirv_real` | 9 968 real ONT SIRV | 26, 12–4595 reads | known transcriptome | real error profile, known answer |
| `droso` | 13 164 real ONT Drosophila | 561, 5–488 reads | genome only | transcriptome scale, real biology |

`sirv_sim_gene` groups reads by **gene**, so isoforms of one gene share a cluster and reads differ
by whole exons — the shape isONform exists to handle — and the largest clusters cross `--max_seqs`
so batching runs for real. Its weakness is the important one: uniform i.i.d. error, no chimeras, no
truncation, which is exactly the profile that misled the isONcorrect port. **Never conclude from it
alone.**

`sirv_real` is the primary accuracy corpus, because SIRV is the one setting where a real ONT error
profile and a trustworthy ground-truth transcriptome coexist. It has no per-read truth, so recall is
computed against all 68 reference transcripts and is a **lower bound**; the harness says so on every
run.

`droso` is the only corpus at transcriptome scale, and has no trustworthy per-isoform truth at all.
It is scored against the **genome** by spliced alignment instead — see below.

Every corpus goes through isONcorrect first. That is not optional preprocessing but the input
contract: run on raw isONclust output, isONform's interval-abundance filter skips essentially every
read and the graph comes out with 2 nodes and 0 edges.

### Scoring against a genome, with no annotation

`bench/accuracy_isoforms_genome.py` aligns each reconstructed isoform with `minimap2 -ax splice` and
reports:

- **error rate** (`NM` / aligned bases) **and aligned fraction, never one alone.** Error rate by
  itself is gamed by aligning less — an isoform whose ends were mangled would soft-clip them away
  and *improve* its error rate. The pair cannot move that way.
- **canonical splice fraction** — of the introns each isoform implies, how many carry `GT..AG` (or
  `CT..AC`). The strongest annotation-free quality signal there is: a correctly reconstructed
  junction lands on a real splice site, an invented one generally does not. It is also **the only
  metric here that emitting fewer isoforms cannot improve**, because it is a rate over junctions
  actually claimed.
- **distinct intron chains** against multi-exon isoform count — the same splice structure reported
  twice is redundancy, measured without any annotation.
- **multi-exon fraction**, because a tool that collapsed everything to single-exon fragments would
  score beautifully on the first two and be useless.

`NM` over a genome includes real biological difference — SNPs against the assembly, sites minimap2
mis-places — so the absolute error rate is an **overestimate**. The inflation is identical across
compared sets, so comparisons hold where the absolute figure does not; do not quote it as "the error
rate". Unaligned isoforms are counted and reported rather than dropped, so a set that emits
sequences minimap2 cannot place does not thereby look clean.

### Scoring against a genome *with* an annotation: SQANTI categories

With `--annotation`, the same alignment is additionally classified into the SQANTI structural
categories, which is a far sharper instrument than canonical-splice fraction — it can distinguish a
correctly reconstructed known transcript from a merely plausible novel one, and it gives a real
recall number. `bench/annotation.py` holds the classifier and its self-tests.

Comparison is on the **splice junction chain** in genomic coordinates. End positions are deliberately
excluded from `FSM`: ONT reads are truncated, and 5'/3' position is not what isoform reconstruction is
being judged on.

| category | meaning | what it says about the tool |
| --- | --- | --- |
| `FSM` | junction chain identical to an annotated transcript's | got it exactly right |
| `ISM` | chain is a *contiguous* sub-chain of an annotated one | right as far as it goes; a fragment |
| `NIC` | every junction, or at least every splice site, is annotated, but the combination is not | real novel isoform, or a chimera of known parts |
| `NNC` | at least one unannotated splice site | the strongest error signal |
| `genic` / `genic_intron` / `antisense` / `fusion` / `intergenic` | the residue | mis-assembly or noise |

`intron_retention` is reported separately rather than as a category, because it cross-cuts `ISM` and
`NIC`. It fires when one of the isoform's exons swallows an annotated intron whole, which for a
graph-based assembler means **a bubble that did not get popped** — a specific, locatable failure.

Two things about this that took measurement rather than reading:

**The recall denominator.** "Fraction of the 34 989 annotated transcripts recovered" is a meaningless
number: 505 isoforms from 13 000 reads cannot cover them, so the percentage measures sequencing depth.
Reported instead are the absolute count of annotated transcripts with at least one `FSM`, and that
count over the transcripts of **only the genes the isoform set actually touched** — which asks "of
the loci we produced anything for, how much of their annotated structure did we recover".

**The fusion test is the easy thing to get wrong, and getting it wrong poisons `NIC`.** The obvious
implementation — a fusion is an isoform whose junctions fall inside more than one gene's span —
reported **13% fusions** on real Drosophila data and left `NIC` at exactly zero, because FlyBase r6.68
contains **1 683 overlapping same-strand gene pairs and 3 305 fully nested ones** (measured), so any
intron spanning a nested gene looks like a fusion. The corrected test requires each contributing gene
to own at least one of the isoform's junctions *as an annotated junction*, and requires the
contributing genes to occupy disjoint spans. Fusions dropped from 13% to 0% and `FSM` rose from 77% to
88%. Both conditions have regression tests in `bench/annotation.py`, including a nested-gene fixture.

Building the annotation file is worth doing once. The FlyBase whole-genome GFF is 6.7 GB of which the
gene models are a rounding error — 19 M `match_part` and 9.3 M `match` records of alignment evidence
against 85 601 exons — so filter it first:

```bash
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=9 && ($3=="gene" || $3=="mRNA" || $3=="exon" \
    || $3=="ncRNA" || $3=="pseudogene" || $3=="pre_miRNA" || $3=="miRNA" || $3=="tRNA" \
    || $3=="snRNA" || $3=="snoRNA" || $3=="rRNA") {print $1,$2,$3,$4,$5,$7,$9}' \
    dmel-all-r6.68.gff > dmel-r6.68-transcripts.tsv
```

That is 139 209 records, 18 MB, and `bench/annotation.py` reads either the 9-column GFF or this
7-column projection.

### What the corpora say about isONform as it stands

Worth recording separately from the determinism fix, because it is the baseline any port has to match
and no number like it existed before. Lexicographic reference, default flags, `--iso_abundance 5`:

| corpus | isoforms | recall | precision | redundancy | identity |
| --- | --- | --- | --- | --- | --- |
| `sirv_sim_gene` | 138 | 89.7% (61/68) | 97.8% | 2.21 | 0.9978 |
| `sirv_real` | 106 | 70.6% (48/68) † | 82.1% | 1.81 | 0.9673 |

† lower bound: real reads carry no source transcript, so the denominator is all 68 references.

`droso`, 505 isoforms against the FlyBase r6.68 annotation:

| | | |
| --- | --- | --- |
| `FSM` | 443 | **88%** |
| `ISM` | 13 | 3% |
| `NIC` | **0** | 0% |
| `NNC` | 12 | 2% |
| `genic` (mostly mono-exonic) | 29 | 6% |
| `genic_intron` / `antisense` / `intergenic` | 7 | 1% |
| `fusion` | 0 | 0% |
| | | |
| annotated transcripts recovered exactly | 409 | |
| of the 966 transcripts in the 427 genes touched | | **42.3%** |
| intron retention | 8 | |
| median error rate / aligned fraction | 0.0134 / 0.978 | |
| canonical introns | 98.2% of 779 | |

Five things to take from this, and the last two are the ones a rewrite should act on:

1. **Simulated data flatters the tool by a lot** — 19 points of recall and 15 of precision against
   real SIRV reads. Any accuracy claim resting on simulated reads alone is not a claim about isONform.
2. **The splice structures are good.** 88% of Drosophila isoforms match an annotated junction chain
   exactly, and 98.2% of implied introns are canonical. Whatever else isONform gets wrong, it is not
   inventing junctions wholesale, and it is not producing fusions.
3. **Redundancy is the weakest number**: 2.21 reconstructed isoforms per recovered transcript on the
   simulated corpus, 1.81 on real SIRV. Roughly half the output duplicates something else in it. On
   Drosophila the same shows up as 359 distinct intron chains from 385 multi-exon isoforms. This is a
   change to what the tool *does* rather than how fast it does it, and it now has a measurement to be
   judged against.
4. **`NIC` is exactly zero, and that is real rather than a classifier bug** — checked directly, after
   one fusion bug had already faked it (above). Not one multi-exon isoform is a novel *combination* of
   annotated junctions. So isONform is not mis-stitching known exons into plausible-but-wrong
   transcripts, which is the failure mode one would most fear from a graph assembler.
5. **Every `NNC` isoform is wrong in exactly one junction.** 12 of 12 for the lexicographic reference,
   11 of 11 for the pre-fix code: each has one unannotated junction and the rest annotated. The
   structural errors are single mis-resolved junctions, not wholesale mis-assembly. That is a narrow,
   locatable target — one wrong edge or one bubble resolved the wrong way — and it is where the graph
   code should be looked at first.

### Running it

```bash
bench/evaluate.sh corpora                            # build all three, once
git archive <commit> | tar -x -C /tmp/before         # the "before" tree
bench/evaluate.sh compare "$PWD" /tmp/before         # after vs before, scored
```

`compare` runs the "before" implementation at three hash seeds rather than one, because with a
seeded hash a single run of it says nothing about what the tool does — the spread *is* the finding.
Extracting the old tree with `git archive` rather than checking it out means the working tree is
never disturbed and both implementations are live at once.

## Method

Carried over from the isONcorrect port. The full version, with the measurements behind each point, is
in that repository's `PORTING.md`; this is the short form.

1. **CLI parity first**, locked by unit tests. Argument names, defaults, validation order, stderr
   text and exit codes. Watch that clap rewrites `field_name` to `--field-name` — every multi-word
   flag needs an explicit `long = "..."`.
2. **Differential oracles, not end-to-end tests.** Wrap the reference without modifying it; dump each
   stage's inputs *and* outputs in a stable line format; replay pure functions from Rust and diff
   stateful ones. End-to-end equivalence tells you *that* something is wrong, never *where*.
3. **Dump from the live driver too, not only from a standalone dump binary.** Oracles replay recorded
   reference inputs, so a port following a different trajectory still passes them. The one real bug in
   the isONcorrect port passed every oracle and was caught only by diffing a dump taken from the
   running driver.
4. **Build a real corpus before trusting anything.** Simulated and spike-in data gave the wrong
   answer twice in isONcorrect — endorsing an alignment band that real data broke, and reversing the
   sign of an accuracy effect. Get real ONT data through the upstream pipeline early.
5. **Profile before optimising, and re-profile after.** Bottlenecks move. Instrument at stage
   granularity, and *remove* sub-stage instrumentation after reading it: on millions of calls the
   timers themselves distorted the table.
6. **Measure, do not reason, about performance.** Recorded null results from isONcorrect: caching an
   edit-distance pattern was slower, reusing the POA engine was worth nothing, 4-bit DP cells were
   slower than 8-bit, hoisting a hash lookup out of a loop was worth nothing.
7. **Set up CI on day one, on Linux *and* macOS, x86_64 *and* arm64.** Three defects in isONcorrect
   existed only because every local check ran on one machine: an aarch64-only build that did not
   compile on x86_64 at all, a C dependency hidden behind a default cargo feature, and
   glibc-2.34-only release binaries that exclude CentOS 7, RHEL 8 and Ubuntu 20.04.
8. **Fix reference bugs upstream once measured, rather than reproducing them.** The largest accuracy
   win in the isONcorrect port was a one-character fix to the *Python*, worth more than both
   deliberate divergences combined. Each such fix is its own commit, with goldens re-recorded.

## The reference is expected to change, and that changes the method

Stated by the owner on 2026-08-23, after seeing the determinism finding: **isONform is unstable and
probably buggy, and it would not be surprising if functions or methods get changed along the way in
this rewrite.** That is a design input, not a caveat, and it shifts the emphasis of everything below.

The consequence is a reordering of what counts as evidence:

- **Byte-identity stays a tool, and drops further from being a goal.** It is how a *port* of an
  unchanged function is proven faithful. It says nothing about whether the function should exist in
  that form, and against a function that is about to be replaced it is worth nothing at all.
- **The accuracy harness becomes the primary instrument, not the secondary one.** When a method
  changes, `bench/equivalence.sh` has nothing to say by construction — the outputs *should* differ.
  What decides whether the change is good is `bench/evaluate.sh` over the three corpora, and the
  isONform-specific point that recall, precision, redundancy and identity must be read together.
  This is why the corpora and the two `accuracy_*` scripts were built before any graph code.
- **Invariants matter more than goldens.** A golden recorded against a function that is about to be
  rewritten becomes a liability — it will fail, and the failure will be uninformative. An invariant
  (`bench/equivalence.sh invariants`) survives the rewrite, because it constrains what any correct
  implementation must do rather than what this one happens to do. Prefer adding invariants over
  adding goldens for anything in the graph stages.
- **Every intentional divergence gets its own commit and an entry below, with a measurement.** The
  determinism fix is the template: defect measured, fix isolated, accuracy compared against the old
  behaviour on three corpora, and the null result (runtime) recorded as a null result.

The redundancy number is the first place this is likely to bite. `sirv_sim_gene` reports 2.21
reconstructed isoforms per recovered transcript and `sirv_real` 1.81 — roughly half the output
duplicates something else in it. Fixing that is a change to what the tool *does*, not to how fast it
does it, and it now has a measurement to be judged against.

## Working agreements

- **Never run `git commit`, `git push`, or any history-changing git command.** Leave finished work in
  the tree, say what changed and why, propose the commit message, and let a human commit it. This
  held throughout the isONcorrect port and holds here — doubly so, since this repository belongs to
  someone else.
- Don't "improve" the algorithm while porting. A behaviour change and the port must not land in the
  same commit; an intentional divergence needs its own commit and a note here.
- **"Behaviour" means observable output, not internal representation.** Different containers, dropping
  provably-dead entries, arena allocation — none of that is a behaviour change if the emitted bytes
  are identical. What is not free: iteration order where it reaches results, tie-breaking, arithmetic
  and rounding.
- When you spot a possible improvement, write it into a *Deferred improvements* section and move on.

## First steps, in order

1. ~~Get it running: install the Python reference in a pinned env (`networkx==2.8.4`,
   `parasail==1.2`, the `spoa` binary), and run it on isONcorrect's output for a small cluster.
   Record what it emits.~~ **Done.** `bench/setup_reference_env.sh` builds `isonform-ref`;
   `bench/corpus/sirv_small` is the input; both entry points run. Note `parasail==1.2` was not
   usable — see *Reconnaissance corrections* 3 for why 1.3.4 was pinned instead. Note also that the
   reference does nothing useful on *uncorrected* reads: on raw isONclust output all 100 reads are
   skipped for low interval abundance and the graph comes out 2 nodes, 0 edges. The corpus has to
   come through isONcorrect first.
2. ~~Check `--disable_numpy` and `--compression` actually work, and check for `PYTHONHASHSEED`
   sensitivity, before recording any goldens.~~ **Done, and it stopped the line.**
   `--disable_numpy` is inert, `--compression` is live, and the tool is non-deterministic run to
   run — findings 3 and 1. `bench/check_seed_sensitivity.py` is the harness; it exited 1 when
   written, and exits 0 now that step 3 has landed.
3. ~~**Make minimizer selection deterministic in the Python.**~~ **Done**, lexicographically, as
   isONcorrect does it — the owner's call. Its own commit, ahead of the port, with the accuracy cost
   measured on three corpora rather than assumed. See *The determinism fix, measured*.
4. Write `bench/equivalence.sh` (**done** — the `cli` layer passes 37 cases; the `stages` layer's
   case matrix is defined and unblocked) and record goldens for both entry points across the flags that
   change behaviour — `--compression`, `--slow`, `--set_w_dynamically`, `--k`/`--w`/`--xmin`/
   `--xmax`, `--delta_len`, `--delta_iso_len_3`/`_5`, `--max_seqs`, `--max_seqs_to_spoa`, and for
   the driver `--split_wrt_batches`, `--iso_abundance`, `--write_fastq`, `--t`, `--keep_old`. Skip
   the four inert flags (finding 3) beyond one test each proving they stay inert. `--record` checks
   the reference agrees with itself across two seeds before writing anything.
5. Port the CLI (**in progress** — `rust/src/cli.rs`, tested against `bench/golden/cli/`), then
   fastq I/O, then work outward from the leaves — the graph construction before the simplification,
   the simplification before isoform generation. The front half of `main` comes across from
   isONcorrect rather than being written (*Reconnaissance corrections* 4).
6. CI on Linux and macOS, x86_64 and arm64, on day one — method point 7, and the reason three
   isONcorrect defects existed at all.

## Deferred improvements

Spotted while porting, deliberately not acted on. Each needs its own commit, and the ones that
change output need a measurement first.

- **`petgraph` is not needed** (already argued above): integer node indices, CSR adjacency in flat
  `Vec`s, attributes in parallel arrays. Worth proposing to the Python side too — networkx's
  per-node object plus per-node and per-edge attribute dicts is a large part of why the reference is
  slow, and the representation change is measurable before the port exists.
- **Don't port the dead branches.** `previously_corrected_regions` is never assigned, so `main:438–
  461` and `main:467–476` are unreachable (finding 3). Not porting them is a representation choice,
  not a behaviour change — nothing observable moves.
- **`main` hard-codes `delta = 0.15`** at `main:560` while `isONform_parallel --delta` defaults to
  `0.1` and only reaches batch merging. Either `main` should take `--delta` and the driver should
  forward it, or the driver's flag should say it applies to merging only. Behaviour change.
- **Flags that cannot reach the work they name**: `isONform_parallel --max_seqs_to_spoa` and
  `--delta` never reach the per-cluster child (*Reconnaissance corrections* 5). Forwarding them is a
  behaviour change; documenting the limit is not.
- **`--exact` is passed on every child invocation and does nothing**, and the parallel driver's
  `--exact_instance_limit` default of 50 therefore also does nothing. Once `--exact` means something
  again — if it should — these defaults need revisiting.
- **Stray debug prints** (finding 4) and the unused `from ast import Param` (finding 7): one-line
  deletions, but the prints are on stdout so goldens move with them.
- **`["python", ...]`** → `[sys.executable, ...]` in `isONform_parallel:79`. Pure bug fix; in the
  Rust port the equivalent is to re-exec the running binary by its own path, not by name on `PATH`.
- **Missing validation**: `--xmin > --xmax` crashes with a traceback, a nonexistent `--fastq` crashes
  with `FileNotFoundError`, and the window-size message says "smaller than 100" where the check
  admits 100. All three want diagnostics.
- **`restructure_isoncorrect_output` mutating its input** (finding 5). At minimum it should be
  opt-in, or write to a new directory.
- ~~**A second corpus from genuinely distinct clusters.**~~ **Done** — three of them, two from real
  ONT data through the upstream pipeline. See *Evaluating isONform* and `bench/corpus/README.md`.
  `sirv_small` stays as the smoke test and for the identical-cluster invariant.
- **Internal low-complexity is untested.** The poly-A measurement behind the minimizer fix was taken
  on SIRV, which has none to speak of. A corpus with genuine internal low-complexity would test the
  one place lexicographic selection could plausibly be worse than hashing.
- **`--slow`, and the flags no corpus exercises.** `bench/evaluate.sh` runs default settings only.
  `--slow`, `--compression` and the `--delta_iso_len_*` cutoffs change output and have no accuracy
  measurement behind them at all.
