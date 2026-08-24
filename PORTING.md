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

- **Graph construction ported and verified against the reference.** `rust/src/graph.rs` and
  `rust/src/graph_build.rs` port `generateGraphfromIntervals`; `bench/dump_reference.py` records the
  stage's inputs and outputs from the live driver and `rust/tests/graph_oracle.rs` replays them and
  diffs node set, `end_mini_seq`, read maps, edges, lengths, edge support and the **exact
  `nx.topological_sort` order**. Passes on `bench/corpus/sirv_small` over 324 reference nodes and
  488 edges, and it is a gate in CI. Two more reference defects came out of it, findings 8 and 9.
- **Bubble detection ported** (`rust/src/simplify.rs`) — the pure, order-sensitive front of
  `SimplifyGraph.py`: candidate starts and ends, the start × end combination enumeration, and the
  two distance helpers. 69 unit tests across the crate.
- **The oracle passes on real data**: 72 recorded graph-construction calls across `sirv_small`,
  `sirv_real` and 40 `droso` clusters — **47 963 reference nodes and 59 320 edges** — matching node
  for node, edge for edge, attribute for attribute. It found one genuine reference defect on the way,
  finding 11, which no amount of reading had turned up.
- **The simplification stage is recorded too** — `bench/dump_reference.py --stage simplify` writes the
  graph before and after each `simplifyGraph` call, since that stage mutates in place and returns
  nothing. 16 records over medium `droso` clusters, with real popping (e.g. 165 edges down to 126).
  A useful invariant fell out: **simplification never changes the node set**, only edges — true in
  all 16.
- **A second reference fix, and this one bought accuracy**: `find_paths` allocated reads to bubble
  paths in CPython `set.pop()` order, which decided which path survives a bubble. Replacing it with a
  defined order raised `sirv_real` recall from 70.6% to 72.1% and cut Drosophila `NNC` isoforms from
  12 to 10. Finding 12.
- **Bubble popping ported, including the driver** — `find_paths`, `remove_edges`,
  `prepare_adding_edges`, `linearize_bubble` and `pop_bubbles` (the
  `new_bubble_popping_routine` state machine), in `rust/src/simplify.rs`. The reference's destructive
  `pop(0)`, its dead branches, its stale-variable hazards, its operator-precedence bug, its
  undercounted pop total and the fact that relinked edges carry **no `length` attribute** are all
  reproduced and pinned by tests. 101 unit tests across the crate.
- **The poppability decision ported too** — `align_bubble_nodes`,
  `generate_consensus_path`, `get_consensus_positions` and `parse_cigar_diversity`, with the
  megabubble memo. 114 unit tests. spoa and parasail sit behind a `Consensus` trait with two
  methods, so everything around them is tested now and wiring them in is mechanical: the spoa
  invocation is the one isONcorrect's `poa.rs` already reproduces, and the alignment is the
  `sg_trace_scan_16` its `parasail.rs` already reproduces.
- **The two-supporting-read path skips spoa entirely** — it takes whichever read's span is longer,
  verbatim. That is where a simplification oracle can start without reproducing spoa at all, and it
  is fully ported and tested.

- **`parasail.rs`, `poa.rs` and `align.rs` carried across from isONcorrect**, so `SpoaParasail`
  provides the real consensus and alignment. 147 unit tests. What came across, and what was left
  behind, is recorded below.

- **The simplification oracle exists** (`rust/tests/simplify_oracle.rs`). It rebuilds the graph the
  reference had on entry, runs the port's `simplify_graph`, and diffs the result. Passing on
  `sirv_small`’s 2 cases (CI), 16 real Drosophila graphs that pop tens of edges each, a
  56-cluster size-stratified Drosophila sample, and a second **disjoint** 56-cluster holdout. The
  holdout is what kept finding bugs the first corpus could not — see below.

- **spoa equivalence is verified on isONform's own inputs. Caveat closed.** 5 368 of 5 368 recorded
  calls across four corpora give a consensus identical to the `spoa` binary, the large majority from
  the bubble-path call site the simplification oracle depends on. That removed the oracle's
  "disagreed but called spoa, so not attributable" escape hatch, so every disagreement now fails.

- **Finding 24, found by that tightening and fixed:** the port's parasail scored a non-`ACGT`
  character as a match against itself, where the real library scores it **0**. All-`X` placeholder
  consensuses are an ordinary occurrence, so the port was popping bubbles the reference refuses. One
  fix closed both cases the tightening had exposed — which had looked like two unrelated bugs.

- **Finding 14 reclassified from latent to live, and the port now diverges deliberately.** The
  holdout's one-edge failure turned out to be a `PYTHONHASHSEED`-dependent pick in
  `prepare_adding_edges` — the *reference* gives two different graphs across eight seeds. The port
  orders the candidates lexicographically, matching finding 1's precedent. This also means finding 1's
  "fixed" was a claim about one cause, not about the tool.

- **parasail is now verified on isONform's inputs too, and that found finding 25.**
  `--record-parasail` writes the cases; **56 549** unique calls across two corpora, **0 score and 0
  CIGAR mismatches**. It started at 12 score and 136 CIGAR errors: the port's semi-global scan
  admitted an "align nothing" end cell parasail excludes, and resolved the remaining end-cell ties by
  a rule parasail does not use. My earlier note that "the parasail side needs no such caveat" was
  wrong on both counts, and the CIGAR half is what `parse_cigar_diversity` actually reads.

- **The front half of `main` is ported and verified.** `minimizers`, `anchors`, `wis`, `intervals`
  and `reads` cover everything from a raw read to the intervals `generateGraphfromIntervals`
  consumes. A fourth oracle (`rust/tests/intervals_oracle.rs`) drives `main` itself and diffs the
  chosen intervals, the `graph_id` assignment and each interval's full instance list: **121 cases,
  3 682 reads, 0 disagreements** across `sirv_small`, both disjoint Drosophila corpora and real SIRV.
  This is the first stage whose input is a read rather than a recorded intermediate, and it closes
  the gap where a wrong minimizer would have gone unnoticed — the graph oracle replays the
  *reference's* intervals, so nothing upstream of them was covered.

- **Isoform generation is ported.** `isoforms` covers the live path of `IsoformGeneration.py` —
  grouping reads by their route through the simplified graph, then merging groups whose consensuses
  are similar. A sixth oracle diffs both halves: **0 wrong partitions and 0 merging failures on all
  114 real cases**. One deliberate divergence, decided and measured: the port orders group members
  ascending where the reference uses CPython set-iteration order, which names every isoform and
  orders every spoa call (finding 28). It fires on 28 of 114 cases and leaves 93% of isoform
  sequences byte-identical there.

- **Batch merging is ported, and it turned out to be a no-op.** `batch_merge` covers the last
  stage. Its merging step has never executed — see finding 31 — so the port reproduces that, and
  what it does port is `write_final_output`'s destination, id and support decisions. A seventh
  oracle, driven from `isONform_parallel` because `main` never reaches this stage: **114/114 cases,
  2 102 output records, 0 disagreements**.

**Every stage is ported, and every recorded corpus agrees except finding 28's ordering.** Intervals
121/121, graph construction on all five corpora, simplification 56/56 on the disjoint holdout and on
the first Drosophila sample, 16/16 on `dsim_mid`, 2/2 on `sirv_small`; isoform partitioning and
merging 114/114; batch merging 114/114. What is left is not another stage but the two things that
need the stages joined up: wiring the driver so the port runs end to end, and the accuracy scoring
`bench/evaluate.sh` was built for.

### The simplification port did not converge, and the oracle is what found it

The oracle's first run against real bubble-popping graphs: **16 of 16 disagreed, and every one hit
the port's iteration cap.** The reference terminates on these graphs in 5 to 9 iterations; the port
ran until the cap stopped it.

That was attributable with no further work, and it is worth being precise about why. The cap fired
because `this_it_pops >= pop_threshold` every iteration — the port kept finding poppable bubbles
where the reference ran out. **No consensus difference can explain that**: spoa can only change
*which* bubbles are judged poppable, never turn a terminating loop into a non-terminating one. So
this was a port bug, not the unverified spoa half, and the oracle classifies a cap hit that way
regardless of how many times spoa was called (its first version filed these under "not attributable",
which was wrong and would have hidden them).

The second symptom pointed the same way: the port's nodes ended up with **far fewer reads** than the
reference's — 6 against 53 on one node, 2 against 67 on another. Reads are added to nodes by
`additional_node_support` during a pop, so fewer reads with more pops is contradictory unless
something is dropping them.

Ruled out along the way, before finding the actual cause:

- **spoa was not returning nothing.** `poa::consensus` on five noisy 116-base variants returns a
  116-base consensus, and parasail on two near-identical sequences returns a sensible single-run
  CIGAR. An empty consensus would have made every bubble look poppable (two short sequences pop
  unconditionally), which was the obvious candidate.
- **The rebuild was correct.** All 16 cases reproduce the reference's exact `nx.topological_sort`
  order before the stage runs, which the oracle asserts as a precondition. A full node/edge-content
  diff at the end of the first iteration also matched byte for byte between the two, which is what
  ruled out graph *construction* and pointed at something order- or state-dependent inside popping
  itself.
- **One real divergence found and fixed on the way**, though it was not the cause: the multi-path
  branch cloned the two paths before linearising them, where the reference passes the same tuples
  that stay in `all_paths` — so `remove_edges`' `pop(0)` is visible to later pairs in the same
  iteration there and was not here. Fixed by mutating in place.

**The actual cause**, found by instrumenting both the port and a directly-replayed reference run
(loading the same recorded "before" graph into a real `networkx.DiGraph` and tracing which edge each
side walked) side by side, iteration by iteration, until a single bubble's alignment decision
diverged: `remove_edges` (`SimplifyGraph.py:279-304`) records an edge for deletion whenever it
*exists*, unconditionally. The port's `remove_edges` instead only lifted an edge when
`Graph::edge_length` returned `Some` — but Finding 16 is exactly the reason that is the wrong test:
an edge `prepare_adding_edges` had re-added on an *earlier* pop carries no `length` at all
(`upsert_edge_support` creates it without one), and the port's check silently treated "exists with no
length" as "does not exist" and skipped removing it. The edge survived every subsequent pop attempt
that touched its endpoints, so the bubble it belonged to kept reappearing as poppable — indefinitely,
since the graph a later iteration inspects is corrupted in exactly the way needed to look poppable
again. The fix: gate lifting on `Graph::has_edge` (existence) and carry the length through as the
`Option<i64>` it already is; `LiftedEdge::length` is now `Option<i64>` rather than `i64` to make that
representable. Regression test:
`a_lengthless_edge_inside_the_bubble_is_still_lifted_and_removed`.

One existing test's premise turned out to rest on the same bug: `an_aligner_that_never_refuses`
documented "the reference's loop has no termination guarantee" using a synthetic three-path bubble
that, with an always-yes aligner, hit the cap before the fix — but after it, the same fixture
converges in three iterations (two pops, six edges down to four), because the re-added edges from the
first pop are now correctly removed on the second. Renamed to
`an_aligner_that_never_refuses_still_converges_by_merging_paths_away` and rewritten to assert the
(correct) convergence. The structural point it meant to document — `while has_combinations` has no
formal termination proof, so the iteration cap stays as a safety net — is still true in principle;
this fixture just was not evidence of it, and no fixture we know of is now.

Verified after the fix: all 16 real Drosophila graphs agree (spoa called thousands of times across
them, e.g. 8 825 times on one), none hit the iteration cap, and every case's rebuild reproduces the
reference's topological order exactly. A further sample of 56 clusters spanning the corpus's size
distribution (1.2 KB to 174 KB of reads, stratified by file size, largest 5 of 561 excluded to bound
runtime) agrees on 54; the other 2 disagree while calling spoa 75 and 178 times respectively, so they
land in the pending-spoa-verification bucket the oracle already had rather than the cap-hit bucket
this fix targets — none of the 56 hit the iteration cap either. All 148 unit tests pass, both oracles
pass (simplify's build gate is "no attributable disagreement", which holds), and
`cargo clippy --all-targets -- -D warnings` is clean.

Bubble popping is now ported and verified against real data. Everything downstream of it (isoform
generation, batch merging) is not started, and the front half of `main` is still owed — which is why
the binaries stop at argument handling and every oracle replays recorded inputs rather than running
the port end to end. That last point is the weakest link in the current evidence, not a formality:
in isONcorrect the one genuine bug that survived every replayed-input oracle was found only by
diffing a dump taken from the *running* driver, because a port that follows a different trajectory
still satisfies replayed inputs.

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

### 4. Much of the front half of `main` is isONcorrect's — but "function for function" was too strong.

The recon framed the reuse as "`cli.rs`, `params.rs` and `validate.rs` are a starting point". It is
larger than that: `main` lines 81–340 are the interval machinery, and most of it does transfer. This
entry originally said the fourteen functions matched **"function for function"**, which was written
from a reading rather than a diff. Diffed, normalising whitespace and comments:

| | functions |
| --- | --- |
| **identical** | `solve_WIS` |
| **cosmetic only** (whitespace, comments, `itertools.` qualification) | `rindex`, `get_intervals_to_correct`, `minimizers_comb_iterator`, `grouper`, `batch`, `get_kmer_minimizers` |
| **substantively different** | `fill_p2`, `get_minimizer_combinations_database`, `get_minimizers_and_positions`, `find_most_supported_span` |
| **isONform only** | `get_kmer_maximizers`, `get_minimizers_and_positions_compressed`, `remove_read_polyA_ends` |

Seven of fourteen transfer essentially verbatim. Four differ in ways that reach the output, and one of
those is not a variation at all:

* **`find_most_supported_span` shares nothing but the name.** isONcorrect's compares the actual
  subsequences, tracks quality values and edit distances, and memoises across reads — sixty lines.
  isONform's keeps a supporting read whenever the two spans differ in *length* by less than
  `delta_len`, and never looks at a base — fifteen. Copying isONcorrect's would have been a
  different algorithm.
* **`fill_p2` is the `solve_WIS` suboptimality, still unfixed here.** See finding 26.
* `get_minimizer_combinations_database` uses a different abundance filter (`> 3 * len(reads)`
  unconditionally, against isONcorrect's poly-A-and-10x rule) and a different argument order.
* `get_minimizers_and_positions` dispatches on `"lex"` where isONcorrect dispatches on `"random"`,
  because the determinism fix of finding 1 moved which branch is the default.

So `minimizers.rs` and `wis.rs` carry across almost intact, `anchors.rs` carries with its filter
replaced, and the interval builder is new. What isONform does differently is still mostly *after* the
intervals are chosen — **the ~2 400 lines of genuinely new algorithm are `GraphGeneration` →
`SimplifyGraph` → `IsoformGeneration` → `batch_merging_parallel`** — but the front half is not the
free ride this entry claimed.

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

### Finding 1 — isONform was non-deterministic run to run. Every output file, every run. **Minimizer selection fixed; see finding 14 for a second source.**

**This blocked recording goldens and blocked the port, so it came first.** The fix is described in
*The determinism fix, measured* below; what follows is the defect as found.

> **Scope correction.** This heading said simply "**Fixed**". It is fixed for the cause described
> here — `hash()` on k-mer strings during minimizer selection, which was the dominant one and the one
> that made *every* output file vary. It is **not** the only source: finding 14 documents an
> independent seed dependency in bubble linearisation that survives this fix and still produces two
> distinct outputs across eight seeds on real data. It is far rarer (one occurrence in 19 831 calls),
> which is why it was not visible in the 24-seed check below — that corpus never reaches it. "Fixed"
> was a claim about a cause, written as though it were a claim about the tool.

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

### Finding 8 — `known_old_node_action` cannot run: it crashes networkx

`modules/GraphGeneration.py:159`:

```python
if not DG.has_edge(previous_node, name):
    DG.add_edge(previous_node, name, this_len)
```

`DiGraph.add_edge` takes `(u, v, **attr)`. A third *positional* argument raises
`TypeError: DiGraph.add_edge() takes 3 positional arguments but 4 were given` — verified directly
against networkx 2.8.4. Every other `add_edge` in the file passes `length=...` by keyword; this one
does not, so it would also have left the edge without the `length` attribute that downstream code
reads.

**It never fires, and it takes some work to see why.** The branch is guarded by `alternatives_filtered`
being non-empty, and `alternatives_filtered` comes from `alternative_nodes[old_node]` — a list that
the *first* visit to `old_node` initialises to `[]` and never appends to (`GraphGeneration.py:427`).
Only `no_node_to_add_to_action`, reached on the second visit, appends. So the crash needs a **third**
visit to the same `old_node`, with a matching predecessor and a length within `delta_len`.

Measured by instrumenting a copy of the reference and running all three corpora:

| | |
| --- | --- |
| times `alternatives_filtered` was computed | **4 541** (2 857 `sirv_real`, 1 684 `droso`) |
| times it was non-empty | **0** |

So it is a latent crash with a narrow trigger, not an active one. The port returns
`BuildError::ReferenceWouldCrash` on that path rather than inventing behaviour the reference does not
have; if a real dataset ever reaches it, the port says so instead of silently differing.

### Finding 9 — a node attribute is computed from the wrong read's sequence

`GraphGeneration.py:431` reads `seq` on a path that never binds it:

```python
if not (old_node in alternative_nodes):
    ...
    end_mini_seq = seq[inter[1]:inter[1] + k]      # `seq` is not assigned on this path
```

`seq` is function-scoped in Python, so it holds whatever the *previous loop iteration* left there.
Within one read that is harmless — same value. Across reads it is **the previous read's sequence**,
and the k-mer taken from it becomes the node's `end_mini_seq` attribute. On the first read to reach
this path before any other branch has run, it would be an outright `UnboundLocalError`.

Measured, by instrumenting the reference to compare the leaked `seq` against `all_reads[r_id][1]`:

| corpus | times the path ran | used a stale sequence |
| --- | --- | --- |
| `sirv_real` | 255 | **252 (98.8%)** |
| `droso` | 941 | **913 (97.0%)** |

`end_mini_seq` is consumed in exactly one place, `SimplifyGraph.new_distance_to_start`, as
`consensus.find(node_seq)` — so a k-mer from the wrong read is simply not found, returning `-1`, and
that bogus position feeds node distances during bubble linearisation.

**The decisive test was to fix it and re-measure**, rather than reason about the consequence:

| | reference | with `seq` bound correctly |
| --- | --- | --- |
| stale-sequence occurrences (`sirv_real`) | 252 | **0** |
| `new_distance_to_start` returning −1 (`sirv_real`) | 2 233 | 2 154 |
| `new_distance_to_start` returning −1 (`droso`) | 1 803 | 1 668 |
| `transcriptome.fasta` (`sirv_real`) | 106 isoforms | **byte-identical** |
| `transcriptome.fasta` (`droso`) | 505 isoforms | 505, **7 sequences differ** |
| every SQANTI category (`droso`) | FSM 443, ISM 13, NIC 0, NNC 12 | **identical** |

So: a real bug, and a small one. Only 3.5% of the `-1` results were caused by it — most are
legitimate. It changes 7 of 505 Drosophila isoform sequences and moves no structural category at
all, and on real SIRV it changes nothing observable. Worth fixing for correctness and to stop it
confusing the next reader, as its own commit; **not** worth describing as an accuracy win. The port
reproduces it by default and fixes it under `BuildOpts::fix_stale_seq`, so the comparison can be
re-run on any future corpus.

The fix is one line — bind `seq = all_reads[r_id][1]` on that path, as every sibling branch does.

### Finding 11 — the cycle-avoidance node is a dead end, and the read's path is severed

**Found by the oracle, not by reading**, which is the point of having one. On real Drosophila data
one edge in one of 40 recorded graphs disagreed: the reference drew it from `"104, 142, 133"` where
the port drew it from `"727, 765, 324"`.

When adding an edge would close a cycle, `GraphGeneration.py:391` calls `cycle_added`, which removes
the offending edge and creates a read-private node to route the *incoming* edge through. It then
rebinds its own local `name` to that new node — and `name` is a Python string, passed by value, so
**the caller never sees it**. `previous_node = name` at `:458` therefore continues the read's path
from the node whose incoming edge was just removed.

The result, visible in the dump:

```text
E 817, 841, 97 -> 727, 765, 324   39   324      # the new node: in-edge...
                                                 # ...and no out-edge at all
E 104, 142, 133 -> 1000, 1050, 2  80   324      # the path continues from the OLD node
N 104, 142, 133  ...  133:104:142:1,324:727:765:1
```

So read 324's path through the graph is in **two disconnected pieces**, and the node created to fix
the cycle contributes nothing but a node count. The reference's own comment admits discomfort here —
"TODO: add alt_cyc_node instead of adding a new node to reduce overall nr nodes in our graph" — but
not this.

Measured before deciding it mattered:

| | |
| --- | --- |
| graphs examined | 70 (`sirv_real` 30, `droso` 40) |
| nodes examined | 47 963 |
| dead-end nodes (in-degree > 0, out-degree 0, not the sink) | **1** |
| fully isolated nodes | 0 |
| `transcriptome.fasta` with the fix, `droso` | **byte-identical**, 505 isoforms |

So: a genuine structural defect that fires **once in ~48 000 nodes** and changes no output on the
corpora we have. Worth fixing, worth not overselling. The one-line fix is to `return name` from
`cycle_added` and assign it at the call site; the port reproduces the reference by default and fixes
it under `BuildOpts::fix_cycle_continuation`, with unit tests for both.

It is also the clearest argument so far for the oracle existing at all. Nothing about the emitted
isoforms would have revealed this, on any corpus.

#### One reference behaviour that looks like a defect and is not reachable

`convert_array_to_hash` drops the interval's own `(r_id, p1, p2)` triple before hashing, so an
interval supported by **no other read** hashes as the empty tuple — and any two such intervals in one
read collide, forcing the second to be judged `is_repetitive` and sending it down a different branch.
That cost a debugging round on a hand-built fixture.

It cannot happen through the real pipeline: `find_most_supported_span` only emits an interval when
`len(relevant_reads) // 3 >= 3`, so every interval reaching the graph carries at least three
occurrences and a non-empty tail. Measured across 188 454 real intervals: **zero** with an empty
tail. Recorded so the next person does not chase it.

### Finding 12 — CPython's `set.pop()` order was choosing how bubbles get linearised. **Fixed.**

`find_paths` (`SimplifyGraph.py:140`) allocated reads to bubble paths by
`node_support_left.pop()` — popping from a Python `set` of read ids. That order is deterministic for
a given set of ints (int hashing is not seed-randomised, which is why finding 1's fix left this in
place and the seed check still passes) but it is **arbitrary**: it is whichever slot CPython's set
table happens to hold first, given the table size and probing.

And it reaches results. Paths are appended to `all_paths` in pop order, and the popping routine then
treats them asymmetrically — `path1 = all_paths_filtered[0][0]`, `path2 = all_paths_filtered[1][0]`,
with different roles in `prepare_adding_edges`. So an implementation detail of CPython's set was
deciding which path survives a bubble.

Measured by replacing the pop with `min()`, on both real corpora:

| `sirv_real`, strict | set-pop | sorted |
| --- | --- | --- |
| isoforms | 106 | 108 |
| recall | 70.6% (48/68) | **72.1% (49/68)** |
| precision | 82.1% | **83.3%** |
| F1 | 0.759 | **0.773** |
| identity | 0.9673 | **0.9697** |
| recall, lenient | 79.4% | **82.4%** |
| F1, lenient | 0.870 | **0.892** |

| `droso` | set-pop | sorted |
| --- | --- | --- |
| isoforms | 505 | 504 |
| `FSM` | 443 | 443 |
| `NNC` | 12 | **10** |
| canonical introns | 0.982 | **0.983** |
| distinct intron chains | 359 | 357 |

A defined order is better on both, and the gain is on real data rather than simulated. Redundancy is
the one metric that moves the wrong way, and marginally (1.81 → 1.84 on `sirv_real`).

**Fixed**, on the owner's decision, by taking the lowest read id. Two reasons beyond the accuracy:
the reference no longer depends on an unspecified property of the interpreter, and a future CPython
that reorganised its set implementation would silently have changed isONform's output. That is the
same class of hazard as finding 1 — a stable-today artifact that nothing pins — and it would have
been much harder to diagnose, because unlike the hash seed it does not vary between runs.

The port therefore targets the fixed behaviour. Reproducing the old one would have meant
reimplementing CPython's set table, probing and pop finger in Rust, and then deleting it.

### Finding 13 — `remove_edges` never advances `prev_node`, so one branch is unreachable

`remove_edges` (`SimplifyGraph.py:279`) sets `prev_node = bubble_start` once per path and never
reassigns it. Inside the loop over the path's nodes:

```python
dist_to_prev = get_dist_to_prev(DG, prev_node, path_node)
if prev_node == bubble_start:          # always true
    prev_to_start_dist = 0
else:
    prev_to_start_dist = node_distances[prev_node]   # unreachable
dist = prev_to_start_dist + dist_to_prev
```

So every node's fallback distance is measured **directly from the bubble start**, never accumulated
hop by hop, and the `else` branch cannot execute. The comment beside it — "we found the distance to
the previous node, however we are still missing the distance of the previous node to s" — describes
chaining that does not happen. The variable, the comment and the dead branch together are strong
evidence the omission is unintended.

Measured by advancing `prev_node` and comparing:

| `sirv_real`, strict | from the start (current) | chained |
| --- | --- | --- |
| isoforms | 108 | 109 |
| recall | 72.1% (49/68) | 72.1% (49/68) |
| precision | 83.3% | 83.5% |
| F1 | 0.773 | 0.774 |
| redundancy | **1.84** | 1.86 |

| `droso` | from the start | chained |
| --- | --- | --- |
| isoforms | 504 | 503 |
| `FSM` | 443 | 443 |
| canonical | 0.983 | 0.983 |
| distinct chains | 357 | 357 |

A wash. Recall identical on SIRV, every SQANTI figure identical on Drosophila, and redundancy
marginally worse chained. **So the recommendation is not to change the behaviour** — measuring from
the bubble start accumulates no error and is arguably the better rule — but to delete the dead branch
and fix the comment, which is a pure no-op cleanup. The port reproduces the current behaviour and
says so in `simplify.rs`.

Worth noting what this run also settled: the two rules coincide whenever every read supports every
node on the path (the distances telescope), and diverge only where the read sets differ. That is why
the effect is small rather than absent.

### Finding 14 — a **live** seed dependency in bubble linearisation. Port diverges deliberately.

`find_connecting_edges` returns a **`set` of `(node, node)` tuples**, and node names are *strings*, so
its iteration order depends on `PYTHONHASHSEED`. `test_conn_end` filters that set into a list and
`prepare_adding_edges` then takes `conn_list[0]` — so if two connecting edges ever end at the same
node, which one is picked is seed-dependent, and finding 1's determinism fix does not cover it.

**This entry previously said "latent, not live". That was wrong, and the correction is the point.**
The original measurement found 0 of 343 non-empty `test_conn_end` results with more than one entry and
concluded the pick was unambiguous on real data. A holdout corpus reached it on the first try.
Re-measured by replaying recorded `simplifyGraph` calls:

| corpus | `test_conn_end` calls | non-empty | **more than one** |
| --- | --- | --- | --- |
| `sirv_small` (2 clusters) | 4 | 0 | 0 |
| Drosophila `dsim_mid` (16) | 3 125 | 2 | 0 |
| Drosophila sample (56) | 8 582 | 81 | 0 |
| Drosophila **holdout** (56) | 8 120 | 57 | **1** |
| total | **19 831** | **140** | **1** |

So it is rare — one occurrence in 19 831 calls, one in 140 non-empty results — and it is real. Rare
enough that the original conclusion was a reasonable read of the evidence available, and wrong enough
to have silently mis-specified the port. The lesson is the one finding 14 already carried in its
second half, turned on itself: **a negative result on a corpus is a statement about the corpus.**
"Zero occurrences in 343 samples" bounds the rate near 1-in-343; it does not establish zero, and it
was written up as though it did.

**And the reference is non-deterministic where it fires.** Replaying the one graph that reaches it
under eight `PYTHONHASHSEED` values produces **two distinct outputs** — seeds 0, 4, 6, 7 give one
graph, seeds 1, 2, 3, 5 another, differing in one edge's support (`18` against `6`). So finding 1's
"isONform was non-deterministic run to run — **fixed**" is true of minimizer selection and not true of
the tool: this is an independent second source, unfixed upstream, and it can change the output.

**What the port does.** Orders the candidates lexicographically by the source node's *name* and takes
the first — the same choice already approved for minimizer selection, and the order the reference's
own string names sort by. On the single observed multi-candidate case this agrees with
`PYTHONHASHSEED=0`, hence with every recorded dump, so the simplification oracle passes; but that is
a sample of one and is **not** evidence that lexicographic equals seed 0 in general. No env var
restores the reference behaviour, because the behaviour is a coin flip rather than a defined order —
there is nothing to restore. Pinned by
`the_connecting_edge_pick_is_the_lexicographically_smallest_source`, whose fixture is built so that
lexicographic and numeric order disagree, and which fails if the ordering is removed.

**The other half of the original finding still stands: the determinism check had been running on a
corpus that barely exercises this stage.** `bench/corpus/sirv_small` pops a single edge (244 → 243),
so "deterministic across 24 seeds" was a statement about graph construction and almost nothing about
simplification. Re-run on 16 medium Drosophila clusters that pop tens of edges each it was still
clean — but note that those 16 register **2** non-empty results and **0** multi-candidate ones, so
they could not have caught this either. **A determinism check is only as good as the paths its corpus
reaches.**

### Finding 16 — simplification strips the `length` attribute from 42% of edges

`prepare_adding_edges` re-adds every bubble edge with `DG.add_edge(u, v, edge_supp=value)`, naming
only the support. networkx merges attribute dicts, so an edge that survived keeps its `length` — but
the bubble edges were *removed* first, so they come back as new edges with **no `length` at all**.

Measured on one recorded Drosophila graph: **97 of 231 edges after simplification carry no `length`**,
against 0 of 380 before. So the graph's attribute schema changes halfway through the pipeline.

Harmless, as it happens: `GraphGeneration.py:402` is the only reader of `length` anywhere, and it runs
during construction. But it is observable, the dump records it as `NA`, and the port has to model it —
which is why `Graph`'s edge length is an `Option`, with `upsert_edge_support` reproducing networkx's
create-without-it / preserve-if-present behaviour exactly.

Not worth fixing on its own. Worth knowing before anyone writes code downstream that assumes every
edge has a length.

### Finding 17 — the printed pop count undercounts, and only on graphs big enough to notice

`new_bubble_popping_routine` adds the previous iteration's pops at the **top** of
its loop (`SimplifyGraph.py:927`) and prints the total after it
(`:1164`), so the final iteration's pops are never counted.

It only shows when the loop exits by the pop threshold rather than by running out
of candidates, and the threshold is `max(initial_edges / 100, 1)` — so below 200
initial edges the threshold is 1, an early exit means zero pops that iteration,
and nothing is lost. Above that it undercounts by the last iteration's total,
bounded by `threshold - 1`.

Measured on eight medium Drosophila clusters:

| cluster | edges at start | threshold | pops per iteration | true total | printed |
| --- | --- | --- | --- | --- | --- |
| 10 | 265 | 2 | 11, 6, 5, 4, 2, 3, 2, 1 | 34 | **33** |
| 11 | 95 | 1 | 3, 2, 3, 3, 1, 1, 0 | 13 | 13 |
| 12 | 219 | 2 | 30, 15, 8, 5, 4, 2, 2, 2, 1 | 69 | **68** |
| 13 | 291 | 2 | 15, 9, 8, 5, 5, 2, 1 | 45 | **44** |
| 14 | 290 | 2 | 7, 5, 3, 2, 1 | 18 | **17** |
| 15 | 245 | 2 | 10, 9, 5, 4, 3, 2, 1 | 34 | **33** |
| 16 | 219 | 2 | 13, 11, 6, 2, 1 | 33 | **32** |
| 17 | 148 | 1 | 12, 6, 4, 4, 3, 2, 1, 1, 0 | 33 | 33 |

Six of eight, and exactly the two with a threshold of 1 are correct — which is the
prediction, so the mechanism is understood rather than just observed.

Cosmetic: the number is only printed, never used. Worth fixing because it is one
line and because a progress figure that is quietly wrong is worse than none.
`PopStats` reports both the true count and the reference's.

### Finding 18 — an operator-precedence bug in the multi-path branch

`SimplifyGraph.py:1100`:

```python
if (not p1[0]) or (not p2[0]) and directpath_marked:
```

`and` binds tighter than `or`, so this reads `(not p1[0]) or ((not p2[0]) and
directpath_marked)`. The comment below it ("Marked") and the symmetry with
`directpath_marked` both say the intent was `((not p1[0]) or (not p2[0])) and
directpath_marked`.

The effect: a pair whose **first** path is a direct start-to-end edge is skipped
unconditionally, whether or not a direct path has already been popped this
iteration; a pair whose *second* path is direct is skipped only when one has.
Reproduced in the port with a comment, not fixed — measuring it needs a corpus
that reaches the multi-path branch with a direct path in the first position, and
that measurement has not been made.

### Finding 19 — the popping loop has no termination guarantee

`while has_combinations` exits on one of two conditions: no candidate survives
filtering, or an iteration pops fewer than the threshold. Neither is guaranteed.
If linearisation reaches a fixed point while every candidate is still *accepted*,
the loop keeps "popping" without changing the graph, records nothing as
not-viable, and never ends.

Found by construction rather than on real data: a three-path bubble plus an
aligner that accepts everything reaches it in three iterations. What prevents it
in practice is that the real aligner refuses most bubbles, and a refusal is what
writes to `not_viable_global` / `not_viable_multibubble` — so progress comes from
*failure*, not success. Nothing in the loop's structure enforces that.

The port carries an iteration cap that the reference does not have, reported in
`PopStats::hit_iteration_cap`, so the failure mode is a flag rather than a hang.
It is never reached on real input. Worth knowing before anyone changes the
aligner's acceptance criteria, because a more permissive aligner is exactly what
would expose this.

### Finding 21 — `generate_consensus_path` files a span under the wrong read, harmlessly

`SimplifyGraph.py:582` and `:585`, the two-supporting-read branch:

```python
if fdist > edist:
    consensus = reads[f_id][1][fstart: (fend + k_size)]
    seq_infos[f_id] = (fstart, fend + k_size, consensus)
else:
    consensus = reads[e_id][1][estart: (eend + k_size)]
    seq_infos[f_id] = (estart, eend + k_size, consensus)   # f_id, not e_id
```

When the second read's span wins, its coordinates are recorded against the
*first* read's id and the second read gets no entry at all.

**It cannot matter, because `seq_infos` is never read.** That took checking, and
the checking turned up more than the bug: `align_bubble_nodes` returns
`(is_poppable, cigar, seq_infos, consensus_log, spoa_count)` and the caller uses
**two of the five**.

| return value | read by the caller? |
| --- | --- |
| `is_poppable` | yes |
| `consensus_log` | yes |
| `cigar` | **no** — unpacked at both call sites, never referenced |
| `seq_infos` | **no** — same |
| `spoa_count` | **no** — threaded back in, incremented, never printed or compared |

The `multi_consensuses` memo stores `seq_infos` too, and its read path takes
elements 0 and 2 only. So the mis-filed entry, the CIGAR, and the spoa counter are
all dead.

Not worth fixing the assignment on its own; worth deleting three of the five
return values, which would have made the bug impossible to write. The port's
`AlignVerdict` carries only the two live ones.

### Finding 23 — read positions go negative after a bubble is popped

`additional_node_support` invents positions for reads that reach a node via the other branch of a
bubble: `newend = prev_end + relative_dist`, then `newstart = newend - avg_len`. Both can come out
**negative**, and do.

Measured across 16 recorded Drosophila simplifications: **57 nodes carry a negative position
afterwards, and none did before.** So a position here is a virtual coordinate that may precede the
read's own start, and anything consuming one has to tolerate that.

Not a defect — the invented coordinates are marked `original_support=False`, which is exactly the
flag for "this read was never really here" — but it has two consequences the port had to be corrected
for:

1. **`ReadInfo`'s positions must be signed.** The port had `u32` and clamped the negatives to zero.
   That was a divergence introduced without measuring it, and the simplification oracle caught it on
   its first run against real data.
2. **Python's negative-index slicing becomes reachable.** `generate_consensus_path` clamps `pos1` to
   zero in its spoa branch (`if pos1 < 0: pos1 = 0`) but **not** in its two-read branch, and
   `align_bubble_nodes` does not clamp its single-read branch either. So a negative position there
   slices from the *end* of the read rather than the start. The port now reproduces Python's slice
   semantics exactly where the reference does not clamp, and clamps where it does.

Both are the kind of thing that would have produced a small, hard-to-localise sequence difference
much later.

### Finding 24 — parasail scores a non-`ACGT` character as 0, even against itself. **Port bug, fixed.**

Not a reference defect — a defect in the port, and the most consequential one found so far. Recorded
here because the reason it was wrong is a reference behaviour nobody would guess.

`consensus.parasail_alignment` builds its substitution matrix with
`parasail.matrix_create("ACGT", match, mismatch)`. That builds a **5×5** matrix — one row and column
per base, plus a single catch-all — and fills the catch-all row and column with **zero**. So any
character outside `ACGT` scores 0 against everything, *including an identical copy of itself*. Read
off the library rather than inferred:

```text
m = parasail.matrix_create("ACGT", 2, -8)
sg("A", "A") ->  2      sg("X", "X") ->  0
sg("A", "C") -> -8      sg("N", "N") ->  0
sg("A", "X") ->  0      sg("a", "A") ->  2   (the lookup is case-insensitive)
```

The port's aligner scored by byte equality: `X == X` ⇒ `match_score`. Which is wrong, and reachable —
`generate_consensus_path` returns `"X" * max_len` whenever **every** path span in a bubble is shorter
than `k`, so two all-`X` placeholders being compared is an ordinary event, not a corner.

The consequence is not subtle. For a 17-`X` against an 18-`X` placeholder:

| | CIGAR | first non-match run | `parse_cigar_diversity` |
| --- | --- | --- | --- |
| parasail (all pairs score 0) | `16I1=17D` | 16 | **False** — 16 > `delta_len` 5 |
| port (X-vs-X scored as a match) | `17=1I` | 1 | **True** |

So the port popped bubbles the reference refuses. Fixed by giving `parasail.rs` a `subst` function
that models the real matrix — class-per-base plus a catch-all, case-insensitive, zero for anything
involving the catch-all — and using it at all five scoring sites. The CIGAR *operation* still comes
from byte equality, which is not an inconsistency: parasail emits `=` for two identical characters
even when the matrix scores them 0, so `X` against `X` is a zero-scoring `=`.

**A residual, measured rather than waved off — and then found to be nothing of the kind.** With every
cell scoring 0 there is no unique optimal alignment, and the port's traceback did not reproduce
parasail's choice: `17I18D` against `16I1=17D`. No setting of `TieBreak` reproduced it (all 24 swept),
so this was written up as **structural**. It was not: **finding 25 closed it**, and the missing piece
was not in `TieBreak` at all but in the end-cell scan admitting cells that consume none of one
sequence. The all-`X` comparison now matches parasail exactly.

The lesson is the same one finding 14 taught: "swept the parameter space and nothing fits" bounds the
*parameter space*, not the problem. Concluding "structural" from it was a category error, and it sat
in this file for one commit.

Worth carrying elsewhere: **isONcorrect's port has the same aligner and may have the same bug.** It
would only fire on non-`ACGT` input, so whether it is reachable there depends on whether anything
upstream can emit one; that has not been checked.

### Finding 25 — parasail's end-cell rule has two non-obvious parts. **Port bug, fixed; CIGAR now exact.**

A second port defect in `parasail.rs`, found the same way as finding 24 and with the same shape: a
reference behaviour nobody would guess, in a function this file had described as needing no
verification.

`parasail_sg` scans the last row and last column for the best score — trailing gaps are free — and
tracebacks from there. Two things about *which* cell it picks are not guessable. Both were measured,
not reasoned about: `result.end_query` / `result.end_ref` report parasail's choice directly, which is
what made them measurable at all.

**1. The scan starts at 1, not 0.** The port scanned from 0, admitting end cell `(n, 0)` — reached by
consuming all of `s1` as a free leading gap and none of `s2` — and its mirror `(0, m)`. Both score 0.
parasail excludes them, so it insists on at least one diagonal step even when every real alignment
scores below zero:

```text
m = parasail.matrix_create("ACGT", 2, -2)
sg("A",    "C")    -> -2  "1X"          not 0
sg("AAAA", "CCCC") -> -2  "3I1X3D"      not 0
sg("AC",   "CA")   ->  2  "1I1=1D"      one pair does match
```

**2. The corner is excluded from both ranges and considered last.** All-A against all-C ties
everywhere, which isolates the rule. The cell parasail picks, 1-based:

```text
  n\m     1      2      3      4      5      6
    1   (1,1)  (1,1)  (1,1)  (1,1)  (1,1)  (1,1)
    2   (1,1)  (2,1)  (2,1)  (2,1)  (2,1)  (2,1)
    3   (1,1)  (3,1)  (3,1)  (3,1)  (3,1)  (3,1)
    4   (1,1)  (4,1)  (4,1)  (4,1)  (4,1)  (4,1)
```

The tie is always between `(n, 1)` on the last row and `(1, m)` on the last column. For `m >= 2`
parasail takes the row one; for `m == 1` it takes `(1, 1)`, which is on the column. **No
row-before-column or column-before-row rule produces both** — which is exactly why sweeping the
`TieBreak` parameter space found nothing and why "structural" looked like a fair conclusion. Excluding
the corner produces both: with `m == 1` the row range `1..m` is *empty*, so the choice falls to the
column and lands on its first maximum.

**Measured, on 56 549 unique recorded calls** — 54 884 from 56 Drosophila clusters, 1 665 from
`sirv_small`:

| | before | after part 1 | after part 2 |
| --- | --- | --- | --- |
| score mismatches (of 54 884) | **12** | 0 | 0 |
| CIGAR mismatches (of 54 884) | 136 | 112 | **0** |
| CIGAR mismatches (of 1 665) | — | 35 | **0** |

Both are now gated exactly. This also closed what finding 24 recorded as a *structural* residual —
the all-`X` placeholder comparison now gives parasail's `16I1=17D` — and it is what closes
`simplify_0051`, the last failing case in the holdout corpus: the reference computed
`div = 0.2055` against a `0.2000` threshold and refused the bubble, while the port's differing tail
lowered the ratio and popped it. Same consensus sequences, same score, different CIGAR, opposite
verdict, because `parse_cigar_diversity` reads the CIGAR and not the score.

**Two method notes, both about how nearly this went wrong.**

*Sweeping bounds the parameter space, not the problem.* "All 48 `TieBreak` settings tried, none fits,
therefore structural" was written into this file and survived one commit. The missing piece was not a
`TieBreak` value; it was outside the parameterisation entirely. A negative result over a search space
says something about the search space.

*Do not tune on the cases you are trying to fix.* A sweep restricted to the 112 failing cases picks
`column_first = true` at a clean 112 / 112 — and takes the full corpus from 54 772 to 54 522. That
setting was applied on the strength of the subset result and the full-corpus oracle caught it
immediately. Finding the setting that fixes the failures is not the same as finding the right one.

### Finding 26 — `fill_p2` is off by one, so `solve_WIS` is measurably suboptimal. **Fix built, off by default.**

`fill_p2` builds `p[j]`, the largest index whose interval finishes at or before interval `j` starts,
and `solve_WIS` then evaluates `OPT[p[j]]` against the **1-based** `OPT` array. `fill_p2` stores a
**0-based** index:

```python
stop_to_max_j = {stop: j for j, (start, stop, w, _) in enumerate(...) if start < stop}
```

Two consequences. Every predecessor is shifted down by one, so compatible earlier intervals are
treated as incompatible; and `j_max` starts at 0 and means both "interval 0" and "nothing precedes
this", so an interval with no predecessor is credited with interval 0's optimum. Both make the
selection *conservative* — fewer intervals than the optimum, never an overlapping set — which is
exactly why it has gone unnoticed.

**This is not a new discovery.** It was found, proven by exhaustive search and fixed in the
isONcorrect port, on both the Rust and Python sides there. isONform still carries the original.
Re-measured here on isONform's own code path against brute-force maximum-weight independent set over
3 000 random instances: the fix is optimal on **3 000 of 3 000**, the reference is worse on **2 085**
and better on **none**. (isONcorrect measured 2 040 on its own generator — the same story from a
different draw.)

`WisOpts::fix_p2` implements it and is **off by default**, so the reference's behaviour is what runs
and every recorded golden still matches. Turning it on changes which intervals reach the graph, so it
is an accuracy question rather than a correctness one, and it wants measuring on the SIRV and
Drosophila corpora before anyone flips it. That measurement has not been done.

### Finding 27 — a third of the interval loop is unreachable, inherited from isONcorrect

`main`'s per-read loop maintains `previously_corrected_regions`, a `defaultdict(list)` that is
created empty, read in four places, and `del`-ed. **Nothing ever appends to it** (`main:423-495`),
and `prev_visited_intervals` beside it is created empty and never appended to either. So on every
path through the reference:

* `read_previously_considered_positions` is always the empty set, which makes
  `not_prev_corrected_spans` equal to `m1_curr_spans` — the filter is a no-op;
* `pos_group` is always empty, so `not_prev_corrected_spans2` is always empty;
* `all_intervals.extend(prev_visited_intervals)` extends by nothing;
* both `if previously_corrected_regions[r_id]:` branches never run.

That is roughly 35 lines including a nested position-grouping loop. In isONcorrect the same machinery
is *live*: correction happens in rounds and later rounds skip what earlier ones already fixed.
isONform runs one pass and never populates it, so the scaffolding came across without the thing that
drives it.

Not a defect — it produces the right answer, just via code that cannot run. Recorded because a reader
would reasonably assume it matters, because the port implements the reachable path only, and because
anyone adding a second pass later needs to know this is already half-built.

### Finding 28 — CPython set-iteration order names every isoform and orders every spoa call. **Deliberate divergence, measured.**

`compute_equal_reads` does:

```python
id = list(current_node_support)[0]
equal_reads[id] = list(current_node_support)
```

Both the representative `id` and the member **order** come from iterating a Python
`set`, and both reach the output: the `id` becomes the **isoform's identifier** in
`mapping.txt` and the consensus file, and the member order is the order sequences
are fed to **spoa**, which is order-sensitive, so it changes the consensus
*sequence*. Measured: `{3, 7, 8, 17}` iterates as `8, 17, 3, 7` and `{1, 5, 18}` as
`1, 18, 5` — slot order in a size-8 table under `x & 7`.

**This is the third time set order has reached output** (finding 12 in `find_paths`,
finding 14 in `prepare_adding_edges`), and it differs from both in the way that
decides the argument: read ids are ints and `hash(int) == int`, so this is **not**
`PYTHONHASHSEED`-dependent. The reference is stable run to run here, so diverging
buys no determinism — only simplicity.

**The decision, and why.** The port uses **ascending order**. Reproducing CPython's
would mean modelling the probing scheme, the `fill * 5 >= mask * 3` resize rule,
*and* insertion order propagated through `set()`, `.intersection()` and `-=`, since
collisions resolve by insertion order — a simplified model (ascending insertion,
linear probing) reproduces only 54 of 64 observed multi-member groups, so the
remainder really does turn on those details. That is a large amount of machinery,
coupled to a CPython implementation detail across versions, to preserve an order
that carries no meaning and that any future fix to this stage would change anyway.
A `BTreeSet` is both simpler and faster. Owner's call, taken on that basis.

**What it costs, measured rather than asserted.** Over 114 recorded real cases the
ordering differs on **28** (16 of 56 Drosophila sample, 12 of 56 holdout); on the
other 86 the two orders coincide and the output is byte-identical. On the 28 where
it does fire:

| | reference | port |
| --- | --- | --- |
| isoforms emitted | 1 066 | 1 073 (+0.7%) |
| byte-identical consensuses | — | **995 (93.3%)** |

So the divergence is confined to a quarter of cases, and within those it leaves
93% of isoform sequences untouched and changes the isoform count by under one
percent. What it does change is *labels*, on every affected group.

**The rest of the stage is verified exactly.** `rust/tests/isoforms_oracle.rs`
separates three outcomes: a wrong **partition** (reads in the wrong groups) fails
the build; a **merging** disagreement fails; a difference that is *only* set order
is reported, counted and costed. Crucially the merge is seeded with the
**reference's own** grouping read from the dump, not the port's — otherwise the 28
affected cases could never be checked at all, since the merge would be judged on
input it was never meant to see. With that: **0 wrong partitions and 0 merging
failures on all 114 cases.**

### Finding 29 — nearly a third of `IsoformGeneration.py` is unreachable

Of its 24 functions, **four are referenced nowhere in the repository** —
`search_last_entries`, `search_first_entries`, `parse_cigar_diversity_isoform_level`
(the older sibling of the `_new` one that is live) and `write_transcriptome_single` —
and two more are reachable only from commented-out call sites,
`generate_isoform_using_spoa` and `generate_isoforms_new`. That is roughly **183 of
631 lines**, and it is a larger version of the same pattern as finding 27.

Two of the four are informative rather than merely dead: `search_last_entries` and
`search_first_entries` take exactly the `delta_len_3`/`delta_len_5` arguments that
`parse_cigar_diversity_isoform_level_new` uses, so they read as helpers that were
inlined when the `_new` variant was written and never removed. And the live
`_new` sits directly above the dead original, which makes it easy to read the wrong
one — worth knowing before editing either.

Only the live path is ported.

### Finding 30 — every interior mismatch counts as `delta_len`, whatever its length

In `parse_cigar_diversity_isoform_level_new`, the accumulator is:

```python
# we want to add up all missmatches to compare to sequence length
miss_match_length += delta_len
```

`+= delta_len`, not `+= cig_len`, immediately under a comment describing the
opposite. Anything longer than `delta_len` has already returned `False` on the line
above, so the total is exactly `delta_len × (number of interior mismatch runs)`: a
one-base mismatch and a `delta_len`-base one contribute identically, and two
one-base mismatches contribute twice as much as one five-base one.

Not obviously wrong — counting *runs* rather than bases is a defensible similarity
measure — but it is not what the comment says and not what a reader would assume,
and the diversity ratio it feeds is compared against a threshold derived from
`delta`. Reproduced, and pinned by
`every_interior_mismatch_counts_as_delta_len_regardless_of_its_length`.

### Finding 31 — batch merging merges nothing. The code has never executed.

`actual_merging_process` is the whole point of the last stage: it compares every
isoform in every batch against every isoform in every *later* batch and folds the
duplicates together. It does nothing, and it never has.

```python
all_infos_list = sorted(all_infos_dict.items(), key=lambda x: x[0], reverse=True)
for b_i, (batchid, id_dict) in enumerate(all_infos_list[:len(all_infos_list) - 1]):
    for b_j, (batchid2, id_dict2) in enumerate(all_infos_list[b_i + 1:]):
        if not batchid2 <= batchid:      # <- never true
```

The list is sorted **descending**, so everything in `all_infos_list[b_i + 1:]` has
a *smaller* batch id than `batchid`. `batchid2 <= batchid` is therefore always
true, `not` makes it always false, and the entire body — the comparison, the
alignment, the merge — is unreachable.

**Measured, not read off.** On data engineered so that every pair is mergeable, the
guard was evaluated 3 times and the body entered **0** times, merging 0 of 5
consensuses. Flipping only the sort to ascending makes it enter immediately, which
isolates the sort/guard combination as the cause. And on the real driver:
`isONform_parallel` over `sirv_small` reports `merging changed anything: False` on
both clusters, and over 112 Drosophila clusters likewise.

**There is a second defect behind the first.** Flip the sort and the body runs
straight into `generate_consensus_path`, which does `all_sequences[id]` where `id`
is a read *accession* and `all_sequences` is `all_reads_dict`, keyed by *batch id*.
It raises `KeyError` on the first read. So the stage cannot be repaired by fixing
the guard: the lookup it needs does not exist in what it is passed. Two independent
defects, stacked, which is itself the evidence that this path has never run — the
first one has been hiding the second.

**What this means for the tool.** Isoforms are never merged across batches. A
cluster with more than `--max_seqs` reads (default 1000) is split into batches,
and an isoform present in two batches is emitted twice, with different ids
(`{cluster}_{batch}_{isoform}`). On corpora where every cluster fits in one batch
nothing is lost, which is why it has gone unnoticed; on larger clusters the
transcriptome carries duplicates.

**The port reproduces the no-op and does not offer a fix flag.** Everywhere else a
reference defect has an opt-in fix (`WisOpts::fix_p2`, `BuildOpts`) there was a
defensible correct behaviour to switch to. Here there is not: making it merge means
inventing the read-sequence lookup `generate_consensus_path` is missing, and
inventing behaviour is the one thing the working agreements forbid. Fixing it
properly is an upstream change, and it needs the owner to decide what
`all_sequences` was meant to be.

The rest of the stage — `write_final_output`, which chooses each isoform's
destination, id and support count — is live, ported, and verified.

### Finding 32 — the low-abundance support file always records 1

In `write_final_output`:

```python
if new_id in all_infos_dict:
    other_support_file.write("{0}: {1}\n".format(new_id, len(all_infos_dict[new_id].reads)))
else:
    other_support_file.write("{0}: {1}\n".format(new_id, 1))
```

`new_id` is a string like `"9_0_1"`; `all_infos_dict` is keyed by integer batch id.
The lookup never succeeds, so the else branch always runs and the support column is
the literal `1` regardless of the isoform's real support. Confirmed against the
reference: an isoform with **four** reads, below an `--iso_abundance` of 5, is
written as `9_0_1: 1`.

Also note `all_infos_dict[new_id].reads` would be wrong even if the lookup
succeeded — the values of that dict are per-batch dicts, not isoforms, so it would
raise `AttributeError`. Cosmetic in practice, since only the low-abundance support
file is affected and its counts are not read back by anything.

One more, recorded and not acted on: `write_final_output` is called **inside** the
per-batch loop in `join_back_via_batch_merging`, and opens its output files with
`'w'` each time. So for a cluster with *n* batches every file is written *n* times,
each rewriting the whole thing. The last write wins and the content is identical,
so it is wasted work rather than a wrong answer.

### Finding 33 — the sink's read lengths belong to the wrong reads, in every batch

`main` builds the length table once, over the whole file:

```python
read_len_dict = get_read_lengths(all_reads)        # main:530, keys 1..N over the file
DG, reads_for_isoforms = GraphGeneration.generateGraphfromIntervals(
    all_intervals_for_graph, k_size, delta_len, read_len_dict, new_all_reads)
```

and `generateGraphfromIntervals` walks it densely and uses the counter as a *read
id*:

```python
for i in range(1, len(read_len_dict) + 1):                      # GraphGeneration:282
    reads_at_end_dict[i] = Read_infos(read_len_dict[i], read_len_dict[i], True)
```

But the graph does not speak in file read ids. It speaks in **graph ids**, which
`main` restarts at 1 for every batch and assigns only to reads that survived the
interval filter. So `i` is a graph id on the left of that assignment and a file read
id on the right, and the two coincide only when the file has exactly one batch and
no read was ever skipped.

Measured on `sirv_small`'s 100-read cluster at `--max_seqs 25`, four batches:

| batch | reads in the graph | entries on the sink | phantom | graph ids whose length is another read's |
|---|---|---|---|---|
| 0 | 17 | 100 | 83 | **17 of 17** |
| 1 | 20 | 100 | 80 | **19 of 20** |
| 2 | 19 | 100 | 81 | **17 of 19** |
| 3 | 22 | 100 | 78 | **21 of 22** |

In all four batches 100 of 100 sink entries equal the length of the *file's* read of
that id, which is what pins the mechanism rather than merely fitting it. Note the
first row: this is not a multi-batch bug that a single-batch corpus escapes. Batch 0
is wrong too — 17 reads had intervals out of 25, so graph id 3 is already not file
read 3. What a single batch escapes is only the *large* version, where the ids
address a different part of the file entirely.

Two separate consequences, and it is worth keeping them apart:

* **phantom entries** — the sink carries a read set larger than the graph's, 100
  against 17–22 here. This was already recorded under finding 22 from the
  single-batch case;
* **wrong lengths** — the entries that *do* correspond to a real graph id mostly
  carry some other read's length. That is new, and it is the more serious half,
  because `start_mini_end`/`end_mini_start` on the sink are positions.

Reproduced, not repaired. Repairing it means choosing what the table was *meant* to
be keyed by, and both candidates change output: keyed by graph id over
`new_all_reads` the phantom entries vanish, keyed by file read id the sink stops
being indexable by graph id at all. That is an upstream design decision.

**This one bit the port.** The first version of the driver derived the lengths from
the batch rather than the file — a reading of "all reads" that is right in English
and wrong here. Every stage oracle still passed, because each replays *recorded*
inputs and none of them owns this wiring; and every corpus passed, because each
cluster fit in one batch, where the port's per-batch table and the reference's
whole-file table agree on the keys that get looked up. It took running both programs
end to end on a cluster forced into four batches to see it: batch 0 agreed, and
batches 1–3 collapsed to roughly one isoform per read (reference 6, 9, 3 against the
port's 16, 16, 17), because graph ids 1..20 were absent from a table keyed 26..50 and
every lookup missed. `driver.rs`'s `read_lengths_span_the_whole_file_not_the_batch`
pins it, and fails if the parameter is ever re-derived per batch.

### Finding 34 — under `--parallel` every batch overwrites the previous batch's output

`main` derives its output filenames from the *input filename* when `--parallel` is
given, and from the loop index otherwise:

```python
if args.parallel:                                   # main:373-378, before the loop
    p_batch_id = args.fastq.split("/")[-1].split("_")[-1].split(".")[0]
    skipfilename = "skip" + p_batch_id + ".fa"

for batch_id, reads in enumerate(batch(all_reads, args.max_seqs)):
    if args.parallel:
        batch_pickle = str(p_batch_id) + "_batch"   # constant across batches
    else:
        batch_pickle = str(batch_id) + "batch"
        skipfilename = "skip" + str(batch_id) + ".fa"
```

and again for the isoform half, `batch_id = p_batch_id` at `main:572`. `p_batch_id`
does not change as the loop runs, so all four output files — `{id}_batch`,
`skip{id}.fa`, `spoa{id}`, `mapping{id}` — are the same four names on every
iteration, opened with `'w'`. A cluster split into several batches keeps only the
last one; the earlier batches are computed in full and then overwritten, silently.

Measured on the same 100-read cluster at `--max_seqs 25`:

```
--parallel True   ->  0_batch  mapping0  skip0.fa  spoa0                     (4 files)
(no --parallel)   ->  0batch..3batch  mapping0..3  skip0.fa..skip3.fa  spoa0..3  (16)
```

Four batches ran in both cases — the reference prints `Working on 25 reads in a
batch` four times either way.

This matters because `--parallel` is not an unusual mode: `isONform_parallel` passes
it on every child invocation. It is invisible on the corpora anyone tests with,
because `--max_seqs` defaults to 1000 and clusters that large are rare, and the loop
then runs once. Reproduced: the port writes the batches in order and lets the last
truncate the rest, which leaves the same files with the same contents.

### Finding 35 — `--iso_abundance` discards silently, because `write_low_abundance` is a local

`isONform_parallel.main` opens with

```python
write_low_abundance = False
```

and never assigns it again. No command-line flag reaches it — there is no
`--write_low_abundance` in the parser at all. It is then passed down to
`join_back_via_batch_merging` and `generate_full_output`, so:

* every isoform with fewer than `--iso_abundance` (default **5**) supporting reads
  is written **nowhere** — not to a low-abundance file, not to a log, nothing;
* the entire low-abundance half of `write_final_output` is unreachable from this
  entry point, including finding 32's always-1 support count. That finding stands
  — it is reachable through `batch_merging_parallel` directly — but it cannot be
  triggered by running the tool.

Measured on `sirv_small`: `main` produces 52 isoforms for cluster 0, and the
default `--iso_abundance 5` leaves exactly **one** in `cluster0_merged.fa`. The
other 51 leave no trace anywhere in the output folder. Whether that is the
intended cutoff behaviour is a product question; that it is unconfigurable is not
obvious from `--help`, which describes `--iso_abundance` as a cutoff without
saying the remainder is dropped.

Reproduced, flag and all.

### Finding 36 — the parallel entry point rewrites its **input** folder in place

`restructure_isoncorrect_output(args.fastq_folder)` is the first thing
`main` calls. If the folder contains no files at the top level, it treats it as
isONcorrect output and:

```python
shutil.move(source_file, target_file)      # subdir/corrected_reads.fastq -> subdir.fastq
...
Parallelization_side_functions.remove_folders(directory)   # rmtree every subdirectory
```

That is a **move**, not a copy, followed by `shutil.rmtree` on every subdirectory
of the input folder — whatever else those subdirectories contained. The decision
rests on a single test: "are there zero files at the top level". A folder of
per-sample subdirectories that happens to have no loose files at the top qualifies,
whether or not it has anything to do with isONcorrect.

This is load-bearing — it is how isONcorrect's output is fed in — so it is
reproduced. It is recorded here because "run the tool on a folder" reading as "the
tool may delete that folder's subdirectories" is not something a `--help` line
prepares anyone for, and because a port that quietly made it non-destructive would
break the pipeline it exists for.

### Finding 37 — seven of `isONform_parallel`'s flags do nothing

The subprocess call is where most of them are lost. `isONform_algorithm_params`
collects eleven settings, and the `subprocess.check_call` that follows passes only
some of them; `--max_seqs` is present but **commented out**:

```python
["python", isONform_exec, "--fastq", read_fastq_file, "--outfolder", outfolder,
 "--exact_instance_limit", ..., #"--max_seqs", str(isONform_algorithm_params["max_seqs"]),
 "--k", ..., "--w", ..., "--xmin", ..., "--xmax", ..., "--delta_len", ...,
 "--exact", "--parallel", "True", "--delta_iso_len_3", ..., "--delta_iso_len_5", ...]
```

Measured rather than read off — each run against the default on `sirv_small`,
comparing every output file:

| flag | effect on output | why |
| --- | --- | --- |
| `--max_seqs 10` | **none** *(unless `--split_wrt_batches`)* | commented out of the subprocess call, so `main` uses its own default of 1000. With `--split_wrt_batches` it *is* live — it controls file splitting: 20 instances against 2. |
| `--set_w_dynamically` | **none** | collected into the params dict, never passed |
| `--delta 0.9` | **none** | only reaches `actual_merging_process`, which is a no-op (finding 31) |
| `--max_seqs_to_spoa 2` | **none** | same — only the no-op merging |
| `--verbose` | **none** | never read by anything on this path |
| `--keep_old` | **none** | see below |
| `--exact_instance_limit 100000` | **none** | passed to `main`, which ignores it (finding 3) |

`--keep_old` deserves its own line, because it is not merely unpassed — it is
checked against a file that does not exist:

```python
candidate_corrected_file = os.path.join(outfolder, "isoforms.fastq")
if os.path.isfile(candidate_corrected_file): ...
```

`isoforms.fastq` appears exactly once in the whole codebase: on that line. Nothing
writes it. So `os.path.isfile` is always false, `compute` stays true, and
`--keep_old` can never skip anything — the resume feature it advertises has never
worked. The intended name was probably `corrected_reads.fastq` (what the isONcorrect
script this was adapted from produced) or one of the `cluster*_merged.*` files.

All seven reproduced, since a port that made them live would change results.

### Finding 38 — output record order is filesystem directory order

Two places read a directory and use the order they get, unsorted:

```python
for batchfile in os.listdir(cl_dir):                        # batch_merging_parallel:235
subfolders = [f.path for f in os.scandir(outfolder) if f.is_dir()]   # side_functions
```

The first decides the insertion order of `all_infos_dict`, which is the order
isoforms appear in `cluster{n}_merged.fa`; the second decides the order clusters
are concatenated into `transcriptome.fasta`. Neither is sorted, so both are
`readdir(3)` order.

Visible as soon as a cluster has more than one batch. On `sirv_small` with
`--split_wrt_batches --max_seqs 25`, cluster 0's four batches come out in the order
`0_2_2, 0_0_2, 0_1_1, 0_3_1` — batch 2 first.

The port uses `std::fs::read_dir`, which is the same `readdir(3)` order, and the
two agree file-for-file on every corpus tested including that split run. They would
not be guaranteed to agree across different filesystems, but neither would two runs
of the reference. Recorded as a property of the tool rather than a port risk:
`transcriptome.fasta`'s record order is not reproducible across machines, and
anything downstream that diffs it needs to sort first.

### Finding 22 — smaller things, recorded and not acted on

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
- **The sink node's read set is larger than the graph.** `read_len_dict` is built in `main` over
  `all_reads` rather than `new_all_reads`, so node `"t"` gets a `reads` entry for every input read
  including those the interval filter skipped: 100 entries against 84 reads with a path, on
  `sirv_small`. The source node is keyed correctly (`graph_id` starts at 1, so there is *no*
  off-by-one there — checked, because there looked like one). Whether the phantom entries matter
  depends on how `SimplifyGraph` uses `DG.nodes[node]['reads']` for the sink, which is the next stage
  to port; reproduced faithfully until then. **Superseded by finding 33**, which shows the phantom
  entries are the lesser half of this: the entries that *do* name a real graph id mostly carry a
  different read's length.
- A nonexistent `--fastq` crashes with a `FileNotFoundError` traceback (exit 1) rather than a
  diagnostic.

## What came across from isONcorrect, and what was left behind

Three modules were copied in rather than rewritten. Recording the trim, because "we reused the
alignment code" is not the same claim as "we reused all of it", and the difference is where the next
person will look for something that is not there.

| module | lines there | lines here | what changed |
| --- | --- | --- | --- |
| `poa.rs` | 179 | 179 | verbatim |
| `parasail.rs` | 1 709 | 1 221 | research scaffolding removed |
| `align.rs` | 604 | 351 | edlib-adjacent code removed |

**Removed from `parasail.rs`:** `global_affine` plus `global_check`/`global_check_report` and their
counters, and `band_check`/`band_check_report`. All of it existed to answer two isONcorrect
questions — could exact affine *global* alignment replace the semi-global one, and was banding safe —
and both were judged by whether `guard::fix_correction` landed on the same corrected sequence.
isONform has no guard, so there is nothing to compare. Both blocks ran only under an `ISONCORRECT_*`
environment variable, so dropping them changes nothing on a normal run. Also removed: the
`blockalign` and `one_case` test modules, which compared against `block-aligner` — a dependency
isONform therefore does not need.

**Removed from `align.rs`:** the `oracle` and `block_option` test modules, and
`speed_and_optimality_at_segment_lengths`. They timed and cross-checked a SIMD edit-distance path
against recorded edlib CIGARs. **isONform never imports edlib** (*Reconnaissance corrections* 3), so
there is no edit-distance path here to accelerate and no decision to make. `Alignment` and `global`
are kept — they are self-contained and cost nothing — but nothing calls them yet.

**Changed:** `Scoring::GUARD` (`match=4, mismatch=-8`) is gone, replaced by isONform's two:

| constant | scoring | call site |
| --- | --- | --- |
| `Scoring::BUBBLE` | `2, -8, 12, 1` | `align_bubble_nodes` — are these two bubble paths the same? |
| `Scoring::MERGE` | `2, -2, 12, 1` | `IsoformGeneration.align_to_merge` — are these two isoforms the same? |

They differ only in the mismatch penalty, and that difference is the whole distinction between the
two questions, so they are separate constants with a test asserting they stay distinct. `-2` is also
`parasail_alignment`'s own default, which `SimplifyGraph` overrides and `IsoformGeneration` does not.

The parasail semantics tests keep the scoring they were **originally verified against** (`match=4`)
rather than being re-derived under `BUBBLE`. Those expected scores were read off the parasail library
itself; recomputing them arithmetically at a new scoring would turn a measurement into an assumption.
What they establish — free end gaps on both sequences, a one-base gap costing `open`, longer gaps
costing `open + (L-1) * ext` — is a property of parasail's semi-global mode and holds at any scoring.

### One verification did not transfer — so it was repeated. **Closed.**

`poa.rs` came with a strong claim: identical consensus to the `spoa` binary on **505 of 505** real
isONcorrect correction intervals. **That number does not carry over.** isONform calls spoa on
different sequences from a different stage — bubble-path consensus, not correction intervals — so the
check had to be repeated here before a simplification oracle could attribute a disagreement to one
side or the other.

Repeated, and it holds: **5 368 of 5 368 recorded isONform spoa calls produce an identical
consensus.** Recording is `bench/dump_reference.py --record-spoa`, which wraps
`IsoformGeneration.run_spoa` and writes the `SPOA_CASES` TSV that `poa.rs`'s `oracle` test already
read — so closing this needed no new test, only isONform's inputs fed to the existing one.

| corpus | unique calls | mismatches |
| --- | --- | --- |
| `sirv_small` (2 clusters) | 32 | 0 |
| Drosophila (56 clusters, size-stratified) | 2 425 | 0 |
| SIRV real (7 clusters) | 820 | 0 |
| Drosophila holdout (56 more, disjoint) | 2 091 | 0 |
| **total** | **5 368** | **0** |

2 to 96 sequences per call, consensus 20 to 1 742 bases, and all three live call sites covered.
On the first three corpora (3 277 cases), the split was `SimplifyGraph.py:570` 3 074 — the bubble-path
consensus the simplification oracle's verdicts depend on — with `IsoformGeneration.py:413` 132 and
`:436` 71; the holdout adds 1 978 / 46 / 67 in the same three. One attribute
patch reaches all of them, because two call sites use `IsoformGeneration.run_spoa(...)` and two a
bare `run_spoa(...)` inside that module, and both forms resolve through the module dict at call time.
`batch_merging_parallel.py:34` is the fourth live site and recorded nothing on these corpora, so it
is covered by construction but not by measurement.

A useful incidental cross-check: `sirv_small` records **zero** calls from `SimplifyGraph.py:570`,
which is independently why the simplification oracle reports "0 spoa calls" on that corpus. Two
instruments built for different purposes agreeing on a fact neither was aimed at.

**What this licenses, precisely.** Given the same input sequences in the same order, `crate::poa`
agrees with `spoa`. It does *not* say the port computes the same input sequences — if span extraction
diverges, spoa is handed different sequences and faithfully returns a different consensus. That is
still a bug on the port's side, which is why the simplification oracle's gate is now unconditional
rather than merely reclassified.

**The parasail side was described here as needing no such caveat — "exact by construction and verified
at the CIGAR level". That was wrong twice over, and finding 25 is the consequence.** The score is
exact by construction; the CIGAR is not, because a semi-global alignment usually has several optimal
paths and which one is reported is a tie-break. And "verified at the CIGAR level" was verified against
*isONcorrect's* recorded calls, which is the same non-transfer this section is about. Recording
isONform's own (`--record-parasail`) found 12 outright **score** errors in 54 884 calls.

### What closing it exposed: one bug wearing two disguises. **Fixed.**

Tightening the simplification oracle to fail on every disagreement turned its two
"reported, not attributable" Drosophila cases into failures. They looked like separate bugs, and I
recorded them as separate bugs:

- **`simplify_0002`** — a *topology* difference. The reference leaves a bubble standing, two parallel
  paths carrying 15 and 3 reads; the port pops it into one chain whose nodes carry 18 = 15 + 3.
- **`simplify_0054`** — topology *agrees exactly*, not one edge added or missing. Only the node
  `reads` maps differ, mostly by one entry each way (36 vs 35, 22 vs 23) with one outlier at 2 vs 31.

That reading was wrong, and the way it was wrong is the useful part. Symptoms that far apart —
one changing the graph's shape, one changing only bookkeeping inside an identical shape — invited
"two different mechanisms", and I wrote that down before checking. **They were the same bug**, and
one fix closed both. Finding 24 below is what it actually was.

The route to it, because the shortcut did not exist: on `0054`, the port's node `reads` differed only
in **synthetic** entries — `original_support == false`, i.e. invented by `additional_node_support` —
which pointed hard at that function, and that function turned out to be faithful. What located the
bug was stepping back out: the port ran the same **15 iterations** but did **94 pops against the
reference's 92**, and a pop-by-pop diff showed the first 81 identical and the 82nd extra. So `0054`
was never a bookkeeping bug at all — it was two surplus pops whose extra `additional_node_support`
calls left synthetic reads behind, while the topology happened to reconverge. The same
poppability-decision bug as `0002`, seen through a graph that hid it.

Worth keeping as a method note: "the diff is only in synthetic reads" was a true observation that
pointed at the wrong function, and the thing that corrected it was a *count* — iterations and pops,
which the reference prints and the port now reports too (`PopStats`, surfaced per case by the
oracle). A cheap aggregate beat a precise-looking local signal.

**The holdout is what kept the job honest.** A second, *disjoint* 56-cluster Drosophila sample — same
stratification, offset so it shares no cluster with the first — failed on 2 of 56 after finding 24's
fix, when the first corpus had gone fully green. Green on the corpus you debugged against is the
weakest evidence there is; both of those turned out to be real, and neither was the mechanism its
symptoms suggested:

- **`simplify_0002`** — **one** edge in the whole graph. Same 21 reads on both sides, same *set*,
  different **multiplicities**. Finding it took extending the oracle to compare edge support as a
  multiset — a set difference reports "no difference" on precisely this case. My first reading, "a
  path-ordering difference in `prepare_adding_edges`", was half right and one step short: the pop
  sequences were **byte-identical** (86 pops), so the divergence was inside a single pop, and it was
  the `PYTHONHASHSEED`-dependent `conn_list[0]` pick of finding 14 — a case where **the reference has
  no single answer to match**.
- **`simplify_0051`** — surplus synthetic reads, the same signature pre-fix `0054` had, and again the
  same *class* of cause but not the same cause. Both sides computed byte-identical consensus
  sequences and reached opposite verdicts, which moved the search out of `simplify.rs` and into the
  aligner: finding 25, parasail's end-cell rule.

Both are fixed and the holdout is 56/56. Three of the four bugs in this stretch were found by a
corpus the previous one could not reach, which is the argument for keeping a holdout rather than
growing the corpus you debug against.

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

### Two order dependencies, one real and one not — and why both are written down

Both were asserted before being checked, and one of the assertions was wrong. Recording the outcome
of the check, not just the conclusion, because "is this order observable" is a question that will
come up again at every stage.

**Not observable: a node's `reads` map insertion order.** Claimed at first to reach results via
`SimplifyGraph.py:84`'s `tuple(DG.nodes[node]['reads'])`. It does not. All five consumers are
order-independent: `:84` and `:103` convert to a `set` and then sort; `get_dist_to_prev` (`:218`) and
`get_avg_interval_length` (`:295`) sum and divide, which commutes; `additional_node_support` (`:311`)
writes one independent value per read id; `:491` merges dicts whose result only feeds the others. The
port still keeps a `Vec` and the dump still records the order unsorted — a linear scan over a
few-entry list beats hashing anyway, and a stricter oracle costs nothing — but the *reason* is
"free and faithful", not "load-bearing".

**Observable: the topological order, and the candidate-bubble enumeration order.** `nx.topological_sort`
returns one of many valid orders, and `SimplifyGraph.py:115` compares the resulting *indices* to
decide which node pairs are candidate bubbles — so a port producing a different valid order finds
different bubbles. networkx 2.8.4's `topological_sort` yields from `topological_generations`, making
it **generation-by-generation Kahn** (not the LIFO variant), seeded in node insertion order and
advancing in parent-order × successor-order. Reproduced in `Graph::topological_sort`, and the oracle
compares the order element by element rather than merely checking it is a valid topological order.

The enumeration order of candidate bubbles matters too, but narrowly: `SimplifyGraph.py:946` sorts
candidates by bubble span, and Python's `sorted` is stable, so the generation order survives only as
the tie-break among equal-span bubbles. That is enough to make a cheaper enumeration a behaviour
change rather than a free win — worth knowing before anyone optimises the quadratic start × end loop.

### Dead code found in the simplification stage

- `filter_combinations` (`SimplifyGraph.py:122`) is **never called**;
  `new_bubble_popping_routine` inlines the same filter as a comprehension at `:939`.
- `find_possible_ends`'s docstring is a copy of `find_possible_starts`'s and says "at least 2 out
  nodes"; the code uses `in_degree`. The code is right.

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
5. ~~Port the CLI~~ (**done** — `rust/src/cli.rs`, 39 differential cases), then work outward from
   the leaves — **graph construction is ported** (`rust/src/graph.rs`, `graph_build.rs`) with its
   dump-based oracle recorded; simplification next, then isoform generation. The front half of `main`
   comes across from isONcorrect rather than being written (*Reconnaissance corrections* 4), and is
   still owed: without it the Rust binaries cannot reach the graph stage on their own, which is why
   the oracle replays recorded inputs for now.
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
- **Stray debug prints** (finding 4) and the unused `from ast import Param` (finding 22): one-line
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
