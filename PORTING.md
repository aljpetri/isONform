# isONform — Rust rewrite

## Status: ported, scored, and 2.2--3.7x faster than the port as first written

Both entry points run end to end, the differential oracles pass, and the defaults
now carry three measured changes on top of the faithful port. **Start here**, then
read findings 40--49; findings 1--39 and the reconnaissance are in
`PORTING-HISTORY.md`.

### Where it stands, measured

| | python | port as first written | **defaults now** |
|---|---|---|---|
| droso | 49.0s | 24.9s | **11.4s** |
| sirv_sim_gene | 354.8s | 107.8s | **30.0s** |
| sirv_real | 241.0s | 220.1s | **59.6s** |
| deep-12 within-batch | --- | 1 403s | **179s** |
| sirv_sim F1 / redundancy | 0.947 / 2.19 | 0.940 / 1.38 | **0.946 / 1.25** |
| sirv_real strict / lenient F1 | 0.773 / 0.892 | 0.734 / 0.873 | 0.722 / 0.862 |
| droso FSM / ISM / NIC | 443 / 14 / 29 | 470 / 12 / 30 | 466 / **9** / **26** |

On by default, each backing out alone: **WFA2** at the two alignment sites
(`ISONFORM_WFA2=0`), the **consensus rebuilt once per group at the end** rather
than per merge (`ISONFORM_MERGE_REBUILD_MAX=50 ISONFORM_FINAL_CONSENSUS=0`), and
**node identity by minimizer pair** (`ISONFORM_PAIR_NODES=0`). The one regression
is `sirv_real` strict F1, 0.734 -> 0.722.

### What is worth doing next, in order

1. **Confirm the recursive merge's batch floor on `droso_deep`.** Finding 52: on
   `sirv_real_deep` a floor of 10 batches beats the pairwise merge on every metric
   and is 1.39x faster. `droso_deep` is measured only at floor 2, where it already
   gains. One run settles whether 10 is the value for both.
2. **Bound the recursive stage's concurrency.** It costs 3.1 -> 9.7 GB, which is
   eight threads times a handful of large clusters --- p99 of the instance size is
   104 sequences and the max is 3 466.
3. **Verify the anchor-chain merge filter on a second deep cluster.**
   `ISONFORM_MERGE_MINSHARE=30` is sound on the two measured so far, and is still
   opt-in.
4. **Find where WFA2 actually loses.** Finding 50 closed the two candidates the
   record named --- the dangle DP and the free-end budget --- and neither was it.
   The remaining `sirv_real` cost is somewhere in WFA2's *core*, unlocated.

~~**Merge in the graph.**~~ **Done** --- finding 52, and it is the largest win in
this round: 1.2--1.4x faster than the pairwise cross-batch merge and better on
every `sirv_real_deep` metric at a batch floor of 10.

~~**The WFA2 dangle DP.**~~ **Done, and negative** --- it costs 0.029 strict F1 on
`sirv_real`. Finding 50.

### A warning about method, earned repeatedly

Every hypothesis in this file that was stated before it was measured turned out
wrong --- twenty-odd of them, including several that looked certain. Comparisons are
only valid at fixed N, `nodes-per-read` is dominated by N at small sizes, and a
binary must be re-frozen after every rebuild or a sweep measures the previous one.
Measure first, then conclude.

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

## Where the rest of the record lives

Findings **1--39**, the reconnaissance notes, and the isONcorrect-transfer analysis
are in **`PORTING-HISTORY.md`** --- settled reference bugs, corrected
reconnaissance, and superseded analyses. Code comments cite findings by number, so
look there for anything below 40. Findings 40 and up, the contract, the evaluation
harness and the method are here.

Finding 45's heading is kept as written even though finding 47 answers it, because
the sequence matters: the prefix sweep looked like a scale effect and was not.

## Goal

Port isONform from Python to Rust. **Identical CLI**, faster, lower memory. Same objectives as the
isONcorrect port, including the one that mattered most there: **the specification is accuracy, not
byte-identity**. Byte-identity is the tool used to prove the port faithful stage by stage; it is not
the thing being optimised, and where the reference is wrong or slow-by-compromise, the port may
diverge — measured, documented, and behind an env var that restores the reference path.

Two entry points must keep their names, flags and defaults:

- `main` — the single-instance corrector/assembler
- `isONform_parallel` — fans out over a folder of per-cluster fastqs

## The bug-fix policy changed on 2026-08-25

Until this point the rule was: reference defects are reproduced by default, fixes
are opt-in and measured. That produced a verified baseline and 38 findings, and
then stopped being the right rule. The rule now is **known bugs are fixed by
default** --- a bug-free implementation is what methodological improvement can be
built on, and a bug-compatible one with better scores is the worse artefact.

This deliberately accepts lower scores where fixing a bug costs accuracy.
`wis_p2` does exactly that (finding 26) and is on anyway.

The mechanism is `ISONFORM_BUG_COMPAT` --- a comma-separated list of bug names, or
`all` --- which puts the named bugs **back**:

| name | reproduces | finding |
| --- | --- | --- |
| `wis_p2` | `fill_p2` off-by-one, a measurably suboptimal WIS | 26 |
| `stale_seq` | a node attribute from the wrong read's sequence | 9 |
| `cycle_continuation` | the cycle-avoidance node severing a read's path | 11 |
| `batch_merge` | cross-batch merging doing nothing at all | 31 |

It has two legitimate callers and no others: the differential oracles, which
replay recorded reference behaviour and therefore need the reference's bugs, and
bisecting an accuracy change to one bug. An unknown name is an error, not a silent
no-op.

`WisOpts::default()` and `BuildOpts::default()` are now the *correct* behaviour and
`::reference()` is the opt-in. The inversion is explicit in the type system rather
than a flipped boolean because several call sites used `default()` to mean "the
reference", and silently changing what that meant would have been a trap --- the
kind that makes an oracle pass while testing nothing.

`bench/compare_end_to_end.py` sets `ISONFORM_BUG_COMPAT=all` on every port
invocation, which is what keeps it an equivalence gate rather than a diff of two
programs meant to differ; `--no-bug-compat` shows what the fixes changed. Getting
that wrong is not hypothetical: the first attempt plumbed the variable through the
parallel path but not the main-entry loop, and the gate immediately reported two
disagreeing clusters.

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

### What the port scores

Run on 2026-08-24, `--t 8 --split_wrt_batches --iso_abundance 5`, reference at
`PYTHONHASHSEED=0`. `bench/evaluate.sh run` takes either implementation as an
impl directory — the repository root for the reference, `rust/target/release`
for the port.

| corpus | metric | reference | port |
| --- | --- | --- | --- |
| `sirv_sim_gene` | isoforms | 138 | 139 |
| | strict recall | 91.2% (62/68) | **91.2% (62/68)** |
| | strict precision | 98.6% | **98.6%** |
| | strict F1 | 0.947 | **0.947** |
| | identity | 0.9976 | **0.9976** |
| | redundancy | 2.19 | 2.21 |
| | wall clock | 377.7 s | **124.0 s** |
| `sirv_real` | isoforms | 108 | 110 |
| | strict F1 | 0.773 | **0.746** |
| | strict recall | 72.1% (49/68) | **70.6% (48/68)** |
| | lenient F1 | 0.892 | **0.892** |
| | lenient recall | 82.4% | **82.4%** |
| | identity (strict) | 0.9697 | **0.9717** |
| | wall clock | 238.7 s | **187.2 s** |
| `droso` | isoforms | 504 | 506 |
| | FSM | 443 (88%) | **446 (88%)** |
| | canonical splice | 0.983 | **0.983** |
| | distinct intron chains | 357 | **357** |
| | median error | 0.0134 | **0.0132** |
| | FSM / tx in touched genes | 42.6% | **42.6%** |
| | intron retention | 9 | 8 |
| | wall clock | 48.7 s | **25.6 s** |

Two of three corpora are a dead heat, and the port is 1.9–3.0× faster. `sirv_real`
is the one that needs explaining, and it is worth doing properly because the first
reading of it is wrong.

**It is not seed noise.** The obvious hypothesis was that finding 14's live seed
dependency was moving the reference around and the port had landed on one side of
it. Measured instead of assumed: the reference is **bit-identical at seeds 0, 1
and 2** on both SIRV corpora — same 108 isoforms, same 0.773, same 49 of 68. So
the gap is reproducible and had to be looked at.

**It is three transcripts, and the port wins one of them.**

| transcript | reference identity | port identity | |
| --- | --- | --- | --- |
| SIRV206 | 0.9571 (pass) | 0.9485 (fail) | port worse by 0.009 |
| SIRV305 | 0.9536 (pass) | 0.9432 (fail) | port worse by 0.010 |
| SIRV615 | 0.9421 (fail) | **0.9824 (pass)** | port better by 0.040 |

All six are inside the length tolerance; every one of them turns on identity
alone, and two of the three sit within 0.01 of the 0.95 line. Aggregate identity
is *higher* for the port (0.9717 against 0.9697) and at the lenient threshold the
two sets are identical to three decimals on recall, F1 and identity. This is the
case the two thresholds exist to distinguish: sequences moved **across a line**,
not sequences that got worse. Quoting the strict F1 alone here would report a
regression that the lenient row shows is not there.

**Where it comes from, checked rather than assumed.** 21 of 26 clusters differ
end to end on this corpus — a far higher rate than the Drosophila corpora's 16 of
56 and 12 of 56, high enough that "probably finding 28" was not good enough. All
four stage oracles were run against reference dumps of this corpus:

| oracle | result |
| --- | --- |
| intervals | **30 of 30 agree** (9 968 reads in, 9 967 kept) |
| graph | **30 cases clean** — 43 143 nodes and 50 596 edges reproduced exactly |
| simplification | **30 of 30 pass** |
| isoforms | **0 wrong partition, 0 merging**, 25 set-order only |

So the whole divergence is finding 28 and nothing else: no wrong partition, no
merging disagreement, and the front half of the pipeline exact on an order of
magnitude more graph than any previous corpus. The higher rate is what finding 28
predicts — its effect scales with how many groups reach the merge scan, and
`sirv_real`'s clusters are much larger (cluster 0 has 4 595 reads and ~75
isoforms, against 5–488 reads per cluster in the Drosophila sets).

**What this licenses saying.** The port reproduces the reference's accuracy on
simulated data and at transcriptome scale, and on real SIRV differs only by which
side of a threshold three marginal consensuses land on, in both directions, as a
consequence of one documented and deliberate divergence. It does not license
"the port is more accurate" — SIRV615 is one transcript, and the aggregate
identity difference is 0.002.

**A caveat on the timings.** They are single runs on a shared laptop and the
seed-1 and seed-2 reference runs overlapped with `droso` scoring, so treat the
ratios as "roughly 2–3×" rather than as benchmarks. The three reference runs of
`sirv_real` came in at 238.7, 239.3 and 240.5 s, which at least says the
measurement is not wildly noisy.

### The three opt-in fixes, measured: all three are dead ends

`ISONFORM_FIX` (comma-separated, or `all`) turns on the fixes built behind
`WisOpts` and `BuildOpts`. An environment variable rather than a flag because both
entry points' argument surfaces are asserted to match the reference's and
`bench/equivalence.sh` drives both binaries — the same convention as
`ISONFORM_TRACE_POPS`. An unknown name exits 1 rather than no-opping, because a
typo would silently measure the baseline twice and read as "the fix changed
nothing".

Baseline is the port with no fixes. Twelve runs, four configurations across all
three corpora:

| fix | `sirv_sim_gene` | `sirv_real` | `droso` |
| --- | --- | --- | --- |
| `cycle_continuation` (finding 11) | **byte-identical** | **byte-identical** | **byte-identical** |
| `stale_seq` (finding 9) | F1 0.947 → 0.939, recall 62 → 61 | F1 0.746 → 0.745 | FSM 446 → 446, chains 357 → 357 |
| `wis_p2` (finding 26) | F1 0.947 → 0.943, recall 62 → 61 | F1 0.746 → **0.764**, recall 48 → **49** | FSM 446 → 446, chains 357 → 357 |
| `all` | F1 0.947 → 0.943, recall 62 → 61 | F1 0.746 → **0.765** | FSM 446 → 447 |

**`cycle_continuation` never fires.** Identical output on all three corpora,
including 561 Drosophila clusters. Finding 11's branch is real in the source and
unreached in practice, so there is no evidence for or against fixing it and no
reason to carry it as a live option. Consistent with the original measurement: one
occurrence in 19 831 calls.

**`stale_seq` changes output and improves nothing.** It moves every corpus and
helps none: worse on simulated, flat on real, flat on Drosophila.

**`wis_p2` is net negative, and the strict numbers alone would have said the
opposite.** On `sirv_real` it gains a transcript at the strict threshold, which
looks like a win. But its **lenient** recall is unchanged at 82.4%, so that gain
is a threshold crossing, not a recovered transcript. On `sirv_sim_gene` the
lenient recall drops 91.2% → 89.7% — a transcript lost at *every* threshold, which
is a real loss. One genuine loss against one threshold crossing.

**Finding 26 needs re-reading in this light.** `fix_p2` is provably the more
correct maximum-weight independent set — brute force says the reference is worse on
2 085 of 3 000 random cases — and it makes isoform reconstruction slightly *worse*.
A better interval selection is not a better transcriptome. The natural reading of
finding 26 was that the off-by-one was costing accuracy; measured, it is not.
Interval-selection optimality is not this tool's bottleneck.

**What did not move at all: redundancy.** 2.21 → 2.26 on simulated, 1.81 → 1.78 on
real, 357 distinct intron chains on Drosophila in every configuration. So none of
these three touches the thing most worth fixing, and the decomposition below still
points where it pointed.

### Where the redundancy actually comes from

Measured on `sirv_sim_gene`, strict matches, port baseline: 137 matching isoforms
for 62 transcripts, redundancy **2.21**. Decomposing the 75 excess isoforms by
whether the duplicate sits in the same batch or a different one:

| | excess isoforms | share |
| --- | --- | --- |
| same transcript, **different batches of one cluster** | 53 | **71%** |
| same transcript, same batch | 22 | 29% |

SIRV301 is emitted seven times — three isoforms from batch 0 and four from batch 1
of cluster 2. Cross-batch duplication is exactly what finding 31 predicts: batch
merging is a no-op, so an isoform present in two batches is emitted twice under
different ids and nothing ever collapses them.

Making cross-batch merging work would take redundancy **2.21 → 1.35** with no
change to recall or precision, which is a larger effect than anything else measured
so far and the reason finding 31 is the next thing to do. It needs a decision that
is not the porting project's to make: `generate_consensus_path` looks up read
accessions in `all_sequences`, which is keyed by batch id, so the lookup it needs
does not exist and inventing one is a design change.

The remaining 29% is intra-batch — `merge_consensuses` declining to merge two
consensuses of the same transcript within one batch — and is a separate question
about the `delta` / `delta_len` / `delta_iso_len_*` thresholds, none of which has
any accuracy measurement behind it (see *Deferred improvements*).

### The four remaining unambiguous bugs, fixed and measured

Fixed 2026-08-25: findings 32, 34, 37 and 30. Each is reversible through
`ISONFORM_BUG_COMPAT` (`batch_names`, `cigar_diversity`, and the two that need no
switch because they change no decision).

**Finding 32 --- low-abundance support.** Writes the isoform's real read count
instead of the literal 1. Unreachable from `isONform_parallel` (finding 35), so it
changes no measured output.

**Finding 37 --- `--keep_old`.** It looked for `isoforms.fastq`, a name that
appears exactly once in the codebase, on the line looking for it. It now checks
that the batch file holds one record per input read. Note the reference's test
could never have passed even with the right filename: it compares `wc -l` of a
fastq (four lines per read) against `wc -l` of the candidate.

**Finding 34 --- the filename collision.** The batch key genuinely needs both
components: `p_batch_id` distinguishes one `main` invocation from another (several
write into one cluster folder under `--split_wrt_batches`), the loop index
distinguishes batches within an invocation, and neither alone is unique. So the
key is `p_batch_id` for a single-batch invocation and `{p}_{b}` for several,
which leaves the common case byte-identical --- including the three-component
`{cluster}_{batch}_{isoform}` id --- and only adds a component where the
alternative was losing a batch. Measured on the four-batch case: **22 isoforms
kept where 3 survived before**, so 86% of the work was being discarded.

That made the batch id a *string*, which is why
[`crate::batch_merge::batch_sort_key`] exists: lexicographically `"10" < "3"`,
which would reorder batches and change which of two mergeable isoforms survives.

#### Finding 30 was the judgement call, and it is the largest win of the four

It was recorded as *not* obviously wrong --- counting mismatch *runs* rather than
bases is a coherent similarity measure --- so fixing it changes the method, not
just correctness. Fixed anyway under the new policy, and it dominates everything
else here.

The direction is the opposite of the obvious guess. Base-counting yields a
*smaller* accumulator (a run longer than `delta_len` has already returned `False`,
so a run contributes at most `delta_len`), and the test is
`miss_match_length <= max_bp_diff`. Smaller accumulator, more permissive, more
merging. On 56 Drosophila clusters:

| `--iso_abundance` | finding 30 present | fixed |
| --- | --- | --- |
| 1 (no filter) | 930 isoforms | **265** |
| 5 (the default) | 16 isoforms | **42** |

Total reads represented is 1 192 either way, so nothing is lost or double-counted;
maximum support rises 50 -> 136. The raw isoform count falls **3.5x** while the
number surviving `--iso_abundance` rises **2.6x**: bigger groups carry enough
support to clear the cutoff, so the count going *up* means fewer isoforms are now
discarded for lack of support.

On the full `droso` corpus, and this is the check that says the extra isoforms are
correct rather than merely more numerous:

| | isoforms | FSM | FSM rate | canonical | chains | genes | median error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| four bugs still in | 507 | 447 | 88% | 0.983 | 357 | 428 | 0.0132 |
| **all six fixed** | **533** | **470** | **88%** | **0.983** | **374** | **450** | 0.0136 |

**+23 full splice matches across +22 genes, at an unchanged 88% FSM rate and an
unchanged 0.983 canonical splice fraction**, with distinct intron chains rising
357 -> 374. So the 26 extra isoforms are overwhelmingly correctly-reconstructed
known transcripts with genuinely distinct splice structures, not admitted noise.
Median error rises 0.0132 -> 0.0136, which is the one number that moves the wrong
way and is 0.4% relative.

SIRV barely notices it: `sirv_sim_gene` is unchanged and `sirv_real` moves by one
isoform. That is consistent rather than puzzling --- SIRV consensuses carry few
interior mismatches, and run-versus-base counting only diverges where there are
many, which is what real cDNA has.

The lesson for the method: this was the finding most likely to be left alone, on
the correct observation that its behaviour was defensible. Leaving it alone would
have cost 23 recovered transcripts. "Defensible" is not the same as "as good as
the alternative", and the only way to tell them apart is to measure.


### The merge threshold is already optimal, and both directions are worse

`--delta_len` (default 5) sets the internal-indel run length above which a pair of
consensuses is "structural" and refuses to merge. Profiling said 82% of all
rejected pairs are rejected by that test, so it looked like the knob controlling
redundancy. Swept properly, it is not a knob to turn --- it is already at its
optimum.

**First, the confound had to go.** `--delta_len` drives three separate things: the
graph's edge-length tolerance, simplification, and the merge test. Sweeping the
flag measures all three at once. `ISONFORM_MERGE_DELTA_LEN` isolates the third; it
defaults to `--delta_len`, so nothing changes unless it is set. Verified to
isolate: on one 1 000-read batch, `--delta_len 5` gives 64 isoforms, `--delta_len
20` gives 28, and `--delta_len 5` with the merge-only threshold at 20 gives **25**
--- distinct from both, and it shows the merge test accounts for essentially all
of the isoform reduction.

**Then the sweep, on `sirv_real` where 68 transcripts are known:**

| merge threshold | isoforms | strict recall | strict F1 | redundancy | lenient recall | lenient F1 |
| --- | --- | --- | --- | --- | --- | --- |
| 2 | 153 | 66.2% | 0.734 | 2.80 | 75.0% | 0.852 |
| 3 | 106 | 66.2% | 0.725 | 1.89 | 77.9% | 0.869 |
| **5 (default)** | **85** | **70.6%** | **0.734** | **1.35** | **82.4%** | **0.873** |
| 10 | 73 | 61.8% | 0.656 | 1.21 | 79.4% | 0.869 |
| 20 | 67 | 58.8% | 0.640 | 1.18 | 76.5% | 0.831 |
| 50 | 63 | 51.5% | 0.575 | 1.17 | 72.1% | 0.802 |

Recall peaks at 5 on both thresholds, lenient F1 peaks at 5, and 5 dominates every
alternative. The two sides fail for different reasons, which is what makes this a
real optimum rather than a coincidence:

* **above 5** --- distinct transcripts are collapsed. Recall falls 48 -> 35
  recovered transcripts by 50, and lenient recall falls with it, so these are
  genuine losses rather than threshold crossings;
* **below 5** --- merging is too timid. Redundancy blows up to 2.80 and recall
  *also* falls, to 66.2%. Counter-intuitive until the mechanism is clear: merging
  pools reads, so a merged consensus is more accurate than its fragments.
  Under-merging leaves poorer consensuses that fail the identity test outright.

So merging buys accuracy until it starts destroying distinctions, and 5 is where
those two curves cross. The shipped default is right, and now has evidence behind
it rather than assumption.

#### Two mistakes this sweep caught, both mine

**Reading a count as a metric.** On a Drosophila batch, raising the threshold took
isoform count 64 -> 25 and I called it a redundancy improvement. Drosophila has no
per-isoform truth, so "collapsed duplicates" and "destroyed distinct transcripts"
look identical there --- and with truth available it is mostly the latter. That is
exactly what method point 4 exists to prevent, walked into by reading a number
with no denominator.

**Reading a distribution as a conclusion.** The run-length profile is sound data:
of 5 017 structural rejects, 2 085 have runs of 6-10 bases, 2 265 of 11-20, 654 of
21-50, and only 13 above 50, with an empty band at 101-200. I read "artefact
scale, therefore safe to merge". The sweep says otherwise: those short indels do
separate genuine isoforms. The histogram was right; the inference from it was not.

#### What this means for the merge cost

Raising the threshold made the cross-batch merge 3.8x faster with 19x fewer
alignments --- by merging more. Since merging more costs recall, that speed is not
available. The cliff from finding 31 (8 hours on one 26-batch cluster) has to be
solved *without* changing which pairs merge, which leaves the vectorised DP and a
candidate filter with a **bounded** false-negative rate. A minhash sketch was
built and then rejected on exactly that ground: it gives no window guarantee, and
a filter that decides which pairs are never examined cannot rest on a rate
measured on one corpus.


### The aligners, measured: parasail vs block-aligner vs WFA2

Both hot alignment sites call `parasail.sg_trace_scan_16` --- `align_to_merge`
(`IsoformGeneration.py:379`, mismatch −2) and `align_bubble_nodes`
(`SimplifyGraph.py:657`, mismatch −8). On a normal 1 000-read batch the merge is
**92% of the wall clock** (`merge=54.77s` of `total=59.62s`), so this is where
runtime lives.

Timed over the 54 884 recorded real calls, all three producing traceback, bucketed
by the identity of the recorded parasail CIGAR:

| bucket | cases | parasail | block-aligner | WFA2 | blk | wfa |
|---|---|---|---|---|---|---|
| <60% | 5 093 | 3.01s | 0.79s | 2.30s | 3.8x | 1.3x |
| 60--80% | 13 995 | 9.19s | 3.02s | 2.28s | 3.0x | 4.0x |
| 80--90% | 32 318 | 18.99s | 6.59s | 2.60s | 2.9x | 7.3x |
| 90--95% | 3 208 | 1.54s | 0.55s | 0.13s | 2.8x | 12.0x |
| >=95% | 270 | 0.03s | 0.01s | 0.00s | 3.3x | 12.9x |
| **total** | | **32.75s** | **10.95s** | **7.31s** | **3.0x** | **4.5x** |

Block-aligner is a uniform 3.0x. WFA2 is 4.5x overall, with the divergence
dependence its complexity predicts (O(n·s) in the alignment penalty): 1.3x on the
most divergent pairs, 12.9x on the most similar, and never slower than parasail.

Two corrections to things asserted earlier in this file's exploratory sections.
The merge's cost is **not** dominated by low-identity pairs --- the 80--90% bucket
holds 59% of the calls and 58% of the time, while `<60%` is 9%. And the `>=95%`
bucket that the block-aligner fidelity work concentrated on is 270 calls and
**0.1% of runtime**. Block-aligner's speed was never measured during that work,
which was the wrong order: the 3.0x was sitting there unrecorded while the effort
went into scores that could not move the total.

### Finding 40 --- WFA2's ends-free mode is not parasail's semi-global

parasail `sg` **maximises a score** with a positive match reward; WFA2 **minimises
a penalty** with `match = 0`. With all four ends free the empty alignment costs
nothing, so it is always optimal --- on two *identical* 20bp sequences WFA2
returns `DDDD…IIII`, aligning nothing. Bounding the free ends removes that, but
WFA2 still pays nothing for terminal *matches* and so declines them: on two
identical 200bp sequences with a 50-base budget it shaved 48 bases off the ends.

`src/wfa.rs` reconciles this in three steps --- bound the free ends, greedily
extend the aligned core outward over matching base pairs, then score the resulting
columns with **parasail's own rules** rather than transforming WFA2's penalty back.
That last step is what makes the arithmetic unfalsifiable-by-drift: WFA2 only
chooses the columns, `Scoring` prices them.

Two errors this shape caught, both mine:

* Scoring *every* gap run outside the aligned core as free made the layer score
  **above the exact DP on 51 of the recorded calls**. parasail's `sg` zeroes row 0
  and column 0 only, so a path starts at `(i,0)` *or* `(0,j)`: one sequence gets a
  free prefix, not both. A left dangle of `10I 5D` pays for the `5D`. The
  "must never beat exact" assertion is what surfaced it.
* A crafted test used a period-12 fixture, where a 48-base offset alignment is
  still all-matches and zero-penalty. That tested the fixture, not the aligner.
  Non-periodic fixtures now, with the periodic case kept deliberately as
  `a_periodic_sequence_admits_a_co_optimal_offset` --- greedy extension cannot
  undo a whole-alignment offset, because such a dangle is all `I` with no `D` to
  pair against.

### Finding 41 --- gate an aligner swap on merge verdicts, not on scores

`align_to_merge` returns a **bool**, and that bool is what collapses two isoforms.
Two co-optimal alignments can score identically and still decide differently,
while two differently-scored alignments can agree. So the gate is
`wfa::oracle::verdicts_match_parasail`, which drives the real `align_to_merge`
through both aligners and compares the decisions. Over the 54 884 recorded calls:

```text
agree 53 684   disagree 1 200   (2.19% of calls, 5.09% of parasail's 23 586 merges)
  of the disagreements WFA2 scored:  worse 902   equal 298   better 0
  worst delta -2457, mean -60.68
```

75% of the disagreements are WFA2 finding a genuinely **worse** alignment, not a
tie-break difference, and the mechanism follows from the objectives: a terminal
stretch of ten columns holding one mismatch earns parasail `9x2 - 2 = +16`, so it
aligns it, while WFA2 sees penalty 4 against 0 for leaving it free and declines.
Greedy extension then halts at that mismatch. The fix is an exact DP over the
`<=50`-base dangles, which is cheap. **Written and measured since, and this
prediction was wrong: it makes `sirv_real` worse. See finding 50.**

WFA2 is wired at both hot sites behind `ISONFORM_WFA2=1`, with a parasail
fallback whenever the layer declines a pair.

#### What the flipped verdicts cost end to end

| corpus | parasail | WFA2 | speedup |
|---|---|---|---|
| droso | 24.9s | 16.8s | **1.48x** |
| sirv_sim_gene | 107.8s | 101.2s | 1.07x |
| sirv_real | 220.1s | 108.7s | **2.03x** |

| metric | parasail | WFA2 | python |
|---|---|---|---|
| droso FSM | 470 | 465 | 443 |
| droso genes / canonical | 450 / 0.983 | 447 / 0.980 | 426 / 0.983 |
| sirv_sim recall / lenient F1 | 89.7% / 0.946 | 89.7% / 0.946 | 91.2% / 0.954 |
| sirv_sim redundancy | 1.38 | 1.36 | 2.19 |
| sirv_real recall (lenient) | 82.4% | 80.9% | 82.4% |
| sirv_real strict F1 | 0.734 | 0.730 | 0.773 |

So 5.09% of merge verdicts flipping costs **5 FSM on droso, about one transcript
on sirv_real, and nothing measurable on sirv_sim_gene**, for between 1.07x and
2.03x wall clock. Both corpora remain well ahead of the reference.

The end-to-end gain is far below the 4.5x on the alignment calls, and unevenly
so. Only `sirv_real` and `droso` are alignment-bound; `sirv_sim_gene` moves 6%.
The "merge is 92% of wall clock" figure came from a single 1 000-read batch and
plainly does not generalise --- where `sirv_sim_gene` spends its time is **not yet
measured**, and may matter more than the aligner does.

Left **off by default**: the precedent is `ISONFORM_BUBBLE_GLOBAL`, rejected for
costing 2 sirv_real transcripts. This costs 1 but buys 2x, so the trade is not the
same one --- and the 902 suboptimal alignments have a known, unwritten fix (an
exact DP over the `<=50`-base dangles) that would likely recover most of the loss.
Turning it on by default is a decision to take after that, not before.

**Superseded on both counts.** WFA2 went on by default anyway (finding 49), and
the DP was written and scored: it recovers none of the loss and costs more
(finding 50). "Would likely recover most of the loss" is one more entry for the
warning at the top of this file.

### Finding 42 --- the bottleneck is POA, not pairwise alignment

Profiling `sirv_sim_gene` (point 3 of the plan) because a 4.5x aligner bought
only 6% there. The corpus has 7 clusters and the driver gives each its own thread
(`n_threads = min(--t, instances)`), so **wall clock is the slowest cluster, not
the sum**. Cluster 0 is the largest (2 242 reads, 3 sequential 1 000-read
batches), and its first batch is the critical path.

The per-cluster `stderr.txt` that would have carried `ISONFORM_STAGE_TIMES` is
deleted by `remove_folders` at the end of a parallel run, so this was measured by
splitting cluster 0 exactly as the driver does (`max_seqs = 1000`, file order) and
running `main` on each batch directly.

```text
WFA2 off  intervals=2.30 graph=0.25 simplify=1.88 grouping=0.05 merge=69.59 total=74.10
          merge-split(poa=57.30s/278 calls  align=12.28s/2119 calls  other=0.01s)
WFA2 on   intervals=2.31 graph=0.27 simplify=3.20 grouping=0.05 merge=55.16 total=61.02
          merge-split(poa=49.99s/252 calls  align= 5.17s/1854 calls  other=0.01s)
```

The merge stage is 94% of the batch --- so "the merge dominates" was right. But
inside it, **POA is 82% of the stage and 77% of the whole batch**, at ~206ms per
call, while the pairwise alignment WFA2 replaces is 16.6% at ~5.8ms per call.

That is a hard Amdahl ceiling: an *infinitely fast* pairwise aligner removes
12.28s of 74.10s, **16.6%**. WFA2 takes 2.4x on that slice (12.28s -> 5.17s),
worth 9.6% of the batch. The remaining apparent gain is not a speedup --- POA fell
57.30s -> 49.99s only because 26 fewer POA calls happened (278 -> 252) as a side
effect of the flipped merge verdicts.

#### Where 278 POA calls come from, for 19 groups

`generate_all_consensuses` runs spoa once per multi-read group. Then **every
successful merge** calls `generate_new_full_consensus`, which rebuilds a full POA
from the union of both groups' raw reads (up to `max_seqs_to_spoa = 200`), guarded
only by `if groups[pos2].1.len() > 50` --- a check on the target group's read
count, not on the size of the union. So roughly 259 of the 278 calls are
merge-triggered rebuilds.

`poa::consensus` already delegates to `spoars`' `SimdEngine`, so the lever is the
**number of calls**, not the cost of one.

#### Per-call POA speed is a closed door --- isONcorrect already proved it

Do not reopen the library question. isONcorrect hit the same bottleneck and
settled it on measurement (`isONcorrect/PORTING.md`, "candidate" table): `spoars`
2.83s and byte-identical on 2 682 real intervals; `abpoa-rs` **calls `exit(1)`
from inside the C library** in both local and global mode; `poa-consensus` 32.7s
(**11.6x slower**) and agreeing on only 62%; `poasta` has no consensus API at all.
That file's conclusion is verbatim that `run_spoa`'s 20% is "not addressable by
swapping libraries", and that a native POA would mean beating both spoa and abPOA
at their own game. Hoisting the engine into a thread-local was also measured and
worth nothing: the time is in the DP, and `spoars`' DP is already vectorised with
real NEON kernels.

#### But the two ports have differently shaped POA problems

| | isONcorrect | isONform |
|---|---|---|
| call sites | 1 (`corrections.rs:112`) | 2 --- per group, **plus one per merge** |
| calls | 11 775 | 278, for 19 groups |
| time per call | **1.80ms** | **206ms** |
| POA input | ~107 seqs x **61bp** intervals | `<=200` seqs x **full-length reads** |

isONcorrect POAs short correction intervals, so its per-call cost was already at
the floor and the library was the only conceivable lever --- which is why that
investigation ended where it did. isONform POAs whole reads, 114x dearer per call,
and calls POA 14x more often than it has groups.

So isONform has two levers isONcorrect never had, and **both are existing
reference parameters rather than method changes**:

1. the `> 50` threshold, which decides how many of the ~259 merge-triggered
   rebuilds happen at all;
2. `max_seqs_to_spoa = 200`, since POA cost grows with the sequence count.

Both sweep like `delta_len` did --- runtime against accuracy, scored rather than
guessed --- and both attack the 77% instead of the 16.6%.

This also corrects the reasoning that sent the last stretch of work at the
aligner. The "merge is 92% of wall clock" figure was right; the inference that the
merge's time was therefore *pairwise alignment* was never checked, and it was
wrong by a factor of five.

### Finding 43 --- skipping the merge rebuild is 1.8x, and costs consensus accuracy not structure

Sweeping `ISONFORM_MERGE_REBUILD_MAX` on `sirv_sim_gene`'s critical-path batch:

| threshold | total | poa | poa calls | align |
|---|---|---|---|---|
| 50 (reference) | 73.31s | 56.66s | 278 | 12.20s |
| 25 | 71.73s | 54.92s | 272 | 12.28s |
| 10 | 69.05s | 52.29s | 261 | 12.20s |
| **0 (never rebuild)** | **20.93s** | **3.53s** | **80** | 12.88s |

The knob is effectively **binary**. 50 -> 25 -> 10 barely moves anything, because
the reference tests the *target group's* read count and most groups are small, so
the skip almost never fires until the threshold reaches 0. The 80 residual calls
are `generate_all_consensuses`: this batch has ~80 groups before merging and the
19 that `stage-times` reports are what survive it.

End to end, with parasail (WFA2 off):

| corpus | parasail | rebuild=0 | speedup | accuracy |
|---|---|---|---|---|
| droso | 24.9s | 18.4s | 1.35x | FSM 470 -> **472**, canonical 0.983 -> **0.985** |
| sirv_sim_gene | 107.8s | **58.5s** | **1.84x** | F1 **identical** (0.940 / 0.946), redundancy 1.38 -> **1.31** |
| sirv_real | 220.1s | 128.3s | **1.72x** | strict F1 0.734 -> **0.660** (recall 70.6% -> 63.2%), lenient F1 0.873 -> 0.871 |

Two corpora come out neutral or better --- droso's FSM 472 is the best figure any
configuration has produced, against the reference's 443. But `sirv_real` loses 5
of its 68 transcripts at **strict** while holding at **lenient** (0.873 -> 0.871),
and that split is the diagnosis rather than noise: the loss is **consensus base
accuracy, not isoform structure**. `sirv_sim_gene`'s identity column shows the
same thing directly, 0.9968 -> 0.9927. Skipping the rebuild leaves a merged group
carrying a consensus built from one group's reads instead of the union of both, so
on noisy real ONT data a few transcripts fall below strict identity while staying
recognisable at lenient.

So `rebuild=0` is **not shippable on its own**, and the interesting consequence is
what it implies about the next knob: what the rebuild buys is base accuracy, so
the thing to try is a rebuild that is *cheaper* rather than *absent* ---
`max_seqs_to_spoa`, which scales the POA's input count while still rebuilding.
Note the reference's own paper parameters for isONcorrect use
`--max_seqs_to_spoa 3` or `5` against isONform's default of **200**, though the
two are POAing different things (61bp correction intervals against full-length
reads) so that is a hint, not a target.

Default unchanged at 50.

### Finding 44 --- rebuild the consensus once per group at the end, not once per merge

Finding 43 left a puzzle: switching the per-merge rebuild off was 1.8x and kept
the isoform *structure*, but cost *base accuracy*. The resolution is that the
rebuild was serving two unrelated purposes --- supplying the consensus that later
merge decisions compare against, and supplying the accuracy of the emitted
sequence. Only the second one needs the union of both read sets, and it only needs
it **once**.

So: merge using the consensuses already in hand, then POA each surviving group
once over the union of its reads. That is `generate_all_consensuses` a second time
on the final membership --- `ISONFORM_FINAL_CONSENSUS=1`, which is only worth
anything alongside `ISONFORM_MERGE_REBUILD_MAX=0`.

On `sirv_sim_gene`'s critical-path batch:

| config | total | poa | poa calls | align |
|---|---|---|---|---|
| reference (rebuild=50) | 72.67s | 55.85s | 278 | 12.20s |
| rebuild=0 | 20.81s | 3.50s | 80 | 12.76s |
| **rebuild=0 + final pass** | **25.19s** | 7.89s | 93 | 12.78s |
| rebuild=0 + final + WFA2 | **18.82s** | 8.22s | 92 | 4.64s |

13 extra POA calls, 4.4s, and 2.9x against the reference. End to end:

| corpus | parasail | final pass | speedup | accuracy |
|---|---|---|---|---|
| droso | 24.9s | 18.6s | 1.34x | FSM 470 -> **476**, ISM 12 -> **9** |
| sirv_sim_gene | 107.8s | 62.3s | 1.73x | strict F1 0.940 -> **0.946**, redundancy 1.38 -> **1.26**, identity restored 0.9927 -> **0.9966** |
| sirv_real | 220.1s | 141.7s | 1.55x | strict F1 0.734 -> **0.735**, lenient F1 0.873 -> **0.884**, recall unchanged at 70.6% / 82.4% |

**Better on all three corpora and faster on all three.** droso's FSM 476 is the
highest any configuration has produced (reference: 443), `sirv_sim_gene`'s
redundancy 1.26 the lowest, and the 5 sirv_real transcripts that `rebuild=0` lost
are all back. This is the one change in this file that costs nothing to buy.

Stacked with WFA2 it reaches **2.35x / 3.62x / 3.60x**, but WFA2 still costs
sirv_real accuracy (strict F1 0.735 -> 0.711, lenient 0.884 -> 0.865) --- its own
verdict flips from finding 41, not the final pass. Note also that the aligner only
becomes worth having *after* POA is fixed: alignment is 25% of the FC batch
against 17% of the reference one.

#### On by default, and it is a deviation

Both defaults are flipped on: `MergeOpts::merge_rebuild_max` is **0** and
`MergeOpts::final_consensus_pass` is **true**. Verified end to end --- a run with
no environment set reproduces the measured configuration exactly, metric for
metric, on all three corpora. `ISONFORM_MERGE_REBUILD_MAX=50
ISONFORM_FINAL_CONSENSUS=0` restores the reference's schedule.

Exactly what differs: the reference's per-merge rebuild concatenates **id1's
reads then id2's**, while the final pass POAs the group's member list in its own
order (id2's originals, then id1's appended). spoa is order-sensitive, so the
emitted consensus differs even on an identical read set. Merge *decisions* also
differ, because they now compare un-rebuilt consensuses.

`ISONFORM_BUG_COMPAT=all` deliberately does **not** cover this --- it is a method
deviation, not a bug fix. Anything comparing the port against Python must pin the
reference's schedule itself, and three places now do: `tests/isoforms_oracle.rs`,
`bench/compare_end_to_end.py`'s `port_env`, and the `isoforms.rs` unit tests.

The scheduling lives in [`MergeOpts`] rather than in environment globals so that
tests and oracles are not at the mercy of the environment. One trap worth
recording: `#[derive(Default)]` would have produced `merge_rebuild_max: 0` **and**
`final_consensus_pass: false` --- the rebuild off with nothing replacing it, which
is precisely finding 43's accuracy-losing configuration --- and the isoform oracle
builds its options from `Default`, so it would have silently run that one. The
impl is written out by hand for this reason.

### Finding 45 --- the degenerate instance is poisoned by SPECIFIC READS, not by scale

**This heading was wrong until it was tested properly, and the correction matters.**
The scale sweep below walks a *prefix* --- reads 1..N --- which conflates "more
reads" with "reads 241+ happen to contain whichever reads cause this". Running
disjoint 250-read blocks instead separates the two:

| block of 250 | nodes | hit rate |
|---|---|---|
| 1--250 | 443 | 94.0% |
| **251--500** | **2 155** | **70.7%** |
| 501--750 | **209** | **97.2%** |
| 751--1000 | **241** | 97.1% |

Blocks 3 and 4 are the cleanest graphs measured anywhere in this file --- 209 and
241 nodes at 97% --- from the same reads that collectively produce 20 421. So the
degeneracy is **not** a function of batch size alone: at a fixed N=250 the
composition of the block decides the outcome, from 209 nodes to 2 155.

#### Read 426: poisons a 250-read batch, irrelevant at 1 000 --- and I got this wrong twice

A **fixed-N** bisect (125 clean base + suspect subset + clean filler, always 250
reads) converges monotonically on a single read:

```text
376-406  204   | 407-437  2959
407-422  198   | 423-437  3327
423-430  3326  | 431-437   201
423-426  3327  | 427-430   199
423-424  197   | 425-426  3383
425-425  202   | 426-426  3410     <- read 426 alone
```

Read 426 added to 249 clean reads: **3 410 nodes against a control of 201**. It is
a pathological read --- 1 592bp against a 1 335 median, **56% AT where the cluster
is 83%**, a **57-base homopolymer**, and one 20-mer present **38 times** where
typical reads have 0--6 duplicated 20-mers at all.

**And removing it from the full 1 000 reads changes nothing: 20 440 nodes at 30.6%,
against 20 421 at 30.7% with it.**

So read 426 is a genuine poisoner of a 250-read population and *causally irrelevant*
at 1 000. Both of the confident claims made about it in this session were wrong in
opposite directions: first "the bisect found it, but the single-read control at
N=101 shows nothing, so it is an artifact" (the control was simply too small ---
the effect needs N >~ 250), then "one read causes the degeneracy" (asserted before
running the full-scale removal, which refutes it). The lesson is the one this file
keeps recording: state the mechanism only after the measurement that would falsify
it.

#### The effect is SUPER-ADDITIVE, and that is where the investigation stops

| population | reads | nodes | hit rate |
|---|---|---|---|
| 501--1000 | 500 | **327** | 98.1% |
| **1--500** | 500 | **7 054** | 52.2% |
| 501--1000 + 251--500 | 750 | 2 382 | 89.5% |
| full file | 1 000 | **20 421** | 30.7% |

500 reads build a 327-node graph; the same file's other 500 build 7 054; together
all 1 000 build 20 421. The whole is far worse than any part, and no part explains
it. Mixing a degenerate population into a clean one degrades it *partially* (2 382
at 750 reads) rather than either dominating or averaging.

And the halves are not two different transcripts: **999 of 1 000 reads are
ambiguous** between k-mer profiles built from each half. The *healthy* set's halves
are much more distinct from one another (profile Jaccard 0.132 against the
fragmented set's 0.607) and still build a clean graph --- inverted, like every other
discriminator tried here.

#### Leave-one-out at fixed N: the poison IS localised

Baseline block 251--500 is 2 155 nodes at 67.6%. Removing one 25-read chunk at a
time, every run at N=225 so they are comparable to each other:

| removed | nodes | hit rate |
|---|---|---|
| 251--275 ... 401--425 (7 chunks) | 2 130--2 145 | 67.6% |
| **426--450** | **182** | **97.4%** |
| 451--475 | 1 509 | 77.2% |
| 476--500 | 1 487 | 77.6% |

Removing 25 reads takes the block from 2 155 nodes to **182** --- cleaner than any
250-read block measured anywhere. Seven of the ten chunks do nothing at all.

Narrowing inside that chunk lands on one read:

| removed from block 251--500 | nodes | hit rate |
|---|---|---|
| nothing (baseline) | **2 155** | 67.6% |
| reads 426--437 (12) | 184 | 97.5% |
| reads 426--430 (5) | 187 | 97.6% |
| **read 426 only (1)** | **191** | **97.6%** |
| reads 438--450 (13) | 1 800 | 74.2% |

**A single read inflates a 250-read graph 11x**, 2 155 nodes to 191, and the hit
rate goes 67.6% -> 97.6%. Independently confirmed from the other direction: read
426 added *to* a clean 249-read population gives 3 410 nodes against a control of
201.

#### And at full scale the cause is concentrated too --- but in different reads

Removing one 50-read chunk at a time from the whole file, N=950 throughout:

| removed | nodes | hit rate |
|---|---|---|
| 1--50 ... 151--200 (4 chunks) | 20 367--20 401 | 27.1% |
| **201--250** | **15 586** | **44.4%** |
| 251--300 ... 451--500 (5 chunks) | 19 066--19 141 | ~31.8% |

One chunk removes 4 835 nodes and lifts the hit rate from 30.7% to 44.4%; five
give a small uniform ~1 300-node improvement each; four do nothing at all.

Narrowing that chunk, still at full scale:

| full 1 000 minus | nodes | hit rate |
|---|---|---|
| 201--225 | 20 407 | 28.9% |
| **226--250** | **15 594** | **45.9%** |
| 201--212 | 20 414 | 29.9% |
| 213--225 | 20 416 | 29.8% |

So the effect is concentrated at both scales --- but **in different reads**. Read
426 dominates at N=250 and is irrelevant at N=1 000, where reads 201--250 dominate
instead. Which reads poison depends on which other reads are present, which is
what makes this a property of the population rather than of a read, and why the
super-additivity above is not a paradox: poisoning is real, findable and localised,
but not additive.

#### The practical consequence: a diagnostic, not a fix

Leave-one-out at fixed N reliably finds poisoners. That is genuinely useful for
*analysis*. It cannot become a prefilter, for two independent reasons:

1. identifying a poisoner requires building the graph --- the expensive thing a
   filter would exist to protect; and
2. the culprits change with batch composition, so a list learned on one batch does
   not transfer.

Together with the composition-screen result below, that leaves a scale-invariant
node assignment (finding 46's chaining, or something like it) as the only
intervention the evidence supports.

#### A composition prefilter would not work

Read 426's profile suggests screening reads on base composition, homopolymer
length and duplicated k-mer count. Scanning all 1 000 reads for that profile
(`AT<75%`, or `maxhp>=30`, or `>=20` duplicated 20-mers) flags **127 reads, split
61/66 between the two halves** --- and one of those halves builds a 327-node graph.

The single worst read by duplicated 20-mers is **read 981: 182 duplicates, five
times read 426's 37 --- and it sits in the clean half**. Reads 134 and 459 have
profiles nearly identical to 426 (max homopolymer 57 and 56, 37 and 36 duplicated
20-mers, top 20-mer at 38 and 37 copies) and cause no trouble at all.

So read 426's pathological composition is not what makes it a poisoner, and a
composition screen would discard 127 reads to no benefit. Poisoners can only be
identified by measurement --- leave-one-out at fixed N.

**The root cause is not identified, and this is the honest end state.** What is
established is a long list of what it is not --- scale, any single read, repeats,
cycles, pair disagreement, tiling phase, transcript mixture, the anchor filters,
and every input-side statistic measured --- plus one sharp positive fact: the failure
is **all-or-nothing per read**, with misses uniformly offset by a few bases along
the whole read. A read joins the graph completely or duplicates it completely.

The practical consequence is unchanged and does not depend on knowing why:
node assignment by exact-coordinate dict lookup is fragile in a way that
population-independent assignment would not be, chaining recovers 6x of it
(finding 46), and the `sirv_real` strict metric is the gate any such change must
pass.

#### A bisect that looked conclusive and was not --- worth recording as a method trap

Bisecting *inside* the bad block appeared to localise the cause to reads 424--427,
and read 426 is a striking outlier: length 1 592 against a ~1 335 median, **56% AT
where the cluster is 83%**, a **57-base homopolymer run**, and one 20-mer present
**38 times**. A tidy story --- and false.

Adding read 426, or any of 424--427, to a clean 100-read baseline changes nothing:
130 nodes -> 130, hit rate 95.5% either way.

The bisect ranked halves by nodes-per-read, and that statistic is dominated by N at
small sizes --- *healthy* reads at N=5 give 17.6 nodes/read. As the bisect shrank
the blocks every subset looked bad, so it followed noise. **Comparisons here are
only valid at fixed N.** The equal-size block table above is sound; the bisect
inside it was not, and read 426's odd composition is real but causally irrelevant.

That also explains the earlier all-or-nothing observation without needing a
population-convergence story: reads 241--250 **on their own** build 41 nodes at an
86.5% hit rate, so they are perfectly capable of sharing with each other. They stop
only in the presence of the wider batch.

The original framing is kept below because the measurements are sound and the
prefix sweep is still the cleanest demonstration that the *reads* are not
individually defective --- it is what that sweep was taken to prove about scale
that was wrong.

### Finding 46 --- co-linear chaining, and what it does to the degenerate instance

Finding 45 established that the degenerate graph is built entirely from the
first-sight branch: no cycles, no within-read repeats, reads simply never agreeing
on interval coordinates while 94% of the misses have a registered coordinate
within 20 bases. The exact lookup asks *"did somebody register this precise
coordinate"*; what it should ask is *"does this read's sequence of anchors follow a
path already in the graph"*.

#### The formulation is one-dimensional, which is what makes it cheap

No graph topology is needed. Registered coordinates for read `r` are positions
**in read `r`**, and so are `r`'s own interval starts, so this is chaining on a
line. Candidates for interval `i` are nodes whose registered coordinate for `r`
lies within `tol`; a chain takes at most one per interval with registered
coordinates strictly increasing in interval order. The DP maximises the number of
chained intervals and, among those, minimises total offset --- roughly 30 intervals
with a handful of candidates each, so the quadratic scan costs nothing.

It consults the chain **only where the exact lookup missed**. Exact hits are
zero-offset candidates and so cannot be lost, which is precisely how finding 39's
tolerant lookup went wrong: it collapsed 20 421 -> 1 754 nodes but exact hits
*fell*, 9 027 -> 8 106, because a greedy nearest match attaches a read to whichever
copy is closest with no regard for order.

`ISONFORM_CHAIN=<tol>`, default 0 (off).

| config | nodes | exact hits | chained hits | first-sight |
|---|---|---|---|---|
| baseline | **20 421** | 9 027 | --- | 20 376 |
| chain=5 | 6 425 | 8 991 | 14 176 | 6 236 |
| chain=10 | **3 427** | 8 823 | 17 818 | 2 762 |

**6x fewer nodes, and the exact hits hold** --- 8 823 against tolerant lookup's
8 106 for a comparable collapse, so the monotonicity constraint is doing real work
rather than merely loosening the match. Still well above the healthy graph's 548,
so this is a large improvement and not a cure. Merge time is unchanged (~145s),
because that stage is POA-bound rather than node-bound (finding 42).

#### Scored on the corpora, which is the only validation that counts

Finding 39 said this can be judged only against known truth, not node counts.

| | default | chain=5 | chain=10 |
|---|---|---|---|
| degenerate instance nodes | 20 421 | 6 425 | **3 427** |
| droso FSM | 476 | 475 | **480** |
| droso genes / canonical | 443 / 0.983 | 446 / 0.982 | 446 / **0.985** |
| sirv_sim_gene F1 | 0.946 | 0.946 | 0.946 |
| sirv_sim redundancy | 1.26 | 1.28 | **1.25** |
| sirv_real **strict** F1 | **0.735** | 0.701 | 0.718 |
| sirv_real **lenient** F1 | 0.884 | **0.889** | 0.880 |
| runtime droso/sim/real | 19.4/63.2/128.8s | 18.8/**58.4**/**119.5**s | 18.4/70.4/132.5s |

A tradeoff, not a win: the degeneracy improves 6x and droso gains (+4 FSM, best
canonical at tol=10), `sirv_sim_gene` is unmoved, and **`sirv_real` strict F1
declines** by 0.017 (tol=10) to 0.034 (tol=5) while its lenient F1 is flat to
better. **Off by default.**

`sirv_real` strict is now the metric that three independent changes have all paid
into --- WFA2 (finding 41), `rebuild=0` (finding 43) and chaining. It is the only
corpus with real ONT reads scored against known transcripts by near-exact
identity, so it is the strictest gate available and worth treating as the binding
constraint on anything that moves node assignment or consensus content.

#### Spacing: inert as a cost, effective as a bound

Adding gap-difference to the DP's cost --- the obvious next increment --- turned out
to be **exactly inert**. Identical chains at weights 0, 1, 4 and 8, at every
tolerance from 10 to 80: same `chained_hits`, same `first_sight`, to the unit.

Two arithmetic reasons, both of which should have been foreseen before writing it:

* within a tolerance the skew `|d_i - d_j|` is bounded by `|d_i| + |d_j|`, the
  per-candidate offset cost the DP already minimises, so it is near-redundant;
* the objective maximises **count** first, so cost only ever breaks ties between
  chains of the same length.

As a **hard bound** it works, because rejecting a transition changes which chains
exist rather than which of the equal-length ones wins. `ISONFORM_CHAIN_SPACING` is
now the largest tolerated skew in bases, 0 meaning unbounded.

| tol | max_skew | first-sight nodes | chained | cycles |
|---|---|---|---|---|
| 40 | unbounded | **467** | 21 104 | 114 |
| 40 | 10 | 971 | 20 001 | 54 |
| 40 | 5 | 1 935 | 18 926 | 139 |
| 40 | 2 | 4 029 | 16 796 | 62 |
| 10 | 2 | 5 673 | 14 796 | 7 |

And widening the tolerance alone almost removes the degeneracy: **467 nodes at
tol=40 against the healthy graph's 394**, from a baseline of 20 376. But tol=40 is
a full interval width, so an interval may attach to the *adjacent tile* --- exactly
the mis-assignment the skew bound exists to prevent, and the reason accuracy
rather than node count has to pick the operating point.

### Finding 47 --- two reads cause it, but their COMPOSITION is not why

The heading first read "two composition-outlier reads cause the whole degeneracy".
Half right. Two reads do cause it, and they are composition outliers --- but the
composition is a **correlate, not the cause**, and a filter built on it would have
been the tailored heuristic it looked like.

The refutations, both measured after the filter was written:

* **The anchor database barely changes.** With the two reads: 12 714 surviving
  anchors, 1 699 264 occurrences. Without: 11 276 and 1 696 388. The 1 438 extra
  anchors carry exactly 2.0 occurrences each, so every one is dropped by
  `find_most_supported_span`'s `relevant.len() < 3` cut --- inert. The other 998
  reads' support values are untouched, which kills the "corrupted support weights"
  mechanism.
* **Position is what matters.** The same 1 000 reads, only the two moved:

| reads 240 + 426 placed | nodes | hit rate |
|---|---|---|
| at the **front** | **26 507** | 9.9% |
| in place (240, 426) | 20 421 | 30.7% |
| at the **end** | **500** | 98.7% |

Identical reads, identical composition, a 53x spread in node count. Earlier is
worse, which is the signature of `prior.entry(key).or_insert(node)` --- **first
writer wins**. A read early in the file registers coordinate keys for every later
read and claims them for its own nodes; the same read last claims nothing.

So the fragility is in **node identity**, not in the reads. See finding 48.

### Finding 47a --- what the two reads are, kept for the record

Removing **two reads of 1 000** takes the degenerate instance from 20 421 nodes at
a 30.7% node-sharing rate to **500 nodes at 98.7%** --- the healthy set is 548 at
98.7%. The instance was never degenerate; it was two reads.

| | length | AT% | effect |
|---|---|---|---|
| **read 240** | 1 628 | **55%** | removing it alone: 20 421 -> 15 602 nodes |
| **read 426** | 1 592 | **56%** | at N=250: 2 155 -> 191 nodes (11x) |
| the other 998 | ~1 340 | **83--84%** | --- |
| **both removed** | | | **20 421 -> 500, hit 30.7% -> 98.7%** |

Reads with `AT < 75%` in this file: **exactly those two.** The distribution is
bimodal with nothing between (min 55%, p5 83%, median 83%, max 84%), and both sit
in the degenerate half, none in the clean half. Each was found independently by
leave-one-out bisection with no reference to composition, so the two methods
converge.

#### Why every population statistic missed it

Two reads in a thousand move no aggregate: identity, indel rate, 5' offsets,
intervals per read, support, candidate counts, tiling density, coordinate
concentration and pair agreement are all averages over ~29 400 intervals, and two
outliers are invisible in them. That is why ~20 hypotheses died here and why the
answer only came from removal experiments. **When an aggregate cannot see the
cause, only intervention can.**

It also explains the puzzles: the "scale effect" (a prefix only includes read 240
from N=240 and read 426 from N=426), the super-additivity, and why which read
matters depends on batch composition --- at N=250 read 426 dominates its block,
at N=1 000 read 240 does.

#### Chaining is NOT stable against these reads --- the filter is

| on the degenerate file | nodes | hit rate | sirv_real strict F1 |
|---|---|---|---|
| as-is | 20 421 | 30.7% | 0.735 |
| chaining, `tol=10` | 3 427 | 44.6% | 0.718 |
| **drop reads 240 + 426** | **500** | **98.7%** | **0.735** (untouched) |

Chaining absorbs about 83% of the inflation and leaves **7x** more nodes than
simply dropping the two reads, while costing 0.017 of `sirv_real` strict F1. On an
already-clean population it is neutral (447 nodes at 97.5% against 500 at 98.7%),
so it does no harm --- but it does not make node assignment stable against a
drifting read, it only dilutes the consequence.

That is the answer to "which approach is stable": **neither chaining nor any
aligner change is**, because the cause is upstream of node assignment. Two reads
whose composition differs from their cluster shift the batch-wide support values
that every read's WIS tiling depends on, and a more robust *assignment* cannot undo
a corrupted *weighting*. Removing the outlier is the only intervention measured
that restores the graph completely, and it is also the cheapest.

#### The filter is composition relative to the cluster, not absolute

The healthy set has **29** reads below 75% AT (51--52%, 2 093--2 430bp) and is
clean, so absolute AT does not predict poisoning. Removing them makes it *worse*:
548 -> **1 005 nodes** (hit 98.7% -> 99.0%), i.e. they were contributing real
structure. They are a coherent sub-population --- 29 reads long enough to build
their own path --- whereas the fragmented file's two are isolated singletons among
998 reads of a different composition.

So a viable screen is **relative**: flag a read whose composition is a sharp
outlier against its own cluster's distribution *and* which has too few companions
to form a group. Here that flags exactly 2 of 1 000. It is cheap (base counting),
needs no graph, and would have prevented a 40x node inflation and the 472s
`simplify` stage that started this investigation. Not implemented; the threshold
and the minimum-group-size rule both need scoring across corpora first.

### Finding 48a --- WHY: the coordinate key is not unique, and `or_insert` drops the losers

Findings 45--47 established *what* (two reads) and left *why* open. This is why.

`add_prior_read_infos` claims keys with `prior.entry(key).or_insert(node)`, so the
**first** node to claim a `(read, start, stop)` triple keeps it and every later
claim by a *different* node is discarded in silence. Counting the discards:

| | claims discarded (different node) | first-sight nodes |
|---|---|---|
| **fragmented** | **7 134 107** | 20 376 |
| minus reads 240 + 426 | **3 292** | 370 |
| healthy | 1 904 | 394 |

**A 2 166x difference.** The coordinate key is not unique: one triple is wanted by
many distinct nodes, and only the first gets it. A read that resolves such a key
later lands on whichever node claimed it first --- not on the node for its own
interval --- so it cannot attach to the good structure earlier reads built, and
creates a fresh node instead.

That answers the objection that this ought to be impossible: later reads *should*
attach to earlier good nodes, and the reason they cannot is that the key they would
use has already been taken by something else. It is a genuine graph-construction
defect, not a property of the data.

It also explains the position law exactly. Claims are made in read order, so two
reads early in the file mis-claim keys for all 998 downstream reads (front: 26 507
nodes), the same two in the middle mis-claim for the reads after them (20 421), and
placed last they claim nothing anyone still needs (500).

And it explains why content keys fix it: a pair hash cannot be mis-claimed, because
the key means "this pair" rather than "this read's coordinates".

### Finding 48 --- identify nodes by CONTENT, and the instance stops being degenerate

`NodeKey::Interval { start, end, r_id }` is per-read coordinates, so it needs a
registration step: each new node writes an entry for every supporting read at that
read's own coordinates, via `or_insert`. That is three fragilities at once ---
coordinate drift, read order, and file position (finding 47's 53x spread).

`ISONFORM_PAIR_NODES=1` keys nodes on a hash of the interval's flanking minimizer
pair instead. The interval starts at `p1 + k` and stops at `p2`, so the pair is
recoverable from the read itself with no plumbing. No registration, no tolerance,
no first-writer race, and order- and position-independent by construction. The
floor is real: the degenerate file has only **332 distinct pairs** across 29 403
selected intervals.

| | coordinate-keyed | **pair-keyed** | chaining `tol=10` | drop the 2 reads |
|---|---|---|---|---|
| degenerate instance | 20 421 | **463** | 3 427 | 500 |
| healthy instance | 548 | 524 | --- | --- |
| droso FSM | 476 | **478** | 480 | --- |
| droso runtime | 19.4s | **17.7s** | 18.4s | --- |
| sirv_sim F1 / redund | 0.946 / 1.26 | 0.946 / **1.23** | 0.946 / 1.25 | --- |
| sirv_real strict F1 | **0.735** | 0.729 | 0.718 | --- |
| sirv_real lenient F1 | **0.884** | 0.880 | 0.880 | --- |
| sirv_real identity | 0.9719 | **0.9746** | --- | --- |
| sirv_real runtime | 128.8s | **112.8s** | 132.5s | --- |

**44x fewer nodes on the degenerate instance, better than dropping the two reads
that cause it**, droso marginally better, `sirv_sim_gene` identical with lower
redundancy, and faster on all three corpora --- for 0.006 of `sirv_real` strict F1.
Chaining costs 0.017 and is slower; the composition filter is a heuristic built on
a correlate. This is the one that should be adopted.

#### The span was tried in the key and measured worse

Pair-only keying merges intervals sharing a pair but differing in span, and the
reference does draw that distinction: **332 distinct pairs against 566 distinct
`(pair, span)`**, 32% of pairs appearing at more than one span, the widest range
21 bases --- above `delta_len = 5`. So `(m1, m2, span / delta_len)` looked like the
more faithful key. It is worse:

| | coordinate | **pair only** | pair + span |
|---|---|---|---|
| degenerate instance nodes | 20 421 | **463** | 544 |
| healthy instance nodes | 548 | 524 | --- |
| droso FSM | 476 | **478** | 477 |
| sirv_sim F1 / redundancy | 0.946 / 1.26 | 0.946 / **1.23** | 0.946 / 1.30 |
| sirv_real strict F1 | **0.735** | 0.729 | 0.732 |
| sirv_real lenient F1 | 0.884 | 0.880 | **0.889** |
| sirv_real identity | 0.9719 | **0.9746** | 0.9721 |
| runtime droso/sim/real | 19.4/63.2/128.8s | **17.7/61.3/112.8s** | 20.2/68.9/152.9s |

More nodes, higher redundancy, and **20--36% slower than the coordinate key**. The
span spread is nearly all drift the reference already tolerates (median 0, p90
4--5), so bucketing splits nodes that ought to merge; it bought only `sirv_real`
lenient F1. Dropped.

Exon structure is not at risk either way: an interval is a 40--80bp window, so a
skipped exon changes *which pairs a read selects* --- its path through the graph ---
not one pair's span. Node identity does not carry exon differences.

#### The reference has this bug too

Node counts before simplification, same inputs and parameters:

| | Python reference | Rust port |
|---|---|---|
| fragmented | **18 771** nodes, 19 851 edges | 20 421 |
| healthy | **503** nodes, 850 edges | 548 |
| ratio | **37x** | 37x |

Under `ISONFORM_BUG_COMPAT=all` (plus the finding-44 schedule pinned to the
reference's) the port is **exact**, node for node and edge for edge:

| | Python | port, bug-compat |
|---|---|---|
| fragmented | 18 771 / 19 851 | **18 771 / 19 851** |
| healthy | 503 / 850 | **503 / 850** |

So the degeneracy is the **reference's**, and the 8% gap in the default
configuration is the port's own bug fixes, not a divergence.
Fixing it by changing node identity is therefore a deliberate method change, not a
port correction, and belongs under the same rule as findings 44 and 46. (The
reference run then dies in `run_spoa`, which shells out to a `spoa` binary absent
from the scratch copy --- an environment artifact at the simplify stage, after the
count is taken.)

Off by default pending a decision, for two reasons worth checking first:

* **the graph oracles diff node *names* against reference dumps**, and these names
  are `p{hash}` rather than `{start}, {end}, {r_id}`. They construct their own
  options so they are unaffected while this is opt-in, but making it default needs
  them switched to the coordinate key explicitly, exactly as finding 44's schedule
  was handled;
* a repeat within a read now maps to the same node, so `cycle_added` fires where
  `is_repetitive` used to catch it. On these corpora that is 0--53 events and
  harmless, but it is a behaviour change, not a no-op.

The composition filter (`ISONFORM_OUTLIER_Z`, default 0 = off) stays in the tree
unused. It works on the one instance it was built from and rests on a correlate;
finding 47 is why it should not be trusted, and this finding is why it is not
needed.

### Finding 49 --- the current defaults, measured together

Findings 41, 44 and 48 are all on by default: WFA2 at the two hot alignment sites,
the consensus rebuilt once per group at the end rather than once per merge, and
node identity by minimizer pair. A plain run with no environment set is this
configuration.

| | Python | port as first written | **defaults now** |
|---|---|---|---|
| droso | 49.0s | 24.9s | **11.4s** |
| sirv_sim_gene | 354.8s | 107.8s | **30.0s** |
| sirv_real | 241.0s | 220.1s | **59.6s** |
| **deep-12 within-batch** | --- | **1 403s** | **179s** |
| sirv_sim F1 / redundancy | 0.947 / 2.19 | 0.940 / 1.38 | **0.946 / 1.25** |
| sirv_real strict F1 | 0.773 | 0.734 | 0.722 |
| sirv_real lenient F1 | 0.892 | 0.873 | 0.862 |
| sirv_real identity | 0.9697 | 0.9716 | **0.9742** |
| droso FSM | 443 | 470 | 466 (**90%**) |
| droso ISM / NIC | 14 / 29 | 12 / 30 | **9 / 26** |

**2.2--3.7x faster than the port, up to 11.8x faster than Python, and 7.8x on the
deep clusters** --- the 12 deepest `droso_deep` clusters (94 603 reads, 102 batches,
deepest 25 493) go from 1 403s to 179s of within-batch work.

The cost is `sirv_real`: strict F1 0.734 -> 0.722, lenient 0.873 -> 0.862. Down from
0.711 before the short-sequence fix below, so roughly half of the earlier
regression was a bug rather than a trade. What remains is WFA2's verdict flips,
whose fix --- an exact DP over the `<=50`-base dangles (finding 41) --- is identified
and unwritten. Everything else improves or holds: `sirv_sim_gene` beats the port on
both F1 and redundancy, droso holds 90% FSM with fewer fragments and fewer NIC, and
base identity on real SIRV is the best of the three.

Each piece backs out on its own: `ISONFORM_WFA2=0`, `ISONFORM_PAIR_NODES=0`,
`ISONFORM_MERGE_REBUILD_MAX=50 ISONFORM_FINAL_CONSENSUS=0`.

#### A bug that only appeared when WFA2 became the default

Below `2 * FREE_ENDS`, and on near-periodic sequence at any length, the bounded
ends-free mode returns a **shifted** alignment: two identical 20bp sequences gave
`4I16=4D` instead of `20=`. A shift costs nothing when the sequence repeats, and
greedy end-extension cannot undo a whole-alignment offset --- the shifted dangle is
all `I` with no `D` to pair against. Short bubbles were therefore aligning
degenerately in every WFA2 run before this, including the first measurement of
these defaults (`sirv_real` strict 0.711).

Guarded by coverage: if the aligned core covers under 90% of the shorter sequence
the layer **declines**, and the call site falls back to parasail. Declining is the
right move --- a wrong alignment is worse than a slower one --- and it is what keeps
`simplify`'s smoke test, whose fixture is `ACGTACGT...` and so exactly periodic,
answering `20=`.

### Finding 50 --- the dangle DP is written, scored, and worse; and 50 free-end bases is right after all

Finding 41 named one unwritten fix for WFA2's verdict flips and finding 49 carried
it forward as the thing that would recover `sirv_real`. It is written
(`anchored_dp` in `src/wfa.rs`), and it does the opposite.

The measurement is one binary with one variable, `ISONFORM_WFA2_DANGLE_DP`:

| corpus | metric | DP on | DP off (greedy) |
|---|---|---|---|
| `sirv_real` | **strict F1** | **0.693** | **0.722** |
| | transcripts recovered | 45 / 68 | **47 / 68** |
| | lenient F1 | **0.874** | 0.862 |
| | isoforms / redundancy | 81 / 1.31 | 86 / 1.38 |
| | identity | **0.9758** | 0.9742 |
| `sirv_sim_gene` | F1 / redundancy | 0.946 / 1.25 | 0.946 / 1.25 |
| `droso` | FSM / NNC / canonical | 466 / 12 / 0.982 | 466 / 12 / 0.982 |

It over-merges: five fewer isoforms, less redundancy, and two real SIRV
transcripts collapsed. Default is now **off**; `ISONFORM_WFA2_DANGLE_DP=1` opts in.
The DP itself is correct --- it never scores above the exact DP, and a new fixture
for finding 41's stated mechanism passes.

#### Why: there is nothing there to repair

Two diagnostics added to `wfa::oracle::verdicts_match_parasail`, over 48 486 droso
and 49 508 `sirv_real` calls sampled from the merge call site
(`IsoformGeneration.py:381`) of a full reference run:

```text
                          droso              sirv_real
                       DP on   DP off      DP on   DP off
verdict flips           1189    1137         347     350
worse-scoring cases    17413   17253       21704   21793
fell back to parasail   4416    4716        7378    7425
two-sided dangles          0     436           0     460
```

A dangle the DP can act on needs **both** sequences to contribute bases. Only
~2% of the losing cases have one --- 436 of 17 253 on droso, 460 of 21 793 on
`sirv_real` --- and the DP consumes all of them. The rest is one-sided overhang:
one sequence running past the other, nothing to pair against, and free for
parasail too.

And the repair backfires through the coverage guard. A better-covered alignment
trips it less, so WFA2 declines less (4 716 -> 4 416 on droso) and answers more
pairs itself instead of handing them to exact parasail --- pushing the flip rate
**up**, 2.35% -> 2.45%. **Declining is worth more than repairing.** That is the
lever this investigation actually found, and it is untuned.

#### And the free-end budget, which was argued rather than measured

The losing cases carry far more unaligned end material than `FREE_ENDS = 50`
admits --- median 56--77 bases, p90 140--195, max 4 395. parasail's semi-global
gives an unbounded free end; WFA2 is capped and pays for the remainder. This file
argued 50 was safe because `delta_iso_len_5` refuses such pairs anyway, so a wider
budget "could not change a verdict". `ISONFORM_WFA2_FREE_ENDS` makes that a
measurement:

| | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| droso flip rate | **2.35%** | 2.58% | 3.51% | 4.50% |
| `sirv_real` flip rate | **0.71%** | 1.23% | 1.06% | 0.96% |
| droso worse-scoring | 17 253 | 15 643 | 14 695 | **14 441** |
| `sirv_real` worse-scoring | 21 793 | 19 966 | 17 579 | **14 411** |
| droso parasail fallbacks | **4 716** | 6 276 | 7 046 | 7 396 |
| `sirv_real` parasail fallbacks | **7 425** | 10 953 | 13 436 | 16 604 |

Widening buys back the scores exactly as the overhang picture predicts, makes the
**verdicts worse anyway**, and costs speed through the fallback rate. 50 is the
best value on both corpora --- the right answer, reached for a reason that was not
the stated one. Finding 41's lesson a second time: score and verdict are different
objectives, and here they move in opposite directions.

#### What this leaves

The residual `sirv_real` cost of WFA2-by-default is not in the ends. It is in the
core WFA2 chooses, and no end repair reaches it. That is now the open question,
and unlike the previous two candidates it has no proposed fix behind it.

### Finding 51 --- `--max_seqs` never reached the work it names, and bigger batches cost memory

`isONform_parallel --max_seqs` only ever changed how cluster *files* are split.
The reference comments the flag out of the child argument list, so the per-batch
size stayed at `main`'s own default of 1 000 whatever was passed --- finding 37's
class, found by sweeping it and getting byte-identical merge statistics at 1 000,
2 000 and 4 000 on three clusters. `ISONFORM_FORWARD_MAX_SEQS=1` forwards it.

The premise behind wanting it is sound. Isoforms per batch grow **sublinearly**
with reads, so fewer, larger batches leave less for the quadratic cross-batch
stage. `droso_deep` cluster 3 (9 261 reads), driving `main` directly:

| max_seqs | batches | isoforms | per batch | cross-batch pairs | within-batch |
|---|---|---|---|---|---|
| 1 000 | 10 | 386 | 38.6 | 65 700 | 28.6s |
| 2 000 | 5 | 343 | 68.6 | 45 805 | 43.2s |
| 4 000 | 3 | 286 | 95.3 | 24 880 | 63.8s |
| 10 000 | 1 | 215 | 215.0 | **0** | 186.6s |

Four times the batch is 2.6x fewer cross-batch pairs for 2.2x more within-batch
work. But end to end on the deep corpora the trade does not pay:

| | | 1 000 | 10 000 |
|---|---|---|---|
| `sirv_real_deep` | wall / peak RSS | 1 859s / **1.8 GB** | 2 792s / **20.7 GB** |
| | recall / redundancy | **82.4%** / 3.48 | 79.4% / **3.06** |
| `droso_deep` | wall / peak RSS | 5 678s / **3.1 GB** | **4 380s** / 13.6 GB |
| | FSM / NNC | 9 339 / 934 | 9 098 / **920** |

Less redundancy and fewer fragments, slightly fewer real transcripts, and **4--11x
the memory**. Wall clock improves only where cross-batch dominated. 10 000 is not
a viable default; the flag is worth forwarding, the value is not worth raising.

### Finding 52 --- merge in the graph, not pairwise: 1.2--1.4x faster and better

The cross-batch merge proves non-mergeability one alignment at a time. Running
isONform's own front end over a cluster's consensuses instead lets the graph do
it, and is not quadratic in the batch count. `ISONFORM_RECURSIVE_MERGE=<passes>`,
opt-in --- it is a different algorithm, not a repair of the reference's.

**What the redundancy actually is.** Of the 61 multi-member groups the first pass
forms on `droso_deep` cluster 3, **57 draw from several source batches and only 4
from within one**. So the initial bubble popping is not failing; what it leaves is
a *batching artifact*, the same isoform independently reconstructed in different
batches.

**Weights are the precondition.** Support is a count of read accessions
everywhere it matters, so on a recursive pass it becomes a count of *consensuses*
--- one to three --- and `--iso_abundance 5` discards nearly everything. Each
consensus therefore enters carrying `:w=<reads behind it>` ([`crate::weights`]),
which is read by interval selection, and each isoform that comes out inherits its
members' **original read accessions**, so the abundance gate and the mapping files
keep their meaning. Weight is conserved: 9 163 in, 9 157 after three passes.

**One pass, and a batch floor.** More passes cost 1.5--2.3x wall clock and move
nothing. The floor is the larger finding, on `sirv_real_deep`:

| min batches | strict F1 | transcripts | recall | precision | redundancy | wall |
|---|---|---|---|---|---|---|
| pairwise | 0.681 | 56 | 82.4% | 58.0% | 3.48 | 1 894s |
| 2 | 0.689 | 55 | 80.9% | 60.0% | 3.38 | 1 397s |
| 5 | 0.703 | 57 | 83.8% | 60.6% | 3.32 | 1 376s |
| **10** | **0.713** | **58** | **85.3%** | **61.2%** | **3.29** | **1 364s** |
| 20 | 0.710 | 58 | 85.3% | 60.9% | 3.38 | 1 370s |

At floor 10 it beats the pairwise merge on **every** metric and is 1.39x faster.
At floor 2 it costs a transcript --- which is why the floor exists, and why "batch
count is not the discriminating variable" (written down before the sweep) was
wrong. `droso_deep` at floor 2 already gains: 9 356 FSM against 9 339, 1 330
off-target genic against 1 459, 1.20x faster.

The cost is memory, 3.1 -> 9.7 GB on `droso_deep`. The recursive instance is
unbatched by construction, and its size is the cluster's isoform count: median
**4**, p95 25, p99 104, max **3 466**. So the peak is a handful of clusters on
eight threads, and the lever is bounding the stage's concurrency rather than
batching it --- batching would reintroduce the cross-batch problem the pass exists
to remove.

#### Two things measured and dropped, and one bug worth recording

* **Weighting the graph's bubble-candidacy gate** (`inter.len() >= 2`, the only
  magnitude decision in simplification) is **byte-identical** on all three shallow
  corpora. The hypothesis --- thinly supported consensuses outvoting well supported
  ones --- is disproved. Reverted.
* **`sirv_real` has exactly one multi-batch cluster**, 4 595 reads in 5 batches, so
  its whole 0.034 strict-F1 loss at floor 2 comes from one cluster. Worth knowing
  before generalising from that corpus.
* **The stage was wired as an if/else**, so a cluster the recursive pass declined
  got *no* cross-batch merging rather than the pairwise one. Invisible at floor 2,
  where declining only happens on single-batch clusters and the pairwise merge is
  a no-op anyway; the floor sweep exposed it immediately --- floor 6 on `sirv_real`
  gave 109 isoforms where the baseline gives 86. Declining now falls through.

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

**Report results, not progress.** While anything is running, say nothing --- no
interim findings, no partial tables, no "here is what I am about to do". One work
stint gets one reply, in three parts and no more:

1. **what was measured**, in a sentence or two;
2. **the numbers**, as a table;
3. **what to do next**, as a proposal.

Never restate a conclusion already given. Corrections are one sentence: what is now
true, not a retrospective on the error. Detail goes in this file, not in the reply
--- the reader will ask if they want more.


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

### Cross-batch merging: measured, filtered, and still the bottleneck

`ISONFORM_MERGE_MINSHARE=<bases>` skips a pair when the co-linear chain of shared
minimizer anchors implies an **indel larger than that**. Kristoffer's criterion,
and the right one: what decides a merge is whether a structural difference exists,
not how much sequence is shared. Two dead ends preceded it, both worth keeping ---

* **a minhash sketch** (`src/sketch.rs`): no false-negative bound, only a rate
  measured once. Parked on that basis.
* **counting shared minimizers**: sound at every threshold and worth nothing.
  Raising the count 15x moved the skip rate 10.5% -> 12.4%. Every isoform in a
  cluster is from the same gene, so mergeable and non-mergeable pairs share nearly
  all their minimizers. The count cannot see structure.

Verified on cluster 74 of `droso_deep` (2 466 reads, 3 batches, ~202 isoforms), with
every skipped pair aligned anyway and counted:

| max chain indel | skipped | merges | false negatives |
|---|---|---|---|
| 5bp (`delta_len`) | 94.9% | 24 | **5** |
| **30bp (`delta_iso_len_3`)** | **38.0%** | 27 | **0** |
| 50bp (`delta_iso_len_5`) | 34.0% | 27 | 0 |
| unbounded | 10.5% | 27 | 0 |

X cannot be `delta_len`: minimizer sampling at `w = 10` is sparser than the indels
the merge tolerates, so the measured offset jump overstates the true one and 5
merges are lost. `delta_iso_len_3 = 30` is sound here and is a value the method
already uses.

On cluster 0 (25 493 reads, 26 batches, 3 466 isoforms):

| | pairs seen | skipped | merges | isoforms out | wall |
|---|---|---|---|---|---|
| indel X=30 | 4 468 468 | 43.2% | **605** | 175 | **2 606s** |
| chain-span x=50% | 4 519 493 | 71.4% | 583 | 176 | --- |

The span variant skips more and **loses ~22 merges**; the indel criterion keeps
more for less skipping. At ~1ms per alignment --- WFA2, far cheaper than the 17ms
parasail figure this analysis started from --- unfiltered would be ~75 minutes, so the
filter is worth about **1.8x**.

#### But the real problem is upstream

605 merges for 4.5M examined pairs, and the `iso_abundance >= 5` filter then
removes 2 707 of the 2 883 survivors. The stage costs 43 minutes on one cluster to
remove 17% of its isoforms, while a cheap abundance rule removes 78%.

3 466 isoforms from 25 493 reads is **133 per 1 000-read batch**, and the
cross-batch cost is quadratic in that number. Halving the isoforms per batch
quarters this stage. So the lever is upstream: merge more in the graph, before
isoforms are emitted, which is both cheaper and attacks the redundancy the corpora
already measure. A faster aligner or a better prefilter only trims the constant.

### Superseded note: cross-batch merging needs a different method

Measured on the 12 deepest `droso_deep` clusters --- 94 603 reads, 102 batches,
deepest cluster 25 493 reads across 26 batches:

| phase | elapsed |
|---|---|
| all within-batch work | **1 403s** (23 min) |
| cross-batch merging | **~2.6 hours, unfinished** |

Seven times the entire rest of the pipeline, on twelve clusters, and this was with
parasail --- WFA2 is off by default and was not in the run. It would not rescue it
either: the cost is the **structure**, `O(batch pairs x isoforms per batch^2)`
all-vs-all, and finding 42 showed the merge stage is POA-dominated rather than
alignment-dominated in the first place. A faster aligner moves a minority of a
minority.

Two directions worth trying before optimising anything inside the current loop,
both from Kristoffer:

* **stop doing all-vs-all.** The pairs that can merge are a tiny fraction of the
  pairs examined (82% are rejected as structural), and the loop currently proves
  that one alignment at a time. A method that proposes candidate pairs instead of
  enumerating them changes the exponent, not the constant. `src/sketch.rs` was an
  attempt and was parked for lacking a window guarantee, which is the right bar.
* **merge more in the graph.** Anything merged before isoforms are emitted never
  reaches the cross-batch stage at all. Cheaper than merging it afterwards, and it
  attacks the redundancy the corpora already show.

Not attempted. Recorded here because the measurement says the aligner and the POA
schedule are both largely spent as levers on this stage, and the remaining gain has
to come from doing fewer comparisons rather than faster ones.


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
