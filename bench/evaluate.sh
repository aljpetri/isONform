#!/usr/bin/env bash
#
# Build the evaluation corpora and score isONform on them.
#
# Equivalence testing answers "do two implementations agree". It cannot answer
# "is the output any good", and that second question is the one that matters
# whenever the reference is *changed* — a determinism fix, a deliberate
# divergence in the port, a performance tradeoff. This script is how that gets
# measured.
#
#   bench/evaluate.sh corpora                    # build all three corpora
#   bench/evaluate.sh run   <name> <impl> <tag>  # run isONform, tagged
#                                                # <impl> is a dir holding main +
#                                                # the parallel driver: the repo
#                                                # root, or rust/target/release
#   bench/evaluate.sh score <name> <tag>...      # score the tagged runs
#                                                # (droso picks up
#                                                #  $DROSO_ANNOTATION if present)
#   bench/evaluate.sh compare <impl-a> <impl-b>  # build if needed, run both, score
#
# Work lands in $ISONFORM_WORK (default ~/work/isonform-corpus), never in the
# repository: the corpora are gigabytes of derived fastq.
#
# ---------------------------------------------------------------------------
# The three corpora, and why three
# ---------------------------------------------------------------------------
#
# PORTING.md method point 4 exists because simulated and spike-in data gave the
# wrong answer twice during the isONcorrect port — endorsing an alignment band
# that real data broke, and reversing the sign of an accuracy effect. One corpus
# is not enough, and the three here fail in different directions on purpose:
#
#   sirv_sim_gene   Simulated SIRV, grouped by *gene* so isoforms of one gene
#                   share a cluster and reads differ by whole exons. Truth is
#                   exact and free: the source transcript is in every read
#                   header, so recall has a real denominator. Reproducible from
#                   the fastq alone, with no clustering run in the loop.
#                   Weakness: uniform i.i.d. error, no chimeras, no truncation.
#
#   sirv_real       Real ONT SIRV spike-in reads through isONclust. Real error
#                   profile, real chimeras, and *still* a known ground-truth
#                   transcriptome, which is the combination that makes SIRV
#                   worth its narrowness. Weakness: 68 transcripts of one
#                   synthetic locus set — not a transcriptome-scale test.
#
#   droso           Real ONT Drosophila cDNA through isONclust. Transcriptome
#                   scale, real biology, hundreds of clusters. No trustworthy
#                   per-isoform truth, so it is scored against the *genome* by
#                   spliced alignment — see accuracy_isoforms_genome.py, which
#                   measures splice-junction canonicality precisely because that
#                   needs no annotation. With DROSO_ANNOTATION set it also
#                   reports the SQANTI structural categories, which do need one
#                   and are worth far more: 88% FSM says something canonical
#                   splice fraction cannot.
#
# A change that improves all three is an improvement. A change that improves
# sirv_sim_gene and not the other two is a change that fits simulated error.
#
# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------
#
# Every *Python* isONform invocation here pins PYTHONHASHSEED, and `run` takes
# the seed so the reference can be run at several. Before minimizer selection was
# made lexicographic, that mattered enormously: the same input gave a different
# transcriptome on every run. It is no longer the dominant effect, but it is not
# belt and braces either — PORTING.md finding 14 is a second, independent seed
# dependency that survives the minimizer fix and fires on real data.
#
# So a single Python run is one *sample*, not a baseline. When comparing the Rust
# port — which has no seeded hash to pin and gives the same answer every time —
# run the reference at several seeds and ask whether the port's single result
# sits inside that spread. `compare` does this for you.
#
# The Rust port is passed as an impl directory like any other:
#
#   bench/evaluate.sh run sirv_real "$PWD/rust/target/release" rs
#
# `run_one` decides Python-vs-native by reading the file, not by its name, because
# macOS is case-insensitive and `isonform_parallel` would otherwise resolve to the
# repository's `isONform_parallel`.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${ISONFORM_WORK:-$HOME/work/isonform-corpus}"
ENV_NAME="${ISONFORM_REF_ENV:-isonform-ref}"

# Inputs. Override any of these in the environment.
SIRV_SIM_FASTQ="${SIRV_SIM_FASTQ:-$HOME/data/lrRNA-seq/sirv/reads_10k_err7%.fastq}"
SIRV_REAL_FASTQ="${SIRV_REAL_FASTQ:-$HOME/data/lrRNA-seq/sirv/SIRV_real_10k.fastq}"
DROSO_FASTQ="${DROSO_FASTQ:-$HOME/data/lrRNA-seq/droso/full_length_output_first_20k.fq}"
SIRV_TRANSCRIPTOME="${SIRV_TRANSCRIPTOME:-$HOME/source/isONcorrect/test_data/sirv_transcriptome.fasta}"
DROSO_GENOME="${DROSO_GENOME:-$HOME/data/genomes/fruitfly.fa}"
# Optional. When present, the droso corpus is additionally scored into the SQANTI
# structural categories. Filter the FlyBase whole-genome GFF first --- it is
# 6.7 GB and mostly alignment evidence; see bench/corpus/README.md.
DROSO_ANNOTATION="${DROSO_ANNOTATION:-$HOME/work/dmel-r6.68-transcripts.tsv}"

# Upstream tools. isONcorrect's Rust build is used because it is equivalent to
# the Python and an order of magnitude faster; the corpus is an input to this
# work, not a subject of it.
RUN_ISONCORRECT="${RUN_ISONCORRECT:-$HOME/source/isONcorrect/rust/target/release/run_isoncorrect}"
MAKE_SIRV_CORPUS="${MAKE_SIRV_CORPUS:-$HOME/source/isONcorrect/bench/make_sirv_corpus.py}"

THREADS="${THREADS:-8}"
ISO_ABUNDANCE="${ISO_ABUNDANCE:-5}"

die() { echo "error: $*" >&2; exit 1; }
note() { printf '\033[36m==>\033[0m %s\n' "$*"; }

# --- environment ------------------------------------------------------------

PREFIX="$(conda env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n {print $NF; exit}')"
[[ -n "$PREFIX" && -x "$PREFIX/bin/python" ]] || die "reference env '$ENV_NAME' not found; run bench/setup_reference_env.sh"
PY_BIN="$PREFIX/bin/python"
export PATH="$PREFIX/bin:$PATH"   # spoa

mkdir -p "$WORK"

# --- corpora ----------------------------------------------------------------

# `correct DIR_IN DIR_OUT` — isONcorrect over a folder of cluster fastqs.
# isONform consumes *corrected* reads: run on raw isONclust output the
# interval-abundance filter skips essentially everything and the graph comes out
# empty. This is not optional preprocessing, it is the input contract.
correct() {
  local in="$1" out="$2"
  [[ -x "$RUN_ISONCORRECT" ]] || die "run_isoncorrect not found at $RUN_ISONCORRECT"
  rm -rf "$out"
  "$RUN_ISONCORRECT" --fastq_folder "$in" --outfolder "$out" \
      --t "$THREADS" --split_wrt_batches >"$out.log" 2>&1 \
    || die "isONcorrect failed; see $out.log"
  note "$(find "$out" -name corrected_reads.fastq | wc -l | tr -d ' ') corrected clusters in $out"
}

# `cluster FASTQ DIR` — isONclust, the documented first stage of the pipeline
# (isON_pipeline.sh). `--N $ISO_ABUNDANCE` on write_fastq drops clusters too
# small to support an isoform, matching what the pipeline script does.
cluster() {
  local fastq="$1" dir="$2"
  command -v isONclust >/dev/null || die "isONclust not on PATH (try: conda activate isonclust-tool)"
  rm -rf "$dir"; mkdir -p "$dir"
  isONclust --t "$THREADS" --ont --fastq "$fastq" --outfolder "$dir/clustering" >"$dir/clustering.log" 2>&1 \
    || die "isONclust failed; see $dir/clustering.log"
  isONclust write_fastq --N "$ISO_ABUNDANCE" --clusters "$dir/clustering/final_clusters.tsv" \
      --fastq "$fastq" --outfolder "$dir/clusters" >>"$dir/clustering.log" 2>&1 \
    || die "isONclust write_fastq failed; see $dir/clustering.log"
  note "$(ls "$dir/clusters" | wc -l | tr -d ' ') clusters in $dir/clusters"
}

corpus_sirv_sim_gene() {
  local d="$WORK/sirv_sim_gene"
  [[ -f "$SIRV_SIM_FASTQ" ]] || die "simulated SIRV reads not found: $SIRV_SIM_FASTQ"
  note "sirv_sim_gene: grouping simulated reads by gene (no clustering needed --- the truth is in the headers)"
  rm -rf "$d"; mkdir -p "$d"
  "$PY_BIN" "$MAKE_SIRV_CORPUS" --fastq "$SIRV_SIM_FASTQ" --outdir "$d/clusters" --group gene
  correct "$d/clusters" "$d/corrected"
}

corpus_sirv_real() {
  local d="$WORK/sirv_real"
  [[ -f "$SIRV_REAL_FASTQ" ]] || die "real SIRV reads not found: $SIRV_REAL_FASTQ"
  note "sirv_real: isONclust on real ONT SIRV reads"
  cluster "$SIRV_REAL_FASTQ" "$d"
  correct "$d/clusters" "$d/corrected"
}

corpus_droso() {
  local d="$WORK/droso"
  [[ -f "$DROSO_FASTQ" ]] || die "Drosophila reads not found: $DROSO_FASTQ"
  note "droso: isONclust on real ONT Drosophila cDNA"
  cluster "$DROSO_FASTQ" "$d"
  correct "$d/clusters" "$d/corrected"
}

# --- running ----------------------------------------------------------------

# `run_one CORPUS IMPL_DIR TAG [SEED]`
#
# IMPL_DIR is a directory holding `main` and the parallel driver — the repository
# root for the Python reference, a `git archive` extraction of any commit, or
# `rust/target/release` for the port. That is how a before/after measurement is
# made without touching the working tree.
#
# Which of the two it is decided by *content*, not filename: macOS mounts a
# case-insensitive filesystem by default, so `$impl/isonform_parallel` happily
# resolves to the repository's `isONform_parallel` and a name test would call
# the Python reference a native build. `impl_entry` reads the first two bytes
# and looks for `#!` instead.
impl_entry() {
  local impl="$1" p
  for p in "$impl/isONform_parallel" "$impl/isonform_parallel"; do
    [[ -f "$p" ]] || continue
    if head -c 2 "$p" 2>/dev/null | grep -q '#!'; then
      printf 'python\t%s\n' "$p"
    else
      printf 'native\t%s\n' "$p"
    fi
    return 0
  done
  return 1
}

run_one() {
  local corpus="$1" impl="$2" tag="$3" seed="${4:-0}"
  local src="$WORK/$corpus/corrected"
  local out="$WORK/eval/${corpus}__${tag}"
  [[ -d "$src" ]] || die "corpus '$corpus' not built; run: bench/evaluate.sh corpora"
  local kind entry
  IFS=$'\t' read -r kind entry < <(impl_entry "$impl") \
    || die "no isONform_parallel or isonform_parallel in $impl"

  rm -rf "$out"; mkdir -p "$out"
  # isONform_parallel's restructure_isoncorrect_output() MUTATES its input
  # directory --- it shutil.moves every <cl>/corrected_reads.fastq up a level and
  # then deletes the folders (PORTING.md finding 5). Every run therefore gets a
  # private copy, or the second run would see a different input from the first.
  cp -R "$src" "$out/in"

  note "$corpus / $tag  (impl=$impl kind=$kind seed=$seed)"
  local start end
  start=$("$PY_BIN" -c 'import time;print(time.time())')
  if [[ "$kind" == python ]]; then
    ( cd "$impl" && PYTHONHASHSEED="$seed" "$PY_BIN" "$entry" \
        --fastq_folder "$out/in" --outfolder "$out/out" --t "$THREADS" \
        --split_wrt_batches --iso_abundance "$ISO_ABUNDANCE" ) >"$out/log" 2>&1
  else
    # No PYTHONHASHSEED: the port has no seeded hash to pin, which is the
    # point of finding 28's ascending order. `seed` is still recorded so the
    # two implementations' rows line up in runs.tsv.
    ( cd "$impl" && "$entry" \
        --fastq_folder "$out/in" --outfolder "$out/out" --t "$THREADS" \
        --split_wrt_batches --iso_abundance "$ISO_ABUNDANCE" ) >"$out/log" 2>&1
  fi
  local ec=$?
  end=$("$PY_BIN" -c 'import time;print(time.time())')

  local secs n
  secs=$("$PY_BIN" -c "print(f'{$end-$start:.1f}')")
  n=$(grep -c '^>' "$out/out/transcriptome.fasta" 2>/dev/null || echo 0)
  printf '    exit=%s  %ss  %s isoforms\n' "$ec" "$secs" "$n"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$corpus" "$tag" "$seed" "$ec" "$secs" "$n" \
      >> "$WORK/eval/runs.tsv"
  # The private input copy is large and reproducible; the corpus is the source
  # of truth.
  rm -rf "$out/in"
  [[ $ec -eq 0 ]] || echo "    (failed --- see $out/log)" >&2
}

# --- scoring ----------------------------------------------------------------

score_corpus() {
  local corpus="$1"; shift
  local -a specs=()
  local tag
  for tag in "$@"; do
    local f="$WORK/eval/${corpus}__${tag}/out/transcriptome.fasta"
    [[ -f "$f" ]] && specs+=(--isoforms "$tag=$f") \
                  || echo "    (no output for $corpus/$tag, skipping)" >&2
  done
  [[ ${#specs[@]} -gt 0 ]] || { echo "nothing to score for $corpus" >&2; return 1; }

  echo
  case "$corpus" in
    sirv_sim_gene)
      # Simulated: the read headers name their source transcript, so the recall
      # denominator is exactly the set of transcripts actually present.
      "$PY_BIN" "$REPO_ROOT/bench/accuracy_isoforms.py" \
        --transcriptome "$SIRV_TRANSCRIPTOME" \
        --expressed-from "$SIRV_SIM_FASTQ" "${specs[@]}"
      ;;
    sirv_real)
      # Real SIRV: no per-read truth, so every reference transcript counts and
      # recall is a lower bound. Said out loud by the script itself.
      "$PY_BIN" "$REPO_ROOT/bench/accuracy_isoforms.py" \
        --transcriptome "$SIRV_TRANSCRIPTOME" "${specs[@]}"
      ;;
    droso)
      local -a ann=()
      if [[ -f "$DROSO_ANNOTATION" ]]; then
        ann=(--annotation "$DROSO_ANNOTATION")
      else
        echo "    (no annotation at $DROSO_ANNOTATION --- reporting the" >&2
        echo "     annotation-free metrics only; see bench/corpus/README.md)" >&2
      fi
      "$PY_BIN" "$REPO_ROOT/bench/accuracy_isoforms_genome.py" \
        --genome "$DROSO_GENOME" --threads "$THREADS" "${ann[@]}" "${specs[@]}"
      ;;
    *) die "unknown corpus '$corpus'" ;;
  esac
}

# --- entry points -----------------------------------------------------------

cmd="${1:-}"; shift || true
case "$cmd" in
  corpora)
    mkdir -p "$WORK/eval"
    corpus_sirv_sim_gene
    corpus_sirv_real
    corpus_droso
    note "corpora in $WORK"
    ;;

  run)
    [[ $# -ge 3 ]] || die "usage: evaluate.sh run <corpus> <impl-dir> <tag> [seed]"
    mkdir -p "$WORK/eval"
    run_one "$@"
    ;;

  score)
    [[ $# -ge 2 ]] || die "usage: evaluate.sh score <corpus> <tag>..."
    score_corpus "$@"
    ;;

  compare)
    # The common case: two implementations, all three corpora, scored side by
    # side. The pre-fix implementation is run at three seeds rather than one,
    # because with a seeded hash a single run says nothing about what the tool
    # does — that spread *is* the finding.
    [[ $# -eq 2 ]] || die "usage: evaluate.sh compare <impl-a-dir> <impl-b-dir>"
    a="$1"; b="$2"
    mkdir -p "$WORK/eval"
    for corpus in sirv_sim_gene sirv_real droso; do
      [[ -d "$WORK/$corpus/corrected" ]] || die "corpus '$corpus' not built; run: bench/evaluate.sh corpora"
      run_one "$corpus" "$a" a 0
      run_one "$corpus" "$b" b0 0
      run_one "$corpus" "$b" b1 1
      run_one "$corpus" "$b" b2 2
      score_corpus "$corpus" a b0 b1 b2
    done
    ;;

  *)
    sed -n '2,60p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
    exit 2
    ;;
esac
