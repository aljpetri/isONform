#!/usr/bin/env bash
#
# Differential equivalence between the Python reference and the Rust port.
#
# Two layers, run independently:
#
#   cli        argument handling: exit codes, stderr diagnostics, the stray
#              stdout prints. Compares the port against the reference live.
#              This layer is complete and passing.
#
#   invariants properties the reference must satisfy whatever the implementation.
#              These need neither a golden nor a second implementation to be
#              trusted, which is what makes them the tool of choice when an
#              oracle fires and the question is *which side is wrong*.
#
#   stages     the algorithm: isoform output for both entry points across the
#              flags that change behaviour. Recordable now that the reference is
#              deterministic; nothing to diff against until the port has stages.
#
# Usage:
#   bench/equivalence.sh              # every runnable layer
#   bench/equivalence.sh cli
#   bench/equivalence.sh invariants
#   bench/equivalence.sh stages
#
# ---------------------------------------------------------------------------
# The determinism gate
# ---------------------------------------------------------------------------
#
# `main` used to select minimizers by `hash()` of the k-mer string, and CPython's
# string hash is seeded randomly per process, so the reference produced a
# different transcriptome on every run. Recording a golden against that would
# have recorded one draw from a lottery. Measured on bench/corpus/sirv_small: 8
# seeds gave 8 distinct `mapping0`s and `transcriptome.fasta` took 5 distinct
# values. See PORTING.md finding 1.
#
# That is fixed --- minimizer selection is lexicographic --- so goldens are
# recordable. The gate is kept, but as a *check* rather than a refusal: before
# recording anything, `--record` runs each case at two different hash seeds and
# refuses if the outputs differ. A hardcoded refusal would have gone stale the
# moment the fix landed; this one cannot.
#
# PYTHONHASHSEED is still exported below as belt and braces.
# bench/check_seed_sensitivity.py is what asserts that is no longer
# load-bearing.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${ISONFORM_REF_ENV:-isonform-ref}"
GOLDEN_DIR="$REPO_ROOT/bench/golden"
CORPUS="$REPO_ROOT/bench/corpus/sirv_small"
RUST_DIR="$REPO_ROOT/rust"

# Every reference invocation gets this. See "The determinism gate" above.
export PYTHONHASHSEED=0

RECORD=0
LAYERS=()
for arg in "$@"; do
  case "$arg" in
    --record) RECORD=1 ;;
    cli|invariants|stages) LAYERS+=("$arg") ;;
    -h|--help) sed -n '2,42p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'; exit 0 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done
[[ ${#LAYERS[@]} -eq 0 ]] && LAYERS=(cli invariants stages)

# --- locate the two implementations ----------------------------------------

PREFIX="$(conda env list 2>/dev/null | awk -v n="$ENV_NAME" '$1==n {print $NF; exit}')"
if [[ -z "$PREFIX" || ! -x "$PREFIX/bin/python" ]]; then
  echo "error: reference env '$ENV_NAME' not found. Run bench/setup_reference_env.sh." >&2
  exit 1
fi
PY_BIN="$PREFIX/bin/python"
export PATH="$PREFIX/bin:$PATH"   # spoa

RS_MAIN="$RUST_DIR/target/release/main"
RS_PAR="$RUST_DIR/target/release/isONform_parallel"
# Always build. A stale release binary passing this script is exactly the kind of
# false green that makes an equivalence harness worse than none.
echo "==> building the Rust port"
cargo build --release --manifest-path "$RUST_DIR/Cargo.toml" >/dev/null || exit 1

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0
SKIP=0

ok()   { PASS=$((PASS+1)); printf '  \033[32mok\033[0m   %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  \033[31mFAIL\033[0m %s\n' "$1"; [[ -n "${2:-}" ]] && printf '       %s\n' "$2"; }
skip() { SKIP=$((SKIP+1)); printf '  \033[33mskip\033[0m %s\n' "$1"; [[ -n "${2:-}" ]] && printf '       %s\n' "$2"; }

# ===========================================================================
# Layer: cli
# ===========================================================================
#
# What is compared, and what is not.
#
# Compared:  exit code, and stderr *content* for the reference's own
#            diagnostics (the ones it writes itself, not argparse's).
# Not compared: --help layout and argparse's usage block. clap wraps and orders
#            differently, nothing downstream parses help text, and the
#            isONcorrect port made the same call. The reference's help output is
#            recorded verbatim in bench/golden/cli/ so a flag change is visible
#            on inspection.
#
# Compared for the stray stdout prints, though, because those are not help
# text: `isONform_parallel` writes `len(sys.argv)` as its first line of stdout
# and that is observable contract (PORTING.md finding 4).

# `pre_work_stderr FILE` — the reference's own diagnostics, with the Python
# traceback stripped.
#
# Cases where validation *passes* run into the algorithm, and the reference then
# dies on a missing input with a traceback while the port stops with its
# not-implemented marker. Comparing whole stderr there compares two unrelated
# failures. What is comparable is everything the reference wrote *before* it
# started work, which is exactly the lines ahead of the traceback.
pre_work_stderr() {
  sed '/^Traceback (most recent call last):/,$d' "$1"
}

# `port_stderr FILE` — the port's stderr with its "not yet implemented" notice
# stripped. Every line of that notice is prefixed with a sentinel precisely so
# it can be removed here; when the algorithm lands the notice goes and this
# becomes the identity.
port_stderr() {
  grep -v '^\[port-unimplemented\] ' "$1"
}

# cli_case DESC ENTRY EXPECT [args...]
#
#   EXPECT is one of:
#     usage      argparse rejects the arguments: both exit 2, and the port's
#                message names the offending flag. (Wording is not compared:
#                clap and argparse phrase it differently and nothing parses it.)
#     rejected   the reference's own validation rejects it: both exit 1 and the
#                whole of stderr matches, verbatim.
#     accepted   validation passes on both sides. The reference goes on to the
#                algorithm; the port stops at exit 3. Compared: the reference's
#                pre-work stderr against the port's stderr, and the first line
#                of stdout. Asserted separately: the reference did *not* exit 2,
#                i.e. it really did get past argparse.
cli_case() {
  local desc="$1" entry="$2" expect="$3"; shift 3
  local py_ec rs_ec rs_bin

  "$PY_BIN" "$REPO_ROOT/$entry" "$@" >"$WORK/py.out" 2>"$WORK/py.err"; py_ec=$?
  [[ "$entry" == "main" ]] && rs_bin="$RS_MAIN" || rs_bin="$RS_PAR"
  "$rs_bin" "$@" >"$WORK/rs.out" 2>"$WORK/rs.err"; rs_ec=$?

  case "$expect" in
    usage)
      if [[ "$py_ec" != 2 ]]; then
        bad "$desc" "reference did not reject it (exit $py_ec) — the case is wrong, not the port"
        return
      fi
      if [[ "$rs_ec" != 2 ]]; then
        bad "$desc" "port exited $rs_ec, expected 2"
        return
      fi
      local flag
      for flag in "$@"; do [[ "$flag" == --* ]] && break; done
      if ! grep -qF -- "$flag" "$WORK/rs.err"; then
        bad "$desc" "port's error does not name $flag"
        return
      fi
      ;;

    rejected)
      if [[ "$py_ec" != "$rs_ec" ]]; then
        bad "$desc" "exit code: reference $py_ec, port $rs_ec"
        return
      fi
      if ! diff -q "$WORK/py.err" "$WORK/rs.err" >/dev/null; then
        bad "$desc" "stderr differs: '$(head -1 "$WORK/py.err")' vs '$(head -1 "$WORK/rs.err")'"
        return
      fi
      ;;

    accepted)
      if [[ "$py_ec" == 2 ]]; then
        bad "$desc" "reference rejected it at parse time — the case expects it to be accepted"
        return
      fi
      if [[ "$rs_ec" != 3 ]]; then
        # 3 is the port's "argument handling done, algorithm not ported" marker.
        # Anything else means it disagreed with the reference about validity.
        bad "$desc" "port exited $rs_ec, expected 3 (accepted-but-unimplemented)"
        return
      fi
      pre_work_stderr "$WORK/py.err" >"$WORK/py.pre"
      port_stderr     "$WORK/rs.err" >"$WORK/rs.pre"
      if ! diff -q "$WORK/py.pre" "$WORK/rs.pre" >/dev/null; then
        bad "$desc" "pre-work stderr differs: '$(head -1 "$WORK/py.pre")' vs '$(head -1 "$WORK/rs.pre")'"
        return
      fi
      local p r
      p="$(head -1 "$WORK/py.out")"; r="$(head -1 "$WORK/rs.out")"
      if [[ "$p" != "$r" ]]; then
        bad "$desc" "first stdout line: '$p' vs '$r'"
        return
      fi
      ;;

    *) bad "$desc" "unknown expectation '$expect'"; return ;;
  esac

  ok "$desc"
}

# cli_stdout_case DESC ENTRY [args...] — compare only the first stdout line and
# that both sides exit 0. For --version and the no-arguments help path, where
# the rest of stdout is help layout and deliberately not compared.
cli_stdout_case() {
  local desc="$1" entry="$2"; shift 2
  local py_ec rs_ec rs_bin

  "$PY_BIN" "$REPO_ROOT/$entry" "$@" >"$WORK/py.out" 2>"$WORK/py.err"; py_ec=$?
  [[ "$entry" == "main" ]] && rs_bin="$RS_MAIN" || rs_bin="$RS_PAR"
  "$rs_bin" "$@" >"$WORK/rs.out" 2>"$WORK/rs.err"; rs_ec=$?

  if [[ "$py_ec" != "$rs_ec" ]]; then
    bad "$desc" "exit code: reference $py_ec, port $rs_ec"
    return
  fi
  if ! diff -q "$WORK/py.err" "$WORK/rs.err" >/dev/null; then
    bad "$desc" "stderr differs: '$(head -1 "$WORK/py.err")' vs '$(head -1 "$WORK/rs.err")'"
    return
  fi
  ok "$desc"
}

layer_cli() {
  echo "== cli =="

  # --version: argparse prints "%(prog)s 0.3.9"; clap must match exactly.
  # Compared on stdout by cli_case's `accepted`? No — argparse exits 0 here, so
  # this needs its own shape.
  local v_py v_rs
  v_py="$("$PY_BIN" "$REPO_ROOT/main" --version 2>&1)"
  v_rs="$("$RS_MAIN" --version 2>&1)"
  [[ "$v_py" == "$v_rs" ]] && ok "main --version ($v_py)" \
    || bad "main --version" "'$v_py' vs '$v_rs'"
  v_py="$("$PY_BIN" "$REPO_ROOT/isONform_parallel" --version 2>&1)"
  v_rs="$("$RS_PAR" --version 2>&1)"
  [[ "$v_py" == "$v_rs" ]] && ok "isONform_parallel --version ($v_py)" \
    || bad "isONform_parallel --version" "'$v_py' vs '$v_rs'"

  # No arguments: exit 0, help on stdout. `main` writes the xmin adjustment to
  # stderr *first* — the order is the point, and stderr is compared verbatim.
  cli_stdout_case "main, no arguments (xmin line precedes help)" main

  # `isONform_parallel` prints `1` — len(sys.argv) — as its first stdout line
  # before the help text. PORTING.md finding 4.
  local a_py a_rs
  a_py="$("$PY_BIN" "$REPO_ROOT/isONform_parallel" 2>/dev/null | head -1)"
  a_rs="$("$RS_PAR" 2>/dev/null | head -1)"
  [[ "$a_py" == "1" && "$a_rs" == "1" ]] && ok "isONform_parallel, no arguments prints argc 1" \
    || bad "isONform_parallel argc on no arguments" "'$a_py' vs '$a_rs'"

  # The same stray print with arguments present, where it reports the real argc.
  cli_case "isONform_parallel argc print" isONform_parallel accepted --fastq_folder /nonexistent

  # The xmin adjustment, announced verbatim on stderr.
  cli_case "xmin raised to 2k (k=20)" main accepted --fastq /nonexistent.fq
  cli_case "xmin raised to 2k (k=25)" main accepted --k 25 --fastq /nonexistent.fq
  cli_case "xmin above 2k, silent"    main accepted --k 10 --xmin 60 --fastq /nonexistent.fq
  cli_case "xmin exactly 2k, silent"  main accepted --k 20 --xmin 40 --fastq /nonexistent.fq

  # The window check: exit 1 plus the reference's own message, verbatim.
  cli_case "w < k rejected"    main rejected --k 20 --w 10 --fastq /nonexistent.fq
  cli_case "w > 100 rejected"  main rejected --w 101       --fastq /nonexistent.fq
  # And the boundaries the message gets wrong: `100 < w` admits 100.
  cli_case "w == 100 accepted" main accepted --w 100        --fastq /nonexistent.fq
  cli_case "w == k accepted"   main accepted --k 31 --w 31  --fastq /nonexistent.fq
  # The driver has no window check at all, so a window the child would reject is
  # accepted one level up.
  cli_case "driver does not check w" isONform_parallel accepted --k 20 --w 5 --fastq_folder /nonexistent

  # Unrecognised arguments: exit 2, and the flag is named.
  cli_case "main rejects unknown flag"   main              usage --bogus
  cli_case "driver rejects unknown flag" isONform_parallel usage --bogus

  # Flags belonging to only one entry point. `--fastq` is absent from this list
  # on purpose: argparse resolves it to `--fastq_folder` by prefix, and so does
  # the port — see the abbreviation cases below.
  local f
  for f in --keep_old --split_wrt_batches --tmpdir --iso_abundance --write_fastq; do
    cli_case "main rejects $f" main usage "$f" x
  done
  for f in --exact --disable_numpy --compression --slow --parallel; do
    cli_case "driver rejects $f" isONform_parallel usage "$f" x
  done

  # argparse accepts any unambiguous prefix of a long option, and the port now
  # matches. These pass parsing on both sides.
  cli_case "main --comp abbreviates --compression"    main accepted --comp --fastq /nonexistent.fq
  cli_case "main --disable_num abbreviates"           main accepted --disable_num --fastq /nonexistent.fq
  cli_case "main --set_w abbreviates"                 main accepted --set_w --fastq /nonexistent.fq
  cli_case "main --fast abbreviates --fastq"          main accepted --fast /nonexistent.fq
  cli_case "driver --fastq abbreviates --fastq_folder" isONform_parallel accepted --fastq /nonexistent
  # Ambiguous prefixes stay rejected on both sides.
  cli_case "main rejects ambiguous --delta"  main usage --delta 7
  cli_case "main rejects ambiguous --m"      main usage --m 7

  # --parallel's type=bool bug (PORTING.md finding 2): every non-empty value is
  # true, including "False". Parsing is what is compared here; that all four
  # values mean *true* is asserted by cli.rs's unit test and by the measurement
  # recorded in PORTING.md.
  for f in True False 0 no; do
    cli_case "--parallel $f parses" main accepted --fastq /nonexistent.fq --parallel "$f"
  done
}

# ===========================================================================
# Layer: stages
# ===========================================================================

# Flag sets that change observable output, as measured (PORTING.md finding 3).
# The four inert flags get one case each, asserting they stay inert, rather
# than a matrix.
main_cases() {
  cat <<'CASES'
default|
compression|--compression
slow|--slow
dynamic_w|--set_w_dynamically
k15|--k 15 --w 25
k25|--k 25 --w 40
w50|--w 50
xspan_narrow|--xmin 45 --xmax 55
xspan_wide|--xmin 40 --xmax 120
delta_len_20|--delta_len 20
iso_len_tight|--delta_iso_len_3 5 --delta_iso_len_5 5
iso_len_loose|--delta_iso_len_3 200 --delta_iso_len_5 200
spoa_cap_tiny|--max_seqs_to_spoa 3
max_seqs_small|--max_seqs 40
parallel_flag|--parallel True
inert_T|--T 0.9
inert_disable_numpy|--disable_numpy
inert_exact|--exact
inert_exact_limit|--exact_instance_limit 100000
CASES
}

parallel_cases() {
  cat <<'CASES'
default|
split_batches|--split_wrt_batches
write_fastq|--write_fastq
iso_abundance_1|--iso_abundance 1
iso_abundance_high|--iso_abundance 50
delta_high|--delta 0.5
spoa_cap_tiny|--max_seqs_to_spoa 3
two_cores|--t 2
keep_old|--keep_old
CASES
}

# ===========================================================================
# Layer: invariants
# ===========================================================================
#
# Properties the reference must satisfy regardless of how it is implemented.
# PORTING.md's "a mismatch is evidence, not a verdict" needs checks that trust
# neither implementation, and this is where they live: when an oracle fires,
# these are what say which side is wrong.
#
# Graph code is good hunting ground for more of these --- acyclicity,
# connectivity, node and edge counts between stages --- and they should be added
# here as the graph stages get ported.

# The identical-cluster invariant.
#
# bench/corpus/sirv_small's two clusters are byte-identical inputs (see
# bench/corpus/README.md), so they MUST produce the same isoform sequences, and
# accessions differing only in the cluster and batch ids. isONform_parallel runs
# them as two separate child processes, so this catches per-instance state
# leakage, fan-out bugs, and --- before it was fixed --- the seed-dependence of
# finding 1, where the two children each drew their own hash seed and emitted
# different isoforms from the same reads.
#
# Note the id normalisation covers BOTH the cluster id and the batch id
# (`0_0_72` vs `1_1_72`). Stripping only the first field, which is the obvious
# thing to write, makes this check fail on output that is in fact correct.
invariant_identical_clusters() {
  local out="$WORK/inv"
  rm -rf "$out"; mkdir -p "$out/in"
  cp "$CORPUS"/*.fastq "$out/in/" 2>/dev/null || { skip "identical-cluster invariant" "no corpus at $CORPUS"; return; }

  ( cd "$REPO_ROOT" && "$PY_BIN" isONform_parallel --fastq_folder "$out/in" \
        --outfolder "$out/out" --t 1 ) >"$out/log" 2>&1
  if [[ $? -ne 0 ]]; then
    bad "identical-cluster invariant" "reference failed; see $out/log"
    return
  fi

  local a="$out/out/cluster0_merged.fa" b="$out/out/cluster1_merged.fa"
  if [[ ! -s "$a" || ! -s "$b" ]]; then
    bad "identical-cluster invariant" "expected cluster0_merged.fa and cluster1_merged.fa"
    return
  fi

  if ! diff <(grep -v '^>' "$a") <(grep -v '^>' "$b") >/dev/null; then
    bad "identical-cluster invariant" "identical inputs produced different isoform sequences"
    return
  fi
  if ! diff <(sed -E 's/^>[0-9]+_[0-9]+_/>X_X_/' "$a") \
            <(sed -E 's/^>[0-9]+_[0-9]+_/>X_X_/' "$b") >/dev/null; then
    bad "identical-cluster invariant" "sequences agree but accessions differ by more than the ids"
    return
  fi
  ok "identical-cluster invariant ($(grep -c '^>' "$a") isoforms, sequences and ids agree)"
}

# Every read is either used or accounted for as skipped.
#
# `main` writes reads it cannot place to skip<batch>.fa and reports the count on
# stdout. Those two must agree, and neither may exceed the input. A drift here
# means reads are being lost between the interval stage and the output, which no
# golden would necessarily reveal.
invariant_skipped_reads_accounted() {
  local out="$WORK/inv_skip"
  local fq="$CORPUS/0.fastq"
  [[ -f "$fq" ]] || { skip "skipped-read accounting" "no corpus at $CORPUS"; return; }
  rm -rf "$out"; mkdir -p "$out"

  ( cd "$REPO_ROOT" && "$PY_BIN" main --fastq "$fq" --outfolder "$out" ) \
      >"$out/stdout" 2>"$out/stderr"
  if [[ $? -ne 0 ]]; then
    bad "skipped-read accounting" "reference failed; see $out/stderr"
    return
  fi

  local n_in n_reported n_written
  n_in=$(( $(wc -l < "$fq") / 4 ))
  n_reported=$(sed -n 's/^Skipped  \([0-9]*\)  reads.*/\1/p' "$out/stdout" | head -1)
  [[ -z "$n_reported" ]] && n_reported=0
  n_written=$(grep -c '^>' "$out/skip0.fa" 2>/dev/null || echo 0)

  if [[ "$n_reported" != "$n_written" ]]; then
    bad "skipped-read accounting" "reported $n_reported skipped, wrote $n_written to skip0.fa"
    return
  fi
  if (( n_written > n_in )); then
    bad "skipped-read accounting" "wrote $n_written skipped reads from $n_in input reads"
    return
  fi
  ok "skipped-read accounting ($n_written of $n_in reads skipped, count matches skip0.fa)"
}

layer_invariants() {
  echo "== invariants =="
  invariant_identical_clusters
  invariant_skipped_reads_accounted
}

# `determinism_ok` --- does the reference agree with itself across hash seeds?
#
# Two runs of `main` on the same input at different PYTHONHASHSEED values must
# produce identical output. This is cheap (one small cluster, two runs) and it is
# the precondition for a golden meaning anything at all. See the header.
determinism_ok() {
  local a="$WORK/det_a" b="$WORK/det_b" seed
  local fq="$CORPUS/0.fastq"
  [[ -f "$fq" ]] || { echo "  (no corpus, cannot check determinism)" >&2; return 1; }
  for seed in 0 12345; do
    local out="$WORK/det_$seed"
    rm -rf "$out"; mkdir -p "$out"
    PYTHONHASHSEED="$seed" "$PY_BIN" "$REPO_ROOT/main" \
        --fastq "$fq" --outfolder "$out" >/dev/null 2>&1 || return 1
  done
  diff -r "$WORK/det_0" "$WORK/det_12345" >/dev/null
}

layer_stages() {
  echo "== stages =="

  if [[ ! -d "$CORPUS" ]]; then
    skip "all stage cases" "corpus missing: $CORPUS (see bench/corpus/README.md)"
    return
  fi

  # The gate, as a check rather than a refusal: prove the reference agrees with
  # itself across hash seeds before recording anything from it.
  if [[ $RECORD -eq 1 ]]; then
    if ! determinism_ok; then
      echo "  refusing to record stage goldens: the reference disagrees with itself" >&2
      echo "  across hash seeds. Something whose iteration order reaches results has" >&2
      echo "  crept back in --- see PORTING.md finding 1 and" >&2
      echo "  bench/check_seed_sensitivity.py before recording anything." >&2
      exit 1
    fi
    ok "reference is seed-independent (checked at 2 seeds before recording)"
  fi

  local n_main n_par
  n_main=$(main_cases | grep -c .)
  n_par=$(parallel_cases | grep -c .)

  # The port has no algorithm yet, so there is nothing to diff. Report the
  # planned coverage rather than printing nothing: a silent layer reads as
  # "covered".
  skip "$n_main main cases, $n_par isONform_parallel cases" \
       "goldens are recordable now that the reference is deterministic, but the Rust port implements argument handling only --- nothing to diff against yet"
  echo "       main:     $(main_cases | cut -d'|' -f1 | tr '\n' ' ')"
  echo "       parallel: $(parallel_cases | cut -d'|' -f1 | tr '\n' ' ')"
}

# ===========================================================================

for layer in "${LAYERS[@]}"; do
  case "$layer" in
    cli)        layer_cli ;;
    invariants) layer_invariants ;;
    stages)     layer_stages ;;
  esac
  echo
done

printf 'passed %d, failed %d, skipped %d\n' "$PASS" "$FAIL" "$SKIP"
[[ $FAIL -eq 0 ]] || exit 1
