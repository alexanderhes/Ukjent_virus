#!/usr/bin/env bash
# Fresh full pipeline run on the local test samples, compared with a baseline.
#
# Everything sample-related stays local, in tests/local/ (gitignored):
#   tests/local/regression.env   settings (template: scripts/regression.env.example)
#   tests/local/baseline/        baseline results (overview, genome tables, genomes)
#   tests/local/runs/<date>/     one folder per run (launch dir, work dir, results)
#
# Usage:
#   scripts/regression_test.sh                 fresh run, compare with baseline
#   scripts/regression_test.sh --set-baseline  fresh run, then make it the baseline
#   scripts/regression_test.sh --compare <results_dir>   compare only, no run
#
# Run it before every commit that changes pipeline code. The run takes a while
# (SPAdes); start it in screen or with nohup.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
LOCAL="$REPO/tests/local"
ENV_FILE="$LOCAL/regression.env"
SET_BASELINE=false
COMPARE_ONLY=""

while [ $# -gt 0 ]; do
    case "$1" in
        --set-baseline) SET_BASELINE=true ;;
        --compare) COMPARE_ONLY="$2"; shift ;;
        -h|--help) sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    shift
done

[ -f "$ENV_FILE" ] || { echo "Missing $ENV_FILE (copy scripts/regression.env.example)" >&2; exit 2; }
# shellcheck disable=SC1090
source "$ENV_FILE"
: "${SAMPLESHEET:?}" "${HOST_INDEX:?}" "${HOST_BLASTDB:?}" "${ESVIRITU_DB:?}" "${LOCAL_CONFIG:?}"
BASELINE="${BASELINE:-$LOCAL/baseline}"
REL_TOL="${REL_TOL:-0.01}"
ABS_TOL="${ABS_TOL:-2}"

compare() {
    if [ ! -d "$BASELINE" ]; then
        echo "No baseline at $BASELINE; run with --set-baseline first." >&2
        exit 2
    fi
    python3 "$REPO/scripts/regression_compare.py" "$BASELINE" "$1" \
        --rel-tol "$REL_TOL" --abs-tol "$ABS_TOL"
}

if [ -n "$COMPARE_ONLY" ]; then
    compare "$COMPARE_ONLY"
    exit $?
fi

RUN_DIR="$LOCAL/runs/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
echo "Fresh run (no -resume) in $RUN_DIR"
echo "Code: $(git -C "$REPO" rev-parse --short HEAD)$(git -C "$REPO" diff --quiet HEAD || echo ' + uncommitted changes')" | tee run_info.txt

# A new launch folder and work dir, and no -resume: nothing is reused.
nextflow -log "$RUN_DIR/nextflow.log" run "$REPO/main.nf" \
    -c "$LOCAL_CONFIG" \
    -w "$RUN_DIR/work" \
    --samplesheet "$SAMPLESHEET" \
    --host_index "$HOST_INDEX" \
    --host_blastdb "$HOST_BLASTDB" \
    --esviritu_db "$ESVIRITU_DB" \
    --outdir "$RUN_DIR/results" \
    --validate ${EXTRA_ARGS:-}

if $SET_BASELINE; then
    rm -rf "$BASELINE"
    mkdir -p "$BASELINE"
    (cd "$RUN_DIR/results" && cp -r --parents overview validation/genome_summary.tsv \
        validation/*/genomes "$BASELINE"/)
    cp run_info.txt "$BASELINE/"
    echo "Baseline set from $RUN_DIR"
else
    compare "$RUN_DIR/results"
fi
