#!/usr/bin/env bash
#SBATCH --job-name=oracle-diff
#SBATCH --partition=cpu
#SBATCH --time=23:00:00          # under the 24h HARD cap; shards target 12-16h
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=%x-%A_%a.out
#
# One array element = one row of work/shards_diff.tsv = one `bin/validate.py`
# invocation (diff ONE combo's already-annotated vepyr VCF against ground
# truth). Array shape (--array=...) is set by slurm/orchestrate.sh from
# `wc -l work/shards_diff.tsv` -- NEVER hard-coded here.
#
# Usage: sbatch --array=1-N%K --dependency=aftercorr:<annotate array id> diff.sbatch <115|116>
#
# `aftercorr` is load-bearing: it lets array element i of THIS array start as
# soon as array element i of the annotate array finishes, instead of waiting
# for the whole annotate array to complete -- diff shards are a 1:1 mirror of
# annotate shards by construction (oracle.shards plans them from the same
# chrom_counts), so this is safe and materially shortens the critical path.
#
# The diff step is NEVER L2 (oracle.shards enforces this: summary.merge()
# requires chromosome-atomic shards) -- so, unlike annotate.sbatch, this never
# needs region/chunk handling.
#
# Resume-safe: bin/validate.py itself exits 0 immediately if its summary JSON
# already exists and is non-empty (see that script's docstring).
set -euo pipefail

VERSION="${1:?usage: sbatch --array=1-N diff.sbatch <115|116>}"
case "$VERSION" in
    115|116) ;;
    *) echo "[FATAL] version must be 115 or 116, got: $VERSION" >&2; exit 2 ;;
esac

WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"

# shellcheck source=lib_common.sh
source "$WORK/slurm/lib_common.sh"

SHARDS="$WORK/shards_diff.tsv"
read_shard_row "$SHARDS" "${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"

if [ "$STEP" != "diff" ]; then
    echo "[FATAL] $SHARDS line $((SLURM_ARRAY_TASK_ID + 1)) has step=$STEP, not" \
         "'diff' -- wrong shards file passed to diff.sbatch?" >&2
    exit 2
fi
if [ "$LEVEL" = "L2" ]; then
    echo "[FATAL] $SHARDS line $((SLURM_ARRAY_TASK_ID + 1)) is L2 -- the diff step" \
         "must never be L2 (oracle.shards.plan_combo() already refuses to emit" \
         "this; a shards_diff.tsv containing L2 rows was not produced by this" \
         "pipeline's plan_shards.py)." >&2
    exit 2
fi

RESULTSDIR="$WORK/results_wgs_$VERSION"

# MUST match the naming annotate.sbatch used for this exact (combo, chrom) --
# diff is never L2, so there is never a chunk suffix to account for here.
SFX=""
[ "$CHROM" != "all" ] && SFX="_$CHROM"
VEPYR_OUT="$RESULTSDIR/vepyr_${COMBO}${SFX}.vcf"

if [ ! -s "$VEPYR_OUT" ]; then
    echo "[FATAL] $VEPYR_OUT missing or empty -- the upstream annotate shard for" \
         "combo=$COMBO chrom=$CHROM has not completed. This diff array should be" \
         "submitted with --dependency=aftercorr:<annotate array id>; if it was and" \
         "you are still seeing this, the annotate element itself failed or was" \
         "preempted -- check its job, not this one." >&2
    exit 2
fi

ARGS=(--combo "$COMBO" --version "$VERSION" --vepyr "$VEPYR_OUT" --outdir "$RESULTSDIR")
[ "$CHROM" != "all" ] && ARGS+=(--chrom "$CHROM")

echo "[diff] combo=$COMBO version=$VERSION level=$LEVEL chrom=$CHROM est_hours=$EST_HOURS"
exec "$PY" "$WORK/bin/validate.py" "${ARGS[@]}"
