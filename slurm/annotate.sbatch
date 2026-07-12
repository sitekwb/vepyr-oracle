#!/usr/bin/env bash
#SBATCH --job-name=oracle-annotate
#SBATCH --partition=cpu
#SBATCH --time=23:00:00          # under the 24h HARD cap; shards target 12-16h
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x-%A_%a.out
#
# One array element = one row of work/shards_annotate.tsv = one
# `bin/run_wgs.py` invocation (vepyr annotate, ONE combo, ONE version, ONE
# chrom/region). Array shape (--array=...) is set by slurm/orchestrate.sh from
# `wc -l work/shards_annotate.tsv` -- NEVER hard-coded here.
#
# Usage: sbatch --array=1-N%K annotate.sbatch <115|116>
#
# The SAME shards_annotate.tsv drives both --version 115 and 116 runs (vepyr's
# per-variant cost is combo-driven, not cache-version-driven) -- only the
# --version flag and the output directory (results_wgs_<version>) change.
#
# Resume-safe: bin/run_wgs.py itself exits 0 immediately if --out already
# exists and is non-empty (see that script's docstring) -- a resubmit after
# preemption or the 24h cap costs nothing beyond re-reading one shard row.
set -euo pipefail

VERSION="${1:?usage: sbatch --array=1-N annotate.sbatch <115|116>}"
case "$VERSION" in
    115|116) ;;
    *) echo "[FATAL] version must be 115 or 116, got: $VERSION" >&2; exit 2 ;;
esac

module load apptainer/1.5.0

WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"
WORKERS="${SLURM_CPUS_PER_TASK:-8}"

# shellcheck source=lib_common.sh
source "$WORK/slurm/lib_common.sh"

SHARDS="$WORK/shards_annotate.tsv"
read_shard_row "$SHARDS" "${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"

if [ "$STEP" != "annotate" ]; then
    echo "[FATAL] $SHARDS line $((SLURM_ARRAY_TASK_ID + 1)) has step=$STEP, not" \
         "'annotate' -- wrong shards file passed to annotate.sbatch?" >&2
    exit 2
fi

RESULTSDIR="$WORK/results_wgs_$VERSION"
mkdir -p "$RESULTSDIR"

SFX=""
[ "$CHROM" != "all" ] && SFX="_$CHROM"
[ "$LEVEL" = "L2" ] && SFX="${SFX}_c${CHUNK}"
OUT="$RESULTSDIR/vepyr_${COMBO}${SFX}.vcf"

ARGS=(--combo "$COMBO" --version "$VERSION" --out "$OUT" --workers "$WORKERS")

case "$LEVEL" in
    L0)
        ;;
    L1)
        ARGS+=(--chrom "$CHROM")
        ;;
    L2)
        # NOTE: vepyr WGS annotation is projected to land at L0/L1 (vepyr is far
        # faster than real VEP -- see the design doc's "current projections"),
        # so this path is expected to be rarely exercised in practice. It is
        # still handled correctly: NUM_CHUNKS is THIS COMBO's own L2 chunk
        # count for THIS chromosome in shards_annotate.tsv. If a DIFFERENT
        # combo of the annotate step disagrees with this one on how many L2
        # chunks chromosome $CHROM needs, ensure_regions_file/make_regions.py
        # FAILS LOUDLY rather than silently reusing a wrong-resolution regions
        # file (see bin/make_regions.py's docstring) -- reconcile the shard
        # plan (e.g. re-measure sec_per_variant so combos agree) if that fires.
        NUM_CHUNKS=$(awk -F'\t' -v c="$COMBO" -v ch="$CHROM" \
            'NR>1 && $1==c && $4==ch && $3=="L2" {n++} END{print n+0}' "$SHARDS")
        ensure_regions_file "$WORK" "$CHROM" "$NUM_CHUNKS" "$PY"
        REGION=$(region_for_chunk "$WORK" "$CHROM" "$CHUNK")
        ensure_region_input "$WORK" "$CHROM" "$REGION" "$BCFTOOLS_SIF"
        ARGS+=(--region "$REGION")
        ;;
    *)
        echo "[FATAL] unknown level $LEVEL in $SHARDS" >&2
        exit 2
        ;;
esac

echo "[annotate] combo=$COMBO version=$VERSION level=$LEVEL chrom=$CHROM chunk=$CHUNK" \
     "est_hours=$EST_HOURS -> $OUT"
exec "$PY" "$WORK/bin/run_wgs.py" "${ARGS[@]}"
