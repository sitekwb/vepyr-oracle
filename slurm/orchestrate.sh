#!/usr/bin/env bash
# Submit the full DAG for one version. Array shapes are READ FROM shards.tsv --
# NEVER hard-coded (see slurm/lib_common.sh's header comment for the
# `shards_<step>.tsv` convention this depends on: shards_annotate.tsv and
# shards_diff.tsv drive both versions; shards_gt.tsv exists only for 116).
#
# Usage: orchestrate.sh <115|116>
#
#   115: annotate -> diff (aftercorr) -> report (merge_summaries + make_report)
#
#   116: (gt -> merge_gt116) IN PARALLEL WITH annotate; then diff (which needs
#        BOTH the vepyr annotate output, per-element via aftercorr, AND the
#        fully-merged 116 ground truth, via afterok on the whole merge array)
#        -> report
#
# This script ONLY SUBMITS -- it never polls or waits, and it never runs a job
# inline. Use bin/status.py to check progress after submitting. Every job id
# this script submits is printed, in submission order.
set -euo pipefail

VERSION="${1:?usage: orchestrate.sh <115|116>}"
case "$VERSION" in
    115|116) ;;
    *) echo "[FATAL] version must be 115 or 116, got: $VERSION" >&2; exit 2 ;;
esac

WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
#: How many array elements of one step run concurrently. Overridable per run
#: (cluster fair-share/queue pressure varies); NOT a substitute for the shard
#: COUNT itself, which always comes from the shards.tsv line count below.
THROTTLE="${VEPYR_ARRAY_THROTTLE:-8}"

# data-row count of a TSV (line count minus its header) -- the array size.
nlines() {
    local path="$1" n
    if [ ! -f "$path" ]; then
        echo "[FATAL] missing $path -- run bin/plan_shards.py first" >&2
        exit 2
    fi
    n=$(($(wc -l < "$path") - 1))
    if [ "$n" -lt 1 ]; then
        echo "[FATAL] $path has no data rows" >&2
        exit 2
    fi
    echo "$n"
}

ANNOTATE_SHARDS="$WORK/shards_annotate.tsv"
DIFF_SHARDS="$WORK/shards_diff.tsv"
N_ANNOTATE=$(nlines "$ANNOTATE_SHARDS")
N_DIFF=$(nlines "$DIFF_SHARDS")
echo "[orchestrate] version=$VERSION annotate=$N_ANNOTATE diff=$N_DIFF shards" \
     "(throttle=$THROTTLE)"

ANN=$(sbatch --parsable --array=1-"$N_ANNOTATE"%"$THROTTLE" \
    "$WORK/slurm/annotate.sbatch" "$VERSION")
echo "[orchestrate] submitted annotate array: $ANN (1-$N_ANNOTATE%$THROTTLE)"

# aftercorr: array element i of the diff array may start as soon as element i
# of the annotate array finishes, rather than waiting for the whole annotate
# array -- safe because plan_shards.py plans diff shards 1:1 against the same
# chrom_counts.tsv as annotate shards (see oracle/shards.py).
DIFF_DEPENDENCY="aftercorr:$ANN"

if [ "$VERSION" = "116" ]; then
    GT_SHARDS="$WORK/shards_gt.tsv"
    MATRIX="$WORK/matrix.tsv"
    N_GT=$(nlines "$GT_SHARDS")
    N_COMBOS=$(nlines "$MATRIX")

    GT=$(sbatch --parsable --array=1-"$N_GT"%"$THROTTLE" "$WORK/slurm/gen_gt116.sh")
    echo "[orchestrate] submitted gen_gt116 array: $GT (1-$N_GT%$THROTTLE, runs" \
         "IN PARALLEL with the annotate array above)"

    MERGE=$(sbatch --parsable --dependency=aftercorr:"$GT" \
        --array=1-"$N_COMBOS" "$WORK/slurm/merge_gt116.sh")
    echo "[orchestrate] submitted merge_gt116 array: $MERGE (1-$N_COMBOS)"

    # diff needs BOTH: its own upstream annotate element (aftercorr:$ANN, as in
    # the 115 case) AND the FULLY merged 116 ground truth -- afterok on the
    # whole merge array, because unlike annotate->diff, diff shards are not
    # 1:1 with merge_gt116 elements (merge is per-COMBO, diff can be per-combo
    # x chromosome), so aftercorr would not be a safe per-element pairing here.
    DIFF_DEPENDENCY="aftercorr:$ANN,afterok:$MERGE"
fi

DIFF=$(sbatch --parsable --dependency="$DIFF_DEPENDENCY" \
    --array=1-"$N_DIFF"%"$THROTTLE" "$WORK/slurm/diff.sbatch" "$VERSION")
echo "[orchestrate] submitted diff array: $DIFF (1-$N_DIFF%$THROTTLE," \
     "dependency=$DIFF_DEPENDENCY)"

REPORT=$(sbatch --parsable --dependency=afterok:"$DIFF" "$WORK/slurm/report.sbatch" "$VERSION")
echo "[orchestrate] submitted report: $REPORT (dependency=afterok:$DIFF)"

echo "[orchestrate] DAG submitted for version $VERSION:"
echo "[orchestrate]   annotate=$ANN"
if [ "$VERSION" = "116" ]; then
    echo "[orchestrate]   gen_gt116=$GT"
    echo "[orchestrate]   merge_gt116=$MERGE"
fi
echo "[orchestrate]   diff=$DIFF"
echo "[orchestrate]   report=$REPORT"
