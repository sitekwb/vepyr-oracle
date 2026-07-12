#!/usr/bin/env bash
#SBATCH --job-name=oracle-report
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# Final DAG stage for one version: merge the per-(combo,chrom) diff summaries
# into one whole-genome summary.json (bin/merge_summaries.py), then render the
# validation PDF from it (bin/make_report.py). Not an array job -- one element,
# submitted by slurm/orchestrate.sh with --dependency=afterok:<diff array id>
# so it only runs once EVERY diff shard has actually succeeded (afterok, not
# aftercorr: a partial diff array must not silently produce a "whole-genome"
# report -- bin/merge_summaries.py enforces this independently too, via its
# expected_chroms completeness contract, but afterok stops a doomed report job
# from even starting).
#
# Usage: sbatch --dependency=afterok:<diff array id> report.sbatch <115|116>
set -euo pipefail

VERSION="${1:?usage: sbatch report.sbatch <115|116>}"
case "$VERSION" in
    115|116) ;;
    *) echo "[FATAL] version must be 115 or 116, got: $VERSION" >&2; exit 2 ;;
esac

WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"

RESULTSDIR="$WORK/results_wgs_$VERSION"
REPORTDIR="$RESULTSDIR/reports"
mkdir -p "$REPORTDIR"

SUMMARY="$RESULTSDIR/summary.json"
PDF="$REPORTDIR/validation_${VERSION}_wgs.pdf"

echo "[report] merging per-shard summaries in $RESULTSDIR"
"$PY" "$WORK/bin/merge_summaries.py" "$RESULTSDIR" --chroms 1-22 -o "$SUMMARY"

echo "[report] rendering $PDF"
"$PY" "$WORK/bin/make_report.py" "$SUMMARY" "$PDF" "(vepyr vs VEP $VERSION, whole WGS)"

echo "[report] done -> $PDF"
