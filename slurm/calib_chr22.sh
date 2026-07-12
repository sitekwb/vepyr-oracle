#!/bin/bash
#SBATCH --job-name=oracle-calib22
#SBATCH --partition=cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x-%j.out
#
# Re-validate chr22 through the FIXED pipeline.
# Purpose (three birds):
#   1. measure sec_per_variant PER COMBO -> feeds bin/plan_shards.py (which refuses to guess)
#   2. prove AlphaMissense works on the am combo
#   3. answer whether the chr22 "--pick_allele drift" was a REAL vepyr bug or just our
#      missing --pick_order (the ground truth used biotype,rank,mane_select,... ; we were
#      passing vepyr's default until now)
set -euo pipefail
# NOTE: SLURM COPIES the batch script to /var/spool/slurmd/jobNNN/, so BASH_SOURCE
# points at the spool dir, not work/slurm/. Must fall back to the real path.
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"

WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_chr22_fixed"
mkdir -p "$OUT"

COMBOS=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv")

for c in $COMBOS; do
  echo "=== $c : annotate chr22 ==="
  "$PY" "$WORK/bin/run_wgs.py" --combo "$c" --version 115 --chrom 22 \
        --out "$OUT/vepyr_${c}_22.vcf" --workers 8
  echo "=== $c : diff vs GT ==="
  "$PY" "$WORK/bin/validate.py" --combo "$c" --version 115 --chrom 22 \
        --vepyr "$OUT/vepyr_${c}_22.vcf" --outdir "$OUT"
done

echo "=== merge + report ==="
"$PY" "$WORK/bin/merge_summaries.py" --resultsdir "$OUT" --chroms 22 || true
"$PY" "$WORK/bin/make_report.py" "$OUT/summary.json" "$OUT/validation_chr22_FIXED.pdf" "(chr22, fixed pipeline)" || true
echo "[calib] DONE"
