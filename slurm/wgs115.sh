#!/bin/bash
#SBATCH --job-name=wgs115
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --array=1-8
#SBATCH --output=%x-%A_%a.out
set -euo pipefail
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_wgs_115"; mkdir -p "$OUT"

COMBO=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "[wgs115] combo=$COMBO (array element $SLURM_ARRAY_TASK_ID)"

"$PY" "$WORK/bin/run_wgs.py"  --combo "$COMBO" --version 115 --out "$OUT/vepyr_${COMBO}.vcf" --workers 8
"$PY" "$WORK/bin/validate.py" --combo "$COMBO" --version 115 --vepyr "$OUT/vepyr_${COMBO}.vcf" --outdir "$OUT"
echo "[wgs115] $COMBO DONE"
