#!/bin/bash
#SBATCH --job-name=wgs116-ann
#SBATCH --partition=gpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-8%3
#SBATCH --output=%x-%A_%a.out
# vepyr annotate against the 116 cache. Independent of the ground truth.
set -euo pipefail
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"; PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_wgs_116"; mkdir -p "$OUT"
COMBO=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
"$PY" "$WORK/bin/run_wgs.py" --combo "$COMBO" --version 116 --out "$OUT/vepyr_${COMBO}.vcf" --workers 8
