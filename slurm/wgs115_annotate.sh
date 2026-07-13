#!/bin/bash
#SBATCH --job-name=wgs115-ann
#SBATCH --partition=gpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=1-8%3
#SBATCH --output=%x-%A_%a.out
# Annotate only. ~9.5 min/combo whole-WGS (measured: 567s, 0.000139 s/variant).
# Memory-hungry (vepyr loads the cache) but SHORT. Resumable: run_wgs.py skips if done.
set -euo pipefail
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"; PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_wgs_115"; mkdir -p "$OUT"
COMBO=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "[ann] $COMBO"
"$PY" "$WORK/bin/run_wgs.py" --combo "$COMBO" --version 115 --out "$OUT/vepyr_${COMBO}.vcf" --workers 8
