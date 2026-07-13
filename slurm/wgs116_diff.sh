#!/bin/bash
#SBATCH --job-name=wgs116-diff
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-8
#SBATCH --output=%x-%A_%a.out
# Diff vepyr@116 against the freshly-minted real-VEP-116 ground truth.
# Streams only -> low memory, so all 8 run concurrently.
set -euo pipefail
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"; PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_wgs_116"
COMBO=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
"$PY" "$WORK/bin/validate.py" --combo "$COMBO" --version 116 --vepyr "$OUT/vepyr_${COMBO}.vcf" --outdir "$OUT"
