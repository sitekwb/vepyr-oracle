#!/bin/bash
#SBATCH --job-name=wgs115-diff
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=1-8
#SBATCH --output=%x-%A_%a.out
# Diff only. SLOW (~1h+: streams a 25GB vepyr VCF against a 15-25GB ground truth over
# NFS) but it only STREAMS -- low memory. So all 8 can run CONCURRENTLY at 12G each
# (96G total) instead of serializing behind the annotate's 48G ask.
set -euo pipefail
source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"; PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
OUT="$WORK/results_wgs_115"
COMBO=$(awk -F'\t' 'NR>1 {print $1}' "$WORK/matrix.tsv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "[diff] $COMBO"
"$PY" "$WORK/bin/validate.py" --combo "$COMBO" --version 115 --vepyr "$OUT/vepyr_${COMBO}.vcf" --outdir "$OUT"
