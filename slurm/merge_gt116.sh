#!/usr/bin/env bash
#SBATCH --job-name=gt116-merge
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%A_%a.out
#
# T13 -- one array element = one combo. `bcftools concat` (via bcftools.sif)
# that combo's per-chrom/per-chunk shards (slurm/gen_gt116.sh's output) into
# $DATA/ground_truth_vep_116/<combo>.vcf.
#
# Usage: sbatch --array=1-8 --dependency=aftercorr:<gen_gt116 array id> merge_gt116.sh
# (8 = the number of matrix.tsv data rows == len(oracle.matrix.COMBOS); read
# from `wc -l work/matrix.tsv` by slurm/orchestrate.sh, never hard-coded here.)
#
# WHY THE SHARD LIST COMES FROM shards_gt.tsv, NOT A FILENAME GLOB:
# `${COMBO}_*.vcf` would ALSO match a DIFFERENT combo whose name happens to be
# prefixed by this one -- e.g. combo "hgvs_merged" globbing for
# "hgvs_merged_*.vcf" would swallow "hgvs_merged_pick_allele_21.vcf" too (this
# is the exact prefix-collision bug bin/merge_summaries.py's docstring
# documents for the diff-summary side of this pipeline). Reading the shard
# list back out of shards_gt.tsv -- which declares each row's combo
# explicitly -- is unambiguous by construction, the same fix
# bin/merge_summaries.py applies by dispatching on the JSON's declared `name`
# rather than its filename.
#
# INTEGRITY CHECK: merged record count MUST equal the SUM of the shard record
# counts. Fails loudly (and does not write a final `<combo>.vcf`) otherwise --
# a silent record-count drift here would corrupt the 116 ground truth every
# downstream diff treats as authoritative.
#
# Resume-safe: skips if the merged output already exists, is non-empty, AND
# already passes the integrity check (a merged file that exists but FAILS the
# check is not trusted just because it is present -- see below).
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh" 2>/dev/null \
  || source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
load_apptainer
DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"

MATRIX="$WORK/matrix.tsv"
SHARDS="$WORK/shards_gt.tsv"
SHARDDIR="$DATA/ground_truth_vep_116/shards"
OUTDIR="$DATA/ground_truth_vep_116"
mkdir -p "$OUTDIR"

TASK_ID="${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"
COMBO=$(sed -n "$((TASK_ID + 1))p" "$MATRIX" | cut -f1)
if [ -z "$COMBO" ]; then
    echo "[FATAL] $MATRIX: no data row at index $TASK_ID (array task id out of range" \
         "for this matrix -- check --array against \`wc -l < $MATRIX\`)" >&2
    exit 3
fi

OUT="$OUTDIR/${COMBO}.vcf"

count_records() {
    grep -vc '^#' "$1" || true
}

# Build the ordered shard-file list for THIS combo from shards_gt.tsv (never a
# filename glob -- see header comment). Ordered by (chrom numeric, chunk
# numeric) so bcftools concat receives its inputs in genome-coordinate order;
# chrom is always numeric here (1..22 -- the "all" L0 case is a single-shard
# combo and sorts trivially), and chunk is empty for L0/L1 (sorts as 0, which
# never collides because a chrom only ever has ONE L0/L1 row).
SHARD_FILES=()
while IFS= read -r f; do
    SHARD_FILES+=("$f")
done < <(
    awk -F'\t' -v c="$COMBO" 'NR>1 && $1==c {print $4"\t"$5}' "$SHARDS" \
        | sort -t$'\t' -k1,1n -k2,2n \
        | while IFS=$'\t' read -r chrom chunk; do
            sfx=""
            [ "$chrom" != "all" ] && sfx="_$chrom"
            [ -n "$chunk" ] && sfx="${sfx}_c${chunk}"
            echo "$SHARDDIR/${COMBO}${sfx}.vcf"
          done
)

if [ "${#SHARD_FILES[@]}" -eq 0 ]; then
    echo "[FATAL] no rows for combo=$COMBO in $SHARDS -- nothing to merge" >&2
    exit 2
fi

SHARD_SUM=0
for f in "${SHARD_FILES[@]}"; do
    if [ ! -s "$f" ]; then
        echo "[FATAL] shard $f is missing or empty -- gen_gt116.sh has not finished" \
             "combo=$COMBO yet (this merge array should be submitted with" \
             "--dependency=aftercorr:<gen_gt116 array id>)." >&2
        exit 2
    fi
    n=$(count_records "$f")
    SHARD_SUM=$((SHARD_SUM + n))
done

if [ -s "$OUT" ]; then
    GOT=$(count_records "$OUT")
    if [ "$GOT" = "$SHARD_SUM" ]; then
        echo "[merge_gt116] $OUT already exists and matches the shard sum" \
             "($GOT records) -- resuming (skip)"
        exit 0
    fi
    echo "[merge_gt116] existing $OUT has $GOT records, shards now sum to" \
         "$SHARD_SUM -- stale (shard set changed since it was written); re-merging" >&2
fi

LIST="$OUTDIR/.${COMBO}.filelist"
printf '%s\n' "${SHARD_FILES[@]}" > "$LIST"

OUT_TMP="$OUT.tmp"
apptainer exec -B "$DATA:$DATA" "$BCFTOOLS_SIF" \
    bcftools concat -f "$LIST" -o "$OUT_TMP" -O v
rm -f "$LIST"

GOT=$(count_records "$OUT_TMP")
if [ "$GOT" != "$SHARD_SUM" ]; then
    echo "[FATAL] merge integrity check failed for combo=$COMBO: shards sum to" \
         "$SHARD_SUM records, merged output has $GOT. Refusing to publish" \
         "$OUT -- something dropped or duplicated records (a non-abutting" \
         "region cut, an out-of-order concat input, or a partial shard write)." >&2
    rm -f "$OUT_TMP"
    exit 1
fi

mv -f "$OUT_TMP" "$OUT"                 # atomic on the same filesystem
echo "[merge_gt116] $COMBO: ${#SHARD_FILES[@]} shards, $GOT records -> $OUT (integrity OK)"
