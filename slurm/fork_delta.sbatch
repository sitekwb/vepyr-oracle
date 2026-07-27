#!/usr/bin/env bash
#SBATCH --job-name=forkdelta
#SBATCH --partition=cpu
#SBATCH --time=03:00:00          # far under the 24h HARD cap; a chr1 pair runs ~10 min
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=%x-%A_%a.out
#
# One array element = one (combo x chrom) row of work/fork_delta_plan.tsv = one
# UNFORKED-vs-FORKED diff of the 116 ground truth. See bin/fork_delta.py's
# docstring for the experiment and for which arm lands in which diff slot.
#
# Usage:
#   sbatch --array=1-N%K fork_delta.sbatch
# where N = `wc -l < work/fork_delta_plan.tsv` minus the header -- NEVER
# hard-coded here (the same rule slurm/gen_gt116.sh follows for shards_gt.tsv).
#
# The plan file is a 2-column TSV (`combo<TAB>chrom`) with a header line, built by
# whoever scopes the comparison -- the comparable-combo INTERSECTION between the
# two arms is an input to this job, not something it may guess. The unforked arm
# has 10 combos and the forked arm 8 (the matrix grew by two `everything*` combos
# between the runs); a plan naming a combo the forked arm never had makes
# fork_delta.py exit 2 rather than invent a comparison.
#
# INTEGRITY, RECORDED PER ELEMENT: the raw non-header record count of BOTH shards
# goes to work/fork_delta/counts/<combo>_<chrom>.tsv BEFORE the diff. A silent
# generation gap (a shard VEP truncated, or one that never covered the whole
# chromosome) would otherwise masquerade as a fork-delta in only_gt/only_vepyr.
# The counts are cross-checkable against shards_gt.tsv's `n_variants` column,
# which is what the input actually holds for that chromosome.
set -euo pipefail

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"

UNFORKED_ROOT="${VEPYR_GT116_UNFORKED_DIR:-$DATA/ground_truth_vep_116_unforked}"
FORKED_ROOT="${VEPYR_GT116_FORKED_DIR:-$DATA/ground_truth_vep_116}"
OUTDIR="${VEPYR_FORK_DELTA_DIR:-$WORK/fork_delta}"

PLAN="$WORK/fork_delta_plan.tsv"
TASK_ID="${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"
LINE=$(sed -n "$((TASK_ID + 1))p" "$PLAN")
if [ -z "$LINE" ]; then
    echo "[FATAL] $PLAN: no data row at index $TASK_ID (--array is wider than the plan)" >&2
    exit 3
fi
IFS=$'\t' read -r COMBO CHROM <<< "$LINE"

mkdir -p "$OUTDIR/counts"
COUNTS="$OUTDIR/counts/${COMBO}_${CHROM}.tsv"

U="$UNFORKED_ROOT/shards/${COMBO}_${CHROM}.vcf"
F="$FORKED_ROOT/shards/${COMBO}_${CHROM}.vcf"

# Structural integrity, cheap and BEFORE the expensive diff: a shard with no CSQ
# header never annotated, and a shard whose last line is not a complete record is
# a truncated write that a resubmit would happily treat as done.
for f in "$U" "$F"; do
    [ -s "$f" ] || { echo "[FATAL] missing/empty shard $f" >&2; exit 2; }
    grep -q '^##INFO=<ID=CSQ' "$f" || { echo "[FATAL] $f has no CSQ header" >&2; exit 2; }
done

n_u=$(grep -vc '^#' "$U" || true)
n_f=$(grep -vc '^#' "$F" || true)
# `tail -c1` is 0 bytes only for an empty file; a truncated VCF write ends
# mid-line, i.e. WITHOUT the trailing newline.
tail_u=$(tail -c1 "$U" | wc -l)
tail_f=$(tail -c1 "$F" | wc -l)
n_plan=$(awk -F'\t' -v c="$COMBO" -v ch="$CHROM" \
    'NR>1 && $1==c && $4==ch {print $6}' "$WORK/shards_gt.tsv")
printf 'combo\tchrom\tn_unforked\tn_forked\tn_plan\tnewline_terminated_unforked\tnewline_terminated_forked\n' > "$COUNTS"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$COMBO" "$CHROM" "$n_u" "$n_f" "${n_plan:-NA}" "$tail_u" "$tail_f" >> "$COUNTS"
echo "[fork_delta] $COMBO chr$CHROM records: unforked=$n_u forked=$n_f plan=${n_plan:-NA}" \
     "newline_terminated=($tail_u,$tail_f)"

"$PY" "$WORK/bin/fork_delta.py" --combo "$COMBO" --chrom "$CHROM" \
    --outdir "$OUTDIR" \
    --unforked-root "$UNFORKED_ROOT" --forked-root "$FORKED_ROOT"
