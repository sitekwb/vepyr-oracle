#!/usr/bin/env bash
#SBATCH --job-name=gt116
#SBATCH --partition=cpu
#SBATCH --time=23:00:00          # under the 24h HARD cap; shards target 12-16h
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=%x-%A_%a.out
#
# T13 -- one array element = one (combo x chrom [x chunk]) row of
# work/shards_gt.tsv = one real-VEP 116 invocation, minting one shard of the
# 116 ground truth. Array shape (--array=...) is set by slurm/orchestrate.sh
# from `wc -l work/shards_gt.tsv` -- NEVER hard-coded here.
#
# Usage: sbatch --array=1-N%K gen_gt116.sh
#
# FLAG FIDELITY: the VEP CLI flags come from matrix.tsv's `vep_flags` column,
# which bin/seed_matrix.py extracted VERBATIM from each combo's EXISTING 115
# ground-truth VCF's `##VEP-command-line=` header -- so the 116 GT is provably
# the SAME invocation as the 115 GT (same --everything/--hgvs/--pick/... flags),
# modulo cache/input/output/fasta/fork, which THIS script re-points at the 116
# cache and our own input/output. Cache/input/output/fasta/vcf/force_overwrite/
# fork are stripped from the recovered flags (via a small inline Python
# filter -- more robust than a regex against an arbitrary real VEP command
# line) and re-supplied explicitly below, so there is never a duplicate or
# conflicting occurrence of any of them.
#
# Output shards carry NO "vepyr_" prefix (they are real-VEP output, not
# vepyr's) -- see oracle/status.py::annotate_out_path's docstring, which is the
# canonical Python-side definition of this naming convention.
#
# Resume-safe: skips if the output VCF already exists and is non-empty.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh" 2>/dev/null \
  || source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
load_apptainer
DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"
VEP_SIF="${VEPYR_VEP116_SIF:-$HOME/bvp/vep116.sif}"
FASTA="$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# CRITICAL: the 116 GT must be built from the SAME input the 115 GT was, or
# the two GTs (and vepyr's output, joined against both) are not comparable to
# each other. Every GT VCF's `##VEP-command-line=` header says
# `--input_file .../HG002_normalized.vcf` -- the raw HG002 benchmark AFTER
# `bcftools norm -m -any` split its 47,781 multi-allelic sites into one row
# per ALT (see slurm/prep_input.sh's header for the record-count proof and
# its GATE 1/GATE 2 verification). The join key downstream is
# (chrom,pos,ref,alt); feeding real VEP the raw benchmark here instead would
# silently drop every multi-allelic site from the comparison, the same way it
# would for vepyr -- see bin/run_wgs.py's WHOLE_WGS_INPUT comment for the full
# explanation. tests/test_input_provenance.py guards against this regressing.
WHOLE_WGS_INPUT="$DATA/HG002_normalized.vcf.gz"

# shellcheck source=lib_common.sh
source "$WORK/slurm/lib_common.sh"

SHARDS="$WORK/shards_gt.tsv"
read_shard_row "$SHARDS" "${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"

if [ "$STEP" != "gt" ]; then
    echo "[FATAL] $SHARDS line $((SLURM_ARRAY_TASK_ID + 1)) has step=$STEP, not" \
         "'gt' -- wrong shards file passed to gen_gt116.sh?" >&2
    exit 2
fi

# matrix.tsv columns (oracle.matrix.MATRIX_HEADER):
#   1=name 2=cache_flavor 3=cache115 4=cache116 5=vepyr_kwargs 6=vep_flags 7=gt115 8=gt116
MATRIX="$WORK/matrix.tsv"
CACHE116=$(awk -F'\t' -v c="$COMBO" 'NR>1 && $1==c {print $4}' "$MATRIX")
FLAGS=$(awk -F'\t' -v c="$COMBO" 'NR>1 && $1==c {print $6}' "$MATRIX")
if [ -z "$CACHE116" ] || [ -z "$FLAGS" ]; then
    echo "[FATAL] combo=$COMBO not found (or has an empty cache116/vep_flags) in" \
         "$MATRIX -- run bin/seed_matrix.py first" >&2
    exit 2
fi
CACHE="$DATA/$CACHE116"

OUTDIR="$DATA/ground_truth_vep_116/shards"
mkdir -p "$OUTDIR"

SFX=""
[ "$CHROM" != "all" ] && SFX="_$CHROM"
[ "$LEVEL" = "L2" ] && SFX="${SFX}_c${CHUNK}"
OUT="$OUTDIR/${COMBO}${SFX}.vcf"

if [ -s "$OUT" ]; then
    echo "[gen_gt116] $OUT already exists -- resuming (skip)"
    exit 0
fi

case "$LEVEL" in
    L0)
        INPUT="$WHOLE_WGS_INPUT"
        ;;
    L1)
        INPUT="$WORK/input/input_${CHROM}.vcf.gz"
        ;;
    L2)
        # See annotate.sbatch's identical L2 branch for the NUM_CHUNKS /
        # cross-combo-mismatch discussion -- this pipeline's real-VEP shards
        # are the ones actually PROJECTED to need L2 (VEP is far slower than
        # vepyr; see the design doc), so this path is the common case here,
        # unlike in annotate.sbatch.
        NUM_CHUNKS=$(awk -F'\t' -v c="$COMBO" -v ch="$CHROM" \
            'NR>1 && $1==c && $4==ch && $3=="L2" {n++} END{print n+0}' "$SHARDS")
        ensure_regions_file "$WORK" "$CHROM" "$NUM_CHUNKS" "$PY"
        REGION=$(region_for_chunk "$WORK" "$CHROM" "$CHUNK")
        ensure_region_input "$WORK" "$CHROM" "$REGION" "$BCFTOOLS_SIF"
        START="${REGION#*:}"; START="${START%-*}"
        END="${REGION##*-}"
        INPUT="$WORK/input/input_${CHROM}_${START}_${END}.vcf.gz"
        ;;
    *)
        echo "[FATAL] unknown level $LEVEL in $SHARDS" >&2
        exit 2
        ;;
esac
if [ ! -s "$INPUT" ]; then
    echo "[FATAL] input $INPUT missing -- run slurm/prep_input.sh first" >&2
    exit 2
fi

# Strip cache/offline/dir_cache/input/output/fasta/vcf/force_overwrite/fork from
# the recovered 115 command line -- all eight are re-supplied explicitly below,
# so a stripped-but-still-present duplicate could never silently override (or
# conflict with) the values THIS script re-points at the 116 cache/our
# input/our output. Everything else (--everything, --hgvs, --pick, --plugin
# AlphaMissense,..., ...) survives verbatim -- that fidelity is the whole point.
# The ground-truth headers were templated: paths appear as the literal string
# "[PATH]". --fasta/--input_file/--output_file are stripped and re-pointed below,
# but --plugin survives VERBATIM (that fidelity is the point), so ITS [PATH] must be
# substituted or VEP is handed a nonexistent AlphaMissense file and the plugin
# silently contributes nothing.
AM_DIR="$DATA/plugin_input/alphamissense"
FLAGS=${FLAGS//\[PATH\]/$AM_DIR}

CLEAN=$("$PY" -c '
import shlex
import sys

STRIP = {"vep", "--cache", "--offline", "--dir_cache", "--dir", "--input_file", "-i",
        "--output_file", "-o", "--fasta", "--fa", "--vcf", "--force_overwrite",
        "--force", "--fork"}
VALUE_FLAGS = {"--dir_cache", "--dir", "--input_file", "-i", "--output_file", "-o",
              "--fasta", "--fa", "--fork"}

toks = shlex.split(sys.argv[1])
out, skip = [], False
for i, t in enumerate(toks):
    if skip:
        skip = False
        continue
    if t in STRIP:
        if t in VALUE_FLAGS and i + 1 < len(toks) and not toks[i + 1].startswith("-"):
            skip = True
        continue
    out.append(t)
print(shlex.join(out))
' "$FLAGS")

# AlphaMissense (hgvs_merged_am) keeps its --plugin flag (CLEAN never strips
# --plugin) but needs its data file mounted -- explicit, even though it already
# lives under $DATA (which is bound wholesale below), to make the dependency
# undeniable rather than accidental.
PLUGIN_BIND=()
if [ "$COMBO" = "hgvs_merged_am" ]; then
    PLUGIN_BIND=(-B "$DATA/plugin_cache:$DATA/plugin_cache")
    echo "[gen_gt116] $COMBO: mounting $DATA/plugin_cache for AlphaMissense"
fi

OUT_TMP="$OUT.tmp"
echo "[gen_gt116] combo=$COMBO level=$LEVEL chrom=$CHROM chunk=$CHUNK" \
     "est_hours=$EST_HOURS cache=$CACHE input=$INPUT -> $OUT"
echo "[gen_gt116] flags (cleaned): $CLEAN"

# shellcheck disable=SC2086  # $CLEAN is a shlex.join()-quoted flag list -- word
# splitting it is the point, and each flag's own quoting survives the split.
apptainer exec -B "$DATA:$DATA" -B "$WORK:$WORK" "${PLUGIN_BIND[@]}" "$VEP_SIF" \
    vep --cache --offline --dir_cache "$CACHE" --fasta "$FASTA" \
        --input_file "$INPUT" --output_file "$OUT_TMP" \
        --vcf --force_overwrite --fork 16 $CLEAN

# Post-run sanity (design doc's verification gate): a valid GT shard must
# carry a CSQ header and at least one non-header record -- an empty or
# truncated output must never be mistaken for a completed shard by a resubmit,
# nor silently feed a wrong-shaped file into merge_gt116.sh's concat.
if ! grep -q '^##INFO=<ID=CSQ' "$OUT_TMP"; then
    echo "[FATAL] $OUT_TMP has no ##INFO=<ID=CSQ header -- VEP did not annotate;" \
         "check the flags above and the job log for errors." >&2
    exit 1
fi
N_RECORDS=$(grep -vc '^#' "$OUT_TMP" || true)
if [ "$N_RECORDS" -eq 0 ]; then
    echo "[FATAL] $OUT_TMP has 0 non-header records -- refusing to accept an empty" \
         "GT shard as complete." >&2
    exit 1
fi

mv -f "$OUT_TMP" "$OUT"                 # atomic on the same filesystem
echo "[gen_gt116] $COMBO chrom=$CHROM chunk=$CHUNK: $N_RECORDS records -> $OUT"
