#!/usr/bin/env bash
# Shared helpers, `source`d by the array .sbatch/.sh scripts. Not standalone
# runnable (it defines functions and does nothing on its own) -- `bash -n` still
# syntax-checks it fine.
#
# Centralising these here is what keeps annotate.sbatch, diff.sbatch and
# gen_gt116.sh (T13) from each re-implementing "read row N of a shards TSV" and
# "materialise an L2 region's input slice" slightly differently -- a drift that
# would otherwise be very easy to introduce silently across three separate bash
# files.
#
# CONVENTION THIS FILE ENCODES (read this before touching any array script):
#   work/shards_<step>.tsv   -- one shards.tsv PER STEP (gt / annotate / diff),
#                                each the output of one bin/plan_shards.py
#                                invocation. Shard GEOMETRY (which combo/chrom/
#                                chunk breakdown) is not version-specific for
#                                gt/annotate/diff, so the SAME shards_annotate.tsv
#                                and shards_diff.tsv drive both --version 115 and
#                                --version 116 runs; only the --version flag
#                                passed to run_wgs.py/validate.py and the output
#                                directory change. shards_gt.tsv exists only for
#                                the 116 half (115's ground truth is static).
#   work/regions_<chrom>.tsv -- one L2 region cut PER CHROMOSOME, shared by every
#                                combo/step that needs L2 chunks of that
#                                chromosome (bin/make_regions.py refuses to
#                                silently reuse a file cut at a different
#                                --num-chunks -- see that script's docstring).

# read_shard_row FILE INDEX
#
# Reads line INDEX+1 of a shards.tsv (INDEX is 1-based, i.e. pass
# $SLURM_ARRAY_TASK_ID directly -- line 1 is the header, so array element 1
# reads data line 2). Exports COMBO, STEP, LEVEL, CHROM, CHUNK, N_VARIANTS,
# EST_HOURS, EST_HOURS_RAW, SAFETY_FACTOR, in oracle.shards.SHARD_HEADER's exact
# column order -- if that header ever changes, this must change with it.
read_shard_row() {
    local file="$1" index="$2"
    local line
    line=$(sed -n "$((index + 1))p" "$file")
    if [ -z "$line" ]; then
        echo "[FATAL] $file: no data row at index $index (array task id out of range" \
             "for this shard file -- check --array against \`wc -l < $file\`)" >&2
        exit 3
    fi
    IFS=$'\t' read -r COMBO STEP LEVEL CHROM CHUNK N_VARIANTS EST_HOURS EST_HOURS_RAW \
        SAFETY_FACTOR <<< "$line"
    export COMBO STEP LEVEL CHROM CHUNK N_VARIANTS EST_HOURS EST_HOURS_RAW SAFETY_FACTOR
}

# ensure_regions_file WORK CHROM NUM_CHUNKS PYTHON
#
# Idempotently (re)builds work/regions_<chrom>.tsv at exactly NUM_CHUNKS chunks
# via bin/make_regions.py, which itself: (a) skips work if a matching file
# already exists, (b) FAILS LOUDLY if an existing file was cut at a different
# chunk count (two combos of the same step disagreeing on how many L2 chunks a
# chromosome needs -- see that script's docstring for why this must not be
# silently papered over).
ensure_regions_file() {
    local work="$1" chrom="$2" num_chunks="$3" py="$4"
    "$py" "$work/bin/make_regions.py" --chrom "$chrom" --num-chunks "$num_chunks"
}

# region_for_chunk WORK CHROM CHUNK
#
# Echoes the "chrom:start-end" spec for CHUNK from work/regions_<chrom>.tsv
# (must already exist -- call ensure_regions_file first).
region_for_chunk() {
    local work="$1" chrom="$2" chunk="$3"
    local regions="$work/regions_$chrom.tsv"
    local spec
    spec=$(awk -F'\t' -v c="$chunk" 'NR>1 && $1==c {print $2}' "$regions")
    if [ -z "$spec" ]; then
        echo "[FATAL] no region for chunk=$chunk in $regions" >&2
        exit 3
    fi
    echo "$spec"
}

# ensure_region_input WORK CHROM REGION BCFTOOLS_SIF
#
# Materialises the pre-sliced L2 input bin/run_wgs.py (and slurm/gen_gt116.sh)
# expect to just READ: work/input/input_<chrom>_<start>_<end>.vcf.gz (+.tbi),
# sliced with `bcftools view -r` from the already-split whole-chromosome input
# (work/input/input_<chrom>.vcf.gz, produced by slurm/prep_input.sh).
#
# Idempotent: skips the slice if the output already exists and is non-empty.
# Staged through a .tmp path so a killed job never leaves a partial slice for
# the next array element (or a resubmit of this one) to read as complete.
ensure_region_input() {
    local work="$1" chrom="$2" region="$3" bcftools_sif="$4"
    local start end out tmp chrom_input
    start="${region#*:}"; start="${start%-*}"
    end="${region##*-}"
    out="$work/input/input_${chrom}_${start}_${end}.vcf.gz"
    chrom_input="$work/input/input_${chrom}.vcf.gz"

    if [ -s "$out" ] && [ -s "$out.tbi" ]; then
        echo "[lib_common] $out already sliced -- resuming (skip)"
        return 0
    fi
    if [ ! -s "$chrom_input" ]; then
        echo "[FATAL] $chrom_input missing -- run slurm/prep_input.sh first" >&2
        exit 2
    fi

    tmp="$out.slice.tmp"
    apptainer exec -B "$work:$work" "$bcftools_sif" \
        bcftools view -r "$region" -O z -o "$tmp" "$chrom_input"
    mv -f "$tmp" "$out"                          # atomic on the same filesystem
    apptainer exec -B "$work:$work" "$bcftools_sif" bcftools index -t -f "$out"
    echo "[lib_common] sliced $chrom_input -r $region -> $out"
}
