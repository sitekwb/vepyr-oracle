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

# load_apptainer
#
# Makes `apptainer` usable inside a SLURM batch job. Call this instead of a bare
# `module load apptainer`. VERIFIED on compute node h52.
#
# THREE separate traps here, all of which bit us for real (jobs 63397/63398):
#
# 1. `module` is a bash FUNCTION, not a binary. A SLURM batch script is
#    non-login/non-interactive and does NOT source any profile, so `module` is
#    simply undefined -> "module: command not found", exit 127, dead in 1 second.
#    The other scripts on this cluster only work by accident: they were submitted
#    from an INTERACTIVE shell, and --export=ALL carried the exported function into
#    the job. Submit from automation (`ssh host bash -s`) and that vanishes.
#    => source the Lmod init EXPLICITLY. Note the compute nodes have Lmod at
#       /opt/apps/lmod/lmod/init/bash; /etc/profile.d/modules.sh exists ONLY on the
#       login node. Do not trust the login node's layout.
#
# 2. MODULEPATH is UNSET on the compute nodes, so even with `module` defined, Lmod
#    finds nothing ("Lmod Warning: MODULEPATH is undefined"). The apptainer
#    modulefile lives on the node-local SSD: /local/ssd/apps/lmod/lmod/modulefiles
#
# 3. NEVER pipe `module load` (e.g. `module load apptainer | head`). A pipeline runs
#    it in a SUBSHELL, so its PATH mutation is discarded and apptainer silently
#    stays missing. This one wasted a debugging cycle -- the module was loading fine
#    and the test was throwing the result away.
#
# Finally: the modulefile SETS APPTAINER_CACHEDIR=/local/ssd/cache/apptainer, whose
# parent is NOT writable. We override to $SCRATCH UNCONDITIONALLY -- a ${VAR:-default}
# would keep the module's unwritable value.
load_apptainer() {
    if ! declare -F module >/dev/null 2>&1; then
        local init
        # The Lmod init script references unset vars (FPATH); `set -u` would abort.
        set +u
        for init in /opt/apps/lmod/lmod/init/bash \
                    /local/ssd/apps/lmod/lmod/init/bash \
                    /etc/profile.d/modules.sh; do
            # shellcheck disable=SC1090
            [ -f "$init" ] && { source "$init"; break; }
        done
        set -u
    fi

    local mp
    for mp in /local/ssd/apps/lmod/lmod/modulefiles /opt/apps/lmod/lmod/modulefiles; do
        [ -d "$mp" ] && export MODULEPATH="$mp${MODULEPATH:+:$MODULEPATH}"
    done

    # NO pipe, NO subshell -- must mutate THIS shell's PATH (trap 3 above).
    module load apptainer/1.5.0 || module load apptainer || true

    : "${SCRATCH:=/scratch/$USER}"
    export SCRATCH
    export APPTAINER_CACHEDIR="$SCRATCH/.apptainer/cache"
    # PER-JOB tmpdir. apptainer has no squashfuse here, so it converts the .sif into a
    # temp SANDBOX on every exec. With a shared tmpdir, concurrent array elements race
    # extracting the same image into the same tree and die with
    #   "Unable to access rootfs path .../rootfs-NNNN/..." (exit 127).
    # Observed once in 120 concurrent GT elements. A per-job tmpdir removes the race.
    export APPTAINER_TMPDIR="$SCRATCH/.apptainer/tmp/${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
    mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
    trap 'rm -rf "$APPTAINER_TMPDIR" 2>/dev/null || true' EXIT

    if ! command -v apptainer >/dev/null 2>&1; then
        echo "FATAL: apptainer not on PATH after load_apptainer." >&2
        echo "  host=$(hostname) MODULEPATH=${MODULEPATH:-<unset>}" >&2
        echo "  Expected the binary at /local/ssd/apps/apptainer/1.5.0/bin/apptainer" >&2
        echo "  and the modulefile at /local/ssd/apps/lmod/lmod/modulefiles/apptainer/" >&2
        exit 1
    fi
}
