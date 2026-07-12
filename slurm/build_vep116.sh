#!/usr/bin/env bash
#SBATCH --job-name=build-vep116
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# T12 -- mint $HOME/bvp/vep116.sif from the official ensembl-vep 116.0 image.
#
# HARD VERSION GATE: a cache/version mismatch here would silently corrupt the
# ENTIRE 116 ground truth this pipeline exists to build -- every downstream
# gen_gt116.sh element would run against the wrong VEP binary and nobody would
# notice until the vepyr-vs-VEP diff turned up drift that is actually just a
# version mismatch, not a real finding. So this script does not just pull the
# image and move on: it runs `vep --help` INSIDE the freshly-built container and
# asserts the version banner literally contains "116". If it does not, this
# exits non-zero and refuses to declare the .sif ready -- fix the pull tag, do
# not patch around this check.
#
# Idempotent: if $SIF already exists AND already passes the version gate, this
# is a no-op. (An existing .sif that FAILS the gate is not silently accepted
# either -- see below: a stale/wrong image from a previous bad pull must not
# be able to hide behind "the file already exists".)
set -euo pipefail

module load apptainer/1.5.0

SIF="${VEPYR_VEP116_SIF:-$HOME/bvp/vep116.sif}"
IMAGE="docker://ensemblorg/ensembl-vep:release_116.0"

mkdir -p "$(dirname "$SIF")"

check_version() {
    local out
    out=$(apptainer exec "$SIF" vep --help 2>&1) || {
        echo "[FATAL] 'vep --help' failed inside $SIF -- the image is not usable" >&2
        return 1
    }
    if ! echo "$out" | grep -q '116'; then
        echo "[FATAL] $SIF does not report VEP 116 in its --help banner:" >&2
        echo "$out" | grep -i -m5 'ensembl-vep\|version' >&2 || true
        return 1
    fi
    echo "$out" | grep -i -m2 'ensembl-vep\|version'
    return 0
}

if [ -f "$SIF" ]; then
    echo "[build_vep116] $SIF already exists -- checking its version before reusing it"
    if check_version; then
        echo "[build_vep116] $SIF already OK (VEP 116 confirmed) -- resuming (skip pull)"
        exit 0
    fi
    echo "[build_vep116] existing $SIF FAILED the version gate -- rebuilding it" >&2
    rm -f "$SIF"
fi

echo "[build_vep116] pulling $IMAGE -> $SIF"
apptainer pull "$SIF" "$IMAGE"

echo "=== version check (must report VEP 116) ==="
if ! check_version; then
    echo "[FATAL] freshly-pulled $SIF is NOT VEP 116 -- wrong tag in $IMAGE?" \
         "Refusing to leave a wrongly-versioned .sif in place: a stale/mismatched" \
         "116 ground truth is worse than no ground truth, because it looks correct" \
         "until someone notices the version banner by hand." >&2
    rm -f "$SIF"
    exit 1
fi

echo "[build_vep116] $SIF OK"
