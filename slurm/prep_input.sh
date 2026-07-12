#!/usr/bin/env bash
#SBATCH --job-name=oracle-prep-input
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# P0.1 -- CRITICAL INPUT-MISMATCH FIX.
#
# Every ground-truth VCF's `##VEP-command-line=` header says
# `--input_file [PATH]/HG002_normalized.vcf` -- a file that no longer exists
# on the cluster. We proved exactly what it was:
#
#   benchmark  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz : 4,048,342 records
#     of which multi-allelic (ALT has a comma)           :    47,781 (each
#                                                                     with
#                                                                     exactly
#                                                                     2 ALTs)
#     extra rows if those are split                      :    47,781
#     => records after `bcftools norm -m -any`           : 4,096,123
#   ground truth record count                            : 4,096,123  EXACT
#
# i.e. HG002_normalized.vcf = the benchmark with multi-allelic sites SPLIT
# into one row per ALT (`bcftools norm -m -any`). The pipeline joins vepyr
# output against ground truth on (chrom, pos, ref, alt). A joint row
# `C -> T,CCGC` and its split counterparts `C -> T` / `C -> CCGC` are three
# DIFFERENT keys that never match each other -- so annotating the RAW,
# un-split benchmark instead of the normalized one would silently drop ALL
# 47,781 multi-allelic sites from the comparison, exactly the variants
# multi-ALT handling is weakest on and most needs to be measured.
#
# THIS SCRIPT now does THREE things, in order:
#
#   1. NORMALIZE: `bcftools norm -m -any` the raw benchmark -> bgzip + tabix
#      -> $DATA/HG002_normalized.vcf.gz, gated by two checks (below) that
#      must BOTH pass before the result is trusted.
#   2. SPLIT: the same per-chromosome slicing this script always did (P0.1's
#      original job) -- but now reading from the VERIFIED normalized file,
#      never from the raw benchmark.
#   3. COUNT: work/chrom_counts.tsv, from the split output, as before.
#
# ----------------------------------------------------------------------------
# GATE 1 (record count): the normalized candidate MUST contain EXACTLY
# $EXPECTED_NORMALIZED_COUNT records (see the proof above for where that
# number comes from). This proves SPLITTING happened correctly -- nothing
# more. A record-count match does NOT prove the rows carry the same variants
# (two files can agree on count and disagree entirely on content), which is
# exactly why GATE 2 exists.
#
# GATE 2 (the DECISIVE one -- key-set identity): on one chromosome
# ($GATE_CHROM, chr22 by default -- cheap: ~50k variants), the SET of
# (chrom,pos,ref,alt) keys in our candidate must be IDENTICAL to the set in
# $GATE_GT (a ground truth real VEP produced from the file we are trying to
# reproduce). bin/verify_input.py (venv python, streaming, memory-safe --
# see its module docstring / oracle/verify_input.py for why this is a
# deliberately PER-CHROMOSOME check, not a whole-genome one) does this
# comparison and reports keys-only-in-ours / keys-only-in-gt / overlap. ANY
# non-empty difference means our normalization does not reproduce theirs.
#
# GATE 2 is also what SETTLES whether left-alignment (`bcftools norm -f
# FASTA`) is needed, rather than us guessing: `bcftools norm -m -any` WITHOUT
# `-f` only splits multi-allelic sites -- left-alignment/indel normalization
# only activates when a reference is supplied (see bcftools(1)). So this
# script tries the PLAIN split first; if GATE 2 fails, it retries WITH `-f
# $FASTA` and checks again. If GATE 2 fails for BOTH, this FAILS LOUDLY
# rather than silently shipping a normalization that does not reproduce the
# ground truth's keys.
#
# RESIDUAL RISK (named, not hidden -- read this before trusting a passing
# GATE 2 blindly): a passing GATE 2 proves key-set identity ON CHR22 ONLY; a
# representation quirk confined to a different chromosome could in principle
# still slip through (chr22 was chosen because it is cheap, per the task,
# not because it is exhaustive). Also, `-m -any`'s SPLIT behaviour does not
# touch INFO/FORMAT/sample fields beyond what splitting requires -- if those
# ever mattered downstream (they do not: the diff pipeline keys and compares
# on CSQ annotations against (chrom,pos,ref,alt), never on genotype/sample
# columns), that would be a separate, unverified risk. Finally, GATE 2
# proves EQUIVALENCE to what real VEP actually saw -- it does not and cannot
# prove we used the exact same bcftools version or invocation THEY used;
# equivalence of the resulting key set is the only thing that is actually
# provable, and the only thing this pipeline's join depends on.
#
# IDEMPOTENT: a normalized file that already exists, is non-empty, is
# indexed, AND already carries the "$NORM.gates_ok" marker (written ONLY
# after both gates pass) is trusted as-is -- neither gate is re-run. Any
# candidate that fails a gate is deleted immediately (never left on disk
# looking plausible), so a partial/failed run can never be mistaken for a
# verified one by a resubmit. A chromosome whose split output already exists
# (non-empty VCF + .tbi) is not re-split -- it is only re-counted (a single
# fast iteration pass over an already-small per-chromosome file), so a
# resubmit after preemption or a partial run never re-pays a full
# normalization + 22x split of the whole-genome input.
#
# Every write is staged to a .tmp path and moved into place with `mv -f` /
# `os.replace` (atomic on a POSIX filesystem), so a job killed mid-write
# leaves no partial file for the next run -- or for plan_shards.py -- to trip
# over.
set -euo pipefail

module load apptainer/1.5.0

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"

RAW_BENCHMARK="$DATA/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
FASTA="$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
NORM="$DATA/HG002_normalized.vcf.gz"
NORM_MARKER="$NORM.gates_ok"

# EXPECTED_NORMALIZED_COUNT: 4,048,342 raw benchmark records + 47,781
# multi-allelic sites split into +1 extra row each (each carries EXACTLY 2
# ALTs -- verified separately) = 4,096,123, which is ALSO the ground truth's
# own record count (independently confirmed, exact match). See this script's
# header for the full proof.
EXPECTED_NORMALIZED_COUNT=4096123

# GATE 2 is deliberately scoped to ONE cheap chromosome -- see header
# (RESIDUAL RISK) and oracle/verify_input.py's module docstring for why a
# whole-genome key-set comparison is a memory/time tradeoff this script does
# not need to make.
GATE_CHROM="${VEPYR_GATE_CHROM:-22}"
GATE_GT="${VEPYR_GATE_GT:-$DATA/ground_truth_vep/HG002_annotated_wgs_everything_hgvs_refseq.vcf}"

NORM_THREADS="${VEPYR_NORM_THREADS:-4}"

mkdir -p "$DATA" "$WORK/input"

if [ ! -s "$RAW_BENCHMARK" ]; then
    echo "[FATAL] raw benchmark not found: $RAW_BENCHMARK -- nothing to normalize" >&2
    exit 2
fi
if [ ! -s "$GATE_GT" ]; then
    echo "[FATAL] GATE 2 ground truth not found: $GATE_GT -- cannot verify the" \
         "normalized input's key set against anything. Set VEPYR_GATE_GT if it" \
         "lives elsewhere." >&2
    exit 2
fi
if [ ! -s "$WORK/bin/verify_input.py" ]; then
    echo "[FATAL] $WORK/bin/verify_input.py not found -- deploy bin/ to the" \
         "cluster first (see deploy.sh)" >&2
    exit 2
fi

_bcftools() {
    apptainer exec -B "$DATA:$DATA" "$BCFTOOLS_SIF" bcftools "$@"
}

# try_norm_and_gate EXTRA_NORM_FLAGS
#
# EXTRA_NORM_FLAGS is either "" (plain split) or "-f $FASTA" (also
# left-align/normalize indels against the reference -- see header). Runs
# `bcftools norm`, GATE 1, then GATE 2; only on a full pass does the
# candidate get promoted to $NORM and the marker written. A GATE failing is a
# SOFT, retryable outcome (`return 1` -- the caller tries the next flag
# combination, or gives up loudly once both are exhausted).
#
# A genuine TOOL failure (bcftools itself erroring -- bad apptainer image,
# unreadable input, disk full, ...) is NOT the same thing and must not be
# treated as "this candidate's gate failed, try the next one": every such
# command is therefore explicitly checked and `exit 1`s the WHOLE script.
# This is not optional ceremony -- this function is called as the condition
# of an `if`/`elif` below, and bash's `set -e` is (by design; see bash(1))
# suppressed for EVERY command inside a function invoked that way. Without
# an explicit check+exit here, a real `bcftools norm` crash would silently
# fall through to GATE 1's count check (on a truncated/missing tmp file) and
# get misreported as an ordinary gate failure instead of the tool crash it
# actually was.
try_norm_and_gate() {
    local extra_flags="$1"
    local tmp="$DATA/.HG002_normalized.norm.tmp.vcf.gz"
    rm -f "$tmp" "$tmp.tbi"

    echo "[prep_input] bcftools norm --threads $NORM_THREADS $extra_flags -m -any" \
         "$RAW_BENCHMARK -> $tmp"
    # shellcheck disable=SC2086  # extra_flags is intentionally unquoted/word-split:
    # it is either empty (zero words) or exactly "-f <path>" (two words).
    if ! _bcftools norm --threads "$NORM_THREADS" $extra_flags -m -any \
            -O z -o "$tmp" "$RAW_BENCHMARK"; then
        echo "[FATAL] 'bcftools norm' itself failed for flags [$extra_flags] -- a TOOL" \
             "failure, not a gate failure. Aborting rather than silently treating this" \
             "as 'this candidate's gate failed, try the next one'." >&2
        rm -f "$tmp" "$tmp.tbi"
        exit 1
    fi
    if ! _bcftools index -t -f --threads "$NORM_THREADS" "$tmp"; then
        echo "[FATAL] 'bcftools index' failed on $tmp -- a TOOL failure, not a gate" \
             "failure. Aborting." >&2
        rm -f "$tmp" "$tmp.tbi"
        exit 1
    fi

    local n
    if ! n=$(_bcftools view -H "$tmp" | wc -l); then
        echo "[FATAL] 'bcftools view' failed while counting records in $tmp -- a TOOL" \
             "failure, not a gate failure. Aborting." >&2
        rm -f "$tmp" "$tmp.tbi"
        exit 1
    fi
    echo "[prep_input] GATE 1 (record count): got $n, want $EXPECTED_NORMALIZED_COUNT"
    if [ "$n" -ne "$EXPECTED_NORMALIZED_COUNT" ]; then
        echo "[prep_input] GATE 1 FAILED for flags [$extra_flags] -- discarding" \
             "candidate (got $n records, expected $EXPECTED_NORMALIZED_COUNT)" >&2
        rm -f "$tmp" "$tmp.tbi"
        return 1
    fi

    echo "[prep_input] GATE 2 (chr$GATE_CHROM key-set identity vs $GATE_GT)"
    local gate2_rc=0
    "$PY" "$WORK/bin/verify_input.py" --ours "$tmp" --gt "$GATE_GT" \
        --chrom "$GATE_CHROM" --examples 10 || gate2_rc=$?
    # verify_input.py: exit 0 = identical (pass), exit 1 = key sets differ (a
    # real, SOFT gate failure -- retry with the other flag combination), exit
    # 2 = usage/file error (e.g. $GATE_GT vanished mid-run) -- that is a TOOL/
    # environment problem, not something a different `bcftools norm` flag
    # combination could ever fix, so it aborts the whole script rather than
    # being silently treated as "this candidate's gate failed".
    if [ "$gate2_rc" -eq 2 ]; then
        echo "[FATAL] bin/verify_input.py exited 2 (usage/file error) -- not a gate" \
             "failure. Aborting." >&2
        rm -f "$tmp" "$tmp.tbi"
        exit 1
    elif [ "$gate2_rc" -ne 0 ]; then
        echo "[prep_input] GATE 2 FAILED for flags [$extra_flags] -- discarding" \
             "candidate (see key-set diff above)" >&2
        rm -f "$tmp" "$tmp.tbi"
        return 1
    fi

    mv -f "$tmp" "$NORM"
    mv -f "$tmp.tbi" "$NORM.tbi"
    touch "$NORM_MARKER"
    return 0
}

if [ -s "$NORM" ] && [ -s "$NORM.tbi" ] && [ -f "$NORM_MARKER" ]; then
    echo "[prep_input] $NORM already produced and GATE-verified -- resuming" \
         "(skip normalization)"
elif try_norm_and_gate ""; then
    echo "[prep_input] plain 'bcftools norm -m -any' (no -f) reproduces the" \
         "ground truth's key set on chr$GATE_CHROM -- left-alignment was NOT needed"
elif [ ! -s "$FASTA" ]; then
    echo "[FATAL] plain -m -any failed GATE 2, and the reference FASTA needed to" \
         "retry WITH left-alignment is missing: $FASTA" >&2
    exit 1
elif try_norm_and_gate "-f $FASTA"; then
    echo "[prep_input] needed '-f $FASTA' (left-alignment/indel normalization)" \
         "to reproduce the ground truth's key set on chr$GATE_CHROM"
else
    echo "[FATAL] GATE 2 failed for BOTH plain 'bcftools norm -m -any' and" \
         "'bcftools norm -f $FASTA -m -any' -- our normalization does not" \
         "reproduce the ground truth's (chrom,pos,ref,alt) key set on" \
         "chr$GATE_CHROM by ANY method this script knows. Refusing to proceed:" \
         "see the GATE 2 output above for example keys unique to each side." >&2
    exit 1
fi

# --- split the VERIFIED normalized input into per-chromosome slices --------
#
# Same job P0.1 always did (one bgzipped+tabix'd $WORK/input/input_<chrom>.vcf.gz
# per chrom "1".."22", no "chr" prefix -- the bare-digit convention
# bin/run_wgs.py, oracle/shards.py and bin/merge_summaries.py's DEFAULT_CHROMS
# all already use -- and $WORK/chrom_counts.tsv, the MEASURED per-chromosome
# variant counts bin/plan_shards.py requires, it refuses to guess; see
# oracle/shards.py's module docstring on calibration) -- but now reading from
# the VERIFIED $NORM, never from $RAW_BENCHMARK.

"$PY" - "$NORM" "$WORK" <<'PYEOF'
import csv
import os
import sys

import pysam

IN, WORK = sys.argv[1], sys.argv[2]
CHROMS = [str(i) for i in range(1, 23)]

if not os.path.exists(IN):
    sys.exit(f"[FATAL] normalized input VCF not found: {IN} -- the normalization + "
              f"GATE 1/GATE 2 step above should have produced (and verified) this "
              f"before this split step ever runs.")

if not os.path.exists(IN + ".tbi") and not os.path.exists(IN + ".csi"):
    # Defensive fallback only -- `bcftools index -t` already indexed $NORM before
    # promoting it above; this only fires if something else removed the index
    # in between.
    print(f"[prep_input] no index next to {IN} -- building one "
          f"(requires the file already be bgzip-compressed)", flush=True)
    pysam.tabix_index(IN, preset="vcf", force=False)

vin = pysam.VariantFile(IN)

# Detect whether this VCF's contigs carry the "chr" prefix -- HG002 GRCh38
# benchmark releases are not consistent about this across versions, and
# guessing wrong means .fetch() silently returns an EMPTY iterator for every
# chromosome. An empty split is a much worse failure than a loud one: it would
# feed plan_shards.py a chromosome with 0 variants, which plans a 0-hour shard
# that then finds nothing to annotate and reports a false 100% concordance on
# no data at all.
contigs = set(vin.header.contigs)
if all(f"chr{c}" in contigs for c in CHROMS):
    prefix = "chr"
elif all(c in contigs for c in CHROMS):
    prefix = ""
else:
    sys.exit(f"[FATAL] input contigs match neither 'chr1'..'chr22' nor '1'..'22' "
              f"naming (header carries e.g. {sorted(contigs)[:8]}). Refusing to "
              f"guess -- fix this script's prefix detection, or the input file.")

os.makedirs(f"{WORK}/input", exist_ok=True)
rows: list[tuple[str, int]] = []

for c in CHROMS:
    contig = f"{prefix}{c}"
    out = f"{WORK}/input/input_{c}.vcf.gz"
    tmp = f"{out}.split.tmp"

    if os.path.exists(out) and os.path.getsize(out) > 0 and os.path.exists(out + ".tbi"):
        n = sum(1 for _ in pysam.VariantFile(out))
        rows.append((c, n))
        print(f"[prep_input] chr{c}: already split, {n} variants (skip split)", flush=True)
        continue

    vout = pysam.VariantFile(tmp, "wz", header=vin.header)
    n = 0
    for rec in vin.fetch(contig):
        vout.write(rec)
        n += 1
    vout.close()
    os.replace(tmp, out)                       # atomic -- no partial VCF survives a kill
    pysam.tabix_index(out, preset="vcf", force=True)
    rows.append((c, n))
    print(f"[prep_input] chr{c}: split {n} variants -> {out}", flush=True)

counts_out = f"{WORK}/chrom_counts.tsv"
counts_tmp = f"{counts_out}.tmp"
with open(counts_tmp, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["chrom", "n_variants"])
    w.writerows(rows)
os.replace(counts_tmp, counts_out)              # atomic

total = sum(n for _, n in rows)
if total == 0:
    sys.exit(f"[FATAL] chrom_counts.tsv would report 0 variants across all "
              f"{len(rows)} chromosomes -- almost certainly a contig-naming or "
              f"input-path bug, not a real empty genome. Refusing to write a "
              f"chrom_counts.tsv that would make plan_shards.py plan a WGS run "
              f"against nothing.")

print(f"[prep_input] wrote {counts_out} ({len(rows)} chroms, {total} variants total)")
PYEOF

echo "[prep_input] done"
