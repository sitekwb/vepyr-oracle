#!/usr/bin/env bash
#SBATCH --job-name=oracle-prep-input
#SBATCH --partition=cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# P0.1 -- split the whole-WGS HG002 benchmark VCF into one bgzipped+tabix'd
# per-chromosome slice: $WORK/input/input_<chrom>.vcf.gz for chrom "1".."22"
# (no "chr" prefix -- the bare-digit convention bin/run_wgs.py, oracle/shards.py
# and bin/merge_summaries.py's DEFAULT_CHROMS all already use), and write
# $WORK/chrom_counts.tsv (chrom, n_variants) -- the MEASURED per-chromosome
# variant counts bin/plan_shards.py requires (it refuses to guess; see
# oracle/shards.py's module docstring on calibration).
#
# IDEMPOTENT: a chromosome whose split output already exists (non-empty VCF +
# .tbi) is not re-split -- it is only re-counted (a single fast iteration pass
# over an already-small per-chromosome file), so a resubmit after preemption or
# a partial run never re-pays a 22x split of the whole-genome input.
#
# Every write is staged to a .tmp path and moved into place with `os.replace`
# (atomic on a POSIX filesystem), so a job killed mid-write leaves no partial
# file for the next run -- or for plan_shards.py -- to trip over.
set -euo pipefail

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"
IN="$DATA/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

mkdir -p "$WORK/input"

"$PY" - "$IN" "$WORK" <<'PYEOF'
import csv
import os
import sys

import pysam

IN, WORK = sys.argv[1], sys.argv[2]
CHROMS = [str(i) for i in range(1, 23)]

if not os.path.exists(IN):
    sys.exit(f"[FATAL] input VCF not found: {IN}")

if not os.path.exists(IN + ".tbi") and not os.path.exists(IN + ".csi"):
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
