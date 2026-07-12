#!/bin/bash
#SBATCH --job-name=fetch-vep116-cache
#SBATCH --partition=gpu
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# Download the NATIVE Ensembl VEP 116 cache.
#
# WHY THIS EXISTS: data/116_GRCh38_{merged,refseq} are vepyr's PARQUET caches
# (exon/ motif/ regulatory/ transcript/ ...). Perl VEP CANNOT read parquet -- it
# needs a native cache laid out as homo_sapiens_merged/116_GRCh38/. There is no
# native cache anywhere on this cluster (the 115 ground truth was produced
# elsewhere). Without this, gen_gt116.sh dies with:
#   "ERROR: Cache directory .../116_GRCh38_merged/homo_sapiens_merged not found"
set -euo pipefail
DEST="${VEPYR_DATA:-$HOME/vepyr/data}/vep_native_cache_116"
mkdir -p "$DEST"
cd "$DEST"
BASE="https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache"

for f in homo_sapiens_merged_vep_116_GRCh38.tar.gz homo_sapiens_refseq_vep_116_GRCh38.tar.gz; do
  marker="${f%.tar.gz}.extracted"
  if [ -f "$marker" ]; then echo "[fetch] $f already extracted -- skip"; continue; fi
  if [ ! -s "$f" ]; then
    echo "[fetch] downloading $f ..."
    curl -fL --retry 5 --retry-delay 20 -C - -o "$f" "$BASE/$f"
  fi
  echo "[fetch] extracting $f ..."
  tar -xzf "$f"
  touch "$marker"
  echo "[fetch] $f done"
done

echo "=== native cache layout (what VEP's --dir_cache expects) ==="
ls -d "$DEST"/homo_sapiens* 2>/dev/null
ls "$DEST"/homo_sapiens_merged/ 2>/dev/null | head -3
echo "[fetch] DONE -> point --dir_cache at $DEST"
