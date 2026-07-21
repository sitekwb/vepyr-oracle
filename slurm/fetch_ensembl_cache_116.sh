#!/usr/bin/env bash
#SBATCH --job-name=fetch-ens116
#SBATCH --partition=cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=%x-%j.out
#
# Natywny (Perl) cache Ensembl VEP 116, wariant ENSEMBL-ONLY -- potrzebny dla combos
# `everything` i `everything_hgvs`. Mamy juz `homo_sapiens_merged` i
# `homo_sapiens_refseq`; ten wariant to trzeci, osobny tarball (~25 GB).
#
# UWAGA: to NIE jest cache parquet vepyr. Perl VEP nie czyta parquetu.
#
# Resume-safe przez sentinel `.extracted`.
set -euo pipefail

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
DEST="$DATA/vep_native_cache_116"
mkdir -p "$DEST"

TARBALL="homo_sapiens_vep_116_GRCh38.tar.gz"
URL="https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/$TARBALL"
SENTINEL="$DEST/homo_sapiens_vep_116_GRCh38.extracted"

if [ -f "$SENTINEL" ]; then
    echo "[fetch-ens116] cache juz rozpakowany, pomijam"
    exit 0
fi

cd "$DEST"
# -C - wznawia przerwane pobranie zamiast zaczynac 25 GB od nowa.
# --retry jak w fetch_vep116_cache.sh: FTP Ensembl potrafi zrywac dlugie transfery.
[ -s "$TARBALL" ] || curl -fSL --retry 5 --retry-delay 20 -C - -o "$TARBALL" "$URL"

echo "[fetch-ens116] rozpakowuje $TARBALL"
tar xzf "$TARBALL"

# GATE: katalog musi istniec i miec >=22 katalogow chromosomow. Pusty katalog po
# udanym `tar` (skrocony tarball) przeszedlby dalej i wywalilby sie dopiero w VEP-ie.
CACHE_DIR="$DEST/homo_sapiens/116_GRCh38"
[ -d "$CACHE_DIR" ] || { echo "[FATAL] brak $CACHE_DIR po rozpakowaniu" >&2; exit 1; }
N_CHROM=$(find "$CACHE_DIR" -maxdepth 1 -type d -name '[0-9]*' | wc -l)
if [ "$N_CHROM" -lt 22 ]; then
    echo "[FATAL] tylko $N_CHROM katalogow chromosomow w $CACHE_DIR (oczekiwano >=22)" >&2
    exit 1
fi
echo "[fetch-ens116] OK: $N_CHROM katalogow chromosomow"
touch "$SENTINEL"
