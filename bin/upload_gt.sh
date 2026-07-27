#!/usr/bin/env bash
# Wyslij gotowa dostawe GT do bucketa GCS i ZWERYFIKUJ ja od strony celu.
#
# Uruchamiany RECZNIE na wezle dostepowym (nie przez sbatch): to transfer sieciowy,
# a wezly obliczeniowe nie sa tu wlasciwym miejscem na ruch wychodzacy.
#
#   bin/upload_gt.sh <katalog-lokalny> <prefiks-w-buckecie>
#   bin/upload_gt.sh ~/vepyr/data/publish_gt116_unforked v2-unforked
#
# ------------------------------------------------------------------------------
# `rclone copy`, NIGDY `rclone sync`.
#
# W buckecie lezy JUZ WYSLANA dostawa (`v1-forked/`). `sync` znaczy "zrob cel
# identycznym ze zrodlem", wiec KASUJE w celu wszystko, czego nie ma w zrodle --
# jedna literowka w prefiksie i `sync` sciera cudza, opublikowana dostawe, po
# ktorej nie ma kosza. `copy` nie ma zadnej sciezki kodu, ktora usuwa obiekt w celu.
# Jesli kiedys bedzie trzeba cos w celu skasowac, zrob to osobna, jawna komenda,
# ktora widac w historii -- nie efektem ubocznym wysylki.
# ------------------------------------------------------------------------------
#
# Wysylamy DOKLADNIE cztery rodzaje plikow. `<combo>.vcf.gz.md5` (per-plikowe sumy
# z jobu publikacyjnego) zostaja lokalnie: odbiorca dostaje jeden `MANIFEST.md5`,
# a dwa zestawy sum w jednym folderze to tylko pytanie "ktory jest wazny".
set -euo pipefail

SRC="${1:?uzycie: upload_gt.sh <katalog-lokalny> <prefiks-w-buckecie>}"
PREFIX="${2:?uzycie: upload_gt.sh <katalog-lokalny> <prefiks-w-buckecie>}"

RCLONE="${RCLONE:-$HOME/bin/rclone}"
REMOTE="${GT_REMOTE:-gcs-vep116}"
BUCKET="${GT_BUCKET:-genomic-benchmarking-vep116-gt}"
DST="$REMOTE:$BUCKET/$PREFIX"

FILTERS=(--include "*.vcf.gz" --include "*.vcf.gz.tbi"
         --include "MANIFEST.md5" --include "README.md")

[ -d "$SRC" ] || { echo "[FATAL] brak katalogu $SRC" >&2; exit 1; }
[ -s "$SRC/MANIFEST.md5" ] || { echo "[FATAL] brak $SRC/MANIFEST.md5 -- nie wysylam dostawy bez manifestu" >&2; exit 1; }
[ -s "$SRC/README.md" ] || { echo "[FATAL] brak $SRC/README.md -- nie wysylam dostawy bez opisu" >&2; exit 1; }

echo "[upload] $SRC -> $DST"
"$RCLONE" copy "$SRC" "$DST" "${FILTERS[@]}" \
    --checksum --transfers 4 --stats 30s --stats-one-line -v

# --checksum, a nie domyslne size+mtime: mtime obiektu w GCS pochodzi z metadanych,
# ktore rclone sam zapisal, wiec porownywanie go z lokalnym mtime sprawdza glownie
# to, czy rclone pamieta wlasny zapis. md5 sprawdza BAJTY.
echo "[upload] weryfikacja rclone check --checksum"
"$RCLONE" check "$SRC" "$DST" "${FILTERS[@]}" --checksum

# Druga, niezalezna weryfikacja: md5 ODCZYTANE Z CELU zestawione z MANIFEST.md5.
# `rclone check` porownuje cel ze ZRODLEM; ten krok porownuje cel z DOKUMENTEM,
# ktory dostanie odbiorca -- czyli sprawdza dokladnie to twierdzenie, ktore
# odbiorca zweryfikuje u siebie przez `md5sum -c`.
echo "[upload] weryfikacja MANIFEST.md5 wzgledem zawartosci bucketa"
"$RCLONE" lsjson --hash "$DST" | python3 -c '
import json, sys
remote = {o["Path"]: (o["Size"], o["Hashes"].get("md5")) for o in json.load(sys.stdin)}
manifest = sys.argv[1]
bad, seen = [], set()
with open(manifest) as fh:
    for line in fh:
        if not line.strip():
            continue
        digest, name = line.rstrip("\n").split("  ", 1)
        seen.add(name)
        if name not in remote:
            bad.append(f"{name}: BRAK w buckecie")
        elif remote[name][1] != digest:
            bad.append(f"{name}: md5 w buckecie {remote[name][1]} != manifest {digest}")
extra = sorted(set(remote) - seen - {"MANIFEST.md5", "README.md"})
if extra:
    bad.append(f"obiekty w buckecie spoza manifestu: {extra}")
for b in bad:
    print(f"[FATAL] {b}", file=sys.stderr)
print(f"[upload] {len(seen)} plikow z manifestu, {len(remote)} obiektow w prefiksie")
sys.exit(1 if bad else 0)
' "$SRC/MANIFEST.md5"

echo "[upload] OK: $DST"
