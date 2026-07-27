#!/usr/bin/env bash
#SBATCH --job-name=verpub116
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=%x-%A_%a.out
#
# Weryfikacja artefaktu publikacyjnego PO fakcie: jeden element arraya = jedno combo.
#
# Bramki w `publish_gt.sh` sprawdzaja plik w momencie, w ktorym powstaje. Ten skrypt
# sprawdza plik taki, jakim LEZY NA DYSKU teraz -- czyli to, co faktycznie poleci do
# odbiorcy. To nie jest ta sama rzecz: miedzy zapisem a wyslaniem stoi NFS, ewentualny
# rerun arraya i czlowiek z `mv`.
#
# TRZY niezalezne bramki:
#   V1  md5sum -c <nazwa>.md5   -- bajty .gz i .tbi zgadzaja sie z suma z jobu
#   V2  tabix -l <nazwa>.gz     -- indeks REALNIE sie otwiera i wylicza kontigi
#   V3  tabix <nazwa>.gz chr:pos-pos -- losowy dostep przez indeks zwraca rekordy
#
# V2 i V3 sa rozdzielone celowo: `tabix -l` czyta sam .tbi, wiec przechodzi nawet dla
# indeksu, ktory nie pasuje juz do pliku danych. V3 dopiero uzywa offsetow z indeksu
# do skoku w .gz, wiec lapie rozjechana pare (.gz z jednego przebiegu, .tbi z innego).
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh" 2>/dev/null \
  || source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
load_apptainer

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"

# Ten sam podzial kontenerow co w publish_gt.sh: htslib z obrazu VEP-a produkowal
# artefakt, bcftools 1.21 jest NIEZALEZNYM czytelnikiem. Tu obraz VEP-a jest uzyty
# tylko do `tabix`, ktorego bcftools.sif nie ma (zweryfikowane, job 65967).
BGZIP_SIF="${VEPYR_VEP116_SIF:-$HOME/bvp/vep116.sif}"
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"

OUT_DIR="${PUBLISH_DIR:-$DATA/publish_gt116_unforked}"

COMBOS=(hgvs_merged hgvs_merged_am hgvs_merged_pick hgvs_merged_pick_allele
        hgvs_merged_pick_allele_gene hgvs_merged_per_gene
        hgvs_merged_flag_pick_allele hgvs_refseq
        everything everything_hgvs)

IDX="${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"
if ! [[ "$IDX" =~ ^[0-9]+$ ]] || [ "$IDX" -lt 1 ] || [ "$IDX" -gt "${#COMBOS[@]}" ]; then
    echo "[FATAL] array index $IDX poza zakresem (1-${#COMBOS[@]})" >&2; exit 1
fi
COMBO="${COMBOS[$((IDX - 1))]}"

NAME=$("$PY" -c \
  "import sys; sys.path.insert(0, '$WORK'); from oracle.publish import drive_filename; print(drive_filename('$COMBO'))")
DST="$OUT_DIR/$NAME"

for f in "$DST" "$DST.tbi" "$DST.md5"; do
    [ -s "$f" ] || { echo "[FATAL] $COMBO: brak lub pusty $f" >&2; exit 1; }
done

echo "[verify] $COMBO -> $NAME"

# --- V1: sumy kontrolne -------------------------------------------------------
( cd "$OUT_DIR" && md5sum -c "$NAME.md5" )
echo "[verify] $COMBO: V1 md5sum -c OK"

# --- V2: indeks sie otwiera ---------------------------------------------------
# `tabix -l` wypisuje kontigi WYLACZNIE z .tbi. Pusta lista = indeks jest niemy;
# publikowalibysmy .tbi, ktory nie indeksuje niczego.
CONTIGS=$(apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" tabix -l "$DST")
N_CONTIGS=$(printf '%s\n' "$CONTIGS" | grep -c . || true)
if [ "$N_CONTIGS" -eq 0 ]; then
    echo "[FATAL] $COMBO: tabix -l nie zwrocil zadnego kontigu" >&2; exit 1
fi
echo "[verify] $COMBO: V2 indeks otwarty, $N_CONTIGS kontigow: $(printf '%s' "$CONTIGS" | tr '\n' ' ')"

# --- V3: losowy dostep przez indeks -------------------------------------------
# Zapytanie o PIERWSZY kontig z indeksu, w pelnym zakresie 1-300Mb (dluzszym niz
# najdluzszy ludzki chromosom), z `head -1`: sprawdzamy, ze skok przez offsety z
# .tbi trafia w realne bajty .gz, nie przeczytawszy 20 GB.
FIRST=$(printf '%s\n' "$CONTIGS" | head -1)
set +o pipefail   # head zamyka pipe -> tabix dostaje SIGPIPE; interesuje nas WYNIK
REC=$(apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" tabix "$DST" "$FIRST:1-300000000" | head -1)
set -o pipefail
if [ -z "$REC" ]; then
    echo "[FATAL] $COMBO: zapytanie tabix $FIRST:1-300000000 nie zwrocilo rekordu -- .tbi nie pasuje do .gz" >&2
    exit 1
fi
echo "[verify] $COMBO: V3 losowy dostep OK, pierwszy rekord: $(printf '%s' "$REC" | cut -f1-5)"

# --- kontrola liczby rekordow, niezaleznym czytelnikiem -----------------------
N_DST=$(apptainer exec -B "$DATA:$DATA" "$BCFTOOLS_SIF" bcftools view -H "$DST" | wc -l)
echo "[verify] $COMBO: $N_DST rekordow (bcftools)"
if [ -n "${EXPECT_RECORDS:-}" ] && [ "$N_DST" -ne "$EXPECT_RECORDS" ]; then
    echo "[FATAL] $COMBO: oczekiwano $EXPECT_RECORDS rekordow, jest $N_DST" >&2; exit 1
fi

echo "[verify] $COMBO: WSZYSTKIE BRAMKI OK"
