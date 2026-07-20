#!/usr/bin/env bash
#SBATCH --job-name=pubgt116
#SBATCH --partition=cpu
#SBATCH --time=04:00:00          # ~30-40 min/combo; zapas na kontencje NFS
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --output=%x-%A_%a.out
#
# Jeden element arraya = jedno combo GT 116: bgzip -> tabix -> md5, pod nazwa
# w konwencji folderu 115.2 (oracle.publish.drive_filename).
#
# --time=04:00:00 celowo NISKIE: wezel h52 backfilluje tylko joby, ktore skoncza
# sie przed startem joba o wyzszym priorytecie. Prosba o 24h zablokowalaby array
# na 2/8 rownoleglych elementow -- patrz reference_ii_hpc_obelix.
#
# Resume-safe: gotowy (.gz + .tbi + .md5, wszystkie niepuste) jest pomijany.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib_common.sh" 2>/dev/null \
  || source "${VEPYR_WORK:-$HOME/vepyr/work}/slurm/lib_common.sh"
load_apptainer

DATA="${VEPYR_DATA:-$HOME/vepyr/data}"
WORK="${VEPYR_WORK:-$HOME/vepyr/work}"
PY="${VEPYR_PYTHON:-$HOME/vepyr/venv/bin/python3}"

# DWA kontenery, celowo -- podzial na "produkuje artefakt" vs "waliduje artefakt".
#
# BGZIP_SIF (vep116.sif, htslib 1.9): bgzip + tabix. bcftools.sif ich NIE MA
# (zweryfikowane, job 65967: `ls /usr/local/bin /usr/bin | grep -E bgzip|tabix`
# zwraca sam `bcftools`). To ten sam obraz, ktorym wygenerowano GT (gen_gt116.sh).
#
# NIE uzywamy `bcftools view -Oz` jako zamiennika bgzip: to nie jest kompresja,
# tylko PARSOWANIE I RESERIALIZACJA VCF-a -- dokleja `##bcftools_viewCommand=` do
# naglowka. Publikujemy to jako ground truth "co wypluil prawdziwy VEP 116";
# naglowek, ktory sami wymyslilismy, jest dokladnie ta roznica, ktora odbiorca
# zobaczy w diffie i bedzie jej szukal przez pol dnia. Md5 ma poswiadczac wyjscie
# VEP-a, a nie nasze przekodowanie tego wyjscia. bgzip rusza wylacznie bajty.
BGZIP_SIF="${VEPYR_VEP116_SIF:-$HOME/bvp/vep116.sif}"

# BCFTOOLS_SIF: tylko GATE 1 i GATE 2, czyli odczyt kontrolny. Walidacja NIEZALEZNA
# implementacja (bcftools 1.21 czyta to, co zapisal bgzip z htslib 1.9) lapie wiecej
# niz sprawdzanie bgzipa samym bgzipem. Ten kontener nie dotyka publikowanych bajtow.
BCFTOOLS_SIF="${VEPYR_BCFTOOLS_SIF:-$HOME/aidiva-onb/bcftools.sif}"

GT_DIR="${GT116_DIR:-$DATA/ground_truth_vep_116}"
OUT_DIR="${PUBLISH_DIR:-$DATA/publish_gt116}"
mkdir -p "$OUT_DIR"

COMBOS=(hgvs_merged hgvs_merged_am hgvs_merged_pick hgvs_merged_pick_allele
        hgvs_merged_pick_allele_gene hgvs_merged_per_gene
        hgvs_merged_flag_pick_allele hgvs_refseq)

IDX="${SLURM_ARRAY_TASK_ID:?must run under a SLURM array}"
# Jawny check granic, NIE ${COMBOS[IDX-1]:?...}: bash 5.2 na wezle (zweryfikowane)
# traktuje ujemny indeks jako offset od KONCA -- IDX=0 -> COMBOS[-1] -> OSTATNI combo,
# cicho, bez bledu. --array=0-8 uruchomiloby wtedy IDX=0 i IDX=8 na TYM SAMYM combo
# rownolegle -> wyscig na tych samych $DST/$DST.tmp/$DST.tbi/$DST.md5. Guard :? nie
# wystrzeliwal (bash: `a=(x y z); echo ${a[-1]}` -> `z`).
if ! [[ "$IDX" =~ ^[0-9]+$ ]] || [ "$IDX" -lt 1 ] || [ "$IDX" -gt "${#COMBOS[@]}" ]; then
    echo "[FATAL] array index $IDX poza zakresem (1-${#COMBOS[@]})" >&2; exit 1
fi
COMBO="${COMBOS[$((IDX - 1))]}"

SRC="$GT_DIR/$COMBO.vcf"
[ -s "$SRC" ] || { echo "[FATAL] brak lub pusty zrodlowy VCF: $SRC" >&2; exit 1; }

NAME=$("$PY" -c \
  "import sys; sys.path.insert(0, '$WORK'); from oracle.publish import drive_filename; print(drive_filename('$COMBO'))")
DST="$OUT_DIR/$NAME"

if [ -s "$DST" ] && [ -s "$DST.tbi" ] && [ -s "$DST.md5" ]; then
    echo "[publish_gt] $COMBO: gotowe, pomijam ($DST)"
    exit 0
fi

echo "[publish_gt] $COMBO: $SRC -> $DST"

# bgzip do .tmp: przerwany job nie moze zostawic pliku, ktory resume uzna za gotowy.
apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" \
    bgzip -c -@ "${SLURM_CPUS_PER_TASK:-16}" "$SRC" > "$DST.tmp"

# GATE 1: naglowek po dekompresji musi byc BAJT W BAJT taki jak zrodlowy.
#
# To nie jest ten sam test co "czy jest CSQ" -- to dowod, ze nie przekodowalismy
# artefaktu. Publikujemy plik opisany jako "wyjscie prawdziwego VEP 116"; gdyby
# cokolwiek dolozylo tu linie (np. `##bcftools_viewCommand=`, gdybysmy kiedys
# podmienili bgzip na `bcftools view -Oz`), odbiorca zobaczylby w diffie roznice,
# ktorej nie ma jak wytlumaczyc inaczej niz czytajac ten skrypt. Gate pilnuje
# tego niezmiennika na stale, nie tylko w smoke runie.
#
# Naglowek zrodla: awk czyta PLIK bezposrednio, wiec wczesny `exit` nie tworzy
# zadnego SIGPIPE.
awk '/^#/{print; next} {exit}' "$SRC" > "$DST.hdr_src"

# Naglowek skompresowanego: tu awk konczy pipe wczesnie, wiec bgzip DOSTAJE
# SIGPIPE (141) i `set -o pipefail` wywrocilby poprawny plik -- ta sama pulapka,
# co przy `grep -q` nizej. Wylaczamy pipefail punktowo i zamiast statusu
# producenta walidujemy WYNIK (niepusty + CSQ + identycznosc), co jest mocniejsze:
# gdyby bgzip naprawde padl, naglowek bylby obciety i cmp by to zlapal.
#
# `< "$DST.tmp"` (STDIN), a NIE `bgzip -dc "$DST.tmp"`: bgzip rozpoznaje wejscie po
# ROZSZERZENIU nazwy i na pliku bez `.gz` konczy sie ostrzezeniem
# "unknown suffix -- ignored", NIC nie wypisujac i wychodzac z kodem 0. Nasza nazwa
# stagingowa `.vcf.gz.tmp` wpada dokladnie w ten przypadek (zlapane w jobie 65971_4).
# Ze stdin nie ma nazwy do sprawdzenia, wiec bgzip po prostu dekompresuje.
set +o pipefail
apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" bgzip -dc < "$DST.tmp" \
    | awk '/^#/{print; next} {exit}' > "$DST.hdr_dst"
set -o pipefail

if [ ! -s "$DST.hdr_dst" ]; then
    echo "[FATAL] $DST.tmp: pusty naglowek po dekompresji" >&2
    rm -f "$DST.tmp" "$DST.hdr_src" "$DST.hdr_dst"; exit 1
fi
# `grep -c` a NIE `grep -q`: przy -q grep konczy sie na pierwszym trafieniu i
# zamyka pipe, przez co producent dostaje SIGPIPE (141). Pod `set -o pipefail`
# caly pipeline jest wtedy niezerowy i gate FALSZYWIE wywalilby poprawny plik.
# Naglowek VEP-a ma setki linii ##contig PO ##INFO, wiec to nie jest teoria.
# Ta sama ostroznosc co przy N_SRC nizej: `grep -c ... || true` polykaloby twardy
# blad grepa jako "0 trafien". Tu ryzyko jest male ($DST.hdr_dst to maly, lokalny
# plik, ktory wlasnie zapisalismy, a cmp nizej i tak jest backstopem), ale wzorzec
# trzymamy jednolity. exit 1 grepa = 0 trafien (poprawne), >=2 = blad.
rc=0
N_CSQ=$(grep -c '^##INFO=<ID=CSQ' "$DST.hdr_dst") || rc=$?
[ "$rc" -le 1 ] || { echo "[FATAL] grep na $DST.hdr_dst padl (exit $rc)" >&2; rm -f "$DST.tmp" "$DST.hdr_src" "$DST.hdr_dst"; exit 1; }
if [ "$N_CSQ" -eq 0 ]; then
    echo "[FATAL] $DST.tmp nie ma naglowka ##INFO=<ID=CSQ" >&2
    rm -f "$DST.tmp" "$DST.hdr_src" "$DST.hdr_dst"; exit 1
fi
if ! cmp -s "$DST.hdr_src" "$DST.hdr_dst"; then
    echo "[FATAL] naglowek NIE jest identyczny ze zrodlowym -- artefakt zostal" \
         "przekodowany, a nie tylko skompresowany. Roznice:" >&2
    diff "$DST.hdr_src" "$DST.hdr_dst" | head -20 >&2
    rm -f "$DST.tmp" "$DST.hdr_src" "$DST.hdr_dst"; exit 1
fi
echo "[publish_gt] $COMBO: naglowek bajt-w-bajt zgodny ($(wc -l < "$DST.hdr_dst") linii)"
rm -f "$DST.hdr_src" "$DST.hdr_dst"

mv "$DST.tmp" "$DST"

apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" tabix -f -p vcf "$DST"

# GATE 2: liczba rekordow po kompresji == przed. bgzip nie powinien gubic linii,
# ale NFS pod obciazeniem potrafi obcinac -- a md5 policzone z obcietego pliku
# jest POPRAWNE dla tego, co lezy, i bezuzyteczne jako dowod kompletnosci.
# grep -vc exit: 0=sa rekordy, 1=zero rekordow (poprawne), >=2=BLAD (plik znikl, NFS
# padl -- dokladnie zagrozenie z komentarza wyzej). Poprzednie `|| true` polykalo
# TAKZE blad -> N_SRC="" -> `[ "" -ne N ]` bleduje ale zwraca nonzero -> if traktuje
# jako false -> GATE 2 (jedyny check dlugosci ciala) POMIJANY na tej samej awarii,
# dla ktorej istnieje. Tolerujemy WYLACZNIE exit 1.
rc=0
N_SRC=$(grep -vc '^#' "$SRC") || rc=$?
[ "$rc" -le 1 ] || { echo "[FATAL] grep na $SRC padl (exit $rc) -- nie licze rekordow z uszkodzonego zrodla" >&2; exit 1; }
N_DST=$(apptainer exec -B "$DATA:$DATA" "$BCFTOOLS_SIF" bcftools view -H "$DST" | wc -l)
if [ "$N_SRC" -ne "$N_DST" ]; then
    echo "[FATAL] liczba rekordow sie rozjechala: zrodlo=$N_SRC skompresowane=$N_DST" >&2
    # Symetrycznie do bramek naglowka wyzej: skasuj kompletnie wygladajacy, ale
    # niezweryfikowany artefakt, zeby nieostrozny downstream (rsync patrzacy tylko
    # na .gz+.tbi) nie wyslal pliku, ktory oblal check liczby rekordow.
    rm -f "$DST" "$DST.tbi"
    exit 1
fi
echo "[publish_gt] $COMBO: $N_DST rekordow, zgadza sie"

# Integralnosc BGZF (w tym blok EOF) TUZ przed md5 -- zamyka okno TOCTOU miedzy
# GATE 2 a md5sum: gdyby NFS obcial $DST po zliczeniu rekordow, md5 opisywalby
# obciete bajty i `md5sum -c` u odbiorcy przeszedlby na uszkodzonym pliku. bgzip -t
# czyta caly plik i wymaga poprawnego markera EOF, ktorego obciety plik nie ma.
apptainer exec -B "$DATA:$DATA" "$BGZIP_SIF" bgzip -t "$DST"

( cd "$OUT_DIR" && md5sum "$NAME" "$NAME.tbi" > "$NAME.md5" )
echo "[publish_gt] $COMBO: gotowe"
cat "$DST.md5"
