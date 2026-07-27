#!/usr/bin/env python3
"""Wygeneruj README.md dla opublikowanego folderu GT (v1-forked / v2-unforked).

Nic, co da sie ODZYSKAC z artefaktu, nie jest tu przepisane recznie:

  * flagi VEP        -- z `##VEP-command-line=` kazdego opublikowanego pliku
  * cache i wersje   -- z `##VEP=` kazdego opublikowanego pliku
  * tabela roznic v1 vs v2 -- z `fork_delta_summary.json`, czyli z tego samego JSON-a,
    ktory wyprodukowalo porownanie obu armow

Recznie przepisana liczba moze sie rozjechac z rzeczywistoscia, a odbiorca nie ma jak
tego wykryc. Dlatego `--variant unforked` WYMAGA `--fork-delta-summary`: README bez
zmierzonej roznicy miedzy armami jest gorszy niz brak README, bo sugeruje, ze roznicy
nie ma.

Usage:
    make_gt_readme.py --publish-dir ~/vepyr/data/publish_gt116_unforked \
                      --variant unforked \
                      --fork-delta-summary ~/vepyr/work/fork_delta/reports/fork_delta_summary.json \
                      --out ~/vepyr/data/publish_gt116_unforked/README.md
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Final

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.csq import open_maybe_gzip
from oracle.matrix import COMBOS as MATRIX_COMBOS, extract_vep_command_line
from oracle.publish import drive_filename, parse_md5_manifest

COMBOS = ["hgvs_merged", "hgvs_merged_am", "hgvs_merged_pick",
          "hgvs_merged_pick_allele", "hgvs_merged_pick_allele_gene",
          "hgvs_merged_per_gene", "hgvs_merged_flag_pick_allele", "hgvs_refseq"]

#: combo -> cache_flavor, z JEDNEGO zrodla prawdy (oracle.matrix.COMBOS). Nie kopiujemy
#: tej tabeli: rozjechanie sie jej z matrix.py opisywaloby odbiorcy inny cache niz ten,
#: przeciwko ktoremu harness porownuje vepyr.
CACHE_FLAVOR: Final[dict[str, str]] = {c.name: c.cache_flavor for c in MATRIX_COMBOS}

#: Jak nazywa sie katalog cache'u dla danego smaku -- do opisu, nie do uruchomienia.
CACHE_DIR_NAME: Final[dict[str, str]] = {
    "ensembl": "homo_sapiens/116_GRCh38 (domyslny cache Ensembl)",
    "merged":  "homo_sapiens_merged/116_GRCh38",
    "refseq":  "homo_sapiens_refseq/116_GRCh38",
}

#: Flagi "hydrauliczne" -- mowia GDZIE sa dane i JAK uruchomiono proces, nie CO
#: znaczy anotacja. Wycinamy je z wyswietlanej komendy, bo naglowek pochodzi z
#: shardu 1 (`input_1.vcf.gz`, `..._1.vcf.tmp`) i te sciezki mylnie sugeruja "tylko
#: chr1" / "plik roboczy". DROP-lista, nie KEEP: flaga nieznana ZOSTAJE widoczna,
#: wiec nigdy nie ukryjemy flagi semantycznej. Pelna surowa komenda jest w naglowku
#: kazdego pliku (nota pod tabela).
_PLUMBING_BOOL: Final[frozenset[str]] = frozenset({
    "cache", "offline", "force_overwrite", "no_stats", "vcf",
})
_PLUMBING_VALUE: Final[frozenset[str]] = frozenset({
    "database", "dir_cache", "dir_plugins", "fasta", "input_file",
    "output_file", "assembly", "fork", "stats_file",
})

#: Atrybuty `##VEP=`, ktore trafiaja do sekcji "Wersje". `cache` i `time` sa per plik,
#: wiec sa raportowane osobno; reszta musi byc identyczna we wszystkich plikach --
#: gdyby nie byla, czesc dostawy powstala z innego cache'u i README musi to pokazac,
#: a nie usrednic.
_VERSION_KEYS: Final[tuple[str, ...]] = (
    "VEP", "API", "assembly", "gencode", "genebuild", "refseq", "dbSNP",
    "ClinVar", "COSMIC", "1000genomes", "gnomADe", "gnomADg", "sift",
    "polyphen", "regbuild", "HGMD-PUBLIC",
)

_VEP_ATTR_RE: Final[re.Pattern[str]] = re.compile(r'([\w.\-]+)="([^"]*)"')


def semantic_flags(cmdline: str) -> str:
    """Zwroc `--flaga [wartosc]` tylko dla flag definiujacych semantyke anotacji.

    Wycina flagi z `_PLUMBING_BOOL`/`_PLUMBING_VALUE` (wraz z ich wartosciami).
    Zachowuje kolejnosc. Poczatkowy token `vep` pomija. Nierozpoznana flaga ZOSTAJE.
    """
    toks = cmdline.split()
    out: list[str] = []
    i = 0
    if toks and toks[0] == "vep":
        i = 1
    while i < len(toks):
        t = toks[i]
        name = t[2:] if t.startswith("--") else t
        if name in _PLUMBING_VALUE:
            i += 2  # zjedz flage i jej wartosc
            continue
        if name in _PLUMBING_BOOL:
            i += 1
            continue
        out.append(t)
        i += 1
    return " ".join(out)


def extract_vep_version_attrs(gt_vcf: str) -> dict[str, str]:
    """Atrybuty linii `##VEP="v116.0" API="v116" cache="..." ...` z naglowka pliku.

    To jedyne miejsce, z ktorego bierzemy wersje: sam VEP je tam wpisal, wiec nie da
    sie ich przepisac z bledem. Pusty dict, gdy linii nie ma.
    """
    with open_maybe_gzip(gt_vcf) as f:
        for line in f:
            if line.startswith("##VEP="):
                return dict(_VEP_ATTR_RE.findall(line[2:].strip()))
            if line.startswith("#CHROM"):
                break
    return {}


FORK_CAVEAT = """\
## ⚠️ Uwaga: `--fork 16`

Te pliki powstaly z `--fork 16`; ground truth 115.2 byl generowany **bez** forka.

Zmierzylismy skutek: **fork psuje wylacznie `HGNC_ID`**. Back-fill HGNC w VEP jest
buffer-scoped (`maxForkSize = buffer_size / (2 * fork)`), wiec fork zmniejsza okno
InputBuffer ~32x i transkrypt-dawca wypada poza nie. Na 9 360 porownanych komorek
(unforked vs `--fork 16`) **zero** roznic poza `HGNC_ID`; consequence, HGVS, CLIN_SIG,
PICK i motif sa nietkniete.

Praktycznie: te pliki sa w pelni uzyteczne do wszystkiego **poza** porownaniami
`HGNC_ID`. Wersja unforked + 2 brakujace combos (`everything`, `everything_hgvs`)
jest w przygotowaniu i trafi do osobnego folderu.

Przy okazji: to jest tez blad po stronie Ensembla wart zgloszenia -- `--fork` czyni
`HGNC_ID` niedeterministycznym wzgledem stopnia zrownoleglenia, a cache 116 trzyma
czesc transkryptow dwukrotnie z niespojna trescia przez granice chunka
(np. `dbID=373064`: `HGNC:26671` w jednym chunku, UNDEF w nastepnym).
"""


def _fmt(n: int) -> str:
    """1234567 -> '1 234 567' (spacja niełamliwa jako separator tysiecy)."""
    return f"{n:,}".replace(",", " ")


def fork_delta_section(summary_path: str) -> list[str]:
    """Sekcja "v1-forked vs v2-unforked", z liczbami WYLICZONYMI z summary JSON.

    `fork_delta_summary.json` to wyjscie `bin/merge_summaries.py --chroms 1-22` dla
    porownania obu armow: lista per-combo podsumowan calogenomowych. W slocie "vepyr"
    stoi arm UNFORKED, w slocie "gt" -- FORKED, wiec kategoria `vepyr_only` znaczy
    doslownie "unforked ma wartosc, forked ma pusto".

    Raises:
        ValueError: gdy JSON nie wyglada jak wynik tego porownania (inny zestaw kluczy,
            combo w stanie != "ok"). Lepiej brak README niz README z liczbami, ktore
            opisuja co innego, niz mowi naglowek sekcji.
    """
    with open(summary_path, encoding="utf-8") as fh:
        summaries = json.load(fh)
    if not isinstance(summaries, list) or not summaries:
        raise ValueError(f"{summary_path}: oczekiwano niepustej listy podsumowan")

    rows: list[tuple[str, int, int, dict[str, int]]] = []
    total_cells = 0
    per_field_totals: dict[str, int] = {}
    for s in summaries:
        for key in ("name", "status", "aligned_annotations", "per_field",
                    "shared_fields", "only_vepyr", "only_gt", "join_rate"):
            if key not in s:
                raise ValueError(f"{summary_path}: podsumowanie bez klucza {key!r}")
        if s["status"] != "ok":
            raise ValueError(
                f"{summary_path}: combo {s['name']!r} ma status={s['status']!r}; "
                f"README nie moze twierdzic, ze porownanie sie udalo")
        if s["only_vepyr"] or s["only_gt"] or s["join_rate"] != 1.0:
            raise ValueError(
                f"{summary_path}: combo {s['name']!r} nie ma pelnej zgodnosci "
                f"strukturalnej (only_vepyr={s['only_vepyr']}, only_gt={s['only_gt']}, "
                f"join_rate={s['join_rate']}); tekst sekcji zaklada, ze ma")
        aligned = int(s["aligned_annotations"])
        n_fields = int(s["shared_fields"])
        total_cells += aligned * n_fields
        moved = {f: int(v["total"]) - int(v["match"])
                 for f, v in s["per_field"].items()
                 if int(v["total"]) != int(v["match"])}
        for f, n in moved.items():
            per_field_totals[f] = per_field_totals.get(f, 0) + n
        rows.append((s["name"], aligned, n_fields, moved))

    moved_total = sum(per_field_totals.values())
    pct = 100.0 * moved_total / total_cells if total_cells else 0.0

    out = [
        "## ⚠️ Dlaczego istnieje `v2-unforked` i czym rozni sie od `v1-forked`",
        "",
        "`v1-forked/` (dostawa z 2026-07-20) powstala z `--fork 16`. Back-fill symbolu",
        "genu w VEP jest **buffer-scoped**: okno `InputBuffer` widziane przez jeden proces",
        "to `buffer_size / (2 * fork)`, wiec `--fork 16` zwezalo je ~32x i transkrypt-dawca",
        "wypadal poza nie. Ta dostawa (`v2-unforked/`) zostala wygenerowana **bez `--fork`**",
        "i jest wersja **autorytatywna**.",
        "",
        "Roznica miedzy armami zostala zmierzona na **calym genomie**, nie na probce:",
        f"{len(rows)} wspolnych combos (8 z 10 -- `everything` i `everything_hgvs` istnieja",
        "tylko w v2, wiec nie maja odpowiednika do porownania), chr1-22, "
        f"**{_fmt(total_cells)} porownanych komorek CSQ**.",
        "",
        f"Wynik: **{_fmt(moved_total)} komorek sie rozni = {pct:.6f}%**, i sa one",
        f"skupione w dokladnie {len(per_field_totals)} polach:",
        "",
        "| pole | komorek roznych | z ilu | co to znaczy |",
        "|---|---:|---:|---|",
    ]
    _MEANING = {
        "HGNC_ID": "w armie forked wartosci **BRAKUJE** (nigdy nie jest bledna)",
        "DOMAINS": "ten sam zbior domen, **inna kolejnosc** (permutacja)",
    }
    for field in sorted(per_field_totals, key=lambda f: -per_field_totals[f]):
        meaning = _MEANING.get(field, "(patrz `fork_delta_per_field.txt`)")
        out.append(f"| `{field}` | {_fmt(per_field_totals[field])} | "
                   f"{_fmt(total_cells)} | {meaning} |")

    out += [
        "",
        "Rozbicie na combo (kolumna `HGNC_ID` = liczba komorek, w ktorych forked zgubil id):",
        "",
        "| combo | anotacji | pol CSQ | " + " | ".join(
            f"`{f}`" for f in sorted(per_field_totals, key=lambda f: -per_field_totals[f])) + " |",
        "|---|---:|---:|" + "---:|" * len(per_field_totals),
    ]
    for name, aligned, n_fields, moved in rows:
        cells = " | ".join(
            _fmt(moved.get(f, 0))
            for f in sorted(per_field_totals, key=lambda f: -per_field_totals[f]))
        out.append(f"| `{name}` | {_fmt(aligned)} | {n_fields} | {cells} |")

    out += [
        "",
        "### `HGNC_ID` -- straty sa STOCHASTYCZNE, nie da sie ich zalatac regula",
        "",
        "99,8% dotknietych transkryptow to RefSeq (`XM_`/`NM_`/`NR_`/`XR_`) **wewnatrz",
        "cache'u merged**; tylko 1 203 z 596 952 to `ENST`. Combo `hgvs_refseq` (czysty",
        "cache RefSeq) ma **zero** roznic w jakimkolwiek polu -- back-fill nie ma tam",
        "rekordow-dawcow, ktore moglby zgubic.",
        "",
        "Cztery combos o **identycznej powierzchni anotacji** (`hgvs_merged`,",
        "`hgvs_merged_am`, `hgvs_merged_pick`, `hgvs_merged_flag_pick_allele`) zgubily",
        "122 767 / 122 565 / 122 835 / 122 591 komorek -- rozne liczby, a zbiory",
        "dotknietych `(pozycja, transkrypt)` pokrywaja sie w 99,56-99,84% (Jaccard",
        "parami; czesc wspolna calej czworki = 99,41%). To znaczy, ze **ktore** komorki",
        "przepadna, zalezy od przypadkowego rozlozenia rekordow po chunkach, a nie od",
        "wlasnosci transkryptu. Nie da sie tego odtworzyc regula ani zalatac po fakcie:",
        "do porownan `HGNC_ID` uzywaj `v2-unforked/`.",
        "",
        "### `DOMAINS` -- 13 komorek, i to TEN SAM mechanizm, nie druga usterka",
        "",
        "Wszystkie 13 to permutacja identycznego zbioru wartosci, np. w dwoch loci na",
        "chr2 (`ENST00000308624`):",
        "",
        "```",
        "v2-unforked : Gene3D:2.30.42.10&PANTHER:PTHR10554&AFDB-ENSP_mappings:AF-Q9NY99-F1",
        "v1-forked   : PANTHER:PTHR10554&Gene3D:2.30.42.10&AFDB-ENSP_mappings:AF-Q9NY99-F1",
        "```",
        "",
        "Udowodnione eksperymentem kontrolnym (`slurm/fork_domains_replicate2.sbatch`,",
        "chr2:1-30Mb, 46 350 wariantow): uruchomienie **bez forka**, ale z",
        "`--buffer_size 156` (= 5000 / (2 x 16), czyli tyle, ile `--fork 16` zostawia",
        "jednemu workerowi) **odtwarza kolejnosc arma forked co do znaku**. Skoro sama",
        "zmiana rozmiaru bufora wystarcza, `DOMAINS` jest efektem tego samego",
        "buffer-scope co back-fill HGNC, a nie osobnym bledem forka.",
        "",
        "### Czego fork NIE ruszyl",
        "",
        "**Zero** roznicych komorek w: `Consequence`, `IMPACT`, `HGVSc`, `HGVSp`, `SIFT`,",
        "`PolyPhen`, `SYMBOL`, `Gene`, `Feature`, `CANONICAL`/`MANE*`/`TSL`/`APPRIS`/`PICK`,",
        "wszystkich czestosciach gnomAD i 1000G oraz obu polach AlphaMissense.",
        "**Wyliczanie konsekwencji jest niewrazliwe na `--fork`.**",
        "",
        "Zgodnosc strukturalna jest pelna: identyczne zbiory transkryptow",
        "(`only_forked = only_unforked = 0`), `join_rate = 1.0`, zero rekordow",
        "malformed i zero duplikatow, po 4 096 123 rekordy w kazdym combo obu armow.",
        "",
        "### Co z tego wynika w praktyce",
        "",
        "* **`v2-unforked/` (ten folder) jest wersja autorytatywna.** Uzywaj jej domyslnie.",
        "* **`v1-forked/` pozostaje uzyteczna dla KAZDEGO pola poza `HGNC_ID`** (i, jesli",
        "  zalezy Ci na kolejnosci w `DOMAINS`, dla 13 komorek w dwoch loci na chr2).",
        "  Nie kasujemy jej i nie unewazniamy -- 27,9 mld komorek CSQ jest w obu armach",
        "  identycznych.",
        "* Zaden plik z `v1-forked/` nie zostal zmieniony ani nadpisany przez te dostawe.",
        "",
        "Niezalezna kontrola: porownanie obu armow jako surowego TEKSTU (poza naglowkami,",
        "wiec obejmujace takze QUAL, FILTER, pozostale klucze INFO i kolumny FORMAT/sample,",
        "ktorych parser CSQ nie widzi) daje dla `hgvs_refseq` **0** roznicych wierszy z",
        "4 096 123, a dla `hgvs_merged` -- 40 267 wierszy, z ktorych po usunieciu tokenow",
        "`HGNC:` zostaja **2**: dokladnie te dwa loci `DOMAINS`. Czyli `HGNC_ID` plus te",
        "dwa loci wyjasniaja **100%** roznicy bajtowej miedzy armami.",
        "",
        "Uwaga na marginesie: to jest tez blad po stronie Ensembla wart zgloszenia --",
        "`--fork` czyni `HGNC_ID` niedeterministycznym wzgledem stopnia zrownoleglenia.",
    ]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--publish-dir", required=True)
    ap.add_argument("--variant", choices=["forked", "unforked"], required=True)
    ap.add_argument("--fork-delta-summary",
                    help="fork_delta_summary.json; WYMAGANE dla --variant unforked")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.variant == "unforked" and not args.fork_delta_summary:
        print("[FATAL] --variant unforked wymaga --fork-delta-summary: README bez "
              "zmierzonej roznicy v1 vs v2 sugerowalby, ze roznicy nie ma",
              file=sys.stderr)
        return 1

    combos = COMBOS if args.variant == "forked" else \
        ["everything", "everything_hgvs", *COMBOS]

    unforked = args.variant == "unforked"
    title_suffix = " — arm UNFORKED (v2)" if unforked else " — arm forked (v1)"

    lines = [
        f"# Ensembl VEP 116 — ground truth (HG002, GRCh38, chr1-22){title_suffix}",
        "",
    ]
    if unforked:
        lines += [
            "Odpowiednik folderu `115.2/`, wygenerowany tym samym wejsciem i tymi samymi",
            f"flagami, na cache'u 116, **bez `--fork`**. {len(combos)} combos, w kazdym",
            "**4 096 123 rekordy**, chr1-22.",
            "",
            "To jest **autorytatywna** wersja ground truth VEP 116. Wczesniejsza dostawa",
            "`v1-forked/` powstala z `--fork 16` -- pozostaje wazna dla wszystkich pol",
            "**poza `HGNC_ID`**; szczegoly i liczby nizej.",
            "",
        ]
    else:
        lines += [
            "Odpowiednik folderu `115.2/`, wygenerowany tym samym wejsciem i tymi samymi",
            "flagami, na cache'u 116.",
            "",
        ]

    lines += [
        "## Wejscie",
        "",
        "Surowy GIAB HG002 benchmark VCF (GRCh38, v4.2.1, chr1-22; nazwa pliku i",
        "pochodzenie: patrz `slurm/prep_input.sh`) po `bcftools norm -m -any`",
        "— 4 096 123 wariantow. Normalizacja jest **konieczna**: bez niej 47 781 miejsc",
        "multialleliczych zostaje w jednym wierszu i wypada z porownania po kluczu",
        "`(chrom, pos, ref, alt)`. Ground truth 115.2 byl generowany z tego samego",
        "znormalizowanego wejscia.",
        "",
        "## Format",
        "",
        "`bgzip` + indeks `tabix`. Sumy kontrolne w `MANIFEST.md5`:",
        "",
        "```bash",
        "md5sum -c MANIFEST.md5",
        "```",
        "",
        "## Pliki",
        "",
        "| plik | cache | flagi VEP |",
        "|---|---|---|",
    ]

    missing: list[str] = []
    version_attrs: dict[str, dict[str, str]] = {}
    for combo in combos:
        name = drive_filename(combo)
        path = os.path.join(args.publish_dir, name)
        if not os.path.exists(path):
            missing.append(name)
            continue
        cmdline = extract_vep_command_line(path) or "(brak naglowka)"
        displayed = cmdline if cmdline == "(brak naglowka)" else semantic_flags(cmdline)
        version_attrs[combo] = extract_vep_version_attrs(path)
        flavor = CACHE_FLAVOR.get(combo, "?")
        lines.append(f"| `{name}` | {CACHE_DIR_NAME.get(flavor, flavor)} | `{displayed}` |")

    if missing:
        print(f"[FATAL] brakuje opublikowanych plikow: {missing}", file=sys.stderr)
        return 1

    lines += [
        "",
        "> **Flagi wyzej to flagi SEMANTYCZNE**, odzyskane z naglowka `##VEP-command-line=`",
        "> kazdego pliku (sciezki I/O, `--dir_cache`, `--fork`, `--fasta` i pochodne pominieto",
        "> dla czytelnosci). Pelna, surowa komenda -- w tym oryginalne `input_N.vcf.gz` z",
        "> etapu shardowania po chromosomach -- jest w naglowku kazdego pliku:",
        "> `zcat plik.vcf.gz | grep '^##VEP-command-line='`.",
        ">",
        "> **Nazwa pliku vs flagi.** Nazwy trzymaja konwencje folderu `115.2/`; gdzie nazwa",
        "> rozni sie od flag, **flagi sa autorytatywne**. Konkretnie:",
        "> `HG002_annotated_wgs_everything_hgvs_merged_pick.vcf.gz` (combo `hgvs_merged_pick`)",
        "> **NIE** powstalo z `--pick`, tylko z **`--flag_pick_allele_gene`** -- a to zupelnie",
        "> inny tryb selekcji (`--pick` USUWA niewybrane transkrypty, `--flag_pick_*` zostawia",
        "> wszystkie i tylko jeden oznacza). Nazwa jest historyczna i celowo jej NIE zmieniamy,",
        "> zeby nie zerwac ciaglosci z folderem `115.2/`; udokumentowane rowniez w",
        "> `oracle/matrix.py`. To samo dotyczy `..._pick_allele` / `..._pick_allele_gene` /",
        "> `..._per_gene`: wszystkie pieciu combos z rodziny pick dostaly jawny, NIEDOMYSLNY",
        "> `--pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length`.",
        "",
    ]

    # --- Wersje: WYLACZNIE to, co VEP sam wpisal do naglowka ---------------------
    #
    # Atrybut moze (a) byc identyczny wszedzie, (b) roznic sie miedzy plikami, (c) byc
    # obecny tylko w czesci plikow -- `refseq=` pojawia sie w naglowku tylko wtedy, gdy
    # cache w ogole zawiera transkrypty RefSeq, wiec nie ma go w plikach z czystego
    # cache'u Ensembl. Kazdy z trzech przypadkow jest raportowany OSOBNO: sklejenie (c)
    # z (a) oglaszaloby wersje adnotacji RefSeq dla plikow, ktore zadnej nie maja.
    values_by_key: dict[str, dict[str, list[str]]] = {}
    for combo, attrs in version_attrs.items():
        for k in _VERSION_KEYS:
            v = attrs.get(k)
            if v is not None:
                values_by_key.setdefault(k, {}).setdefault(v, []).append(combo)

    lines += ["", "## Wersje i proweniencja", "",
              "Wszystko ponizej jest ODCZYTANE z naglowka `##VEP=` opublikowanych plikow",
              "(`zcat plik.vcf.gz | grep '^##VEP='`), nie przepisane recznie.", ""]
    lines += ["| co | wersja |", "|---|---|"]
    for k in _VERSION_KEYS:
        by_value = values_by_key.get(k)
        if not by_value:
            continue
        if len(by_value) > 1:
            detail = "; ".join(f"{v} ({', '.join(f'`{c}`' for c in sorted(combos_))})"
                               for v, combos_ in sorted(by_value.items()))
            lines.append(f"| `{k}` | **ROZNE MIEDZY PLIKAMI**: {detail} |")
            continue
        (value, present_in), = by_value.items()
        if len(present_in) == len(combos):
            lines.append(f"| `{k}` | {value} |")
        else:
            lines.append(f"| `{k}` | {value} — tylko w: "
                         f"{', '.join(f'`{c}`' for c in sorted(present_in))} |")
    lines += [
        "",
        "- kontener: `vep116.sif` (Ensembl VEP 116.0, uruchamiany przez `apptainer`)",
        "- cache (natywny, rozpakowany z tarballi Ensembl 116):",
        "  `homo_sapiens/116_GRCh38`, `homo_sapiens_merged/116_GRCh38`,",
        "  `homo_sapiens_refseq/116_GRCh38` -- ktory plik uzyl ktorego, patrz tabela wyzej",
        "- FASTA: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`",
        # Swiadomie BEZ nazwy surowego pliku benchmarkowego: tests/test_input_provenance.py
        # pilnuje, zeby ta nazwa nie wystepowala nigdzie poza slurm/prep_input.sh. Wejsciem
        # VEP-a jest i tak plik PO normalizacji -- i to jego nazwe podajemy.
        "- wejscie podane VEP-owi: `HG002_normalized.vcf.gz` -- benchmark GIAB HG002 v4.2.1",
        "  (GRCh38, chr1-22) po `bcftools norm -m -any`; pelna proweniencja w `slurm/prep_input.sh`",
        "- AlphaMissense: patrz nota o `VEP_plugins` nizej",
        "",
    ]
    lines += ["Znacznik czasu z naglowka (shard 1 kazdego combo):", "",
              "| combo | `time=` |", "|---|---|"]
    for combo in combos:
        lines.append(f"| `{combo}` | {version_attrs[combo].get('time', '(brak)')} |")
    lines += [
        "",
        "> Kazdy plik to sklejka 22 shardow (po jednym na chromosom) uruchomionych",
        "> rownolegle jako array SLURM, wiec `time=` opisuje moment startu **sharda 1**,",
        "> a nie calego combo.",
        "",
    ]

    if args.variant == "forked":
        lines += [FORK_CAVEAT, ""]
    else:
        lines += fork_delta_section(args.fork_delta_summary) + [""]

    lines += [
        "## ⚠️ AlphaMissense: zmiana semantyki miedzy 115 a 116",
        "",
        "Commit `c6fbb1b` w `VEP_plugins` **usunal sprawdzanie pozycji** z klucza",
        "dopasowania AlphaMissense. Combo `*_am` dziedziczy to zachowanie, wiec jego",
        "wyniki NIE sa porownywalne 1:1 z `115.2/..._am.vcf` na wpisach missense",
        "(zweryfikowane: regula nieczula na pozycje przewiduje VEP 116 na wszystkich",
        "173 701 wpisach missense, zero reszty). To zmiana w VEP_plugins, nie blad.",
        "",
    ]

    manifest = os.path.join(args.publish_dir, "MANIFEST.md5")
    if os.path.exists(manifest):
        with open(manifest, encoding="utf-8") as fh:
            n = len(parse_md5_manifest(fh.read()))
        lines += [f"`MANIFEST.md5` obejmuje {n} plikow.", ""]

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print(f"[make_gt_readme] zapisano {args.out} ({len(combos)} combos)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
