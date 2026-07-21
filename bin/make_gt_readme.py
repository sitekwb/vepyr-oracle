#!/usr/bin/env python3
"""Wygeneruj README.md dla folderu GT na Drive.

Flagi VEP NIE sa przepisywane recznie: kazde combo dostaje swoj `##VEP-command-line=`
wyciagniety wprost z opublikowanego VCF-a. Recznie przepisana flaga moze sie rozjechac
z rzeczywistoscia, a odbiorca nie ma jak tego wykryc.

Usage:
    make_gt_readme.py --publish-dir ~/vepyr/data/publish_gt116 \
                      --variant forked --out ~/vepyr/data/publish_gt116/README.md
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Final

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import extract_vep_command_line
from oracle.publish import drive_filename, parse_md5_manifest

COMBOS = ["hgvs_merged", "hgvs_merged_am", "hgvs_merged_pick",
          "hgvs_merged_pick_allele", "hgvs_merged_pick_allele_gene",
          "hgvs_merged_per_gene", "hgvs_merged_flag_pick_allele", "hgvs_refseq"]

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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--publish-dir", required=True)
    ap.add_argument("--variant", choices=["forked", "unforked"], required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    combos = COMBOS if args.variant == "forked" else \
        ["everything", "everything_hgvs", *COMBOS]

    lines = [
        "# Ensembl VEP 116 — ground truth (HG002, GRCh38, chr1-22)",
        "",
        "Odpowiednik folderu `115.2/`, wygenerowany tym samym wejsciem i tymi samymi",
        "flagami, na cache'u 116.",
        "",
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
        "| plik | flagi VEP |",
        "|---|---|",
    ]

    missing: list[str] = []
    for combo in combos:
        name = drive_filename(combo)
        path = os.path.join(args.publish_dir, name)
        if not os.path.exists(path):
            missing.append(name)
            continue
        cmdline = extract_vep_command_line(path) or "(brak naglowka)"
        displayed = cmdline if cmdline == "(brak naglowka)" else semantic_flags(cmdline)
        lines.append(f"| `{name}` | `{displayed}` |")

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
        "> rozni sie od flag (np. `..._pick` powstalo z `--flag_pick_allele_gene`, nie `--pick`),",
        "> **flagi sa autorytatywne**. To nie pomylka: flagi odzyskano dokladnie z naglowka",
        "> Twojego pliku `115.2/`, wiec z definicji ZGADZAJA sie z 115.2 -- nazwa jest",
        "> historyczna, naglowek jest prawda.",
        "",
    ]

    lines += ["", "## Wersje", "",
              "- Ensembl VEP **116** (kontener `vep116.sif`)",
              "- cache: `homo_sapiens_merged/116_GRCh38`, `homo_sapiens_refseq/116_GRCh38`",
              "- FASTA: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`",
              "- AlphaMissense: patrz commit `VEP_plugins` ponizej", ""]

    if args.variant == "forked":
        lines += [FORK_CAVEAT, ""]

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
