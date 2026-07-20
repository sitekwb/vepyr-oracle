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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import extract_vep_command_line
from oracle.publish import drive_filename, parse_md5_manifest

COMBOS = ["hgvs_merged", "hgvs_merged_am", "hgvs_merged_pick",
          "hgvs_merged_pick_allele", "hgvs_merged_pick_allele_gene",
          "hgvs_merged_per_gene", "hgvs_merged_flag_pick_allele", "hgvs_refseq"]

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
        lines.append(f"| `{name}` | `{cmdline}` |")

    if missing:
        print(f"[FATAL] brakuje opublikowanych plikow: {missing}", file=sys.stderr)
        return 1

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
