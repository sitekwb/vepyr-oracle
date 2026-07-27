#!/usr/bin/env python3
"""Zloz `MANIFEST.md5` z per-plikowych `<nazwa>.md5` wyprodukowanych przez publish_gt.sh.

Manifest NIE liczy md5 od nowa. Sumy powstaja w jobie `slurm/publish_gt.sh`, tuz po
`bgzip -t` (bramka integralnosci BGZF), na tym samym wezle, ktory zapisal plik --
przeliczanie ich pozniej, na loginie, po NFS-ie, dokladalo by tylko drugie zrodlo
prawdy i drugi sposob na ciche rozjechanie sie z rzeczywistoscia. Weryfikacja
`md5sum -c` (czyli: czy to, co lezy na dysku, nadal zgadza sie z manifestem) to
osobny krok -- `slurm/verify_publish.sh`.

Kolejnosc wierszy jest deterministyczna (kolejnosc combos), zeby dwa przebiegi na
tym samym katalogu dawaly bajt w bajt ten sam manifest.

Usage:
    make_manifest.py --publish-dir ~/vepyr/data/publish_gt116_unforked \
                     --variant unforked --out .../MANIFEST.md5
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.publish import drive_filename, format_md5_manifest, parse_md5_manifest

#: Te same listy, co w bin/make_gt_readme.py -- forked (v1) mial 8 combos,
#: unforked (v2) ma dodatkowo dwa na czystym cache'u Ensembl.
COMBOS_FORKED = ["hgvs_merged", "hgvs_merged_am", "hgvs_merged_pick",
                 "hgvs_merged_pick_allele", "hgvs_merged_pick_allele_gene",
                 "hgvs_merged_per_gene", "hgvs_merged_flag_pick_allele",
                 "hgvs_refseq"]
COMBOS_UNFORKED = ["everything", "everything_hgvs", *COMBOS_FORKED]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--publish-dir", required=True)
    ap.add_argument("--variant", choices=["forked", "unforked"], required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    combos = COMBOS_FORKED if args.variant == "forked" else COMBOS_UNFORKED

    entries: list[tuple[str, str]] = []
    problems: list[str] = []
    for combo in combos:
        name = drive_filename(combo)
        md5_path = os.path.join(args.publish_dir, f"{name}.md5")
        if not os.path.exists(md5_path):
            problems.append(f"{combo}: brak {name}.md5 -- job publikacyjny nie skonczyl")
            continue
        with open(md5_path, encoding="utf-8") as fh:
            rows = parse_md5_manifest(fh.read())
        # publish_gt.sh zapisuje DWA wiersze: `.vcf.gz` i `.vcf.gz.tbi`. Jeden wiersz
        # znaczy, ze `md5sum` przerwal w polowie -- manifest z brakujacym .tbi wyglada
        # poprawnie i przechodzi `md5sum -c`, wiec nikt by tego nie zauwazyl.
        got = [n for _, n in rows]
        if got != [name, f"{name}.tbi"]:
            problems.append(f"{combo}: {name}.md5 opisuje {got}, oczekiwano "
                            f"[{name!r}, {name + '.tbi'!r}]")
            continue
        for digest, fname in rows:
            target = os.path.join(args.publish_dir, fname)
            if not os.path.exists(target):
                problems.append(f"{combo}: manifest wskazuje na nieistniejacy {fname}")
                continue
            entries.append((digest, fname))

    if problems:
        for p in problems:
            print(f"[FATAL] {p}", file=sys.stderr)
        return 1

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(format_md5_manifest(entries))
    print(f"[make_manifest] zapisano {args.out}: {len(entries)} plikow "
          f"({len(combos)} combos)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
