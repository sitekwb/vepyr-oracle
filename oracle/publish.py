"""Nazewnictwo i manifest dla publikacji ground truth na Google Drive.

Folder 115.2 Marka uzywa jednej konwencji (`HG002_annotated_wgs_everything_*`),
a shardy na Obelixie drugiej (`hgvs_merged.vcf`). To mapowanie jest JAWNA tabela,
nie regexem: plik z bledna nazwa jest gorszy niz brak pliku, bo Marek nie ma jak
tego wykryc inaczej niz czytajac naglowek ##VEP-command-line= kazdego z osmiu.
"""
from __future__ import annotations

import re
from typing import Final

_PREFIX: Final[str] = "HG002_annotated_wgs_everything"

#: combo (nazwa na Obelixie) -> sufiks w konwencji folderu 115.2.
#: Pusty sufiks = sam prefiks (combo `everything`).
_SUFFIXES: Final[dict[str, str]] = {
    "everything":                   "",
    "everything_hgvs":              "_hgvs",
    "hgvs_merged":                  "_hgvs_merged",
    "hgvs_merged_am":               "_hgvs_merged_am",
    "hgvs_merged_pick":             "_hgvs_merged_pick",
    "hgvs_merged_pick_allele":      "_hgvs_merged_pick_allele",
    "hgvs_merged_pick_allele_gene": "_hgvs_merged_pick_allele_gene",
    "hgvs_merged_per_gene":         "_hgvs_merged_per_gene",
    "hgvs_merged_flag_pick_allele": "_hgvs_merged_flag_pick_allele",
    "hgvs_refseq":                  "_hgvs_refseq",
}

_MD5_RE: Final[re.Pattern[str]] = re.compile(r"^([0-9a-f]{32}) [ *](.+)$")


def drive_filename(combo: str, *, compressed: bool = True) -> str:
    """Nazwa pliku dla `combo` w konwencji folderu 115.2.

    Args:
        combo: nazwa combo tak, jak wystepuje w matrix.tsv i w nazwach shardow.
        compressed: czy dokleic `.gz` (publikujemy bgzip; `False` dla surowych).

    Returns:
        Nazwe pliku, np. `HG002_annotated_wgs_everything_hgvs_merged.vcf.gz`.

    Raises:
        ValueError: gdy `combo` nie jest znane. Fallback nadalby plikowi
            prawdopodobnie wygladajaca, ale bledna nazwe.
    """
    try:
        suffix = _SUFFIXES[combo]
    except KeyError:
        raise ValueError(
            f"unknown combo {combo!r}; znane: {sorted(_SUFFIXES)}"
        ) from None
    return f"{_PREFIX}{suffix}.vcf{'.gz' if compressed else ''}"


def format_md5_manifest(entries: list[tuple[str, str]]) -> str:
    """Zapisz `(md5, nazwa)` w formacie czytelnym dla `md5sum -c`.

    Separator to DWIE spacje (tryb tekstowy). Jedna spacja powoduje
    `no properly formatted checksum lines found` po stronie odbiorcy.
    """
    return "".join(f"{digest}  {name}\n" for digest, name in entries)


def parse_md5_manifest(text: str) -> list[tuple[str, str]]:
    """Odwrotnosc `format_md5_manifest`; toleruje tez format binarny (` *nazwa`).

    Raises:
        ValueError: gdy ktorakolwiek niepusta linia nie jest poprawna linia md5sum.
    """
    out: list[tuple[str, str]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        m = _MD5_RE.match(line)
        if not m:
            raise ValueError(f"malformed md5 line {lineno}: {line!r}")
        out.append((m.group(1), m.group(2)))
    return out
