"""CSQ header + record parsing. Pure: no vepyr, no pysam."""
from __future__ import annotations
import gzip, re


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def norm_chrom(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def csq_format_fields(path: str) -> list[str]:
    """The `Format: a|b|c` field list from the ##INFO=<ID=CSQ ...> header."""
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format:\s*([^\">]+)", line)
                return m.group(1).split("|") if m else []
            if line.startswith("#CHROM"):
                break
    return []


def chrom_ranks(path: str) -> dict[str, int]:
    """Contig order from ##contig headers -> used to merge-join two sorted VCFs."""
    ranks: dict[str, int] = {}
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith("##contig"):
                m = re.search(r"ID=([^,>]+)", line)
                if m:
                    ranks[norm_chrom(m.group(1))] = len(ranks)
            elif line.startswith("#CHROM"):
                break
    return ranks


def parse_record(line: str, fields: list[str], feat_i: int | None):
    """-> ((chrom, pos, ref, alt), {feature: {field: value}}) or None."""
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return None
    m = re.search(r"(?:^|;)CSQ=([^;\t]+)", cols[7])
    if not m:
        return None
    feats: dict[str, dict[str, str]] = {}
    for entry in m.group(1).split(","):
        parts = entry.split("|")
        feat = parts[feat_i] if feat_i is not None and feat_i < len(parts) else str(len(feats))
        feats[feat] = {f: (parts[i] if i < len(parts) else "") for i, f in enumerate(fields)}
    return (norm_chrom(cols[0]), int(cols[1]), cols[3], cols[4]), feats
