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
    """-> ((chrom, pos, ref, alt), {feature: {field: value}}, n_malformed) or None.

    The Feature value is the join key across the two files. An entry too short to
    carry that column (or a CSQ header with no Feature column at all) therefore has
    NO join key: it is excluded from `feats` and counted as malformed instead.

    It deliberately gets no synthesised key. The previous positional fallback
    (`str(len(feats))`) minted "0", "1", ... independently in each file, so two
    unrelated keyless entries collided on the same key and were diffed against each
    other -- reporting field-level `value_diff`s between two different transcripts.
    Counting the entry is honest; guessing its identity is worse than dropping it.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return None
    m = re.search(r"(?:^|;)CSQ=([^;\t]+)", cols[7])
    if not m:
        return None
    feats: dict[str, dict[str, str]] = {}
    malformed = 0
    for entry in m.group(1).split(","):
        parts = entry.split("|")
        if feat_i is None or feat_i >= len(parts):
            malformed += 1
            continue
        feats[parts[feat_i]] = {f: (parts[i] if i < len(parts) else "") for i, f in enumerate(fields)}
    return (norm_chrom(cols[0]), int(cols[1]), cols[3], cols[4]), feats, malformed
