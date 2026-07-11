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


def parse_record(line: str, fields: list[str], feat_i: int | None,
                 allele_i: int | None = None):
    """-> ((chrom, pos, ref, alt), {(allele, feature): {field: value}}, n_malformed).

    None if the line carries no CSQ.

    The join key is **(Allele, Feature)**, not Feature alone. A multi-allelic record
    (`REF=C ALT=T,CCGC`) packs one CSQ entry per *allele x transcript* pair into a
    single VCF record, all sharing the same Feature. Keyed by Feature alone, the
    later allele's entry overwrote the earlier one: half the annotations vanished,
    and whichever survived was diffed against the other file's arbitrary survivor --
    potentially a different allele entirely (an insertion's annotation compared
    against an SNV's). `Allele` is the first CSQ field precisely to disambiguate this.

    An entry that cannot produce a well-formed key -- missing EITHER component, or a
    CSQ header lacking the column outright -- has no identity, so it is excluded from
    the join and counted as malformed. It gets no synthesised key: the old positional
    fallback (`str(len(feats))`) minted "0", "1", ... independently in each file, so
    unrelated keyless entries collided and were diffed against each other. A repeated
    key is ambiguous for the same reason: keep the first, count the surplus. Counting
    is honest; guessing an entry's identity is worse than declining to match it.
    """
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return None
    m = re.search(r"(?:^|;)CSQ=([^;\t]+)", cols[7])
    if not m:
        return None
    feats: dict[tuple[str, str], dict[str, str]] = {}
    malformed = 0
    for entry in m.group(1).split(","):
        parts = entry.split("|")
        if feat_i is None or feat_i >= len(parts) or allele_i is None or allele_i >= len(parts):
            malformed += 1
            continue
        key = (parts[allele_i], parts[feat_i])
        if key in feats:          # ambiguous -- never overwrite, never silently drop
            malformed += 1
            continue
        feats[key] = {f: (parts[i] if i < len(parts) else "") for i, f in enumerate(fields)}
    return (norm_chrom(cols[0]), int(cols[1]), cols[3], cols[4]), feats, malformed
