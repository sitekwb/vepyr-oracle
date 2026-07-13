"""Slice the exhaustive 116 mismatch TSVs. Pure: no vepyr, no pyarrow, no cluster.

`bin/triage.py` is the cluster-side CLI wrapper: it reads the real TSVs (via this
module) and additionally probes the vepyr parquet caches (via pyarrow, which this
module deliberately does NOT depend on -- see the repo-wide rule that oracle/ never
imports vepyr/pysam/polars, and Task 1's rule that oracle/triage.py additionally
never imports pyarrow, so the slicing logic stays testable on a laptop with none of
that installed).

The mismatch TSVs (`results_wgs_116/mismatches/<combo>.tsv`) share ONE fixed header:
`chrom  pos  ref  alt  allele  feature  field  vepyr_val  vep_val  category`
-- see oracle/diff.py / bin/validate.py for where they are written.
"""
from __future__ import annotations
import collections, csv
from enum import StrEnum


class FeatureKind(StrEnum):
    """Which annotation-source family a CSQ `Feature` id belongs to.

    This is the whole question behind triage cluster 2 (HGNC_ID): is a drift
    RefSeq-specific (a RefSeq-cache gap) or does it also hit Ensembl transcripts
    (a genuine vepyr bug)? `field_rollup`'s `by_feature_kind` answers it directly.
    """
    ENSEMBL = "ENSEMBL"   # ENST... / ENSR... / ENSM... (Ensembl-style stable ids)
    REFSEQ = "REFSEQ"     # XM_ / NM_ / XR_ / NR_ / NP_ (RefSeq-style accessions)
    NONE = "NONE"         # empty Feature (intergenic, or a motif with no id)
    OTHER = "OTHER"       # anything else (e.g. a bare regulatory/motif stable id)


#: RefSeq accession prefixes that can appear as a CSQ `Feature` value. Checked
#: against the first 3 characters, so `XM_`/`NM_`/`XR_`/`NR_`/`NP_` all match
#: regardless of the accession's numeric suffix or version (`.3`, `.5`, ...).
_REFSEQ_PREFIXES = ("XM_", "NM_", "XR_", "NR_", "NP_")


def classify_feature(feature: str) -> FeatureKind:
    """Bucket one CSQ `Feature` value into ENSEMBL / REFSEQ / NONE / OTHER.

    Order matters: the empty-string check MUST come first -- `"".startswith("ENS")`
    is False anyway, but checking emptiness explicitly documents that an empty
    Feature is its own category (no feature at all), not merely "unrecognised".
    """
    if not feature:
        return FeatureKind.NONE
    if feature.startswith("ENS"):
        return FeatureKind.ENSEMBL
    if feature[:3] in _REFSEQ_PREFIXES:
        return FeatureKind.REFSEQ
    return FeatureKind.OTHER


def load_mismatches(path: str) -> list[dict]:
    """Load one `mismatches/<combo>.tsv` in full, as a list of column->value dicts.

    The 116 mismatch TSVs top out at ~235k rows (hgvs_merged_am) -- small enough
    that loading a whole combo into memory with the stdlib `csv` module is fine,
    and keeps this module dependency-free (no pyarrow/polars) for Task 1's TDD.
    """
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def field_rollup(rows: list[dict]) -> dict:
    """Roll mismatch rows up per CSQ `field`.

    Returns `field -> {total, by_category, by_feature_kind, samples}`:
      - `total`: row count for this field.
      - `by_category`: counts per `oracle.drift.Category` value seen (as plain
        str keys, since the TSV already carries the category as a string).
      - `by_feature_kind`: counts per `FeatureKind` of the row's `feature` column
        -- e.g. `{"ENSEMBL": 1, "REFSEQ": 1}` tells you a drift is NOT
        RefSeq-specific; `{"REFSEQ": N}` alone tells you it is.
      - `samples`: the first 5 raw rows for this field, for eyeballing.
    """
    out: dict = {}
    for r in rows:
        f = r["field"]
        e = out.setdefault(f, {
            "total": 0,
            "by_category": collections.Counter(),
            "by_feature_kind": collections.Counter(),
            "samples": [],
        })
        e["total"] += 1
        e["by_category"][r["category"]] += 1
        e["by_feature_kind"][str(classify_feature(r["feature"]))] += 1
        if len(e["samples"]) < 5:
            e["samples"].append(r)
    for e in out.values():
        e["by_category"] = dict(e["by_category"])
        e["by_feature_kind"] = dict(e["by_feature_kind"])
    return out
