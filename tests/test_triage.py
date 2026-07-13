"""Tests for oracle.triage -- slicing the exhaustive 116 mismatch TSVs.

Pure: no vepyr, no pyarrow, no cluster. See oracle/triage.py's module docstring.
"""
from __future__ import annotations

from oracle.triage import load_mismatches, field_rollup, classify_feature, FeatureKind

TSV = """chrom\tpos\tref\talt\tallele\tfeature\tfield\tvepyr_val\tvep_val\tcategory
1\t100\tA\tG\tG\tENST1\tHGNC_ID\tHGNC:5\t\tvepyr_only
1\t200\tC\tT\tT\tXM_1\tHGNC_ID\tHGNC:6\t\tvepyr_only
1\t300\tA\tG\tG\tENST2\tam_class\t\tlikely_benign\tvep_only
1\t400\tA\tG\tG\t\tConsequence\tTF_binding_site_variant\tintergenic_variant\tvalue_diff
"""


def test_load_and_rollup(tmp_path):
    p = tmp_path / "m.tsv"; p.write_text(TSV)
    rows = load_mismatches(str(p))
    assert len(rows) == 4
    r = field_rollup(rows)
    assert r["HGNC_ID"]["total"] == 2
    assert r["HGNC_ID"]["by_category"] == {"vepyr_only": 2}
    assert r["am_class"]["by_category"] == {"vep_only": 1}


def test_feature_kind_distinguishes_refseq_from_ensembl():
    assert classify_feature("ENST00000155674") is FeatureKind.ENSEMBL
    assert classify_feature("XM_011540537.3") is FeatureKind.REFSEQ
    assert classify_feature("NM_145172.5")    is FeatureKind.REFSEQ
    assert classify_feature("")               is FeatureKind.NONE


def test_rollup_splits_by_feature_kind(tmp_path):
    p = tmp_path / "m.tsv"; p.write_text(TSV)
    r = field_rollup(load_mismatches(str(p)))
    # Is the HGNC_ID drift RefSeq-specific? That is the whole question for cluster 2.
    assert r["HGNC_ID"]["by_feature_kind"] == {"ENSEMBL": 1, "REFSEQ": 1}
    assert r["Consequence"]["by_feature_kind"] == {"NONE": 1}


def test_load_mismatches_tolerates_crlf(tmp_path):
    """The TSVs on the cluster were written with csv.writer's default CRLF terminator.
    A trailing \r on the last column made `awk '$10=="vep_only"'` match ZERO rows --
    a whole cluster silently looked empty. The reader must not inherit that trap."""
    p = tmp_path / "crlf.tsv"
    p.write_bytes(TSV.replace("\n", "\r\n").encode())
    rows = load_mismatches(str(p))
    assert rows[2]["category"] == "vep_only"          # NOT "vep_only\r"
    assert field_rollup(rows)["am_class"]["by_category"] == {"vep_only": 1}
