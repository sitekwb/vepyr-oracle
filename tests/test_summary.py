import pytest

from oracle.summary import new_accumulator, record_match, record_mismatch, finalize, merge
from oracle.drift import Category

def _acc(field="SIFT"):
    a = new_accumulator([field])
    return a

def test_finalize_computes_per_field_pct_and_by_category():
    a = _acc()
    for _ in range(9):
        record_match(a, "SIFT")
    record_mismatch(a, "SIFT", Category.VEP_ONLY, ["22:100 A>G ENST1", "", "deleterious"])
    s = finalize(a, name="c", cache="115_GRCh38_merged")
    pf = s["per_field"]["SIFT"]
    assert (pf["total"], pf["match"], pf["pct"]) == (10, 9, 90.0)
    assert pf["by_category"] == {"vep_only": 1}
    assert s["overall_pct"] == 90.0
    assert s["samples"]["SIFT"][0][2] == "deleterious"

def test_merge_is_additive_across_shards():
    a, b = _acc(), _acc()
    for _ in range(5):
        record_match(a, "SIFT")
    record_mismatch(a, "SIFT", Category.VEP_ONLY, ["chr1", "", "x"])
    for _ in range(3):
        record_match(b, "SIFT")
    record_mismatch(b, "SIFT", Category.VALUE_DIFF, ["chr2", "p", "q"])

    m = merge([finalize(a, "c", "k"), finalize(b, "c", "k")])
    pf = m["per_field"]["SIFT"]
    assert (pf["total"], pf["match"]) == (10, 8)          # 6+4 total, 5+3 match
    assert pf["by_category"] == {"vep_only": 1, "value_diff": 1}
    assert m["overall_pct"] == 80.0
    assert len(m["samples"]["SIFT"]) == 2                  # samples concatenated

def test_merge_rejects_shards_with_mismatched_field_sets():
    a = new_accumulator(["SIFT"])
    record_match(a, "SIFT")
    b = new_accumulator(["SIFT", "PolyPhen"])
    record_match(b, "SIFT")
    record_match(b, "PolyPhen")
    with pytest.raises(ValueError, match="field set"):
        merge([finalize(a, "c", "k"), finalize(b, "c", "k")])

def test_merge_rejects_shards_from_different_combos():
    a = new_accumulator(["SIFT"]); record_match(a, "SIFT")
    b = new_accumulator(["SIFT"]); record_match(b, "SIFT")
    with pytest.raises(ValueError, match="combo"):
        merge([finalize(a, "combo_x", "k"), finalize(b, "combo_y", "k")])

def test_merge_rejects_shards_from_different_caches():
    a = new_accumulator(["SIFT"]); record_match(a, "SIFT")
    b = new_accumulator(["SIFT"]); record_match(b, "SIFT")
    with pytest.raises(ValueError, match="cache"):
        merge([finalize(a, "hgvs_merged", "115_GRCh38_merged"),
               finalize(b, "hgvs_merged", "116_GRCh38_merged")])
