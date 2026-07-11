import pytest

from oracle.summary import new_accumulator, record_match, record_mismatch, finalize, merge
from oracle.drift import Category


def _acc(field="SIFT"):
    return new_accumulator([field])


def _finalize(acc, chrom, *, name="c", cache="k", only_vepyr=0, only_gt=0,
              malformed_vepyr=0, malformed_gt=0, shared_fields=1,
              vepyr_only_fields=(), gt_only_fields=(), status="ok"):
    """finalize() with EXACTLY the extras diff_files() attaches -- a real shard."""
    seen = acc["aligned"] + only_vepyr + only_gt
    s = finalize(
        acc, name=name, cache=cache,
        only_vepyr=only_vepyr, only_gt=only_gt,
        join_rate=round(acc["aligned"] / seen, 4) if seen else None,
        shared_fields=shared_fields,
        vepyr_only_fields=list(vepyr_only_fields), gt_only_fields=list(gt_only_fields),
        malformed_vepyr=malformed_vepyr, malformed_gt=malformed_gt,
        mismatches_tsv=f"/out/{chrom}.tsv", chrom=chrom,
    )
    s["status"] = status
    return s


def _shard(chrom, *, match=1, mismatch=0, field="SIFT", **over):
    a = new_accumulator([field])
    for _ in range(match):
        record_match(a, field)
    for i in range(mismatch):
        record_mismatch(a, field, Category.VALUE_DIFF, [f"{chrom}:{i}", "p", "q"])
    a["aligned"] = match + mismatch
    return _finalize(a, chrom, **over)


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

    m = merge([_finalize(a, "21"), _finalize(b, "22")])
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
        merge([_finalize(a, "21"), _finalize(b, "22", shared_fields=2)])


def test_merge_rejects_shards_from_different_combos():
    with pytest.raises(ValueError, match="combo"):
        merge([_shard("21", name="combo_x"), _shard("22", name="combo_y")])


def test_merge_rejects_shards_from_different_caches():
    with pytest.raises(ValueError, match="cache"):
        merge([_shard("21", name="hgvs_merged", cache="115_GRCh38_merged"),
               _shard("22", name="hgvs_merged", cache="116_GRCh38_merged")])


# --- CRITICAL 2: merge() must not silently DROP the coverage counters. ----------
# --- The whole-genome summary is the one the PDF reads; if the unjoined-variant --
# --- evidence evaporates at merge, every safety signal we built is decorative. ---

#: keys that are meaningful per-shard only -- merge folds them into plural forms
SHARD_ONLY_KEYS = {"chrom", "mismatches_tsv"}


def test_merge_propagates_every_key_a_shard_carries():
    m = merge([_shard("21"), _shard("22")])
    shard = _shard("21")
    missing = (set(shard) - SHARD_ONLY_KEYS) - set(m)
    assert missing == set(), f"merge() DROPPED shard keys: {sorted(missing)}"
    # ...and the folded forms replace what it dropped
    assert m["chroms"] == ["21", "22"]
    assert m["mismatches_tsvs"] == ["/out/21.tsv", "/out/22.tsv"]


def test_merge_sums_the_coverage_counters():
    m = merge([
        _shard("21", match=3, only_vepyr=2, only_gt=1, malformed_vepyr=4, malformed_gt=5),
        _shard("22", match=7, only_vepyr=8, only_gt=9, malformed_vepyr=1, malformed_gt=2),
    ])
    assert m["aligned_annotations"] == 10
    assert m["only_vepyr"] == 10
    assert m["only_gt"] == 10
    assert m["malformed_vepyr"] == 5
    assert m["malformed_gt"] == 7


def test_merge_recomputes_join_rate_from_the_summed_counters():
    """Averaging per-shard rates is NOT associative -- it must be recomputed."""
    # shard A: 1 aligned / 1 seen = 1.0 ; shard B: 1 aligned / 101 seen = 0.0099
    # naive mean = 0.505 ; truthful = 2 / 102 = 0.0196
    m = merge([_shard("21", match=1), _shard("22", match=1, only_gt=100)])
    assert m["join_rate"] == round(2 / 102, 4)


def test_merge_join_rate_is_none_when_nothing_was_seen():
    m = merge([_shard("21", match=0), _shard("22", match=0)])
    assert m["join_rate"] is None


def test_merge_carries_the_field_asymmetry_lists():
    a = _shard("21", vepyr_only_fields=["VEPYR_ONLY"], gt_only_fields=["GT_ONLY"])
    b = _shard("22", vepyr_only_fields=["VEPYR_ONLY"], gt_only_fields=["GT_ONLY"])
    m = merge([a, b])
    assert m["shared_fields"] == 1
    assert m["vepyr_only_fields"] == ["VEPYR_ONLY"]
    assert m["gt_only_fields"] == ["GT_ONLY"]


def test_merge_rejects_shards_that_disagree_on_the_field_asymmetry():
    with pytest.raises(ValueError, match="vepyr_only_fields"):
        merge([_shard("21"), _shard("22", vepyr_only_fields=["VEPYR_ONLY"])])


def test_merge_refuses_a_shard_key_it_does_not_know_how_to_combine():
    """The regression guard: add a key to finalize() and forget merge() -> boom.

    This is the ONLY thing standing between a future counter and silent deletion
    at exactly the step that produces the published number.
    """
    a, b = _shard("21"), _shard("22")
    a["brand_new_counter"] = 1
    b["brand_new_counter"] = 2
    with pytest.raises(ValueError, match="brand_new_counter"):
        merge([a, b])


def test_merge_refuses_a_half_schema_shard():
    """A stale JSON from an earlier schema would contribute counters we cannot see."""
    a, b = _shard("21"), _shard("22")
    del b["only_gt"]
    with pytest.raises(ValueError, match="only_gt"):
        merge([a, b])
