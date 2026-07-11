import csv, pathlib
from oracle.diff import diff_files

FIX = pathlib.Path(__file__).parent / "fixtures"
COMBO_NO_HGVS = {"everything": True}   # -> HGVSc mismatches are flag_expected

def test_diff_classifies_and_writes_tsv(tmp_path):
    tsv = tmp_path / "m.tsv"
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS, tsv_path=str(tsv))

    assert s["aligned_annotations"] == 2
    # SIFT: deleterious vs tolerated -> a real value_diff
    assert s["per_field"]["SIFT"]["by_category"] == {"value_diff": 1}
    # HGVSc: vepyr computed it, GT empty, combo has no hgvs -> flag_expected (x2)
    assert s["per_field"]["HGVSc"]["by_category"] == {"flag_expected": 2}
    assert s["per_field"]["Feature"]["pct"] == 100.0

    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert len(rows) == 3                                  # 1 SIFT + 2 HGVSc
    sift = [r for r in rows if r["field"] == "SIFT"][0]
    assert (sift["pos"], sift["vepyr_val"], sift["vep_val"], sift["category"]) \
        == ("100", "deleterious", "tolerated", "value_diff")

def test_chrom_filter_restricts_the_diff(tmp_path):
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"), chrom="1")
    assert s["aligned_annotations"] == 0
