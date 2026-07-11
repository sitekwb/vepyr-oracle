import csv, pathlib
import pytest
from oracle.diff import diff_files
from oracle.summary import merge

FIX = pathlib.Path(__file__).parent / "fixtures"
COMBO_NO_HGVS = {"everything": True}                 # -> HGVSc diffs are flag_expected
COMBO_HGVS = {"everything": True, "hgvs": True}      # -> HGVSc compared normally

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
    # 2 annotations x 5 shared fields = 10 comparisons; 1 SIFT + 2 HGVSc mismatch
    assert s["overall_pct"] == 70.0

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


@pytest.mark.parametrize("chrom", ["22", "chr22"])
def test_chrom_filter_KEEPS_the_requested_chromosome(tmp_path, chrom):
    """The positive counterpart to the test above -- which also passes if the filter
    drops EVERYTHING. Chrom filtering is how every WGS shard is produced: a filter
    that quietly matched nothing would make each shard a truthful-looking 0/0 and the
    whole-genome merge a sum of nothings. Both spellings must select the contig."""
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"), chrom=chrom)
    assert s["aligned_annotations"] == 2
    assert s["chrom"] == chrom


# --- Fix A: the merge-join must never silently consume an unsorted stream -------

def test_unsorted_vepyr_stream_raises(tmp_path):
    with pytest.raises(ValueError, match="not sorted"):
        diff_files(str(FIX / "vepyr_unsorted.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))


def test_unsorted_gt_stream_raises(tmp_path):
    with pytest.raises(ValueError, match="not sorted"):
        diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_unsorted.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))


def test_unsorted_error_names_the_stream_and_both_loci(tmp_path):
    with pytest.raises(ValueError) as ei:
        diff_files(str(FIX / "vepyr_unsorted.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))
    msg = str(ei.value)
    assert "vepyr_vcf" in msg          # which stream
    assert "22:100" in msg             # the offending locus
    assert "22:200" in msg             # the locus it came after
    assert "##contig" in msg           # where contig rank comes from -> diagnosable


# --- Fix B: fields present on only one side are dropped from the diff, so the ---
# --- summary must at least SAY so (else a never-implemented field reads 100%) ---

def test_summary_records_fields_present_on_only_one_side(tmp_path):
    s = diff_files(str(FIX / "vepyr_extra_field.vcf"), str(FIX / "gt_extra_field.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))

    assert s["vepyr_only_fields"] == ["VEPYR_ONLY"]
    assert s["gt_only_fields"] == ["GT_ONLY"]
    # the dropped fields really are excluded from the comparison itself
    assert s["shared_fields"] == 4
    assert "VEPYR_ONLY" not in s["per_field"]
    assert "GT_ONLY" not in s["per_field"]


def test_no_field_asymmetry_reports_empty_lists(tmp_path):
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))
    assert s["vepyr_only_fields"] == []
    assert s["gt_only_fields"] == []


# --- Fix C: a CSQ entry too short to carry Feature has NO join key. It must be --
# --- excluded from the join and counted -- never cross-paired against the other -
# --- file's equally-keyless entry (which is how false value_diffs were minted). -

def test_truncated_entries_never_cross_pair_and_are_counted(tmp_path):
    tsv = tmp_path / "m.tsv"
    # Both fixtures carry ENST_A + ENST_B with IDENTICAL values (listed in a
    # DIFFERENT order) plus one truncated, Feature-less entry each. Any value_diff
    # here is therefore false by construction.
    s = diff_files(str(FIX / "vepyr_truncated.vcf"), str(FIX / "gt_truncated.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS, tsv_path=str(tsv))

    assert s["malformed_vepyr"] == 1
    assert s["malformed_gt"] == 1
    # only the two real transcripts joined; the keyless entries did NOT pair up
    assert s["aligned_annotations"] == 2
    produced = {c for pf in s["per_field"].values() for c in pf["by_category"]}
    assert "value_diff" not in produced
    assert s["overall_pct"] == 100.0
    assert list(csv.DictReader(tsv.open(), delimiter="\t")) == []


# --- Fix D: CSQ annotations are keyed by (Allele, Feature). A multi-allelic ------
# --- record carries one entry per (allele, transcript); keying by Feature alone --
# --- made same-transcript entries overwrite each other, so half the annotations --
# --- vanished and the survivors were diffed ACROSS alleles. -----------------------

def test_multiallelic_alleles_do_not_overwrite_each_other(tmp_path):
    tsv = tmp_path / "m.tsv"
    # One multi-allelic record (C -> T,CCGC), 4 CSQ entries = 2 alleles x 2
    # transcripts, IDENTICAL values on both sides but with the ALLELES LISTED IN
    # OPPOSITE ORDER. Any mismatch here is false by construction.
    s = diff_files(str(FIX / "vepyr_multiallelic.vcf"), str(FIX / "gt_multiallelic.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS, tsv_path=str(tsv))

    assert s["aligned_annotations"] == 4          # allele x transcript pairs, not 2
    assert s["overall_pct"] == 100.0
    produced = {c for pf in s["per_field"].values() for c in pf["by_category"]}
    assert produced == set()
    assert s["malformed_vepyr"] == 0
    assert s["malformed_gt"] == 0
    assert list(csv.DictReader(tsv.open(), delimiter="\t")) == []


# --- Fix E: the OUTER record key is (REF, ALT). ALT is a comma-separated list ----
# --- whose ORDER carries no meaning, so it must not decide whether two records ---
# --- join. And the join rate must be honest about what was actually compared. ----

def test_alt_listing_order_does_not_break_the_join(tmp_path):
    tsv = tmp_path / "m.tsv"
    # Same variant, same 4 annotations; vepyr lists ALT=CCGC,T and VEP lists T,CCGC.
    s = diff_files(str(FIX / "vepyr_altorder.vcf"), str(FIX / "gt_altorder.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS, tsv_path=str(tsv))

    assert s["aligned_annotations"] == 4     # not 0 -- the record must still join
    assert s["overall_pct"] == 100.0
    assert s["only_vepyr"] == 0
    assert s["only_gt"] == 0
    assert s["join_rate"] == 1.0
    assert list(csv.DictReader(tsv.open(), delimiter="\t")) == []


def test_tsv_keeps_the_raw_alt_string_actually_present_in_the_file(tmp_path):
    """Normalisation is for the join key only; the human must see the real ALT."""
    tsv = tmp_path / "m.tsv"
    gt = tmp_path / "gt.vcf"   # perturb one allele's SIFT so a row is emitted
    gt.write_text((FIX / "gt_altorder.vcf").read_text()
                  .replace("T|missense|ENST_A|deleterious", "T|missense|ENST_A|tolerated"))
    diff_files(str(FIX / "vepyr_altorder.vcf"), str(gt),
               name="c", cache="k", combo_kwargs=COMBO_HGVS, tsv_path=str(tsv))

    rows = list(csv.DictReader(tsv.open(), delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["alt"] == "CCGC,T"        # vepyr's raw ALT, not a sorted tuple
    assert rows[0]["allele"] == "T"          # the CSQ Allele still disambiguates
    assert rows[0]["feature"] == "ENST_A"


def test_join_rate_exposes_annotations_that_never_got_compared(tmp_path):
    """vepyr emits ENST_A + ENST_B; the GT has only ENST_A.

    Comparing just ENST_A and calling that 100% is the hollowed-out-sample failure:
    an annotation nobody compared must not be able to hide inside a perfect score.
    """
    s = diff_files(str(FIX / "vepyr_partial.vcf"), str(FIX / "gt_partial.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))

    assert s["aligned_annotations"] == 1
    assert s["only_vepyr"] == 1              # ENST_B was NOT compared -- say so
    assert s["only_gt"] == 0
    assert s["join_rate"] == 0.5             # we compared half the annotations
    assert s["overall_pct"] == 100.0         # ...of the half we did compare


def test_join_rate_is_none_when_nothing_was_read(tmp_path):
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"), chrom="1")
    assert s["aligned_annotations"] == 0
    assert s["join_rate"] is None            # no denominator -> not "0%", not "100%"


def test_clean_files_report_zero_malformed(tmp_path):
    s = diff_files(str(FIX / "vepyr_mini.vcf"), str(FIX / "gt_mini.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_NO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))
    assert s["malformed_vepyr"] == 0
    assert s["malformed_gt"] == 0


# --- CRITICAL 4: an unknown contig must NOT default to a shared sentinel rank. ---
# --- Every contig absent from the ##contig header used to rank 9999, and the -----
# --- merge-join's equal-branch fires on RANK equality without ever comparing -----
# --- contig NAMES -- so two DIFFERENT unranked contigs became "the same contig". --

def test_two_unranked_contigs_raise_instead_of_being_joined_together(tmp_path):
    """vepyr has chrUn_A, the GT has chrUn_B, the GT header lists neither.

    Before: both ranked 9999, the equal-branch fired, and chrUn_A's annotations were
    diffed against chrUn_B's -- aligned_annotations=3, join_rate=1.0, and four
    FABRICATED value_diff rows in the TSV. Nothing warned.
    """
    with pytest.raises(ValueError) as ei:
        diff_files(str(FIX / "vepyr_unranked.vcf"), str(FIX / "gt_unranked.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))
    msg = str(ei.value)
    assert "Un_A" in msg          # names the offending contig
    assert "##contig" in msg      # ...and the header it must appear in


def test_contig_missing_from_the_gt_contig_header_raises(tmp_path):
    """The quiet variant: norm_chrom('chrM') == 'M' but the GT header says 'MT'.

    Before: 'M' ranked 9999, sorted last, the sort guard passed -- and the ENTIRE
    MITOCHONDRION dropped into only_vepyr/only_gt while overall_pct read 100.0.
    """
    with pytest.raises(ValueError) as ei:
        diff_files(str(FIX / "vepyr_mito.vcf"), str(FIX / "gt_mito.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS,
                   tsv_path=str(tmp_path / "m.tsv"))
    msg = str(ei.value)
    assert "'M'" in msg
    assert "##contig" in msg


# --- THE INVARIANT THE WHOLE SHARDING STRATEGY RESTS ON --------------------------
# --- Every published number is produced as merge(per-chromosome diffs). If that ---
# --- is not equal to a single-pass diff of the whole file, the sharding itself ----
# --- is the bug -- and nothing asserted it. --------------------------------------

#: per-shard keys that merge() deliberately folds into plural forms
_FOLDED = {"chrom", "mismatches_tsv", "chroms", "mismatches_tsvs"}


def test_merging_per_chromosome_shards_equals_diffing_the_whole_file(tmp_path):
    whole = diff_files(str(FIX / "vepyr_twochrom.vcf"), str(FIX / "gt_twochrom.vcf"),
                       name="c", cache="k", combo_kwargs=COMBO_HGVS,
                       tsv_path=str(tmp_path / "whole.tsv"))
    shards = [
        diff_files(str(FIX / "vepyr_twochrom.vcf"), str(FIX / "gt_twochrom.vcf"),
                   name="c", cache="k", combo_kwargs=COMBO_HGVS,
                   tsv_path=str(tmp_path / f"{c}.tsv"), chrom=c)
        for c in ("21", "22")
    ]
    merged = merge(shards, expected_chroms={"21", "22"})

    # the fixture is not trivial: real mismatches, and unjoined annotations BOTH ways
    assert whole["aligned_annotations"] == 4
    assert (whole["only_vepyr"], whole["only_gt"]) == (1, 2)
    assert whole["overall_pct"] == 90.0
    assert whole["join_rate"] == round(4 / 7, 4)

    for k in set(whole) - _FOLDED:
        assert merged[k] == whole[k], f"sharding changed {k!r}: {merged[k]} != {whole[k]}"
    assert merged["chroms"] == ["21", "22"]


def test_a_shard_set_missing_a_chromosome_cannot_masquerade_as_the_whole_file(tmp_path):
    """The other half of the invariant: a HOLE must not merge into a genome number."""
    only_21 = diff_files(str(FIX / "vepyr_twochrom.vcf"), str(FIX / "gt_twochrom.vcf"),
                         name="c", cache="k", combo_kwargs=COMBO_HGVS,
                         tsv_path=str(tmp_path / "21.tsv"), chrom="21")
    with pytest.raises(ValueError, match="22"):
        merge([only_21], expected_chroms={"21", "22"})
