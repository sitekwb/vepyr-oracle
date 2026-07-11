"""oracle.report -- the pure logic behind the PDF.

The PDF is the deliverable: it is where a wrong number gets PUBLISHED, and every
adversarial finding in this codebase ultimately cashes out here. The single most
dangerous shape a report can take is the one measured on real chr22 data:

    overall_pct = 100.0   join_rate = 0.843

"100.0% concordant, zero mismatches" -- while 16% of annotations were never
compared at all. `overall_pct` speaks ONLY for the annotations that joined; on a
`--pick` combo (one annotation per variant) a vepyr that picks the WRONG transcript
contributes nothing to the numerator OR the denominator, so a `--pick` implementation
that is wrong about everything scores a perfect 100%.

`coverage()` is the guard: it exists to make that shape impossible to publish
without a warning printed next to it, and these tests are what hold it to that.
"""
from __future__ import annotations

import pytest

from oracle.drift import Category
from oracle.report import (coverage, drift_rows, headline_risk, rollup_categories,
                           true_pct)


def _pf(total: int, match: int, **by_category: int) -> dict:
    return {"total": total, "match": match,
            "pct": round(100 * match / total, 3) if total else None,
            "by_category": dict(by_category)}


def _summary(**over) -> dict:
    """A whole-genome summary with the exact key set oracle.summary.merge() emits."""
    s = {
        "name": "hgvs_merged", "cache": "115_GRCh38_merged", "status": "ok",
        "aligned_annotations": 1000,
        "overall_pct": 100.0,
        "join_rate": 1.0,
        "only_vepyr": 0, "only_gt": 0,
        "malformed_vepyr": 0, "malformed_gt": 0,
        "duplicate_records_vepyr": 0, "duplicate_records_gt": 0,
        "per_field": {"Consequence": _pf(1000, 1000)},
        "samples": {},
        "shared_fields": 1, "vepyr_only_fields": [], "gt_only_fields": [],
        "chroms": [str(i) for i in range(1, 23)],
        "mismatches_tsvs": [],
    }
    s.update(over)
    return s


# =============================================================================
# true_pct -- flag-semantics expected-diffs count as AGREEMENT, nothing else does
# =============================================================================

def test_true_pct_is_overall_pct_when_there_are_no_flag_expected_diffs():
    s = _summary(per_field={"Consequence": _pf(100, 90, value_diff=10)})
    assert true_pct(s) == 90.0


def test_true_pct_counts_flag_expected_mismatches_as_agreement():
    """HGVSc on a non-hgvs combo: vepyr computes it, the GT was never asked for it."""
    s = _summary(per_field={
        "Consequence": _pf(100, 100),
        "HGVSc": _pf(100, 0, flag_expected=100),
    })
    assert s["per_field"]["HGVSc"]["pct"] == 0.0     # raw: looks catastrophic
    assert true_pct(s) == 100.0                      # true: not a bug, not a diff


def test_true_pct_does_not_launder_a_real_gap_hiding_in_a_flag_expected_field():
    """A field can be PART flag-expected and part real: only the flag-expected half
    is forgiven. A vepyr regression that stopped emitting HGVSc lands in vep_only,
    and must still drag the true % down."""
    s = _summary(per_field={
        "HGVSc": _pf(100, 0, flag_expected=90, vep_only=10),
    })
    assert true_pct(s) == 90.0


def test_true_pct_ignores_fields_nothing_was_compared_on():
    s = _summary(per_field={"Consequence": _pf(100, 90, value_diff=10),
                            "NEVER_SEEN": _pf(0, 0)})
    assert true_pct(s) == 90.0


def test_true_pct_is_none_when_nothing_was_compared_at_all():
    assert true_pct(_summary(per_field={})) is None
    assert true_pct({"name": "x", "cache": "y", "status": "no_gt"}) is None


# =============================================================================
# coverage -- THE honesty panel
# =============================================================================

def test_coverage_of_a_clean_combo_is_trustworthy_and_silent():
    cov = coverage(_summary())
    assert cov["trustworthy"] is True
    assert cov["warnings"] == []
    assert cov["join_rate"] == 1.0
    assert cov["aligned"] == 1000


def test_coverage_reports_every_counter_the_summary_carries():
    cov = coverage(_summary(aligned_annotations=10, only_vepyr=1, only_gt=2,
                            malformed_vepyr=3, malformed_gt=4,
                            duplicate_records_vepyr=5, duplicate_records_gt=6,
                            join_rate=round(10 / 13, 4)))
    assert cov["aligned"] == 10
    assert cov["only_vepyr"] == 1
    assert cov["only_gt"] == 2
    assert cov["malformed_vepyr"] == 3
    assert cov["malformed_gt"] == 4
    assert cov["duplicate_records_vepyr"] == 5
    assert cov["duplicate_records_gt"] == 6
    assert cov["not_compared"] == 3            # only_vepyr + only_gt


# --- THE regression test: the real chr22 --pick_allele shape -------------------

def test_the_real_chr22_pick_allele_shape_is_NOT_trustworthy():
    """MEASURED on real chr22 data: overall_pct=100.0 with join_rate=0.843.

    A perfect score on the 84.3% of annotations that joined, and total silence
    about the 15.7% that never got compared. If coverage() ever calls this
    trustworthy, the PDF publishes "100.0% concordant" for a --pick implementation
    that could be wrong about every single variant it picked differently.
    """
    s = _summary(name="hgvs_merged_pick_allele", overall_pct=100.0, join_rate=0.843,
                 aligned_annotations=84_300, only_vepyr=7_850, only_gt=7_850,
                 per_field={"Consequence": _pf(84_300, 84_300)})
    cov = coverage(s)

    assert cov["trustworthy"] is False
    joined = " ".join(cov["warnings"])
    assert "15.7%" in joined, f"the un-compared fraction must be NAMED: {cov['warnings']}"
    assert "never compared" in joined
    # ...and it must say WHY that makes the headline number untrustworthy, not just
    # print a statistic next to it.
    assert "conditional on the two tools already agreeing" in joined


def test_a_join_rate_below_the_floor_names_the_uncompared_percentage():
    cov = coverage(_summary(join_rate=0.9, aligned_annotations=90, only_gt=10))
    assert cov["trustworthy"] is False
    assert any("10.0%" in w and "never compared" in w for w in cov["warnings"])


def test_a_join_rate_at_the_floor_is_still_trustworthy():
    """0.99 is the floor, not a value below it -- the boundary must not drift."""
    cov = coverage(_summary(join_rate=0.99, aligned_annotations=99))
    assert not any("never compared" in w for w in cov["warnings"])


def test_a_null_join_rate_means_nothing_was_read_and_is_never_trustworthy():
    cov = coverage(_summary(join_rate=None, aligned_annotations=0, overall_pct=None,
                            per_field={}))
    assert cov["trustworthy"] is False
    assert any("no annotations" in w.lower() for w in cov["warnings"])


def test_only_gt_annotations_warn_as_potential_vepyr_gaps():
    """VEP emitted them, vepyr did not: the vep_only direction is the REAL-GAP
    direction, and these never even reached the classifier."""
    cov = coverage(_summary(only_gt=42))
    assert cov["trustworthy"] is False
    w = " ".join(cov["warnings"])
    assert "42" in w
    assert "potential vepyr gaps" in w
    assert "not compared" in w


def test_only_vepyr_annotations_warn_too():
    cov = coverage(_summary(only_vepyr=7))
    assert cov["trustworthy"] is False
    assert any("7" in w and "vepyr" in w for w in cov["warnings"])


def test_malformed_entries_warn_on_both_sides():
    cov = coverage(_summary(malformed_vepyr=3, malformed_gt=5))
    assert cov["trustworthy"] is False
    w = " ".join(cov["warnings"])
    assert "3" in w and "5" in w
    assert w.count("could not be keyed") == 2      # one warning per side


def test_duplicate_records_warn_on_both_sides():
    cov = coverage(_summary(duplicate_records_vepyr=2, duplicate_records_gt=9))
    assert cov["trustworthy"] is False
    w = " ".join(cov["warnings"])
    assert "2" in w and "9" in w
    assert "duplicate" in w.lower()


def test_a_field_only_one_side_emits_warns_that_a_WHOLE_COLUMN_went_uncompared():
    """A CSQ column vepyr never implements is dropped from the diff entirely --
    leaving a report that is 100% green on the fields it does emit."""
    cov = coverage(_summary(gt_only_fields=["SpliceAI_pred_DS_AG", "am_pathogenicity"],
                            vepyr_only_fields=["VEPYR_INTERNAL"]))
    assert cov["trustworthy"] is False
    w = " ".join(cov["warnings"])
    assert "SpliceAI_pred_DS_AG" in w and "am_pathogenicity" in w
    assert "VEPYR_INTERNAL" in w
    assert w.lower().count("not compared at all") == 2


def test_a_non_ok_combo_is_never_trustworthy():
    cov = coverage({"name": "x", "cache": "y", "status": "no_gt"})
    assert cov["trustworthy"] is False
    assert any("no_gt" in w for w in cov["warnings"])
    assert cov["join_rate"] is None


def test_coverage_warnings_are_plain_renderable_strings():
    """The PDF prints these verbatim -- no objects, no nesting, no empties."""
    cov = coverage(_summary(join_rate=0.5, only_gt=1, malformed_gt=1,
                            duplicate_records_vepyr=1, gt_only_fields=["X"]))
    assert len(cov["warnings"]) >= 5
    assert all(isinstance(w, str) and w.strip() for w in cov["warnings"])


# =============================================================================
# headline_risk -- the ONE line page 1 prints next to an untrustworthy combo
# =============================================================================
#
# Page 1 is the page everyone reads and most people read ONLY. If its one-line
# summary of the risk can itself understate the risk, the trust page behind it is
# decoration. The bug this function exists to fix: page 1 used to render the risk as
# "N% of annotations were never compared", which is a statement about the JOIN RATE
# alone -- so a combo that joined perfectly but never compared two whole CSQ COLUMNS
# was announced on page 1 as "0.0% of annotations NEVER COMPARED (0 annotations)".
# Technically true, and the exact opposite of the warning it was standing in for.

def test_headline_risk_is_none_for_a_trustworthy_combo():
    assert headline_risk(_summary()) is None


def test_headline_risk_names_the_uncompared_fraction_when_the_join_is_the_problem():
    risk = headline_risk(_summary(overall_pct=100.0, join_rate=0.843,
                                  aligned_annotations=84_300,
                                  only_vepyr=7_850, only_gt=7_850))
    assert "15.7%" in risk
    assert "NEVER COMPARED" in risk
    assert "15,700" in risk            # the absolute count, thousands-separated


def test_headline_risk_does_not_report_a_column_gap_as_zero_percent_uncompared():
    """A perfect join rate does NOT mean everything was compared: two CSQ columns VEP
    emits and vepyr does not were dropped from the diff wholesale. Page 1 must say
    THAT, not '0.0% of annotations never compared'."""
    risk = headline_risk(_summary(join_rate=1.0, gt_only_fields=["SpliceAI", "am_path"]))
    assert "0.0%" not in risk
    assert "CSQ column" in risk


def test_headline_risk_speaks_the_leading_problem_when_the_join_is_clean():
    risk = headline_risk(_summary(join_rate=1.0, malformed_gt=17))
    assert "17" in risk and "keyed" in risk


def test_headline_risk_says_how_many_further_warnings_are_waiting_on_page_2():
    risk = headline_risk(_summary(join_rate=0.9, aligned_annotations=90, only_gt=10,
                                  malformed_gt=3, duplicate_records_vepyr=2,
                                  gt_only_fields=["X"]))
    assert "10.0%" in risk
    assert "+4 more" in risk          # only_gt, malformed_gt, duplicates, column


def test_headline_risk_of_an_undiffed_combo_says_it_was_never_diffed():
    risk = headline_risk({"name": "x", "cache": "y", "status": "no_gt"})
    assert "no_gt" in risk


def test_headline_risk_when_nothing_was_read_at_all():
    risk = headline_risk(_summary(join_rate=None, aligned_annotations=0, per_field={}))
    assert "no annotations" in risk.lower() or "nothing" in risk.lower()


# =============================================================================
# rollup_categories -- vep_only FIRST (real gaps), flag_expected LAST
# =============================================================================

def test_rollup_totals_every_category_across_combos_and_fields():
    a = _summary(per_field={"SIFT": _pf(10, 8, vep_only=1, value_diff=1),
                            "HGVSc": _pf(10, 5, flag_expected=5)})
    b = _summary(per_field={"SIFT": _pf(10, 9, vep_only=1)})
    roll = rollup_categories([a, b])
    assert roll["vep_only"] == 2
    assert roll["value_diff"] == 1
    assert roll["flag_expected"] == 5


def test_rollup_puts_the_real_gaps_first_and_the_non_bugs_last():
    roll = rollup_categories([_summary()])
    keys = list(roll)
    assert keys[0] == Category.VEP_ONLY
    assert keys[-1] == Category.FLAG_EXPECTED


def test_rollup_lists_every_category_even_at_zero():
    """A category that dropped to zero must read as an explicit 0, not vanish --
    an absent row is indistinguishable from a row nobody computed."""
    roll = rollup_categories([_summary()])
    assert set(roll) == {str(c) for c in Category}
    assert set(roll.values()) == {0}


def test_rollup_skips_combos_that_were_never_diffed():
    roll = rollup_categories([_summary(per_field={"SIFT": _pf(10, 9, vep_only=1)}),
                              {"name": "x", "cache": "y", "status": "no_gt"}])
    assert roll["vep_only"] == 1


# =============================================================================
# drift_rows -- one row per drifting field, worst first
# =============================================================================

def test_drift_rows_skips_fields_with_no_mismatches():
    rows = drift_rows(_summary(per_field={"Consequence": _pf(100, 100),
                                          "SIFT": _pf(100, 99, value_diff=1)}))
    assert [r["field"] for r in rows] == ["SIFT"]


def test_drift_rows_carries_the_counts_pct_and_dominant_cause():
    rows = drift_rows(_summary(per_field={
        "SIFT": _pf(100, 70, vep_only=25, value_diff=5),
    }))
    (row,) = rows
    assert row["mismatches"] == 30
    assert row["pct"] == 70.0
    assert row["dominant"] == Category.VEP_ONLY
    assert row["by_category"] == {"vep_only": 25, "value_diff": 5}
    assert row["flag_expected"] is False


def test_drift_rows_are_ordered_worst_first():
    rows = drift_rows(_summary(per_field={
        "A": _pf(100, 98, value_diff=2),
        "B": _pf(100, 50, value_diff=50),
        "C": _pf(100, 90, value_diff=10),
    }))
    assert [r["field"] for r in rows] == ["B", "C", "A"]


def test_a_purely_flag_expected_field_is_TAGGED_and_sorts_below_real_drift():
    """HGVSc's 100k expected-diffs must not crowd a 5-mismatch REAL gap off the top
    of the table -- the whole point of the table is to surface real drift."""
    rows = drift_rows(_summary(per_field={
        "HGVSc": _pf(100_000, 0, flag_expected=100_000),
        "SIFT": _pf(100, 95, vep_only=5),
    }))
    assert [r["field"] for r in rows] == ["SIFT", "HGVSc"]
    assert rows[0]["flag_expected"] is False
    assert rows[1]["flag_expected"] is True
    assert rows[1]["real_mismatches"] == 0


def test_a_partly_flag_expected_field_is_NOT_tagged_as_a_non_bug():
    """90 expected + 10 real vep_only: tagging this "flag-expected" would hide a
    real vepyr gap behind a "not a bug" label."""
    rows = drift_rows(_summary(per_field={
        "HGVSc": _pf(100, 0, flag_expected=90, vep_only=10)}))
    (row,) = rows
    assert row["flag_expected"] is False
    assert row["real_mismatches"] == 10
    assert row["mismatches"] == 100


def test_drift_rows_of_a_never_diffed_combo_is_empty():
    assert drift_rows({"name": "x", "cache": "y", "status": "no_gt"}) == []
