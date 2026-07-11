"""Pure report logic: what the PDF says, decided here and tested here.

`bin/make_report.py` is a rendering shell over this module -- it draws boxes and
picks fonts, and it decides NOTHING. Every number, every ordering and every warning
the report publishes is computed by a function below, because the PDF is the
artefact a scientist reads and trusts, and a figure is not made true by being
rendered in a nice font.

THE FAILURE MODE THIS MODULE EXISTS TO PREVENT
---------------------------------------------
`overall_pct` speaks ONLY for the annotations that were actually COMPARED. The diff
is a join: an (allele x transcript) annotation contributes to the concordance figure
only if BOTH tools emitted it. So an annotation vepyr never produced, or produced for
a different transcript than VEP chose, does not count as a mismatch -- it does not
count at all. On a `--pick` combo (one annotation per variant) this is lethal: if
vepyr picks the WRONG transcript, that variant vanishes from both the numerator and
the denominator, and a `--pick` implementation that is wrong about 99% of variants
prints

    100.0% concordant, zero mismatches

Measured on real chr22 data: `overall_pct=100.0` at `join_rate=0.843` -- a perfect
score while 15.7% of annotations were never looked at.

`join_rate` (aligned / (aligned + only_vepyr + only_gt)) is the signal that makes the
headline number readable, and `coverage()` below turns it into warnings the report
prints in full. A concordance % without its join rate is untrustworthy NO MATTER HOW
CORRECT THE COUNTERS ARE, which is why `coverage()` is the heart of this module and
why the two must never be rendered apart.
"""
from __future__ import annotations

from typing import Any

from .drift import Category

#: Below this join rate the headline concordance figure stops being a statement about
#: vepyr and starts being a statement about the subset of annotations that happened to
#: agree well enough to be joined at all. 0.99 is deliberately strict: the join key is
#: (chrom, pos, ref, alt) x (allele, transcript), which the two tools should agree on
#: essentially always -- anything more than a 1% shortfall is a finding, not noise.
JOIN_RATE_FLOOR = 0.99

#: Report order for drift categories: REAL GAPS FIRST, non-bugs last.
#:
#: `vep_only` (VEP populated the field, vepyr left it empty) is a vepyr GAP and the
#: single most actionable thing the run can find, so it leads. `flag_expected` is the
#: one category that is definitionally NOT a bug (the combo's flags never asked VEP for
#: the field, so VEP left it empty while vepyr computed it anyway), so it goes last --
#: it is the only category allowed to be large without being alarming, and putting it
#: anywhere but the bottom lets it dominate the table and bury the rest.
CATEGORY_ORDER: tuple[Category, ...] = (
    Category.VEP_ONLY,      # vepyr GAP -- highest priority
    Category.VALUE_DIFF,    # both populated and different -- the real signal
    Category.VEPYR_ONLY,    # vepyr populated, VEP empty
    Category.ORDER_DIFF,    # same set, different order -- cosmetic
    Category.NUMERIC_TOL,   # agree within tolerance -- rounding
    Category.FLAG_EXPECTED, # not a bug, by construction -- always last
)

_CATEGORY_RANK: dict[str, int] = {str(c): i for i, c in enumerate(CATEGORY_ORDER)}

#: The counters `coverage()` reads out of a summary. A summary that predates one of
#: them (or a `status != "ok"` stub, which carries none of them) reads 0 rather than
#: raising -- but see `coverage()`: a stub can never be `trustworthy`, so a missing
#: counter cannot masquerade as a clean one.
_COUNTERS = ("only_vepyr", "only_gt", "malformed_vepyr", "malformed_gt",
             "duplicate_records_vepyr", "duplicate_records_gt")


def _per_field(summary: dict[str, Any]) -> dict[str, dict]:
    """per_field of a diffed combo; empty for a `no_gt`/failed stub."""
    return summary.get("per_field") or {}


def true_pct(summary: dict[str, Any]) -> float | None:
    """Concordance with flag-semantics expected-diffs counted as AGREEMENT.

    A combo that does not pass `--hgvs` never asks VEP for HGVSc/HGVSp, so VEP leaves
    them empty while vepyr computes them anyway. Those are not disagreements about
    biology, they are disagreements about what was requested, and `overall_pct` (which
    counts them as mismatches) understates vepyr on exactly the combos where it is
    doing MORE work than it was asked to.

    The forgiveness is per-ANNOTATION, never per-field. `oracle.drift.is_flag_expected`
    fires only on the one shape the flag premise actually explains (vepyr non-empty AND
    VEP empty), so a field can be part flag-expected and part real: a vepyr regression
    that stopped emitting HGVSc lands in `vep_only`, a wrong HGVSc lands in
    `value_diff`, and both still drag this number down. The legacy report excluded the
    whole FIELD by name whenever the combo lacked `--hgvs`, which erased those real
    findings from the headline figure entirely.

    Returns None when nothing was compared (an undiffed combo has no concordance, and
    must not be rendered as 0% or 100%).
    """
    total = matched = 0
    for stats in _per_field(summary).values():
        if not stats["total"]:
            continue
        total += stats["total"]
        matched += stats["match"] + stats.get("by_category", {}).get(
            str(Category.FLAG_EXPECTED), 0)
    return round(100 * matched / total, 3) if total else None


def coverage(summary: dict[str, Any]) -> dict[str, Any]:
    """The honesty panel: what did we NOT compare, and can the headline % be trusted?

    Returns the raw coverage counters plus a `trustworthy: bool` and a list of
    plain-string `warnings` the report prints VERBATIM and PROMINENTLY.

    `trustworthy` is not a quality judgement about vepyr -- it is a statement about
    whether `overall_pct` / `true_pct` can be read at face value at all. It is False
    whenever ANY annotation escaped comparison, because every one of the conditions
    below hollows out the denominator the headline % is computed over:

    * join rate below `JOIN_RATE_FLOOR`, or undefined (nothing was read),
    * annotations present on ONE side only -- `only_gt` is the vepyr-GAP direction
      (VEP emitted a transcript annotation vepyr never produced) and `only_vepyr` its
      mirror; on a `--pick` combo the two are the SAME variants, disagreeing about
      which transcript to pick, and neither reaches the classifier,
    * CSQ entries that could not be keyed (`malformed_*`) or repeated
      (chrom,pos,ref,alt) records whose surplus was excluded (`duplicate_records_*`),
    * whole CSQ COLUMNS only one side emits (`vepyr_only_fields` / `gt_only_fields`):
      these are dropped from the diff outright, leaving a report that is 100% green on
      the fields that did get compared,
    * a combo that was never diffed at all (`status != "ok"`).

    Being strict here is the point. A false "trustworthy" is the one error this whole
    pipeline cannot survive; a false alarm costs a sentence of explanation in the PDF.
    """
    counters = {k: summary.get(k, 0) or 0 for k in _COUNTERS}
    join_rate = summary.get("join_rate")
    aligned = summary.get("aligned_annotations", 0) or 0
    status = summary.get("status", "ok")
    vepyr_only_fields = list(summary.get("vepyr_only_fields") or [])
    gt_only_fields = list(summary.get("gt_only_fields") or [])

    warnings: list[str] = []

    if status != "ok":
        warnings.append(
            f"this combo has status {status!r}: it was NOT diffed against ground "
            f"truth, and contributes no evidence about vepyr's correctness."
        )

    if join_rate is None:
        if status == "ok":
            warnings.append(
                "join rate is undefined: NO annotations were read from one or both "
                "VCFs, so nothing whatsoever was compared."
            )
    elif join_rate < JOIN_RATE_FLOOR:
        warnings.append(
            f"{(1 - join_rate) * 100:.1f}% of annotations were never compared "
            f"(join rate {join_rate:.1%}) - the concordance figure below is "
            f"conditional on the two tools already agreeing on which "
            f"(variant, allele, transcript) annotations to emit. An annotation only "
            f"one tool emits is not a mismatch here: it is invisible."
        )

    if counters["only_gt"]:
        warnings.append(
            f"{counters['only_gt']:,} annotations present in VEP but absent from vepyr "
            f"(potential vepyr gaps) were not compared - they are missing from BOTH "
            f"the numerator and the denominator of the concordance figure."
        )
    if counters["only_vepyr"]:
        warnings.append(
            f"{counters['only_vepyr']:,} annotations present in vepyr but absent from "
            f"VEP were not compared - on a --pick combo these are the same variants as "
            f"the VEP-only ones above, with the two tools picking different transcripts."
        )

    for side, label in (("vepyr", "vepyr"), ("gt", "VEP")):
        if counters[f"malformed_{side}"]:
            warnings.append(
                f"{counters[f'malformed_{side}']:,} CSQ entries in {label} could not be "
                f"keyed and were excluded from the join."
            )
        if counters[f"duplicate_records_{side}"]:
            warnings.append(
                f"{counters[f'duplicate_records_{side}']:,} duplicate (chrom,pos,ref,alt) "
                f"records in {label}: the surplus copies were excluded from the join."
            )

    for fields, label in ((gt_only_fields, "VEP"), (vepyr_only_fields, "vepyr")):
        if fields:
            warnings.append(
                f"{len(fields)} CSQ column(s) emitted only by {label} were not compared "
                f"at all (no counterpart column on the other side): "
                f"{', '.join(fields)}."
            )

    return {
        "join_rate": join_rate,
        "aligned": aligned,
        **counters,
        "not_compared": counters["only_vepyr"] + counters["only_gt"],
        "vepyr_only_fields": vepyr_only_fields,
        "gt_only_fields": gt_only_fields,
        "trustworthy": not warnings,
        "warnings": warnings,
    }


def headline_risk(summary: dict[str, Any]) -> str | None:
    """The ONE line page 1 prints beside an untrustworthy combo. None if trustworthy.

    Page 1 is the page everyone reads, and the page most people read ONLY -- so its
    one-line statement of the risk must not be able to understate the risk it stands
    for. Rendering it as "N% of annotations were never compared" (a statement about the
    JOIN RATE alone) meant a combo that joined perfectly but silently dropped two whole
    CSQ COLUMNS from the diff was announced on page 1 as

        0.0% of annotations NEVER COMPARED (0 annotations)

    -- technically true, and the exact inverse of the warning it was standing in for.

    So: lead with the un-compared FRACTION when annotations went uncompared (the
    dominant and most quantitative failure), otherwise lead with the first warning that
    actually fired, and in both cases say how many further warnings are waiting on the
    coverage page. The line always points AT the full warnings; it never replaces them.
    """
    cov = coverage(summary)
    if cov["trustworthy"]:
        return None

    if summary.get("status", "ok") != "ok":
        return f"status {summary['status']!r}: never diffed against ground truth"

    jr, extra = cov["join_rate"], len(cov["warnings"]) - 1
    more = f"  (+{extra} more warning{'s' if extra > 1 else ''} on p.2)" if extra > 0 else ""

    if jr is None:
        return f"no annotations were compared AT ALL (join rate undefined){more}"
    if cov["not_compared"]:
        return (f"{(1 - jr) * 100:.1f}% of annotations NEVER COMPARED "
                f"({cov['not_compared']:,} annotations){more}")
    # A clean join is NOT a clean bill of health -- something else fired.
    lead = cov["warnings"][0].rstrip(".")
    return f"{lead}{more}"


def rollup_categories(summaries: list[dict[str, Any]]) -> dict[str, int]:
    """Drift-category totals across every combo, in `CATEGORY_ORDER`.

    Every category is present even at zero: an absent row is indistinguishable from a
    row nobody computed, and "vep_only: 0" (no vepyr gaps found) is a RESULT worth
    printing, not an absence worth hiding.
    """
    totals = {str(c): 0 for c in CATEGORY_ORDER}
    for summary in summaries:
        for stats in _per_field(summary).values():
            for category, n in stats.get("by_category", {}).items():
                totals[str(category)] = totals.get(str(category), 0) + n
    return totals


def drift_rows(summary: dict[str, Any]) -> list[dict[str, Any]]:
    """One row per field that drifted, WORST FIRST.

    Ordered by REAL (non-flag-expected) mismatches, so HGVSc's 100k expected-diffs on a
    non-hgvs combo cannot crowd a 5-mismatch genuine vepyr gap off the top of the table
    -- surfacing real drift is the only reason the table exists.

    `flag_expected` tags a field whose mismatches are ENTIRELY flag-semantics, i.e. the
    one case a reader may safely dismiss. A field that is 90% expected and 10% real gap
    is NOT tagged: labelling it "not a bug" is precisely how a real gap gets hidden.
    """
    rows = []
    for field, stats in _per_field(summary).items():
        mismatches = stats["total"] - stats["match"]
        if mismatches <= 0:
            continue
        by_category = dict(stats.get("by_category", {}))
        expected = by_category.get(str(Category.FLAG_EXPECTED), 0)
        rows.append({
            "field": field,
            "mismatches": mismatches,
            "real_mismatches": mismatches - expected,
            "pct": stats["pct"],
            "dominant": _dominant(by_category),
            "by_category": by_category,
            "flag_expected": expected == mismatches,
        })
    rows.sort(key=lambda r: (-r["real_mismatches"], -r["mismatches"], r["field"]))
    return rows


def _dominant(by_category: dict[str, int]) -> str | None:
    """The category behind most of a field's mismatches; ties go to the more serious
    one (`CATEGORY_ORDER`), never to whichever happened to be inserted first."""
    if not by_category:
        return None
    return max(by_category.items(),
               key=lambda kv: (kv[1], -_CATEGORY_RANK.get(str(kv[0]), len(CATEGORY_ORDER))))[0]
