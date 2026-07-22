"""Smoke tests for the report CLIs (bin/merge_summaries.py, bin/make_report.py).

The pure logic these two render is tested in tests/test_report.py; what matters here
is that they RUN end-to-end on a realistic summary, and that the merge step cannot
quietly publish a whole-genome number it does not have the shards for.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from oracle.drift import Category
from oracle.matrix import COMBOS
from oracle.summary import finalize, new_accumulator, record_match, record_mismatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MERGE = os.path.join(ROOT, "bin", "merge_summaries.py")
REPORT = os.path.join(ROOT, "bin", "make_report.py")

FIELDS = ["Consequence", "SIFT", "HGVSc"]


def _run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], capture_output=True, text=True,
                          timeout=120)


def _load_make_report():
    """bin/make_report.py as a module (it is a script, not a package member)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_make_report", REPORT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _shard(name: str, cache: str, chrom: str, *, match: int = 100, vep_only: int = 0,
           only_gt: int = 0, only_vepyr: int = 0) -> dict:
    """A shard summary carrying EXACTLY the key set diff_files() emits."""
    acc = new_accumulator(FIELDS)
    for field in FIELDS:
        for _ in range(match):
            record_match(acc, field)
    for i in range(vep_only):
        record_mismatch(acc, "SIFT", Category.VEP_ONLY,
                        [f"{chrom}:{1000 + i} A>G G|ENST{i}", "", "deleterious(0.01)"])
    for i in range(match):
        record_mismatch(acc, "HGVSc", Category.FLAG_EXPECTED,
                        [f"{chrom}:{2000 + i} A>G G|ENST{i}", "ENST1.1:c.10A>G", ""])
    acc["aligned"] = match
    seen = acc["aligned"] + only_vepyr + only_gt
    return finalize(
        acc, name=name, cache=cache,
        only_vepyr=only_vepyr, only_gt=only_gt,
        join_rate=round(acc["aligned"] / seen, 4) if seen else None,
        shared_fields=len(FIELDS), vepyr_only_fields=[], gt_only_fields=[],
        malformed_vepyr=0, malformed_gt=0,
        duplicate_records_vepyr=0, duplicate_records_gt=0,
        mismatches_tsv=f"mismatches/{name}_{chrom}.tsv", chrom=chrom,
    )


def _seed_results(dirpath, chroms=("21", "22"), **per_combo_over) -> None:
    """One shard JSON per (combo, chrom), exactly as bin/validate.py writes them."""
    os.makedirs(dirpath, exist_ok=True)
    for combo in COMBOS:
        for chrom in chroms:
            over = per_combo_over.get(combo.name, {})
            shard = _shard(combo.name, f"115_GRCh38_{combo.cache_flavor}", chrom, **over)
            with open(os.path.join(dirpath, f"summary_{combo.name}_{chrom}.json"), "w") as fh:
                json.dump(shard, fh)


# --- merge_summaries -----------------------------------------------------------

def test_merge_writes_a_list_of_per_combo_whole_genome_summaries(tmp_path):
    results = tmp_path / "results"
    _seed_results(str(results))
    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode == 0, r.stderr

    merged = json.loads((results / "summary.json").read_text())
    assert isinstance(merged, list)
    assert [s["name"] for s in merged] == [c.name for c in COMBOS]
    assert all(s["chroms"] == ["21", "22"] for s in merged)
    assert all(s["aligned_annotations"] == 200 for s in merged)     # 100 x 2 shards
    # the merge step itself must never print a concordance % without its join rate
    assert "overall_pct=" in r.stdout and "join_rate=" in r.stdout


def _seed_subset(dirpath, names, chrom="22") -> None:
    """Seed shard JSONs for ONLY the named combos (a gate running a subset)."""
    os.makedirs(dirpath, exist_ok=True)
    by_name = {c.name: c for c in COMBOS}
    for name in names:
        combo = by_name[name]
        shard = _shard(name, f"115_GRCh38_{combo.cache_flavor}", chrom)
        with open(os.path.join(dirpath, f"summary_{name}_{chrom}.json"), "w") as fh:
            json.dump(shard, fh)


def test_merge_combos_subset_requires_only_the_subset(tmp_path):
    """--combos restricts the completeness contract to the named subset, so a gate
    that ran fewer than every combo does not trip the 'combo omitted' FATAL."""
    results = tmp_path / "results"
    _seed_subset(str(results), ["hgvs_merged", "hgvs_refseq"])
    r = _run([MERGE, str(results), "--chroms", "22",
              "--combos", "hgvs_merged,hgvs_refseq"])
    assert r.returncode == 0, r.stderr
    merged = json.loads((results / "summary.json").read_text())
    assert [s["name"] for s in merged] == ["hgvs_merged", "hgvs_refseq"]


def test_merge_without_combos_still_requires_all(tmp_path):
    """Default (no --combos) keeps the whole-genome contract: a subset dir FATALs."""
    results = tmp_path / "results"
    _seed_subset(str(results), ["hgvs_merged", "hgvs_refseq"])
    r = _run([MERGE, str(results), "--chroms", "22"])
    assert r.returncode == 2
    assert "no shard summaries" in r.stderr


def test_merge_combos_rejects_unknown_name(tmp_path):
    results = tmp_path / "results"
    _seed_subset(str(results), ["hgvs_merged"])
    r = _run([MERGE, str(results), "--chroms", "22", "--combos", "hgvs_merged,nonsuch"])
    assert r.returncode == 2
    assert "not in oracle.matrix.COMBOS" in r.stderr


def test_merge_refuses_to_publish_a_whole_genome_number_from_21_chromosomes(tmp_path):
    """expected_chroms is the whole point: a dead chromosome's missing shard must NOT
    silently yield a 'whole-genome' number computed from the survivors."""
    results = tmp_path / "results"
    _seed_results(str(results), chroms=("21", "22"))
    os.remove(results / f"summary_{COMBOS[0].name}_21.json")

    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode != 0
    assert "21" in (r.stderr + r.stdout)
    assert "missing" in (r.stderr + r.stdout).lower()
    assert not (results / "summary.json").exists()


def test_merge_defaults_to_the_22_autosomes(tmp_path):
    """The default must be the strict one -- a 2-chromosome results dir is NOT a
    whole-genome run, and merging it without saying so is the bug we are preventing."""
    results = tmp_path / "results"
    _seed_results(str(results), chroms=("21", "22"))
    r = _run([MERGE, str(results)])
    assert r.returncode != 0
    assert "missing" in (r.stderr + r.stdout).lower()


def test_merge_does_not_swallow_a_sibling_combos_shards_via_prefix_collision(tmp_path):
    """`glob("summary_hgvs_merged_*.json")` also matches
    `summary_hgvs_merged_pick_allele_21.json`. Dispatching on the name INSIDE the JSON
    is what makes this impossible -- if it ever regresses, the pick combos' shards get
    merged into hgvs_merged and the counters silently double."""
    results = tmp_path / "results"
    _seed_results(str(results))
    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode == 0, r.stderr

    merged = {s["name"]: s for s in json.loads((results / "summary.json").read_text())}
    assert merged["hgvs_merged"]["aligned_annotations"] == 200
    assert merged["hgvs_merged"]["chroms"] == ["21", "22"]


def test_merge_passes_a_no_gt_combo_through_instead_of_dropping_it(tmp_path):
    results = tmp_path / "results"
    _seed_results(str(results))
    dead = COMBOS[1].name
    for chrom in ("21", "22"):
        (results / f"summary_{dead}_{chrom}.json").write_text(json.dumps(
            {"name": dead, "cache": "115_GRCh38_merged", "status": "no_gt"}))

    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode == 0, r.stderr
    merged = {s["name"]: s for s in json.loads((results / "summary.json").read_text())}
    assert merged[dead]["status"] == "no_gt"
    assert len(merged) == len(COMBOS)


def test_merge_refuses_a_combo_that_is_half_diffed_half_no_gt(tmp_path):
    results = tmp_path / "results"
    _seed_results(str(results))
    (results / f"summary_{COMBOS[0].name}_21.json").write_text(json.dumps(
        {"name": COMBOS[0].name, "cache": "115_GRCh38_merged", "status": "no_gt"}))
    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode != 0


def test_merge_fails_when_a_combo_has_no_shards_at_all(tmp_path):
    results = tmp_path / "results"
    _seed_results(str(results))
    for chrom in ("21", "22"):
        os.remove(results / f"summary_{COMBOS[-1].name}_{chrom}.json")
    r = _run([MERGE, str(results), "--chroms", "21,22"])
    assert r.returncode != 0
    assert COMBOS[-1].name in r.stderr


# --- make_report ---------------------------------------------------------------

@pytest.fixture
def mixed_summary(tmp_path):
    """One trustworthy combo and one that publishes 100.0% at join_rate 0.843."""
    results = tmp_path / "results"
    _seed_results(str(results))
    _run([MERGE, str(results), "--chroms", "21,22"])
    merged = json.loads((results / "summary.json").read_text())

    liar = next(s for s in merged if s["name"] == "hgvs_merged_pick_allele")
    liar.update(overall_pct=100.0, join_rate=0.843, aligned_annotations=84_300,
                only_vepyr=7_850, only_gt=7_850)
    for f, stats in liar["per_field"].items():
        stats.update(total=84_300, match=84_300, pct=100.0, by_category={})
    liar["samples"] = {}

    path = tmp_path / "summary.json"
    path.write_text(json.dumps(merged))
    return path


def test_make_report_renders_a_pdf_and_says_which_combos_cannot_be_trusted(
        mixed_summary, tmp_path):
    out = tmp_path / "report.pdf"
    r = _run([REPORT, str(mixed_summary), str(out), "(smoke)"])
    assert r.returncode == 0, r.stderr
    assert out.exists() and out.stat().st_size > 20_000
    assert out.read_bytes().startswith(b"%PDF")
    assert "NOT TRUSTWORTHY" in r.stdout
    assert "hgvs_merged_pick_allele" in r.stdout


def test_the_coverage_page_precedes_the_heatmaps(mixed_summary, tmp_path, monkeypatch):
    """Page order is a safety property, not a layout preference: the trust page must
    not sit behind the pretty pictures where a reader can skip past it."""
    mr = _load_make_report()
    order: list[str] = []
    for fn in ("page_overview", "page_trust", "page_drift_summary", "page_heatmaps",
               "page_combo_drift"):
        monkeypatch.setattr(mr, fn, (lambda name: lambda *a, **k: order.append(name))(fn))

    mr.build_report(json.loads(mixed_summary.read_text()), str(tmp_path / "o.pdf"))

    assert order[:4] == ["page_overview", "page_trust", "page_drift_summary",
                         "page_heatmaps"]
    assert order.index("page_trust") < order.index("page_heatmaps")


def test_the_trust_page_is_handed_the_uncompared_fraction_verbatim(mixed_summary,
                                                                   tmp_path, monkeypatch):
    """The join_rate=0.843 combo's warning must reach the PAGE, not just the summary
    dict -- i.e. the renderer prints coverage() warnings rather than paraphrasing them."""
    mr = _load_make_report()
    printed: list[str] = []
    real = mr._text_pages
    monkeypatch.setattr(mr, "_text_pages",
                        lambda pdf, title, lines, **kw: (
                            printed.extend(t for t, _c, _s in lines) if title.startswith(
                                "COVERAGE") else None,
                            real(pdf, title, lines, **kw))[1])

    mr.build_report(json.loads(mixed_summary.read_text()), str(tmp_path / "o.pdf"))

    # normalise whitespace: the warnings are word-wrapped for the page, so a phrase may
    # legitimately straddle a line break -- the CONTENT is what must survive, not the
    # column it landed in.
    page = " ".join(" ".join(printed).split())
    assert "hgvs_merged_pick_allele" in page
    assert "15.7% of annotations were never compared" in page
    assert "conditional on the two tools already agreeing" in page
    assert "potential vepyr gaps" in page


def test_page_1_cannot_understate_a_risk_it_stands_for(mixed_summary, tmp_path,
                                                       monkeypatch):
    """A reader who gets no further than page 1 must not be misled.

    The regression: page 1 phrased every risk as "N% of annotations were never
    compared" -- a statement about the JOIN RATE alone -- so a combo that joined
    perfectly but dropped two whole CSQ COLUMNS from the diff was announced as
    "0.0% of annotations NEVER COMPARED (0 annotations)". Page 1 now prints
    oracle.report.headline_risk(), which leads with whatever ACTUALLY fired.
    """
    import matplotlib.figure

    summaries = json.loads(mixed_summary.read_text())
    victim = next(s for s in summaries if s["name"] == "hgvs_refseq")
    victim.update(join_rate=1.0, only_vepyr=0, only_gt=0,
                  gt_only_fields=["SpliceAI_pred_DS_AG", "am_pathogenicity"])

    drawn: list[str] = []
    real_text = matplotlib.figure.Figure.text
    monkeypatch.setattr(matplotlib.figure.Figure, "text",
                        lambda self, x, y, s, **kw: (drawn.append(str(s)),
                                                     real_text(self, x, y, s, **kw))[1])

    mr = _load_make_report()
    mr.page_overview(_NullPdf(), summaries, [mr.coverage(s) for s in summaries], "")

    page1 = " ".join(" ".join(drawn).split())
    assert "hgvs_refseq" in page1
    assert "CSQ column" in page1, "page 1 must NAME the columns that went uncompared"
    assert "hgvs_refseq: overall 99.906%, but 0.0%" not in page1
    # ...and the join-rate case still leads with the fraction, as before
    assert "15.7% of annotations NEVER COMPARED" in page1


class _NullPdf:
    """PdfPages stand-in: page_overview() draws, we inspect, nothing is written."""

    def savefig(self, fig) -> None:
        pass


def test_make_report_rejects_a_summary_that_is_not_a_list(tmp_path):
    bad = tmp_path / "summary.json"
    bad.write_text(json.dumps({"name": "hgvs_merged"}))
    r = _run([REPORT, str(bad), str(tmp_path / "x.pdf")])
    assert r.returncode != 0
    assert "merge_summaries" in r.stderr


def test_make_report_refuses_to_emit_an_empty_pdf(tmp_path):
    """An empty PDF renders as a completed validation that found nothing wrong."""
    empty = tmp_path / "summary.json"
    empty.write_text("[]")
    out = tmp_path / "x.pdf"
    r = _run([REPORT, str(empty), str(out)])
    assert r.returncode != 0
    assert "no combos" in r.stderr
    assert not out.exists()


def test_make_report_renders_a_combo_that_was_never_diffed(tmp_path):
    """A `no_gt` combo has no per_field/samples at all -- it must still render (as a
    combo contributing NO evidence), not crash the whole report."""
    path = tmp_path / "summary.json"
    path.write_text(json.dumps([{"name": "hgvs_refseq", "cache": "115_GRCh38_refseq",
                                 "status": "no_gt"}]))
    out = tmp_path / "x.pdf"
    r = _run([REPORT, str(path), str(out)])
    assert r.returncode == 0, r.stderr
    assert out.exists()
    assert "NOT TRUSTWORTHY" in r.stdout
