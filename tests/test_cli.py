"""Smoke tests for the one-combo CLIs (bin/run_wgs.py, bin/validate.py).

No vepyr, no cluster, no network. Two properties matter above all else:

1. Both CLIs must stay importable/runnable (at least for --help and the resume /
   no-gt short-circuits) on a machine with no vepyr installed -- that is what keeps
   the whole oracle/ test suite laptop-runnable, and it is the exact property
   `legacy/validate.py` broke by importing vepyr at module top.
2. `oracle.matrix.resolve_kwargs()` must never let the plugin-cache placeholder
   silently reach `vepyr.annotate()` -- see test_matrix.py for the exhaustive
   sentinel-vs-string coverage; the checks here just re-confirm the CLI-facing
   contract inline with the rest of this file's scenarios.

Paths are made testable via `VEPYR_DATA` / `VEPYR_WORK` env var overrides, which
both CLIs read instead of hardcoding `~/vepyr/data` and `~/vepyr/work` -- see
bin/run_wgs.py and bin/validate.py.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

import pytest

from oracle.matrix import MATRIX_HEADER, NEEDS_PLUGIN_CACHE, resolve_kwargs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_WGS = os.path.join(ROOT, "bin", "run_wgs.py")
VALIDATE = os.path.join(ROOT, "bin", "validate.py")
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")


def _write_matrix(path: str, **overrides) -> None:
    row = {
        "name": "hgvs_merged",
        "cache_flavor": "merged",
        "cache115": "115_GRCh38_merged",
        "cache116": "116_GRCh38_merged",
        "vepyr_kwargs": json.dumps({"everything": True, "hgvs": True}),
        "vep_flags": "vep --everything --hgvs",
        "gt115": "hgvs_merged.vcf",
        "gt116": "hgvs_merged.vcf",
    }
    row.update(overrides)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MATRIX_HEADER, delimiter="\t")
        w.writeheader()
        w.writerow(row)


def _run(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], env=env,
                          capture_output=True, text=True, timeout=30)


# --- resolve_kwargs: the CLI-facing contract ---------------------------------
#
# The exhaustive sentinel-vs-round-tripped-string matrix lives in test_matrix.py;
# these three confirm the shape bin/run_wgs.py actually depends on (it calls
# resolve_kwargs() on kwargs freshly parsed from matrix.tsv JSON, i.e. the STRING
# form, not the live sentinel object).

def test_resolve_kwargs_raises_when_the_plugin_path_is_missing():
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs({"plugin_cache_root": NEEDS_PLUGIN_CACHE})
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs(json.loads(json.dumps({"plugin_cache_root": "@PLUGIN@"})))


def test_resolve_kwargs_substitutes_when_given():
    out = resolve_kwargs({"hgvs": True, "plugin_cache_root": NEEDS_PLUGIN_CACHE},
                         plugin_cache_root="/data/plugin_cache")
    assert out["plugin_cache_root"] == "/data/plugin_cache"
    assert out["hgvs"] is True


def test_resolve_kwargs_leaves_non_plugin_combos_untouched():
    kw = {"everything": True, "hgvs": True, "pick": True}
    assert resolve_kwargs(kw) == kw
    assert resolve_kwargs(kw, plugin_cache_root="/anything") == kw


# --- lazy import: pins the property the whole suite depends on --------------
#
# This environment (the .venv these tests run under) has no vepyr installed --
# confirmed separately, and load-bearing here: if `import vepyr` were ever moved
# to module top in either CLI, THESE tests would fail with ModuleNotFoundError
# surfacing on stderr and a non-zero returncode, instead of a clean --help exit.
# That is what makes this a pin and not just a smoke check.

@pytest.mark.parametrize("script", [RUN_WGS, VALIDATE])
def test_help_does_not_import_vepyr(script):
    result = _run([script, "--help"], env=os.environ.copy())
    assert result.returncode == 0, result.stderr
    assert "ModuleNotFoundError" not in result.stderr
    assert "No module named 'vepyr'" not in result.stderr
    assert "usage:" in result.stdout.lower()


def test_vepyr_is_in_fact_not_installed_here():
    """If this ever fails (vepyr got installed in the test venv), the pin above
    stops meaning anything -- it would pass even with `import vepyr` at module
    top. Fails loudly instead of silently testing nothing."""
    result = subprocess.run([sys.executable, "-c", "import vepyr"],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "ModuleNotFoundError" in result.stderr


# --- resume: both CLIs must be no-ops on an existing, non-empty output ------

def test_run_wgs_resume_skips_without_importing_vepyr(tmp_path):
    out = tmp_path / "out.vcf"
    out.write_text("not a real vcf, just needs to be non-empty\n")
    env = {**os.environ, "VEPYR_DATA": str(tmp_path / "data"),
          "VEPYR_WORK": str(tmp_path / "work")}
    result = _run([RUN_WGS, "--combo", "hgvs_merged", "--version", "115",
                  "--out", str(out)], env=env)
    assert result.returncode == 0, result.stderr
    assert "already exists" in result.stdout
    assert "ModuleNotFoundError" not in result.stderr


def test_validate_resume_skips_an_existing_summary(tmp_path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "summary_hgvs_merged.json").write_text('{"status": "ok"}')
    env = {**os.environ, "VEPYR_DATA": str(tmp_path / "data"),
          "VEPYR_WORK": str(tmp_path / "work")}
    result = _run([VALIDATE, "--combo", "hgvs_merged", "--version", "115",
                  "--vepyr", "/nonexistent.vcf", "--outdir", str(outdir)], env=env)
    assert result.returncode == 0, result.stderr
    assert "already exists" in result.stdout


# --- validate.py: no_gt must exit 0, never crash the array ------------------

def test_validate_writes_no_gt_summary_and_exits_0_when_gt_is_absent(tmp_path):
    data_dir = tmp_path / "data"
    work_dir = tmp_path / "work"
    (data_dir / "ground_truth_vep").mkdir(parents=True)   # dir exists, file does not
    work_dir.mkdir()
    _write_matrix(str(work_dir / "matrix.tsv"))
    outdir = tmp_path / "out"

    env = {**os.environ, "VEPYR_DATA": str(data_dir), "VEPYR_WORK": str(work_dir)}
    result = _run([VALIDATE, "--combo", "hgvs_merged", "--version", "115",
                  "--vepyr", str(tmp_path / "vepyr_out.vcf"),
                  "--outdir", str(outdir)], env=env)

    assert result.returncode == 0, result.stderr
    summary_file = outdir / "summary_hgvs_merged.json"
    assert summary_file.exists()
    summary = json.loads(summary_file.read_text())
    assert summary == {"name": "hgvs_merged", "cache": "115_GRCh38_merged", "status": "no_gt"}
    # no mismatches TSV should have been produced -- diff_files() was never called
    assert not (outdir / "mismatches").exists()


def test_validate_runs_the_real_diff_when_gt_is_present(tmp_path):
    """End-to-end (no vepyr needed: diff only ever reads two already-annotated
    VCFs) -- confirms the CLI wires combo_kwargs / cache / paths correctly
    against oracle.diff.diff_files(), using the same fixtures test_diff.py uses."""
    data_dir = tmp_path / "data"
    work_dir = tmp_path / "work"
    gt_dir = data_dir / "ground_truth_vep"
    gt_dir.mkdir(parents=True)
    work_dir.mkdir()
    (gt_dir / "hgvs_merged.vcf").write_text(
        (open(os.path.join(FIXTURES, "gt_mini.vcf")).read()))
    _write_matrix(str(work_dir / "matrix.tsv"))
    outdir = tmp_path / "out"

    env = {**os.environ, "VEPYR_DATA": str(data_dir), "VEPYR_WORK": str(work_dir)}
    result = _run([VALIDATE, "--combo", "hgvs_merged", "--version", "115",
                  "--vepyr", os.path.join(FIXTURES, "vepyr_mini.vcf"),
                  "--outdir", str(outdir)], env=env)

    assert result.returncode == 0, result.stderr
    summary = json.loads((outdir / "summary_hgvs_merged.json").read_text())
    assert summary["status"] == "ok"
    assert "overall_pct" in summary and "join_rate" in summary
    assert (outdir / "mismatches" / "hgvs_merged.tsv").exists()
    # the printed headline must carry the coverage signals alongside overall_pct
    assert "overall_pct=" in result.stdout and "join_rate=" in result.stdout


# --- an UNSEEDED row must stop both CLIs dead -------------------------------
#
# bin/seed_matrix.py writes EMPTY vepyr_kwargs/vep_flags for a combo whose ground
# truth it could not read (missing header, unmappable flag, cache-flavor mismatch) --
# deliberately, so that "we do not know this combo's semantics" can never be mistaken
# for "this combo needs no flags". Both CLIs must therefore refuse the row outright
# rather than crash on json.loads("") -- or, far worse, treat it as `{}` and annotate
# with vepyr's defaults against a ground truth generated with --everything --hgvs.

@pytest.mark.parametrize("cli", ["run_wgs", "validate"])
def test_cli_refuses_a_combo_whose_kwargs_were_never_seeded(tmp_path, cli):
    data_dir = tmp_path / "data"
    work_dir = tmp_path / "work"
    (data_dir / "ground_truth_vep").mkdir(parents=True)
    work_dir.mkdir()
    _write_matrix(str(work_dir / "matrix.tsv"), vepyr_kwargs="", vep_flags="")

    env = {**os.environ, "VEPYR_DATA": str(data_dir), "VEPYR_WORK": str(work_dir)}
    args = ([RUN_WGS, "--combo", "hgvs_merged", "--version", "115",
             "--chrom", "22", "--out", str(tmp_path / "out.vcf")]
            if cli == "run_wgs" else
            [VALIDATE, "--combo", "hgvs_merged", "--version", "115",
             "--vepyr", "/nonexistent.vcf", "--outdir", str(tmp_path / "out")])
    result = _run(args, env=env)

    assert result.returncode != 0
    assert "hgvs_merged" in result.stderr
    assert "seed_matrix" in result.stderr    # tell the operator what to do about it


def test_validate_unknown_combo_exits_nonzero_with_a_clear_message(tmp_path):
    data_dir = tmp_path / "data"
    work_dir = tmp_path / "work"
    (data_dir / "ground_truth_vep").mkdir(parents=True)
    work_dir.mkdir()
    _write_matrix(str(work_dir / "matrix.tsv"))

    env = {**os.environ, "VEPYR_DATA": str(data_dir), "VEPYR_WORK": str(work_dir)}
    result = _run([VALIDATE, "--combo", "does_not_exist", "--version", "115",
                  "--vepyr", "/nonexistent.vcf", "--outdir", str(tmp_path / "out")], env=env)
    assert result.returncode != 0
    assert "does_not_exist" in result.stderr


# --- GT_DIRS: version -> ground-truth-directory map --------------------------
#
# ground_truth_vep_116 was ground truth minted with `--fork N>1`; deleted 2026-07-30
# because forking narrows the annotation buffer window ~30x and drops HGNC_ID, so it
# was never a valid standard (see bin/validate.py's GT_DIRS comment). The
# replacement, ground_truth_vep_116_unforked, sits beside it under the same
# DATA_DIR. Unlike the CLI-behaviour tests above, this pins the DICT ITSELF (not
# just an observable exit code / stdout line), so a future edit that repoints
# GT_DIRS at a deleted or forked directory fails loudly here instead of silently
# producing status=no_gt summaries for every run, as actually happened from
# 2026-07-30 until this fix.

def _load_validate():
    """bin/validate.py as a module (it is a script, not a package member) -- same
    load-by-path mechanism as _load_make_report() in test_report_cli.py."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_validate", VALIDATE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_gt_dirs_116_points_at_the_unforked_ground_truth(monkeypatch, tmp_path):
    # VEPYR_DATA must be set BEFORE the module executes -- DATA_DIR/GT_DIRS are
    # computed at module-exec time, same as every other CLI in bin/.
    monkeypatch.setenv("VEPYR_DATA", str(tmp_path))
    mod = _load_validate()

    assert mod.GT_DIRS[116].endswith("ground_truth_vep_116_unforked")
    # Positive control: 115 must be untouched -- proves the fix moved exactly one
    # path, not both.
    assert mod.GT_DIRS[115] == os.path.join(mod.DATA_DIR, "ground_truth_vep")
    # DATA_DIR must still be driven by VEPYR_DATA -- the hook a containerised
    # runner uses to point the oracle at its own data directory, unrelated to
    # this fix and must survive it.
    assert mod.DATA_DIR == str(tmp_path)


# --- MATRIX_PATH: $VEPYR_MATRIX > $WORK_DIR/matrix.tsv > the shipped copy ----
#
# matrix.tsv used to exist ONLY on the cluster, at $VEPYR_WORK/matrix.tsv -- so a
# fresh clone of this repo (the container image being built around it, a laptop, a
# CI runner) had no matrix.tsv anywhere and this CLI could not run at all. A copy
# fetched from the cluster is now shipped at <repo root>/matrix.tsv (see
# matrix.tsv.README beside it) as the last link of a 3-way fallback chain -- see
# the comment beside MATRIX_PATH in bin/validate.py for the exact order.
#
# Each test below is built so it FAILS if that one link of the chain is missing or
# out of order -- e.g. test_matrix_path_falls_back_to_the_shipped_copy_when_
# work_dir_has_none would fail against the PRE-FIX code (unconditional
# WORK_DIR/matrix.tsv), because that resolves to a tmp_path under WORK_DIR that
# does not exist, never to the shipped copy. A fallback chain whose failing
# direction is never exercised is not actually tested -- see the task's
# self-review note.

def test_matrix_path_falls_back_to_the_shipped_copy_when_work_dir_has_none(
        monkeypatch, tmp_path):
    monkeypatch.delenv("VEPYR_MATRIX", raising=False)
    monkeypatch.setenv("VEPYR_DATA", str(tmp_path / "data"))
    work = tmp_path / "work"
    work.mkdir()                              # exists, but no matrix.tsv inside
    monkeypatch.setenv("VEPYR_WORK", str(work))

    mod = _load_validate()

    assert mod.MATRIX_PATH == os.path.join(ROOT, "matrix.tsv")
    # Not just the right path -- the shipped copy must actually BE there, or a
    # fresh clone is exactly as unable to run as before this fix.
    assert os.path.exists(mod.MATRIX_PATH)


def test_matrix_path_prefers_the_cluster_copy_when_work_dir_has_one(monkeypatch, tmp_path):
    """The cluster's own freshly-seeded matrix.tsv must keep winning over the
    shipped copy -- this is what pins "existing cluster behaviour is unchanged".
    Against a fix that got the fallback order backwards (shipped copy preferred
    over VEPYR_WORK), this would resolve to <repo root>/matrix.tsv instead and
    fail."""
    monkeypatch.delenv("VEPYR_MATRIX", raising=False)
    monkeypatch.setenv("VEPYR_DATA", str(tmp_path / "data"))
    work = tmp_path / "work"
    work.mkdir()
    cluster_matrix = work / "matrix.tsv"
    cluster_matrix.write_text("name\tcache_flavor\tcache115\tcache116\t"
                              "vepyr_kwargs\tvep_flags\tgt115\tgt116\n")
    monkeypatch.setenv("VEPYR_WORK", str(work))

    mod = _load_validate()

    assert mod.MATRIX_PATH == str(cluster_matrix)
    assert mod.MATRIX_PATH != os.path.join(ROOT, "matrix.tsv")


def test_matrix_path_explicit_override_wins_over_both(monkeypatch, tmp_path):
    """VEPYR_MATRIX must win even when a cluster copy ALSO exists under
    VEPYR_WORK -- the harder case to satisfy by accident (a fix that forgot the
    override entirely would still pass a test where no VEPYR_WORK copy exists).
    Against that wrong fix, this resolves to the VEPYR_WORK copy instead and
    fails."""
    monkeypatch.setenv("VEPYR_DATA", str(tmp_path / "data"))
    work = tmp_path / "work"
    work.mkdir()
    (work / "matrix.tsv").write_text("name\tcache_flavor\tcache115\tcache116\t"
                                     "vepyr_kwargs\tvep_flags\tgt115\tgt116\n")
    explicit = tmp_path / "pinned_matrix.tsv"
    explicit.write_text("name\tcache_flavor\tcache115\tcache116\t"
                        "vepyr_kwargs\tvep_flags\tgt115\tgt116\n")
    monkeypatch.setenv("VEPYR_WORK", str(work))
    monkeypatch.setenv("VEPYR_MATRIX", str(explicit))

    mod = _load_validate()

    assert mod.MATRIX_PATH == str(explicit)
