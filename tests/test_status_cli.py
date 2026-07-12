"""Smoke tests for bin/status.py.

The state/overrun LOGIC is exhaustively tested in tests/test_status.py against
oracle.status directly; these confirm the CLI wiring: which shards.tsv files get
read per --version, path defaults via VEPYR_DATA/VEPYR_WORK, and that a missing
shards.tsv fails loudly rather than silently reporting nothing.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

from oracle.shards import SHARD_HEADER

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATUS = os.path.join(ROOT, "bin", "status.py")


def _run(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, STATUS, *args], env=env,
                          capture_output=True, text=True, timeout=30)


def _write_shards(path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SHARD_HEADER, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _row(combo="hgvs_merged", step="annotate", level="L1", chrom="22", chunk="",
        n_variants=1000, est_hours=10.0) -> dict:
    return {"combo": combo, "step": step, "level": level, "chrom": chrom, "chunk": chunk,
           "n_variants": n_variants, "est_hours": est_hours,
           "est_hours_raw": est_hours / 1.4, "safety_factor": 1.4}


def test_missing_shards_files_fails_loudly(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    env = {**os.environ, "VEPYR_WORK": str(work), "VEPYR_DATA": str(tmp_path / "data")}
    r = _run(["--version", "115"], env=env)
    assert r.returncode != 0
    assert "no shards.tsv" in r.stderr.lower() or "no shards.tsv" in r.stdout.lower()


def test_115_does_not_read_the_gt_step(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    _write_shards(work / "shards_annotate.tsv", [_row()])
    _write_shards(work / "shards_diff.tsv", [_row(step="diff")])
    # a gt shards file exists too, but --version 115 must ignore it (no gt step
    # for 115 -- its ground truth already exists as static files)
    _write_shards(work / "shards_gt.tsv", [_row(step="gt")])

    env = {**os.environ, "VEPYR_WORK": str(work), "VEPYR_DATA": str(tmp_path / "data")}
    r = _run(["--version", "115"], env=env)
    assert r.returncode == 0, r.stderr
    assert "=== annotate " in r.stdout
    assert "=== diff " in r.stdout
    assert "=== gt " not in r.stdout


def test_116_reads_all_three_steps(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    _write_shards(work / "shards_annotate.tsv", [_row()])
    _write_shards(work / "shards_diff.tsv", [_row(step="diff")])
    _write_shards(work / "shards_gt.tsv", [_row(step="gt")])

    env = {**os.environ, "VEPYR_WORK": str(work), "VEPYR_DATA": str(tmp_path / "data")}
    r = _run(["--version", "116"], env=env)
    assert r.returncode == 0, r.stderr
    assert "=== annotate " in r.stdout
    assert "=== diff " in r.stdout
    assert "=== gt " in r.stdout


def test_reports_missing_and_overrun_elements(tmp_path):
    work = tmp_path / "work"
    results = work / "results_wgs_115"
    results.mkdir(parents=True)
    _write_shards(work / "shards_annotate.tsv", [
        _row(combo="hgvs_merged", chrom="21", est_hours=10.0),
        _row(combo="hgvs_merged", chrom="22", est_hours=1.0),
    ])
    _write_shards(work / "shards_diff.tsv", [])

    # chrom 21: overran (actual way over est_hours)
    (results / "vepyr_hgvs_merged_21.vcf").write_text("x")
    (results / "vepyr_hgvs_merged_21.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 10.0 * 3600 * 3}))
    # chrom 22: missing entirely

    env = {**os.environ, "VEPYR_WORK": str(work), "VEPYR_DATA": str(tmp_path / "data")}
    r = _run(["--version", "115"], env=env)
    assert r.returncode == 0, r.stderr
    assert "OVERRAN" in r.stdout
    assert "missing" in r.stdout
    assert "chrom=21" in r.stdout
    assert "chrom=22" in r.stdout
    assert "WARNING: 1 element(s) ran more than" in r.stdout


def test_help_works_without_a_shards_file():
    r = subprocess.run([sys.executable, STATUS, "--help"],
                       capture_output=True, text=True, timeout=10)
    assert r.returncode == 0, r.stderr
    assert "usage:" in r.stdout.lower()
