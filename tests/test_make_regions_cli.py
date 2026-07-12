"""Smoke tests for bin/make_regions.py.

The quantile-cutting logic itself is exhaustively tested in tests/test_regions.py
against oracle.regions directly -- these confirm the CLI wiring: path defaults via
VEPYR_WORK, the resume/mismatch behaviour around an existing regions file, and
that a missing input fails loudly instead of writing a bogus/empty regions file.
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAKE_REGIONS = os.path.join(ROOT, "bin", "make_regions.py")

HEADER = "##fileformat=VCFv4.2\n##contig=<ID=22>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _run(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, MAKE_REGIONS, *args], env=env,
                          capture_output=True, text=True, timeout=30)


def _write_input(work_dir, chrom: str, positions: list[int]) -> None:
    # oracle.csq.open_maybe_gzip only branches on a literal ".gz" suffix, so a
    # plain-text ".vcf" (passed explicitly via --input below) exercises the
    # non-gzip path without needing an actual bgzip dependency in the test.
    input_dir = work_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    lines = [f"22\t{p}\t.\tA\tG\t.\t.\t.\n" for p in positions]
    (input_dir / f"input_{chrom}.vcf").write_text(HEADER + "".join(lines))


def test_writes_a_regions_file_and_prints_a_summary(tmp_path):
    work = tmp_path / "work"
    _write_input(work, "22", list(range(1, 101)))
    env = {**os.environ, "VEPYR_WORK": str(work)}
    out = work / "regions_22.tsv"
    r = _run(["--chrom", "22", "--num-chunks", "4",
             "--input", str(work / "input" / "input_22.vcf"), "--out", str(out)], env=env)
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert len(rows) == 4
    assert sum(int(row["n_variants"]) for row in rows) == 100
    assert "4 regions" in r.stdout


def test_resumes_when_the_existing_file_matches_num_chunks(tmp_path):
    work = tmp_path / "work"
    _write_input(work, "22", list(range(1, 51)))
    env = {**os.environ, "VEPYR_WORK": str(work)}
    out = work / "regions_22.tsv"
    inp = str(work / "input" / "input_22.vcf")
    first = _run(["--chrom", "22", "--num-chunks", "3", "--input", inp, "--out", str(out)],
                env=env)
    assert first.returncode == 0, first.stderr
    written = out.read_text()

    second = _run(["--chrom", "22", "--num-chunks", "3", "--input", inp, "--out", str(out)],
                 env=env)
    assert second.returncode == 0, second.stderr
    assert "resuming" in second.stdout
    assert out.read_text() == written    # untouched, not silently rewritten


def test_refuses_to_reuse_a_regions_file_cut_at_a_different_resolution(tmp_path):
    work = tmp_path / "work"
    _write_input(work, "22", list(range(1, 51)))
    env = {**os.environ, "VEPYR_WORK": str(work)}
    out = work / "regions_22.tsv"
    inp = str(work / "input" / "input_22.vcf")
    first = _run(["--chrom", "22", "--num-chunks", "3", "--input", inp, "--out", str(out)],
                env=env)
    assert first.returncode == 0, first.stderr

    mismatched = _run(["--chrom", "22", "--num-chunks", "5", "--input", inp, "--out", str(out)],
                     env=env)
    assert mismatched.returncode != 0
    assert "disagree" in mismatched.stderr


def test_missing_input_fails_loudly(tmp_path):
    work = tmp_path / "work"
    (work / "input").mkdir(parents=True)
    env = {**os.environ, "VEPYR_WORK": str(work)}
    out = work / "regions_22.tsv"
    r = _run(["--chrom", "22", "--num-chunks", "2",
             "--input", str(work / "input" / "input_22.vcf"), "--out", str(out)], env=env)
    assert r.returncode != 0
    assert "does not exist" in r.stderr
    assert not out.exists()


def test_help_does_not_require_pysam():
    """oracle/regions.py reads VCFs with plain text/gzip, not pysam -- make_regions.py
    must therefore stay runnable (at least --help) with no pysam installed, same
    laptop-testability property bin/run_wgs.py pins for vepyr."""
    r = subprocess.run([sys.executable, MAKE_REGIONS, "--help"],
                       capture_output=True, text=True, timeout=10)
    assert r.returncode == 0, r.stderr
    assert "usage:" in r.stdout.lower()
