"""Smoke tests for bin/verify_input.py -- GATE 2's CLI wrapper.

The key-set comparison logic itself is exhaustively tested in
tests/test_verify_input.py against oracle.verify_input directly -- these
confirm the CLI wiring: exit codes, the printed report, and that a missing
file fails loudly instead of silently comparing against nothing.
"""
from __future__ import annotations

import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VERIFY_INPUT = os.path.join(ROOT, "bin", "verify_input.py")

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _line(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n"


def _run(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, VERIFY_INPUT, *args],
                          capture_output=True, text=True, timeout=30)


def test_identical_key_sets_exit_0(tmp_path):
    ours = tmp_path / "ours.vcf"
    gt = tmp_path / "gt.vcf"
    ours.write_text(HEADER + _line("22", 100, "A", "G") + _line("22", 200, "C", "T"))
    gt.write_text(HEADER + _line("22", 100, "A", "G") + _line("22", 200, "C", "T"))

    r = _run(["--ours", str(ours), "--gt", str(gt), "--chrom", "22"])
    assert r.returncode == 0, r.stderr
    assert "GATE 2 PASSED" in r.stdout


def test_differing_key_sets_exit_1_with_examples(tmp_path):
    ours = tmp_path / "ours.vcf"
    gt = tmp_path / "gt.vcf"
    ours.write_text(HEADER + _line("22", 100, "C", "T,CCGC"))
    gt.write_text(HEADER + _line("22", 100, "C", "T") + _line("22", 100, "C", "CCGC"))

    r = _run(["--ours", str(ours), "--gt", str(gt), "--chrom", "22", "--examples", "5"])
    assert r.returncode == 1
    assert "GATE 2 FAILED" in r.stderr
    assert "only in --ours" in r.stderr
    assert "only in --gt" in r.stderr
    assert "22:100 C>T,CCGC" in r.stderr


def test_chr_prefix_mismatch_between_files_still_compares_correctly(tmp_path):
    ours = tmp_path / "ours.vcf"
    gt = tmp_path / "gt.vcf"
    ours.write_text(HEADER + _line("chr22", 100, "A", "G"))
    gt.write_text(HEADER + _line("22", 100, "A", "G"))

    r = _run(["--ours", str(ours), "--gt", str(gt), "--chrom", "22"])
    assert r.returncode == 0, r.stderr
    assert "GATE 2 PASSED" in r.stdout


def test_missing_ours_file_fails_loudly(tmp_path):
    gt = tmp_path / "gt.vcf"
    gt.write_text(HEADER)
    r = _run(["--ours", str(tmp_path / "nope.vcf"), "--gt", str(gt), "--chrom", "22"])
    assert r.returncode == 2
    assert "--ours" in r.stderr and "does not exist" in r.stderr


def test_missing_gt_file_fails_loudly(tmp_path):
    ours = tmp_path / "ours.vcf"
    ours.write_text(HEADER)
    r = _run(["--ours", str(ours), "--gt", str(tmp_path / "nope.vcf"), "--chrom", "22"])
    assert r.returncode == 2
    assert "--gt" in r.stderr and "does not exist" in r.stderr


def test_help_does_not_require_pysam():
    """oracle/verify_input.py reads VCFs with plain text/gzip, not pysam --
    verify_input.py must therefore stay runnable (at least --help) with no
    pysam installed, same laptop-testability property bin/make_regions.py
    and bin/run_wgs.py pin for their own dependencies."""
    r = subprocess.run([sys.executable, VERIFY_INPUT, "--help"],
                       capture_output=True, text=True, timeout=10)
    assert r.returncode == 0, r.stderr
    assert "usage:" in r.stdout.lower()
