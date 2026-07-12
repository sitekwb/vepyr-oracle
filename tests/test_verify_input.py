"""oracle/verify_input.py: GATE 2 -- key-set identity, the check that actually
proves the (chrom,pos,ref,alt) join will work (a record-count match alone does
not -- see slurm/prep_input.sh's GATE 1/GATE 2 split and the module docstring).
"""
from __future__ import annotations

from oracle.verify_input import compare_key_sets, format_examples, read_keys

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _write(path, lines: list[str]) -> None:
    path.write_text(HEADER + "".join(lines))


def _line(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n"


# --- read_keys -------------------------------------------------------------

def test_read_keys_reads_matching_chrom_only(tmp_path):
    p = tmp_path / "x.vcf"
    _write(p, [_line("21", 100, "A", "G"), _line("22", 200, "C", "T"),
              _line("22", 300, "G", "A")])
    assert read_keys(str(p), "22") == {("22", 200, "C", "T"), ("22", 300, "G", "A")}


def test_read_keys_normalizes_chr_prefix_on_both_sides(tmp_path):
    p = tmp_path / "x.vcf"
    _write(p, [_line("chr22", 200, "C", "T")])
    assert read_keys(str(p), "22") == {("22", 200, "C", "T")}
    assert read_keys(str(p), "chr22") == {("22", 200, "C", "T")}


def test_read_keys_ignores_header_lines(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text(HEADER + "##extra header line, not data\n" + _line("22", 1, "A", "G"))
    assert read_keys(str(p), "22") == {("22", 1, "A", "G")}


def test_read_keys_keeps_a_joint_multiallelic_alt_literal(tmp_path):
    """The whole point: an un-split ALT is ONE key, never exploded on comma --
    see the module docstring on why exploding would hide the exact bug this
    module exists to catch."""
    p = tmp_path / "x.vcf"
    _write(p, [_line("22", 100, "C", "T,CCGC")])
    assert read_keys(str(p), "22") == {("22", 100, "C", "T,CCGC")}


def test_read_keys_on_an_empty_chromosome_returns_empty_set(tmp_path):
    p = tmp_path / "x.vcf"
    _write(p, [_line("21", 1, "A", "G")])
    assert read_keys(str(p), "22") == set()


# --- compare_key_sets: pure set arithmetic ----------------------------------

def test_identical_sets_are_identical():
    ours = {("22", 100, "A", "G"), ("22", 200, "C", "T")}
    comparison = compare_key_sets(ours, set(ours), chrom="22")
    assert comparison.identical
    assert comparison.only_ours == frozenset()
    assert comparison.only_gt == frozenset()
    assert comparison.overlap == 2


def test_split_vs_joint_multiallelic_difference_is_detected():
    """OUR candidate still has the joint row; the GT (built from the already
    split input) has the two split rows. Neither key matches the other --
    exactly the CRITICAL bug this whole gate exists to catch."""
    ours = {("22", 100, "C", "T,CCGC")}
    gt = {("22", 100, "C", "T"), ("22", 100, "C", "CCGC")}
    comparison = compare_key_sets(ours, gt, chrom="22")
    assert not comparison.identical
    assert comparison.only_ours == frozenset({("22", 100, "C", "T,CCGC")})
    assert comparison.only_gt == frozenset({("22", 100, "C", "T"), ("22", 100, "C", "CCGC")})
    assert comparison.overlap == 0


def test_left_alignment_difference_is_detected():
    """Same biological indel, two different (pos,ref,alt) representations --
    an un-left-aligned vs. a left-aligned encoding. Caught as a plain literal
    set difference; no indel-aware special-casing needed."""
    ours = {("22", 1001, "CA", "CAA")}          # not left-aligned
    gt = {("22", 1000, "C", "CA")}               # left-aligned
    comparison = compare_key_sets(ours, gt, chrom="22")
    assert not comparison.identical
    assert comparison.only_ours == frozenset({("22", 1001, "CA", "CAA")})
    assert comparison.only_gt == frozenset({("22", 1000, "C", "CA")})


def test_overlap_counts_shared_keys_only():
    ours = {("22", 1, "A", "G"), ("22", 2, "A", "G"), ("22", 3, "A", "G")}
    gt = {("22", 2, "A", "G"), ("22", 3, "A", "G"), ("22", 4, "A", "G")}
    comparison = compare_key_sets(ours, gt, chrom="22")
    assert comparison.overlap == 2
    assert comparison.only_ours == frozenset({("22", 1, "A", "G")})
    assert comparison.only_gt == frozenset({("22", 4, "A", "G")})
    assert comparison.ours_total == 3
    assert comparison.gt_total == 3


def test_empty_vs_empty_is_identical():
    comparison = compare_key_sets(set(), set(), chrom="22")
    assert comparison.identical
    assert comparison.overlap == 0


# --- format_examples ---------------------------------------------------------

def test_format_examples_is_sorted_and_capped():
    keys = frozenset({("22", 300, "A", "G"), ("22", 100, "A", "G"), ("22", 200, "A", "G")})
    assert format_examples(keys, limit=2) == ["22:100 A>G", "22:200 A>G"]


def test_format_examples_renders_a_still_joint_alt_readably():
    keys = frozenset({("22", 100, "C", "T,CCGC")})
    assert format_examples(keys) == ["22:100 C>T,CCGC"]
