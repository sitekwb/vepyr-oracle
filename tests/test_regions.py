"""oracle/regions.py: the L2 region cutter.

The property under test throughout is the one the module docstring leads with:
cutting at empirical variant-count QUANTILES keeps every chunk's variant count
(and therefore its wall-clock estimate) close to equal, even when the underlying
genomic density is wildly uneven -- which an equal-LENGTH split would not.
"""
from __future__ import annotations

import random

import pytest

from oracle.regions import RegionError, load_regions, quantile_regions, read_positions, write_regions

FIELDS_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##contig=<ID=chr22>\n'
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
)


def _vcf_line(chrom: str, pos: int) -> str:
    return f"{chrom}\t{pos}\t.\tA\tG\t.\t.\t.\n"


# --- read_positions ------------------------------------------------------------

def test_read_positions_reads_in_file_order(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text(FIELDS_HEADER + _vcf_line("chr22", 100) + _vcf_line("chr22", 250)
                 + _vcf_line("chr22", 9000))
    assert read_positions(str(p)) == [100, 250, 9000]


def test_read_positions_filters_by_chrom(tmp_path):
    p = tmp_path / "x.vcf"
    p.write_text(FIELDS_HEADER + _vcf_line("chr21", 5) + _vcf_line("chr22", 100)
                 + _vcf_line("chr21", 6))
    assert read_positions(str(p), chrom="22") == [100]
    assert read_positions(str(p), chrom="chr21") == [5, 6]


# --- quantile_regions: basic correctness ---------------------------------------

def test_even_split_produces_equal_sized_chunks():
    positions = list(range(1, 1001))       # 1000 evenly-spaced variants
    regions = quantile_regions("22", positions, 4)
    assert [r.n_variants for r in regions] == [250, 250, 250, 250]


def test_chunks_abut_exactly_no_gap_no_overlap():
    positions = sorted(random.Random(7).sample(range(1, 2_000_000), 5000))
    regions = quantile_regions("22", positions, 9)
    for prev, cur in zip(regions, regions[1:]):
        assert cur.start == prev.end + 1


def test_self_check_counts_sum_to_the_chromosome_total():
    positions = sorted(random.Random(11).sample(range(1, 500_000), 777))
    regions = quantile_regions("1", positions, 5)
    assert sum(r.n_variants for r in regions) == len(positions)


def test_every_position_falls_inside_its_assigned_region():
    positions = sorted(random.Random(3).sample(range(1, 300_000), 613))
    regions = quantile_regions("2", positions, 6)
    idx = 0
    for r in regions:
        for _ in range(r.n_variants):
            assert r.start <= positions[idx] <= r.end
            idx += 1
    assert idx == len(positions)


def test_first_region_starts_at_the_first_position_last_ends_at_the_last():
    positions = [10, 200, 3000, 40000, 500000]
    regions = quantile_regions("3", positions, 3)
    assert regions[0].start == positions[0]
    assert regions[-1].end == positions[-1]


def test_single_chunk_covers_everything():
    positions = [5, 6, 7, 1000]
    (r,) = quantile_regions("4", positions, 1)
    assert r.n_variants == 4
    assert r.start == 5 and r.end == 1000


# --- quantile_regions: same-position ties ---------------------------------------

def test_a_position_tie_spanning_the_natural_cut_is_never_split_across_regions():
    # 12 positions, split 2-way -> chunk_sizes(12, 2) == [6, 6], so the NATURAL
    # (unsnapped) cut lands exactly at index 6 -- right where positions[5..7]
    # (three records at pos=50, e.g. an SNV + two indels at the same base) tie.
    # An unsnapped cut would slice that tie run in half.
    positions = [1, 2, 3, 4, 5, 50, 50, 50, 60, 70, 80, 90]
    regions = quantile_regions("5", positions, 2)
    assert sum(r.n_variants for r in regions) == len(positions)
    # every position must be inside exactly the region its index falls in --
    # the tie run [50,50,50] must be wholly inside ONE region, not split (a
    # split would put one of the ==50 positions outside its region's bounds,
    # since a boundary cannot sit strictly at pos=50 while the other region
    # also claims pos=50).
    idx = 0
    for r in regions:
        for _ in range(r.n_variants):
            assert r.start <= positions[idx] <= r.end
            idx += 1
    assert len(regions) == 2
    # the cut was pushed past the whole tie run, so chunk 0 grew from the
    # target 6 to 8 -- correctness (never split a tie) over exactness.
    assert regions[0].n_variants == 8
    assert regions[1].n_variants == 4


def test_all_positions_tied_still_yields_a_single_valid_covering_region():
    positions = [42] * 10
    regions = quantile_regions("6", positions, 3)
    assert sum(r.n_variants for r in regions) == 10
    for r in regions:
        assert r.start == 42 and r.end == 42


# --- quantile_regions: refuses bad input -----------------------------------------

def test_raises_on_empty_positions():
    with pytest.raises(RegionError, match="no variant positions"):
        quantile_regions("22", [], 4)


def test_raises_on_unsorted_positions():
    with pytest.raises(RegionError, match="not sorted"):
        quantile_regions("22", [100, 50, 200], 2)


def test_raises_on_num_chunks_below_one():
    with pytest.raises(ValueError, match="num_chunks"):
        quantile_regions("22", [1, 2, 3], 0)


# --- write_regions / load_regions round-trip -------------------------------------

def test_roundtrip(tmp_path):
    positions = list(range(1, 101))
    regions = quantile_regions("22", positions, 4)
    p = tmp_path / "regions_22.tsv"
    write_regions(str(p), regions)
    back = load_regions(str(p))
    assert back == regions


def test_written_region_spec_matches_run_wgs_region_format(tmp_path):
    positions = [10, 20, 30, 40]
    (r,) = quantile_regions("22", positions, 1)
    p = tmp_path / "regions_22.tsv"
    write_regions(str(p), [r])
    rows = p.read_text().splitlines()
    assert rows[0] == "chunk\tregion\tn_variants"
    assert rows[1] == f"0\t22:{r.start}-{r.end}\t4"


# --- THE PROOF: quantile cutting vs equal-length cutting on an uneven chromosome --
#
# Synthetic chromosome: 90% of the variants packed into the first 10% of the
# genomic length, the remaining 10% spread thinly over the other 90% -- modelling
# a gene-dense arm next to a near-empty centromere/telomere stretch. This is the
# exact shape the module docstring names as the failure mode an equal-LENGTH
# split gets badly wrong.

CHROM_LENGTH = 10_000_000
DENSE_FRACTION = 0.10          # first 10% of the length ...
DENSE_VARIANT_SHARE = 0.90     # ... holds 90% of the variants
N_VARIANTS = 9_000
N_CHUNKS = 9


def _synthetic_uneven_positions() -> list[int]:
    rng = random.Random(42)
    dense_end = int(CHROM_LENGTH * DENSE_FRACTION)
    n_dense = int(N_VARIANTS * DENSE_VARIANT_SHARE)
    n_sparse = N_VARIANTS - n_dense
    dense = rng.sample(range(1, dense_end), n_dense)
    sparse = rng.sample(range(dense_end, CHROM_LENGTH), n_sparse)
    return sorted(dense + sparse)


def _equal_length_chunk_counts(positions: list[int], num_chunks: int) -> list[int]:
    """What a LENGTH-based split (chrom_length / num_chunks, the approach this
    module explicitly rejects) would have produced -- built independently of
    oracle.regions, purely to demonstrate the contrast."""
    width = CHROM_LENGTH // num_chunks
    counts = [0] * num_chunks
    for pos in positions:
        i = min(pos // max(width, 1), num_chunks - 1)
        counts[i] += 1
    return counts


def test_equal_length_split_would_have_been_badly_uneven():
    """Sanity check on the synthetic fixture itself: confirms the length-based
    split really is the bad baseline it claims to be, before contrasting it
    with the quantile-based one below."""
    positions = _synthetic_uneven_positions()
    counts = _equal_length_chunk_counts(positions, N_CHUNKS)
    assert sum(counts) == N_VARIANTS
    # the first (dense) length-chunk swallows ~90% of all variants...
    assert counts[0] > 0.85 * N_VARIANTS
    # ...while several of the sparse chunks get next to nothing.
    assert min(counts) < 0.02 * N_VARIANTS
    # the resulting per-chunk wall-clock projections (proportional to count)
    # would then differ by orders of magnitude -- exactly what the 12-16h band
    # cannot tolerate.
    assert max(counts) / max(min(counts), 1) > 50


def test_quantile_split_keeps_chunks_close_to_equal_on_the_same_uneven_chromosome():
    """The module's actual fix, applied to the IDENTICAL synthetic chromosome
    the previous test showed a length-based split mishandling."""
    positions = _synthetic_uneven_positions()
    regions = quantile_regions("uneven", positions, N_CHUNKS)

    sizes = [r.n_variants for r in regions]
    assert sum(sizes) == N_VARIANTS
    target = N_VARIANTS / N_CHUNKS
    # every chunk lands within one variant of perfectly even (chunk_sizes()'s
    # own [base+1]*rem + [base]*(n-rem) guarantee) -- a world away from the
    # length-based split's 50x spread above.
    assert max(sizes) - min(sizes) <= 1
    assert max(sizes) / min(sizes) < 1.1
    for s in sizes:
        assert abs(s - target) <= 1

    # and the self-check the module runs internally holds on this exact input:
    # no gap, no overlap, full coverage.
    for prev, cur in zip(regions, regions[1:]):
        assert cur.start == prev.end + 1
    assert regions[0].start == positions[0]
    assert regions[-1].end == positions[-1]
