"""Genomic region slicing for L2 shards: turn "chromosome X, N chunks" into
concrete `chrom:start-end` regions that abut exactly and hold VARIANT-COUNT-equal
slices of the chromosome -- never equal-LENGTH slices.

THE PROBLEM THIS SOLVES: variant density along a chromosome is wildly uneven
(centromeres are near-empty, gene-dense arms are packed). An equal-LENGTH region
split (chromosome_length / num_chunks) therefore produces wildly UNEQUAL variant
counts per chunk -- and `oracle.shards`'s `est_hours` for an L2 chunk is a
straight `n_variants * sec_per_variant`, so a length-based split silently
invalidates the wall-clock estimate the whole 12-16h plan rests on. A chunk that
happens to land over a gene-dense arm can carry many times the variants its
neighbours do while looking identical on a length-based ruler -- and the operator
would not find out until that one chunk blew the 24h cap.

THE FIX: cut the SORTED variant POSITION STREAM at empirical count quantiles,
using `oracle.shards.chunk_sizes` -- the exact same partition `_chunk_chrom()`
used to decide HOW MANY L2 chunks a chromosome needs. One shared formula means
the region this module emits always realises precisely the `n_variants`
`shards.tsv` already committed an `est_hours` to; two independent
implementations of "the same" partition is exactly how a plan and its execution
drift apart.

NO GAP, NO OVERLAP: chunk i's end is chunk i+1's start - 1, so the emitted
regions partition `[positions[0], positions[-1]]` exactly. A gap silently drops
the variants that fall inside it; an overlap double-counts them when `bcftools
concat` stitches the chunks back into a chromosome. Both are CHECKED before
anything is written, never assumed.

SAME-POSITION TIES: two VCF records can legitimately share one POS (e.g. an SNV
and an indel both starting at the same base, recorded as separate lines). A
region is a genomic *coordinate* range, not a record index, so a cut can never
literally land between two records at the same position -- both would fall
inside `[X, X]` and a `--region` fetch cannot split them. `quantile_regions`
detects this and pushes the cut forward past the whole tied run, so realised
chunk sizes can differ from the exactly-even target by the (small) size of a
tie run -- correctness beats exactness here.

No pysam import at module level: `read_positions` is a plain text/gzip scan
(via `oracle.csq.open_maybe_gzip`), so this module stays laptop-testable like
the rest of `oracle/`.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

from .csq import norm_chrom, open_maybe_gzip
from .shards import chunk_sizes

REGIONS_HEADER = ["chunk", "region", "n_variants"]


class RegionError(RuntimeError):
    """The region cut cannot honour its own contract (sorted input required,
    counts must sum to the chromosome total, chunks must abut with no gap and
    no overlap). Refusing to emit a region file `run_wgs.py`/`gen_gt116.sh`
    would silently drop or double-count variants against."""


@dataclass(frozen=True)
class Region:
    chunk: int
    chrom: str
    start: int
    end: int
    n_variants: int

    @property
    def spec(self) -> str:
        """`"chrom:start-end"` -- what `--region` on `bin/run_wgs.py` expects."""
        return f"{self.chrom}:{self.start}-{self.end}"


def read_positions(vcf_path: str, chrom: str | None = None) -> list[int]:
    """Every variant POSITION in `vcf_path`, in file order.

    The file must already be position-sorted (a per-chromosome slice produced by
    `slurm/prep_input.sh` is, by construction) -- `quantile_regions` checks this
    explicitly rather than trusting it silently.

    `chrom`, if given, restricts to records on that contig (after
    `oracle.csq.norm_chrom` normalisation) -- lets this be pointed at a
    whole-genome file for testing without needing a pre-split per-chromosome one.
    """
    cf = norm_chrom(chrom) if chrom else None
    positions: list[int] = []
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t", 2)
            if len(cols) < 2:
                continue
            if cf is not None and norm_chrom(cols[0]) != cf:
                continue
            positions.append(int(cols[1]))
    return positions


def _monotonic_clip(cuts: list[int], n: int) -> list[int]:
    """Force `cuts` non-decreasing and within `[0, n]`.

    Snapping a cut forward past a same-position tie run (see module docstring)
    can, in principle, push it past where the NEXT cut would naturally have
    landed (an extreme, densely-tied run). Rather than assume that never
    happens, clamp explicitly: a later cut is never allowed to sit before an
    earlier one, which only ever COLLAPSES a chunk to empty (silently dropped
    by the caller, exactly like `_chunk_chrom`'s own `if sz > 0` filter) -- it
    can never produce a gap or an overlap.
    """
    out = []
    prev = 0
    for c in cuts:
        c = max(prev, min(c, n))
        out.append(c)
        prev = c
    return out


def quantile_regions(chrom: str, positions: list[int], num_chunks: int) -> list[Region]:
    """Cut SORTED `positions` into `num_chunks` count-equal, gap/overlap-free regions.

    Uses `oracle.shards.chunk_sizes(len(positions), num_chunks)` for the TARGET
    per-chunk counts (see module docstring for why this must be the same formula
    the shard planner used), then snaps each cut point away from any run of
    variants sharing one POSITION so a region boundary never needs to split them.

    A boundary between chunk i and chunk i+1 is placed at the integer midpoint
    between chunk i's last position and chunk i+1's first position -- both are
    guaranteed strictly increasing by the tie-snap above, so the midpoint always
    satisfies `chunk_i.end >= last position in chunk i` and
    `chunk_{i+1}.start <= first position in chunk i+1`: every real variant stays
    inside the region it was assigned to.
    """
    if num_chunks < 1:
        raise ValueError(f"quantile_regions(): num_chunks must be >= 1, got {num_chunks!r}")
    if not positions:
        raise RegionError(
            f"quantile_regions(): chromosome {chrom!r} has no variant positions to "
            f"split -- refusing to emit a region set for an empty chromosome."
        )
    for a, b in zip(positions, positions[1:]):
        if b < a:
            raise RegionError(
                f"quantile_regions(): positions for chromosome {chrom!r} are not "
                f"sorted ascending ({a} then {b}). The region cut assumes a "
                f"position-sorted per-chromosome input (see slurm/prep_input.sh); "
                f"an unsorted stream would silently produce overlapping regions."
            )

    n = len(positions)
    target_sizes = chunk_sizes(n, num_chunks)
    cuts: list[int] = []
    acc = 0
    for sz in target_sizes[:-1]:
        acc += sz
        # never let a cut fall strictly between two records at the SAME position
        # -- a region is a coordinate range and cannot split a tie (see docstring)
        while 0 < acc < n and positions[acc - 1] == positions[acc]:
            acc += 1
        cuts.append(acc)
    cuts = _monotonic_clip(cuts, n)

    bounds = [0, *cuts, n]
    regions: list[Region] = []
    for lo, hi in zip(bounds, bounds[1:]):
        if lo >= hi:
            continue          # a tie-snap collapsed this chunk to empty -- drop it
        chunk_positions = positions[lo:hi]
        start = chunk_positions[0] if not regions else regions[-1].end + 1
        if hi == n:
            end = chunk_positions[-1]
        else:
            end = (chunk_positions[-1] + positions[hi]) // 2
        regions.append(Region(chunk=len(regions), chrom=chrom, start=start, end=end,
                              n_variants=hi - lo))

    _check_partition(chrom, regions, positions)
    return regions


def _check_partition(chrom: str, regions: list[Region], positions: list[int]) -> None:
    """NO GAP, NO OVERLAP, and the counts sum to the chromosome total -- checked,
    never assumed. A silent gap drops variants; a silent overlap double-counts
    them in `bcftools concat`; a count mismatch invalidates every `est_hours` the
    plan committed to."""
    total = sum(r.n_variants for r in regions)
    if total != len(positions):
        raise RegionError(
            f"quantile_regions(): chromosome {chrom!r} chunk counts sum to {total}, "
            f"not the {len(positions)} variants in the input. Refusing to emit a "
            f"region set that would silently drop or fabricate variants."
        )
    for prev, cur in zip(regions, regions[1:]):
        if cur.start != prev.end + 1:
            raise RegionError(
                f"quantile_regions(): chromosome {chrom!r} chunk {prev.chunk} ends at "
                f"{prev.end} but chunk {cur.chunk} starts at {cur.start} (expected "
                f"{prev.end + 1}). A gap or overlap here would drop or double-count "
                f"variants once the chunks are concatenated back into a chromosome."
            )
    for r in regions:
        if r.start > r.end:
            raise RegionError(
                f"quantile_regions(): chromosome {chrom!r} chunk {r.chunk} has "
                f"start {r.start} > end {r.end} -- refusing to emit an inverted region."
            )


def write_regions(path: str, regions: list[Region]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(REGIONS_HEADER)
        for r in regions:
            w.writerow([r.chunk, r.spec, r.n_variants])


def load_regions(path: str) -> list[Region]:
    with open(path, newline="") as fh:
        out = []
        for row in csv.DictReader(fh, delimiter="\t"):
            chrom, span = row["region"].split(":", 1)
            start_s, end_s = span.split("-", 1)
            out.append(Region(chunk=int(row["chunk"]), chrom=chrom, start=int(start_s),
                              end=int(end_s), n_variants=int(row["n_variants"])))
        return out
