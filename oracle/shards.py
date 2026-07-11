"""Measured shard planner: how finely to split a WGS run across SLURM array elements.

THE PROBLEM THIS SOLVES: the cluster hard-caps every job at 24h wall-clock (jobs are
also preempted, which makes anything close to the cap risky in practice). A
whole-genome real-VEP run does not fit in that. So work is split -- but splitting
too finely wastes scheduling overhead and array-slot budget, and splitting too
coarsely blows the cap. The operator's rule: **target 12-16h per element, 24h is a
hard ceiling and never a target.** Anything projected over ~16h gets split one
level finer.

THREE LEVELS:
  L0 -- whole WGS, one element per combo.
  L1 -- one element per chromosome.
  L2 -- one element per genomic chunk WITHIN a chromosome (chr1/chr2 carry ~8x
        chr22's variants, so an L1 element that is comfortable on chr22 can still
        blow the band on chr1 -- L2 exists for exactly that chromosome).

THE PLANNER IS PER-CHROMOSOME, NOT ALL-OR-NOTHING: whether a run overflows L0 is a
single whole-genome decision, but once it does, each chromosome is judged on its
OWN projected runtime. A run can legitimately come back with chr22 at L1 and chr1
at L2 in the same plan -- that unevenness is the entire reason L2 exists.

THE gt/annotate vs diff ASYMMETRY: gt and annotate shards are VCFs, merged back
together with `bcftools concat` -- a VCF-level, position-associative operation that
does not care how finely the genome was sliced. diff shards are summary JSON
merged by `summary.merge()`, which is hard-coded to require **chromosome-atomic**
shards (see summary.py's `_check_shard_set` -- it raises on anything sub-
chromosomal, precisely because summing two summaries for the *same* chromosome
double-counts, and there is no principled way to tell a legitimate two-shard
chromosome apart from a resubmitted-and-duplicated one without banning the shape
outright). So `Step.DIFF` must never come out of this planner as L2: if diff
overflows the band even at L1, that is a hard planning error the operator must go
fix (re-measure `sec_per_variant`, or accept a >16h diff job) -- never a silent L2
emission that `merge()` would later reject anyway, on a cluster, hours into a run.

CALIBRATION, NOT GUESSWORK: `sec_per_variant` is a MEASURED quantity (time one real
element, divide by its variant count) supplied by the caller. This module never
estimates it -- see bin/plan_shards.py's docstring and the Task 7 report for the
caveats on what such a measurement does and does not capture.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from enum import StrEnum

#: The target band: split anything projected to run longer than BAND_HI_H, aim for
#: TARGET_H once splitting. BAND_LO_H is not a splitting threshold by itself (a
#: run that comes in comfortably under 12h, e.g. chr22 alone, is NOT re-merged with
#: its neighbour just to fill the band) -- it exists so callers can report how many
#: shards actually landed inside the operator's stated sweet spot.
BAND_LO_H = 12.0
BAND_HI_H = 16.0
#: The cluster's hard wall-clock cap. Never a target -- see module docstring.
CAP_H = 24.0
#: Aim point when chunking: comfortably inside the band (not hugging BAND_HI_H),
#: so that the ceil()/remainder rounding that spreads N variants over K chunks
#: cannot itself push a chunk back over BAND_HI_H.
TARGET_H = 14.0


class Level(StrEnum):
    L0 = "L0"
    L1 = "L1"
    L2 = "L2"


class Step(StrEnum):
    GT = "gt"
    ANNOTATE = "annotate"
    DIFF = "diff"


class WalltimeError(RuntimeError):
    """Raised when a plan cannot honour the wall-clock band/cap or the diff
    step's chromosome-atomic requirement -- refuse rather than emit a plan that
    would silently blow the 24h cap or that `summary.merge()` would reject."""


@dataclass(frozen=True)
class Shard:
    combo: str
    step: Step
    level: Level
    chrom: str          # "all" for L0
    chunk: int | None   # chunk index within the chromosome; None unless L2
    n_variants: int
    est_hours: float


def _hours(n_variants: int, sec_per_variant: float) -> float:
    return n_variants * sec_per_variant / 3600.0


def _chrom_sort_key(chrom: str) -> tuple[int, object]:
    """Natural/numeric order: '2' < '10', not '10' < '2'.

    Numeric chromosome names sort first, by integer value; anything else (X, Y,
    MT, chrUn_*, ...) sorts after, lexicographically. This keeps the shard TSV
    (and any log built from it) in the order an operator actually reads it.
    """
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)


def _chunk_chrom(combo: str, step: Step, chrom: str, n_variants: int,
                  sec_per_variant: float, max_chunks_per_chrom: int) -> list[Shard]:
    """Split one overflowing chromosome into L2 chunks, each sized to land near
    TARGET_H, none ever left over BAND_HI_H.

    Chunk count starts from the TARGET_H-sized estimate, then grows (never
    shrinks) until the LARGEST resulting chunk actually fits the band -- integer
    division means an even split is never bigger than ceil(n / chunks), so this
    converges immediately in practice; the loop only exists to make that a
    checked invariant rather than an assumption.
    """
    if sec_per_variant <= 0:
        raise ValueError(f"sec_per_variant must be > 0, got {sec_per_variant!r}")

    variants_per_chunk = max(1, int(TARGET_H * 3600 / sec_per_variant))
    num_chunks = max(1, math.ceil(n_variants / variants_per_chunk))

    while True:
        if num_chunks > max_chunks_per_chrom:
            raise WalltimeError(
                f"plan_combo(): chromosome {chrom!r} of combo {combo!r} needs more "
                f"than max_chunks_per_chrom={max_chunks_per_chrom} L2 chunks to bring "
                f"every chunk under the {BAND_HI_H}h band (n_variants={n_variants}, "
                f"sec_per_variant={sec_per_variant}). Refusing to emit a chunk that "
                f"would overflow the {CAP_H}h hard cluster cap -- re-measure "
                f"sec_per_variant, or raise max_chunks_per_chrom if the estimate is "
                f"simply coarse."
            )
        base, rem = divmod(n_variants, num_chunks)
        sizes = [base + 1 if i < rem else base for i in range(num_chunks)]
        if _hours(max(sizes), sec_per_variant) <= BAND_HI_H:
            break
        num_chunks += 1

    return [
        Shard(combo, step, Level.L2, chrom, i, sz, _hours(sz, sec_per_variant))
        for i, sz in enumerate(sizes)
        if sz > 0
    ]


def plan_combo(combo: str, step: Step, chrom_counts: dict[str, int],
               sec_per_variant: float, max_chunks_per_chrom: int = 64) -> list[Shard]:
    """Plan the shard set for one combo, one step, over the given per-chromosome
    variant counts.

    Escalates L0 -> L1 -> L2 only as far as needed: the whole WGS run is tried
    first (one shard), then each chromosome independently (a chromosome that
    fits stays L1; one that does not is chunked to L2) -- EXCEPT for
    `Step.DIFF`, which raises `WalltimeError` instead of ever producing L2 (see
    module docstring: `summary.merge()` requires chromosome-atomic shards).
    """
    if not chrom_counts:
        raise ValueError("plan_combo(): chrom_counts must not be empty")

    total = sum(chrom_counts.values())
    whole_hours = _hours(total, sec_per_variant)
    if whole_hours <= BAND_HI_H:
        return [Shard(combo, step, Level.L0, "all", None, total, whole_hours)]

    chroms = sorted(chrom_counts, key=_chrom_sort_key)
    per_chrom_hours = {c: _hours(chrom_counts[c], sec_per_variant) for c in chroms}
    overflowing = [c for c in chroms if per_chrom_hours[c] > BAND_HI_H]

    if step is Step.DIFF:
        if overflowing:
            worst = max(overflowing, key=lambda c: per_chrom_hours[c])
            raise WalltimeError(
                f"plan_combo(): diff for combo {combo!r} projects "
                f"{per_chrom_hours[worst]:.1f}h on chromosome {worst!r} alone -- over "
                f"the {BAND_HI_H}h band even at L1 (whole-genome would be "
                f"{whole_hours:.1f}h). The diff step's shards must be "
                f"chromosome-atomic: summary.merge() requires exactly one summary "
                f"JSON per chromosome and rejects sub-chromosomal shards outright, "
                f"so this planner refuses to emit L2 for diff. Fix the wall-clock "
                f"projection (re-measure sec_per_variant, split gt/annotate finer "
                f"and diff a smaller region set, or accept a >16h diff job) instead."
            )
        return [Shard(combo, step, Level.L1, c, None, chrom_counts[c], per_chrom_hours[c])
                for c in chroms]

    shards: list[Shard] = []
    for c in chroms:
        if per_chrom_hours[c] <= BAND_HI_H:
            shards.append(Shard(combo, step, Level.L1, c, None, chrom_counts[c], per_chrom_hours[c]))
        else:
            shards.extend(_chunk_chrom(combo, step, c, chrom_counts[c],
                                        sec_per_variant, max_chunks_per_chrom))
    return shards


def plan_all(combos: list[str], step: Step, chrom_counts: dict[str, int],
             sec_per_variant: float, max_chunks_per_chrom: int = 64) -> list[Shard]:
    """Plan every combo's shard set for one step, concatenated into one list.

    A belt-and-braces check: even though `plan_combo` is built so no shard it
    emits can exceed BAND_HI_H (and BAND_HI_H < CAP_H by construction), this
    re-checks every shard against the hard CAP_H before returning -- so a future
    change to the per-shard logic that quietly breaks that invariant fails loud
    here rather than shipping a plan that blows the 24h cap on the cluster.
    """
    shards: list[Shard] = []
    for combo in combos:
        shards.extend(plan_combo(combo, step, chrom_counts, sec_per_variant,
                                  max_chunks_per_chrom=max_chunks_per_chrom))

    over_cap = [s for s in shards if s.est_hours > CAP_H]
    if over_cap:
        worst = max(over_cap, key=lambda s: s.est_hours)
        raise WalltimeError(
            f"plan_all(): shard {worst.combo}/{worst.step}/{worst.level}/"
            f"{worst.chrom}#{worst.chunk} projects {worst.est_hours:.1f}h, over the "
            f"hard {CAP_H}h cluster cap. This should be unreachable given "
            f"plan_combo()'s own band checks -- treat it as a bug in the planner, "
            f"not just a bad measurement."
        )
    return shards


SHARD_HEADER = ["combo", "step", "level", "chrom", "chunk", "n_variants", "est_hours"]


def write_shards(path: str, shards: list[Shard]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SHARD_HEADER)
        for s in shards:
            w.writerow([s.combo, s.step.value, s.level.value, s.chrom,
                        "" if s.chunk is None else s.chunk,
                        s.n_variants, repr(s.est_hours)])


def load_shards(path: str) -> list[Shard]:
    with open(path, newline="") as fh:
        return [
            Shard(
                combo=row["combo"],
                step=Step(row["step"]),
                level=Level(row["level"]),
                chrom=row["chrom"],
                chunk=int(row["chunk"]) if row["chunk"] else None,
                n_variants=int(row["n_variants"]),
                est_hours=float(row["est_hours"]),
            )
            for row in csv.DictReader(fh, delimiter="\t")
        ]
