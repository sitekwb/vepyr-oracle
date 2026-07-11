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
estimates it. Two things follow from the fact that it is a MODEL and not a law:

  * SAFETY MARGIN. `sec_per_variant` is a proxy: VEP's true per-variant cost tracks
    the number of overlapping TRANSCRIPTS and the indel-vs-SNV mix, which a single
    averaged rate smooths over -- a chunk that happens to be transcript-dense or
    indel-rich runs slower than the model says. The downside is ASYMMETRIC:
    under-provisioning costs a job KILLED at the 24h cap after burning up to 24h of
    compute (repeatedly, on a preemption-prone cluster), while over-provisioning
    costs a few more array elements, which run concurrently anyway. So the planner
    is deliberately biased conservative: every decision is made against
    `sec_per_variant * safety_factor` (default 1.4), never the raw rate. A shard
    planned right at the 16h band ceiling would have to overrun the model by >50%
    before it reached the 24h cap.

    `Shard.est_hours` reports the SAFETY-ADJUSTED estimate -- that is the number the
    operator should reason about -- while `Shard.est_hours_raw` carries the
    unadjusted model output, so the post-hoc calibration loop (compare planned vs.
    ACTUAL wall time; re-measure if a combo systematically diverges) can compare
    like with like instead of against a number the margin already inflated.

  * PER-COMBO RATES. One global rate does not transfer across combos:
    `hgvs_merged_am` runs the AlphaMissense plugin and `hgvs_refseq` does not, so
    they do not cost the same per variant. `plan_all` therefore accepts either a
    single float (measured once, applied to every combo) or a `dict[combo -> rate]`.
    If a dict is given and a combo is MISSING from it, that RAISES -- silently
    falling back to a global average is exactly how a slow combo blows the cap
    unnoticed.
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

#: Conservative bias applied to the MEASURED rate before any planning decision --
#: see the module docstring. 1.4 means a shard planned at the BAND_HI_H ceiling
#: (16h) has to overrun the model by more than 50% before it reaches CAP_H (24h),
#: and one planned at TARGET_H (14h) by more than 70%. The cost of being wrong in
#: this direction is a few extra concurrent array elements; the cost of being wrong
#: in the other direction is a killed job and a lost day.
DEFAULT_SAFETY_FACTOR = 1.4


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
    chrom: str            # "all" for L0
    chunk: int | None     # chunk index within the chromosome; None unless L2
    n_variants: int
    est_hours: float      # SAFETY-ADJUSTED -- what the planner decided against
    est_hours_raw: float  # the unadjusted model -- for planned-vs-ACTUAL calibration
    safety_factor: float  # the margin this plan assumed, carried so it is self-documenting


def _hours(n_variants: int, sec_per_variant: float) -> float:
    return n_variants * sec_per_variant / 3600.0


def _check_rate(sec_per_variant: float, safety_factor: float) -> float:
    """Validate the measured rate + margin, and return the EFFECTIVE rate to plan against."""
    if sec_per_variant <= 0:
        raise ValueError(f"sec_per_variant must be > 0, got {sec_per_variant!r}")
    if safety_factor < 1.0:
        raise ValueError(
            f"safety_factor must be >= 1.0 (it is a conservative margin on a proxy "
            f"model, never a discount), got {safety_factor!r}"
        )
    return sec_per_variant * safety_factor


def _chrom_sort_key(chrom: str) -> tuple[int, object]:
    """Natural/numeric order: '2' < '10', not '10' < '2'.

    Numeric chromosome names sort first, by integer value; anything else (X, Y,
    MT, chrUn_*, ...) sorts after, lexicographically. This keeps the shard TSV
    (and any log built from it) in the order an operator actually reads it.
    """
    return (0, int(chrom)) if chrom.isdigit() else (1, chrom)


def _shard(combo: str, step: Step, level: Level, chrom: str, chunk: int | None,
           n_variants: int, sec_per_variant: float, safety_factor: float) -> Shard:
    """Build a Shard, carrying BOTH the adjusted estimate the planner decided on and
    the raw model output the calibration loop will compare actual wall time against."""
    return Shard(
        combo=combo, step=step, level=level, chrom=chrom, chunk=chunk,
        n_variants=n_variants,
        est_hours=_hours(n_variants, sec_per_variant * safety_factor),
        est_hours_raw=_hours(n_variants, sec_per_variant),
        safety_factor=safety_factor,
    )


def _chunk_chrom(combo: str, step: Step, chrom: str, n_variants: int,
                 sec_per_variant: float, safety_factor: float,
                 max_chunks_per_chrom: int) -> list[Shard]:
    """Split one overflowing chromosome into L2 chunks, each sized to land near
    TARGET_H *on the safety-adjusted rate*, none ever left over BAND_HI_H.

    Chunk count starts from the TARGET_H-sized estimate, then grows (never
    shrinks) until the LARGEST resulting chunk actually fits the band -- integer
    division means an even split is never bigger than ceil(n / chunks), so this
    converges immediately in practice; the loop only exists to make that a
    checked invariant rather than an assumption.
    """
    eff = _check_rate(sec_per_variant, safety_factor)

    variants_per_chunk = max(1, int(TARGET_H * 3600 / eff))
    num_chunks = max(1, math.ceil(n_variants / variants_per_chunk))

    while True:
        if num_chunks > max_chunks_per_chrom:
            raise WalltimeError(
                f"plan_combo(): chromosome {chrom!r} of combo {combo!r} needs more "
                f"than max_chunks_per_chrom={max_chunks_per_chrom} L2 chunks to bring "
                f"every chunk under the {BAND_HI_H}h band (n_variants={n_variants}, "
                f"sec_per_variant={sec_per_variant}, safety_factor={safety_factor} "
                f"-> effective {eff}s/variant). Refusing to emit a chunk that would "
                f"overflow the {CAP_H}h hard cluster cap -- re-measure "
                f"sec_per_variant, or raise max_chunks_per_chrom if the estimate is "
                f"simply coarse."
            )
        base, rem = divmod(n_variants, num_chunks)
        sizes = [base + 1 if i < rem else base for i in range(num_chunks)]
        if _hours(max(sizes), eff) <= BAND_HI_H:
            break
        num_chunks += 1

    return [
        _shard(combo, step, Level.L2, chrom, i, sz, sec_per_variant, safety_factor)
        for i, sz in enumerate(sizes)
        if sz > 0
    ]


def plan_combo(combo: str, step: Step, chrom_counts: dict[str, int],
               sec_per_variant: float, max_chunks_per_chrom: int = 64,
               safety_factor: float = DEFAULT_SAFETY_FACTOR) -> list[Shard]:
    """Plan the shard set for one combo, one step, over the given per-chromosome
    variant counts.

    Every wall-clock decision is made against the SAFETY-ADJUSTED rate
    (`sec_per_variant * safety_factor`), never the raw measurement -- see the module
    docstring on why the planner is deliberately biased conservative.

    Escalates L0 -> L1 -> L2 only as far as needed: the whole WGS run is tried first
    (one shard), then each chromosome independently (a chromosome that fits stays L1;
    one that does not is chunked to L2) -- EXCEPT for `Step.DIFF`, which raises
    `WalltimeError` instead of ever producing L2 (see module docstring:
    `summary.merge()` requires chromosome-atomic shards).
    """
    if not chrom_counts:
        raise ValueError("plan_combo(): chrom_counts must not be empty")
    eff = _check_rate(sec_per_variant, safety_factor)

    def mk(level: Level, chrom: str, chunk: int | None, n: int) -> Shard:
        return _shard(combo, step, level, chrom, chunk, n, sec_per_variant, safety_factor)

    total = sum(chrom_counts.values())
    whole_hours = _hours(total, eff)
    if whole_hours <= BAND_HI_H:
        return [mk(Level.L0, "all", None, total)]

    chroms = sorted(chrom_counts, key=_chrom_sort_key)
    per_chrom_hours = {c: _hours(chrom_counts[c], eff) for c in chroms}
    overflowing = [c for c in chroms if per_chrom_hours[c] > BAND_HI_H]

    if step is Step.DIFF:
        if overflowing:
            worst = max(overflowing, key=lambda c: per_chrom_hours[c])
            raise WalltimeError(
                f"plan_combo(): diff for combo {combo!r} projects "
                f"{per_chrom_hours[worst]:.1f}h on chromosome {worst!r} alone "
                f"(sec_per_variant={sec_per_variant} x safety_factor={safety_factor}) "
                f"-- over the {BAND_HI_H}h band even at L1 (whole-genome would be "
                f"{whole_hours:.1f}h). The diff step's shards must be "
                f"chromosome-atomic: summary.merge() requires exactly one summary "
                f"JSON per chromosome and rejects sub-chromosomal shards outright, "
                f"so this planner refuses to emit L2 for diff. Fix the wall-clock "
                f"projection (re-measure sec_per_variant, split gt/annotate finer "
                f"and diff a smaller region set, or accept a >16h diff job) instead."
            )
        return [mk(Level.L1, c, None, chrom_counts[c]) for c in chroms]

    shards: list[Shard] = []
    for c in chroms:
        if per_chrom_hours[c] <= BAND_HI_H:
            shards.append(mk(Level.L1, c, None, chrom_counts[c]))
        else:
            shards.extend(_chunk_chrom(combo, step, c, chrom_counts[c], sec_per_variant,
                                       safety_factor, max_chunks_per_chrom))
    return shards


def rate_for(sec_per_variant: float | dict[str, float], combo: str) -> float:
    """The measured rate for one combo: a single float applies to all; a dict must
    CONTAIN the combo.

    A missing entry raises rather than falling back to some global average, because
    the fallback is exactly the failure this indirection exists to prevent: combos do
    not cost the same per variant (`hgvs_merged_am` runs the AlphaMissense plugin,
    `hgvs_refseq` does not), so an unmeasured combo silently inheriting a cheap
    combo's rate is how a slow combo gets under-provisioned and killed at the cap.
    """
    if isinstance(sec_per_variant, dict):
        if combo not in sec_per_variant:
            raise ValueError(
                f"combo {combo!r} has no measured sec_per_variant (per-combo rates were "
                f"given for {sorted(sec_per_variant)}). Refusing to fall back to a "
                f"global average: combos do not cost the same per variant, and an "
                f"unmeasured combo inheriting a cheaper combo's rate is precisely how a "
                f"slow combo gets under-provisioned and killed at the {CAP_H}h cap. "
                f"Measure it, or pass a single float to apply one rate to every combo."
            )
        return sec_per_variant[combo]
    return sec_per_variant


def plan_all(combos: list[str], step: Step, chrom_counts: dict[str, int],
             sec_per_variant: float | dict[str, float], max_chunks_per_chrom: int = 64,
             safety_factor: float = DEFAULT_SAFETY_FACTOR) -> list[Shard]:
    """Plan every combo's shard set for one step, concatenated into one list.

    `sec_per_variant` is either a single MEASURED float applied to every combo, or a
    `dict[combo -> measured rate]`. A dict missing a combo raises (see `rate_for`).

    A belt-and-braces check: even though `plan_combo` is built so no shard it emits
    can exceed BAND_HI_H (and BAND_HI_H < CAP_H by construction), this re-checks every
    shard against the hard CAP_H before returning -- so a future change to the
    per-shard logic that quietly breaks that invariant fails loud here rather than
    shipping a plan that blows the 24h cap on the cluster.
    """
    shards: list[Shard] = []
    for combo in combos:
        shards.extend(plan_combo(combo, step, chrom_counts,
                                 rate_for(sec_per_variant, combo),
                                 max_chunks_per_chrom=max_chunks_per_chrom,
                                 safety_factor=safety_factor))

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


#: `est_hours_raw` and `safety_factor` are on every row on purpose: the plan must be
#: self-documenting about the margin it assumed (a shards.tsv read six weeks later
#: must not need the command line that produced it to be interpretable), and the
#: margin is PER-COMBO-rate-independent metadata that the calibration loop needs
#: alongside the raw estimate to compare planned vs. actual honestly.
SHARD_HEADER = ["combo", "step", "level", "chrom", "chunk", "n_variants",
                "est_hours", "est_hours_raw", "safety_factor"]


def write_shards(path: str, shards: list[Shard]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SHARD_HEADER)
        for s in shards:
            w.writerow([s.combo, s.step.value, s.level.value, s.chrom,
                        "" if s.chunk is None else s.chunk, s.n_variants,
                        repr(s.est_hours), repr(s.est_hours_raw), repr(s.safety_factor)])


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
                est_hours_raw=float(row["est_hours_raw"]),
                safety_factor=float(row["safety_factor"]),
            )
            for row in csv.DictReader(fh, delimiter="\t")
        ]
