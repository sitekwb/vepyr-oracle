#!/usr/bin/env python3
"""Plan the SLURM shard set for one pipeline step, across all 8 in-scope combos.

Usage:
  plan_shards.py <step: gt|annotate|diff> <rate> <chrom_counts.tsv> <out shards.tsv>
                 [--safety-factor F]

`<rate>` is EITHER a single measured `sec_per_variant` (a float, applied to every
combo) OR a path to a per-combo rates TSV with `combo` + `sec_per_variant` columns.
Prefer the TSV: one global rate does not transfer across combos -- `hgvs_merged_am`
runs the AlphaMissense plugin and `hgvs_refseq` does not, so they do not cost the
same per variant. A rates TSV that is MISSING a combo is a hard error, never a
silent fallback to some average (see oracle/shards.py: rate_for()).

Either way the rate MUST be MEASURED, never guessed: time one real element on the
cluster (one chromosome or one chunk of one combo, for the step being planned) and
divide by the variant count it processed.

The planner then applies `--safety-factor` (default 1.4) to that measured rate
BEFORE deciding anything, because sec_per_variant is a proxy -- VEP's true cost
tracks transcript density and the indel/SNV mix -- and the downside is asymmetric: an
under-provisioned job is KILLED at the 24h cap after burning a day of compute, while
an over-provisioned one just means a few more array elements, which run concurrently
anyway. `est_hours` in the output is the safety-adjusted number; `est_hours_raw`
carries the unadjusted model so planned-vs-ACTUAL calibration stays honest.

`chrom_counts.tsv` needs a header row with (at least) `chrom` and `n_variants`
columns, tab-separated -- e.g. derived from `bcftools index --stats`.

Exits non-zero if any chromosome cannot be brought under the wall-clock band/cap:
for `diff`, that means a chromosome projects over the band even at L1 (diff cannot
go L2 -- summary.merge() requires chromosome-atomic shards); for `gt`/`annotate`, it
means a chromosome would need more L2 chunks than max_chunks_per_chrom allows.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import COMBOS
from oracle.shards import (BAND_HI_H, BAND_LO_H, DEFAULT_SAFETY_FACTOR, Step, WalltimeError,
                           plan_all, write_shards)


def _read_tsv(path: str, required: set[str]) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no rows (need a header + at least one data row)")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path}: missing column(s) {sorted(missing)} (got {sorted(rows[0])})")
    return rows


def _read_chrom_counts(path: str) -> dict[str, int]:
    rows = _read_tsv(path, {"chrom", "n_variants"})
    return {r["chrom"]: int(r["n_variants"]) for r in rows}


def _read_rate(arg: str) -> float | dict[str, float]:
    """A bare float, or a path to a per-combo rates TSV (`combo`, `sec_per_variant`)."""
    try:
        return float(arg)
    except ValueError:
        pass
    if not os.path.exists(arg):
        raise ValueError(
            f"rate {arg!r} is neither a number nor an existing per-combo rates TSV "
            f"(needs `combo` + `sec_per_variant` columns)"
        )
    rows = _read_tsv(arg, {"combo", "sec_per_variant"})
    return {r["combo"]: float(r["sec_per_variant"]) for r in rows}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="plan_shards.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("step", choices=[s.value for s in Step])
    ap.add_argument("rate", help="measured sec_per_variant (float), OR a per-combo rates TSV")
    ap.add_argument("chrom_counts", help="TSV with `chrom` + `n_variants` columns")
    ap.add_argument("out", help="output shards.tsv")
    ap.add_argument("--safety-factor", type=float, default=DEFAULT_SAFETY_FACTOR,
                    help=f"conservative margin applied to the measured rate before "
                         f"planning (default {DEFAULT_SAFETY_FACTOR}; must be >= 1.0)")
    args = ap.parse_args(argv)

    step = Step(args.step)
    try:
        rate = _read_rate(args.rate)
        chrom_counts = _read_chrom_counts(args.chrom_counts)
    except (OSError, ValueError) as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        return 2

    combos = [c.name for c in COMBOS]

    try:
        shards = plan_all(combos, step, chrom_counts, rate, safety_factor=args.safety_factor)
    except WalltimeError as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        return 1
    except ValueError as e:            # unmeasured combo, safety_factor < 1.0, ...
        print(f"[FATAL] {e}", file=sys.stderr)
        return 2

    write_shards(args.out, shards)

    max_h = max(s.est_hours for s in shards)
    max_raw_h = max(s.est_hours_raw for s in shards)
    in_band = sum(1 for s in shards if BAND_LO_H <= s.est_hours <= BAND_HI_H)
    by_level: dict[str, int] = {}
    for s in shards:
        by_level[s.level.value] = by_level.get(s.level.value, 0) + 1

    rate_desc = (f"per-combo rates for {len(rate)} combos ({args.rate})"
                 if isinstance(rate, dict) else f"sec_per_variant={rate}")

    print(f"[plan_shards] step={step.value} combos={len(combos)} "
          f"chroms={len(chrom_counts)} {rate_desc} safety_factor={args.safety_factor}")
    print(f"[plan_shards] wrote {len(shards)} shards -> {args.out} "
          f"({', '.join(f'{lvl}={n}' for lvl, n in sorted(by_level.items()))})")
    print(f"[plan_shards] max est_hours={max_h:.2f}h "
          f"(raw model {max_raw_h:.2f}h, before the {args.safety_factor}x margin)")
    print(f"[plan_shards] {in_band}/{len(shards)} shards inside the "
          f"{BAND_LO_H:.0f}-{BAND_HI_H:.0f}h target band")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
