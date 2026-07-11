#!/usr/bin/env python3
"""Plan the SLURM shard set for one pipeline step, across all 8 in-scope combos.

Usage: plan_shards.py <step: gt|annotate|diff> <sec_per_variant> <chrom_counts.tsv> <out shards.tsv>

`sec_per_variant` MUST be a MEASURED quantity, never a guess: time one real
element on the cluster (one chromosome or one chunk of one combo, for the step
being planned) and divide by the variant count it processed. This CLI does not
compute it for you and does not sanity-check it against anything -- see the
Task 7 report for why `sec_per_variant` alone is a coarse model of real-VEP
runtime (transcript density and indel-vs-SNV mix both move the true cost) and
what safety margin / mid-run check compensates for that.

`chrom_counts.tsv` needs a header row with (at least) `chrom` and `n_variants`
columns, tab-separated -- e.g. the output of `bcftools stats` per contig.

Exits non-zero (with the WalltimeError's own message, unmodified) if any
chromosome cannot be brought under the wall-clock band/cap: for `diff`, that
means a chromosome projects over the band even at L1 (diff cannot go L2 -- see
oracle/shards.py); for `gt`/`annotate`, it means a chromosome would need more
L2 chunks than the planner's max_chunks_per_chrom limit.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import COMBOS
from oracle.shards import BAND_HI_H, BAND_LO_H, Step, WalltimeError, plan_all, write_shards


def _read_chrom_counts(path: str) -> dict[str, int]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no rows (need a header + at least one chrom)")
    missing = {"chrom", "n_variants"} - set(rows[0])
    if missing:
        raise ValueError(f"{path}: missing column(s) {sorted(missing)} "
                          f"(got {sorted(rows[0])})")
    return {r["chrom"]: int(r["n_variants"]) for r in rows}


def main(argv: list[str]) -> int:
    if len(argv) != 5:
        print(f"Usage: {argv[0]} <step: gt|annotate|diff> <sec_per_variant> "
              f"<chrom_counts.tsv> <out shards.tsv>", file=sys.stderr)
        return 2
    _, step_arg, sec_arg, counts_path, out_path = argv

    try:
        step = Step(step_arg)
    except ValueError:
        print(f"[FATAL] unknown step {step_arg!r}; must be one of "
              f"{[s.value for s in Step]}", file=sys.stderr)
        return 2

    try:
        sec_per_variant = float(sec_arg)
    except ValueError:
        print(f"[FATAL] sec_per_variant must be a number, got {sec_arg!r}", file=sys.stderr)
        return 2

    try:
        chrom_counts = _read_chrom_counts(counts_path)
    except (OSError, ValueError) as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        return 2

    combos = [c.name for c in COMBOS]

    try:
        shards = plan_all(combos, step, chrom_counts, sec_per_variant)
    except WalltimeError as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        return 1

    write_shards(out_path, shards)

    max_h = max(s.est_hours for s in shards)
    in_band = sum(1 for s in shards if BAND_LO_H <= s.est_hours <= BAND_HI_H)
    by_level: dict[str, int] = {}
    for s in shards:
        by_level[s.level.value] = by_level.get(s.level.value, 0) + 1

    print(f"[plan_shards] step={step.value} combos={len(combos)} "
          f"chroms={len(chrom_counts)} sec_per_variant={sec_per_variant}")
    print(f"[plan_shards] wrote {len(shards)} shards -> {out_path} "
          f"({', '.join(f'{lvl}={n}' for lvl, n in sorted(by_level.items()))})")
    print(f"[plan_shards] max est_hours={max_h:.2f}h, "
          f"{in_band}/{len(shards)} shards inside the "
          f"{BAND_LO_H:.0f}-{BAND_HI_H:.0f}h target band")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
