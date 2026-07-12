#!/usr/bin/env python3
"""Turn an L2 chunk index into a genomic region string.

Usage:
  make_regions.py --chrom C --num-chunks N [--input FILE] [--out FILE]

Streams the chromosome's variant positions ONCE (`oracle.regions.read_positions`,
pure stdlib -- no vepyr, no pysam, so this stays laptop-testable) and cuts at
empirical variant-count quantiles (`oracle.regions.quantile_regions`; see that
module's docstring for why an EQUAL-LENGTH split would silently invalidate the
wall-clock estimate the whole 12-16h plan rests on). Writes
`$VEPYR_WORK/regions_<chrom>.tsv`: `chunk`, `region` (`"chrom:start-end"`),
`n_variants`. Regions abut exactly -- no gap, no overlap -- and the self-check
(chunk counts sum to the chromosome total) runs, and would raise, before
anything is written.

`--num-chunks` MUST equal the L2 chunk count `shards.tsv` planned for this
(step, chromosome) -- chunk index N in `shards.tsv` only means the same region
as chunk N in `regions_<chrom>.tsv` if both were cut at the SAME resolution. If
a regions file already exists for this chromosome but holds a DIFFERENT number
of chunks (e.g. two combos of the same step disagreed on how many L2 chunks
this chromosome needs -- see oracle/shards.py's per-combo-rate docstring for why
that can happen), this REFUSES to silently reuse it: a stale-resolution regions
file is exactly how a shard's `--region` ends up pointing at the wrong slice of
the genome. Reconcile the shard plan (or pass a distinct `--out`) instead of
overriding this check.

Resumable: if `--out` already exists and matches `--num-chunks`, this is a
no-op. A large chromosome's position stream can be millions of lines; a shard
must not re-pay that scan on every array element that happens to need it.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.regions import RegionError, load_regions, quantile_regions, read_positions, write_regions

#: Overridable so tests can point this CLI at a tmp dir, without touching
#: ~/vepyr/work at all -- same convention as bin/run_wgs.py / bin/validate.py.
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="make_regions.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--chrom", required=True, help="e.g. '22' -- no 'chr' prefix")
    ap.add_argument("--num-chunks", type=int, required=True)
    ap.add_argument("--input", help="default: $VEPYR_WORK/input/input_<chrom>.vcf.gz "
                                    "(the output of slurm/prep_input.sh)")
    ap.add_argument("--out", help="default: $VEPYR_WORK/regions_<chrom>.tsv")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    inp = args.input or os.path.join(WORK_DIR, "input", f"input_{args.chrom}.vcf.gz")
    out = args.out or os.path.join(WORK_DIR, f"regions_{args.chrom}.tsv")

    if os.path.exists(out) and os.path.getsize(out) > 0:
        existing = load_regions(out)
        if len(existing) == args.num_chunks:
            print(f"[make_regions] {out} already cut into {len(existing)} chunks -- "
                  f"resuming (skip)")
            return 0
        print(f"[FATAL] {out} already exists but holds {len(existing)} chunk(s), not "
              f"the requested {args.num_chunks}. Two shard rows for chromosome "
              f"{args.chrom!r} disagree on the L2 chunk count -- refusing to silently "
              f"reuse a different-resolution regions file (chunk index N would then "
              f"point at the wrong slice of the genome for one of them). Reconcile the "
              f"shard plan, or pass a distinct --out.", file=sys.stderr)
        return 2

    if not os.path.exists(inp):
        print(f"[FATAL] input {inp!r} does not exist (run slurm/prep_input.sh first)",
              file=sys.stderr)
        return 2

    positions = read_positions(inp, chrom=args.chrom)
    try:
        regions = quantile_regions(args.chrom, positions, args.num_chunks)
    except (ValueError, RegionError) as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        return 1

    out_dir = os.path.dirname(os.path.abspath(out))
    os.makedirs(out_dir, exist_ok=True)
    out_tmp = f"{out}.tmp"
    write_regions(out_tmp, regions)
    os.replace(out_tmp, out)          # atomic: no partial regions file survives a kill

    sizes = [r.n_variants for r in regions]
    print(f"[make_regions] chr{args.chrom}: {len(positions)} variants -> "
          f"{len(regions)} regions (sizes {min(sizes)}-{max(sizes)}, "
          f"sum={sum(sizes)}) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
