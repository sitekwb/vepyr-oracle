#!/usr/bin/env python3
"""Annotate ONE combo, ONE vepyr version, in ONE process.

THIS IS THE 24H FIX: `legacy/validate.py` loops over all 10 combos in a single
process, annotate-then-diff, serially. On whole-genome data that is ~10 x 5h = 50h+,
and the cluster hard-caps every job at 24h -- that monolithic loop is why a
whole-genome run was never possible. Splitting annotation into one-combo-per-
invocation CLIs is what lets SLURM run combos as parallel array elements, each
safely inside the cap. `bin/validate.py` is this script's diff-side counterpart;
`oracle/shards.py` / `bin/plan_shards.py` decide how finely a combo itself must be
split (L0 whole-WGS / L1 per-chromosome / L2 per-chunk) to stay in the 12-16h band.

Usage:
  run_wgs.py --combo NAME --version {115,116} --out FILE
             [--chrom C] [--region CHROM:START-END] [--workers N]

Input resolution (--chrom and --region are mutually exclusive):
  --region CHROM:START-END -> the pre-sliced L2 chunk at
      $VEPYR_WORK/input/input_<chrom>_<start>_<end>.vcf.gz
  --chrom C                 -> the pre-sliced L1 chromosome at
      $VEPYR_WORK/input/input_<chrom>.vcf.gz
  neither                    -> the whole-WGS input at
      $VEPYR_DATA/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  This script only READS the slice; an earlier prep step is responsible for
  producing it (bcftools view / tabix, not this script's job).

Resumable: if --out already exists and is non-empty, this exits 0 immediately --
before importing vepyr, before touching matrix.tsv. Preemption and the 24h cap
make re-running a shard that already finished routine, and this must be a no-op.

`import vepyr` is LAZY: it happens only inside main(), after every argument check
and the resume short-circuit. That is what keeps `--help` (and the whole oracle/
test suite) working on a laptop with no vepyr installed -- see oracle/matrix.py
and oracle/csq.py, which import ONLY the stdlib for the same reason.

Emits `<out>.timing.json` with `{combo, version, chrom, n_variants, elapsed_s,
sec_per_variant}` -- the measured per-combo rate `plan_shards.py` requires (it
refuses to guess; see oracle/shards.py's `rate_for()`).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.csq import open_maybe_gzip
from oracle.matrix import load_matrix, resolve_kwargs

#: Overridable so the smoke tests (tests/test_cli.py) can point this CLI at a tmp
#: dir instead of the real cluster paths, without touching ~/vepyr at all.
DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))

MATRIX_PATH = os.path.join(WORK_DIR, "matrix.tsv")
WHOLE_WGS_INPUT = os.path.join(DATA_DIR, "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
REFERENCE_FASTA = os.path.join(DATA_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
PLUGIN_CACHE_ROOT = os.path.join(DATA_DIR, "plugin_cache")


def _parse_region(region: str) -> tuple[str, int, int]:
    """`"CHROM:START-END"` -> `(chrom, start, end)`."""
    try:
        chrom, span = region.split(":", 1)
        start_s, end_s = span.split("-", 1)
        return chrom, int(start_s), int(end_s)
    except ValueError as e:
        raise ValueError(f"--region must be CHROM:START-END, got {region!r}") from e


def resolve_input(args: argparse.Namespace) -> tuple[str, str]:
    """-> (input_vcf_path, chrom label to record in the timing JSON)."""
    if args.region:
        chrom, start, end = _parse_region(args.region)
        path = os.path.join(WORK_DIR, "input", f"input_{chrom}_{start}_{end}.vcf.gz")
        return path, chrom
    if args.chrom:
        return os.path.join(WORK_DIR, "input", f"input_{args.chrom}.vcf.gz"), args.chrom
    return WHOLE_WGS_INPUT, "all"


def count_variants(vcf_path: str) -> int:
    """Cheap line count of non-header records.

    The input is bgzipped; bgzip is a valid gzip stream, so the stdlib `gzip`
    module (via `open_maybe_gzip`) reads it fine -- no pysam dependency needed just
    to count lines. See the module docstring / task report for the tradeoff on a
    whole-WGS (~4.5M line) input.
    """
    n = 0
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="run_wgs.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--combo", required=True, help="combo name, must exist in matrix.tsv")
    ap.add_argument("--version", type=int, required=True, choices=[115, 116])
    ap.add_argument("--out", required=True, help="output annotated VCF path")
    loc = ap.add_mutually_exclusive_group()
    loc.add_argument("--chrom", help="e.g. '22' -- L1 per-chromosome input")
    loc.add_argument("--region", help="CHROM:START-END -- L2 per-chunk input")
    ap.add_argument("--workers", type=int, default=8)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    if os.path.exists(args.out) and os.path.getsize(args.out) > 0:
        print(f"[run_wgs] {args.combo}/{args.version}: {args.out} already exists and "
              f"is non-empty -- resuming (skip)")
        return 0

    matrix = load_matrix(MATRIX_PATH)
    if args.combo not in matrix:
        print(f"[FATAL] combo {args.combo!r} not found in {MATRIX_PATH} "
              f"(known combos: {sorted(matrix)})", file=sys.stderr)
        return 2
    row = matrix[args.combo]

    cache_col = f"cache{args.version}"
    cache_dir = os.path.join(DATA_DIR, row[cache_col])
    kwargs = resolve_kwargs(json.loads(row["vepyr_kwargs"]),
                            plugin_cache_root=PLUGIN_CACHE_ROOT)

    inp, chrom_label = resolve_input(args)
    if not os.path.exists(inp):
        print(f"[FATAL] input {inp!r} does not exist", file=sys.stderr)
        return 2

    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    import vepyr  # LAZY -- see module docstring; oracle/ stays importable without it

    t0 = time.time()
    vepyr.annotate(inp, cache_dir, reference_fasta=REFERENCE_FASTA,
                   output_vcf=args.out, show_progress=False, workers=args.workers,
                   **kwargs)
    elapsed_s = time.time() - t0

    n_variants = count_variants(inp)
    sec_per_variant = elapsed_s / n_variants if n_variants else None

    timing = {
        "combo": args.combo,
        "version": args.version,
        "chrom": chrom_label,
        "n_variants": n_variants,
        "elapsed_s": round(elapsed_s, 1),
        "sec_per_variant": round(sec_per_variant, 6) if sec_per_variant is not None else None,
    }
    with open(f"{args.out}.timing.json", "w") as fh:
        json.dump(timing, fh, indent=2)

    print(f"[run_wgs] {args.combo}/{args.version} chrom={chrom_label}: "
          f"{n_variants} variants in {elapsed_s:.1f}s "
          f"({timing['sec_per_variant']} s/variant) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
