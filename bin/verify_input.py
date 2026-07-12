#!/usr/bin/env python3
"""GATE 2: prove OUR normalized input carries the same variant keys real VEP saw.

Usage:
  verify_input.py --ours FILE --gt FILE --chrom C [--examples N]

`--ours` is our candidate normalized VCF (whole-genome or already per-chrom
sliced -- either works, this only ever looks at the ONE chromosome `--chrom`
names). `--gt` is a ground-truth VCF real VEP produced FROM the input we are
trying to reproduce (e.g. `data/ground_truth_vep/
HG002_annotated_wgs_everything_hgvs_refseq.vcf` on the cluster).

Exits 0 and prints "GATE 2 PASSED" if the (chrom,pos,ref,alt) key sets are
IDENTICAL on `--chrom`. Exits 1 and prints up to `--examples` keys unique to
each side otherwise. Exits 2 on a usage/file error (missing `--ours`/`--gt`).

See oracle/verify_input.py's module docstring for why literal key-set
identity (not just a record-count match) is what actually proves the
downstream (chrom,pos,ref,alt) join will work, and why this is deliberately
scoped to ONE chromosome at a time (memory characteristic: `--chrom` is
required, never "all").

No pysam import at module level -- see oracle/verify_input.py; this stays
runnable (at least --help) with no pysam/vepyr installed, same
laptop-testability property the rest of bin/ pins.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.verify_input import compare_key_sets, format_examples, read_keys


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="verify_input.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--ours", required=True, help="our candidate normalized VCF")
    ap.add_argument("--gt", required=True, help="ground-truth VCF real VEP produced")
    ap.add_argument("--chrom", required=True,
                    help="ONE chromosome to check, e.g. '22' -- see module "
                         "docstring on why this is never 'all'")
    ap.add_argument("--examples", type=int, default=10,
                    help="how many example-only keys to print per side on a mismatch "
                         "(default 10)")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    for label, path in (("--ours", args.ours), ("--gt", args.gt)):
        if not os.path.exists(path):
            print(f"[FATAL] {label} file does not exist: {path}", file=sys.stderr)
            return 2

    ours = read_keys(args.ours, args.chrom)
    gt = read_keys(args.gt, args.chrom)
    comparison = compare_key_sets(ours, gt, chrom=args.chrom)

    print(f"[verify_input] chrom={comparison.chrom} ours_total={comparison.ours_total} "
          f"gt_total={comparison.gt_total} overlap={comparison.overlap} "
          f"only_ours={len(comparison.only_ours)} only_gt={len(comparison.only_gt)}")

    if comparison.identical:
        print(f"[verify_input] GATE 2 PASSED -- chr{comparison.chrom} key sets are "
              f"IDENTICAL ({comparison.ours_total} variants)")
        return 0

    print(f"[verify_input] GATE 2 FAILED -- chr{comparison.chrom} key sets differ",
          file=sys.stderr)
    if comparison.only_ours:
        print(f"[verify_input]   only in --ours ({len(comparison.only_ours)}), examples:",
              file=sys.stderr)
        for ex in format_examples(comparison.only_ours, args.examples):
            print(f"[verify_input]     {ex}", file=sys.stderr)
    if comparison.only_gt:
        print(f"[verify_input]   only in --gt ({len(comparison.only_gt)}), examples:",
              file=sys.stderr)
        for ex in format_examples(comparison.only_gt, args.examples):
            print(f"[verify_input]     {ex}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
