#!/usr/bin/env python3
"""Fork-delta: diff ONE combo x ONE chromosome of the UNFORKED vs FORKED VEP-116 GT.

WHAT THIS ANSWERS
-----------------
The 116 ground truth exists in two arms, generated from the SAME input
(work/input/input_<chrom>.vcf.gz), the SAME native 116 cache and the SAME VEP
container, with command lines that differ in EXACTLY one token -- `--fork 16`:

  forked   (Faza A)  $VEPYR_DATA/ground_truth_vep_116/shards/
  unforked (Faza B)  $VEPYR_DATA/ground_truth_vep_116_unforked/shards/

The standing claim is that `--fork` corrupts ONLY `HGNC_ID`, via VEP's
buffer-scoped gene-symbol back-fill (maxForkSize = buffer_size / (2 * fork), so
`--fork 16` shrinks the InputBuffer ~32x and the transcript that would have
donated the HGNC id falls outside it). That claim was only ever spot-checked.
This script tests it at whole-genome scale, per FIELD, using the SAME comparison
engine (`oracle.diff.diff_files`) that every vepyr-vs-VEP number in this project
is computed with -- so the fork delta is measured on exactly the same footing.

SLOT ASSIGNMENT (read this before reading any output)
-----------------------------------------------------
`diff_files()` is written for vepyr-vs-VEP, so its columns and categories are
named for that. Here BOTH sides are real VEP, and the slots are:

    "vepyr" slot  <- UNFORKED  (Faza B)   -- columns `vepyr_val`, cat VEPYR_ONLY
    "gt"    slot  <- FORKED    (Faza A)   -- columns `vep_val`,  cat VEP_ONLY

So `VEPYR_ONLY` reads "unforked populated it, forked left it EMPTY" -- which is
precisely what the back-fill hypothesis predicts for HGNC_ID -- and `VEP_ONLY`
reads "forked populated it, unforked left it empty" (the hypothesis predicts
NONE of these). `VALUE_DIFF` means both arms emitted a value and they disagree,
which the hypothesis predicts for NO field at all.

The FORKED tree is opened READ-ONLY and is never written to by anything here.

WHY PER-CHROMOSOME SHARDS AND NOT THE MERGED VCFs
-------------------------------------------------
A merged combo VCF is ~28 GB of text; two of them streamed through Python is a
single-threaded hour per combo. The per-chrom shards are already chromosome-
atomic -- exactly `oracle.summary.merge()`'s shard contract -- so 22 shards diff
in parallel as a SLURM array and `bin/merge_summaries.py` folds them into one
whole-genome summary WITH its completeness check (a chromosome whose element
died cannot silently vanish from the denominator).

Usage:
  fork_delta.py --combo NAME --chrom C --outdir DIR
                [--unforked-root DIR] [--forked-root DIR]

Writes:
  <outdir>/summary_<combo>_<chrom>.json   (schema-identical to bin/validate.py's,
                                           so bin/merge_summaries.py just works)
  <outdir>/mismatches/<combo>_<chrom>.tsv

Resumable: an existing non-empty summary JSON short-circuits to exit 0.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.diff import diff_files
from oracle.matrix import kwargs_from_row, load_matrix

DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))
MATRIX_PATH = os.path.join(WORK_DIR, "matrix.tsv")

#: Defaults mirror slurm/gen_gt116.sh's OUTDIR convention: the FORKED tree is the
#: bare `ground_truth_vep_116`, the unforked regeneration lives beside it under
#: the `_unforked` suffix (that is what VEPYR_GT116_DIR selected at generation time).
DEFAULT_UNFORKED_ROOT = os.path.join(DATA_DIR, "ground_truth_vep_116_unforked")
DEFAULT_FORKED_ROOT = os.path.join(DATA_DIR, "ground_truth_vep_116")

_COVERAGE_KEYS = ("join_rate", "only_vepyr", "only_gt", "malformed_vepyr",
                  "malformed_gt", "duplicate_records_vepyr", "duplicate_records_gt")


def shard_path(root: str, combo: str, chrom: str) -> str:
    """The L1 (per-chromosome) shard slurm/gen_gt116.sh wrote for this combo.

    MUST mirror gen_gt116.sh's SFX construction: `SFX="_$CHROM"` for a non-"all"
    chromosome, plus `_c<chunk>` for L2 only. The whole 116 GT was planned L1
    (shards_gt.tsv is 22 per-chromosome rows per combo, no L2), so only the L1
    shape is constructed here -- an L2 shard would need a chunk and there is no
    honest way to guess one.
    """
    return os.path.join(root, "shards", f"{combo}_{chrom}.vcf")


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="fork_delta.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--combo", required=True, help="combo name, must exist in matrix.tsv")
    ap.add_argument("--chrom", required=True, help="chromosome (bare, e.g. '22')")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--unforked-root", default=DEFAULT_UNFORKED_ROOT,
                    help=f"GT tree that goes in the 'vepyr' slot (default: {DEFAULT_UNFORKED_ROOT})")
    ap.add_argument("--forked-root", default=DEFAULT_FORKED_ROOT,
                    help=f"GT tree that goes in the 'gt' slot, READ-ONLY "
                         f"(default: {DEFAULT_FORKED_ROOT})")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    out_summary = os.path.join(args.outdir, f"summary_{args.combo}_{args.chrom}.json")
    if os.path.exists(out_summary) and os.path.getsize(out_summary) > 0:
        print(f"[fork_delta] {out_summary} already exists -- resuming (skip)")
        return 0

    matrix = load_matrix(MATRIX_PATH)
    if args.combo not in matrix:
        print(f"[FATAL] combo {args.combo!r} not found in {MATRIX_PATH} "
              f"(known: {sorted(matrix)})", file=sys.stderr)
        return 2
    row = matrix[args.combo]
    try:
        combo_kwargs = kwargs_from_row(row)
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        return 2

    unforked = shard_path(args.unforked_root, args.combo, args.chrom)
    forked = shard_path(args.forked_root, args.combo, args.chrom)
    # NEVER a no_gt-style soft stub here (unlike bin/validate.py, where a missing
    # GT is an expected mid-run state). Both arms are supposed to exist for every
    # comparable combo; a missing shard means the comparable-combo intersection
    # was computed wrong, and answering "no difference" from a file that is not
    # there is the exact silent-wrongness this pipeline exists to prevent.
    for label, path in (("unforked", unforked), ("forked", forked)):
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"[FATAL] {label} shard missing or empty: {path}. This combo/chrom "
                  f"is NOT part of the comparable intersection -- refusing to report a "
                  f"fork delta computed from one arm.", file=sys.stderr)
            return 2

    out_mismatches = os.path.join(args.outdir, "mismatches",
                                  f"{args.combo}_{args.chrom}.tsv")
    os.makedirs(os.path.dirname(out_mismatches), exist_ok=True)

    summary = diff_files(unforked, forked, name=args.combo, cache=row["cache116"],
                         combo_kwargs=combo_kwargs, tsv_path=out_mismatches,
                         chrom=args.chrom)

    with open(out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    coverage = " ".join(f"{k}={summary[k]}" for k in _COVERAGE_KEYS)
    print(f"[fork_delta] {args.combo} chr{args.chrom}: overall_pct={summary['overall_pct']} "
          f"aligned={summary['aligned_annotations']} {coverage}")
    moved = {f: v for f, v in summary["per_field"].items()
             if v["total"] and v["match"] != v["total"]}
    if moved:
        for f, v in sorted(moved.items(), key=lambda kv: kv[1]["match"] - kv[1]["total"]):
            print(f"[fork_delta]   FIELD MOVED {f}: {v['total'] - v['match']}/{v['total']} "
                  f"cells differ (pct={v['pct']}) {v['by_category']}")
    else:
        print("[fork_delta]   no field moved")
    print(f"[fork_delta] wrote {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
