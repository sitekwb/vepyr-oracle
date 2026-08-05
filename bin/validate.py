#!/usr/bin/env python3
"""Diff ONE combo's already-annotated vepyr VCF against ground truth, ONE process.

This is `bin/run_wgs.py`'s diff-side counterpart -- see that script's docstring for
why one-combo-per-invocation is what fits a whole-genome validation run inside the
cluster's 24h wall-clock cap. `run_wgs.py` produces the annotated VCF; this script
consumes it (`--vepyr FILE`) and never calls vepyr itself, so it has no lazy-import
concern of its own -- but it still avoids importing anything heavier than
`oracle.diff` / `oracle.matrix`, both pure stdlib.

Usage:
  validate.py --combo NAME --version {115,116} --vepyr FILE --outdir DIR [--chrom C]

Writes:
  <outdir>/summary_<combo>[_<chrom>].json
  <outdir>/mismatches/<combo>[_<chrom>].tsv

Resumable: if the summary JSON already exists (and is non-empty), this exits 0
immediately -- re-running a shard's diff after preemption or a resubmit is a no-op.

If the ground-truth VCF for this combo/version is missing, this writes a minimal
`{"name":..., "cache":..., "status":"no_gt"}` summary and exits 0. A missing GT is a
known, expected gap (not every combo necessarily has both 115 and 116 ground truth
at every point in the run) -- it must NEVER crash a SLURM array job that has dozens
of healthy siblings still running.

Otherwise this calls `oracle.diff.diff_files()` and dumps its summary JSON, then
prints the headline number ALONGSIDE the coverage signals that can silently hollow
it out: `overall_pct` and `join_rate` are always printed together (a diff that only
joined 2% of annotations can still show 100% overall_pct on the scraps it managed to
compare), plus `only_vepyr` / `only_gt` / `malformed_*` / `duplicate_records_*` --
see oracle/diff.py and oracle/summary.py for the concrete failure mode each of these
guards against.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
from oracle.diff import diff_files
from oracle.matrix import kwargs_from_row, load_matrix, resolve_matrix_path

#: Overridable so the smoke tests (tests/test_cli.py) can point this CLI at a tmp
#: dir instead of the real cluster paths, without touching ~/vepyr at all.
DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))

# matrix.tsv is the combos matrix (one row per oracle.matrix.COMBOS entry): which
# vepyr kwargs to run each combo with, DERIVED by bin/seed_matrix.py from the
# ground truth's own `##VEP-command-line=` headers (see that script's docstring).
# It used to exist ONLY on the cluster, at $VEPYR_WORK/matrix.tsv -- so a fresh
# clone of this repo (e.g. the container image being built around it) had no
# matrix.tsv anywhere, and this CLI could not run at all. A copy fetched from the
# cluster (2026-07-23) is now versioned at <repo root>/matrix.tsv -- see
# matrix.tsv.README beside it for provenance (sha256 included) -- as a fallback
# for exactly that case. It is not regenerated at clone/build time: regeneration
# needs bin/seed_matrix.py's multi-gigabyte ground truth input, which a laptop or
# CI runner will not have.
#
# The 3-way fallback chain ($VEPYR_MATRIX > $WORK_DIR/matrix.tsv > the shipped
# copy) now lives in oracle.matrix.resolve_matrix_path() -- see its docstring for
# the exact order and reasoning -- so every consumer shares ONE implementation of
# this precedence rule instead of each keeping its own copy that can silently
# drift out of sync (that drift is exactly how bin/run_wgs.py ended up with no
# fallback chain at all: this rule used to be written out inline, here only).
MATRIX_PATH = resolve_matrix_path(WORK_DIR, _REPO_ROOT)

# The ground truth lives in a DIFFERENT directory per VEP version:
#   115 -> data/ground_truth_vep/                (shipped with the project)
#   116 -> data/ground_truth_vep_116_unforked/    (minted by slurm/gen_gt116.sh +
#                                                  merge_gt116.sh with VEPYR_GT116_DIR
#                                                  pointed at the _unforked tree)
# Hardcoding the 115 dir made every --version 116 diff report status=no_gt while a
# perfectly good 116 ground truth sat in the other directory.
#
# 2026-07-30/2026-08-02: the 116 entry used to be plain "ground_truth_vep_116" --
# ground truth minted with `--fork N>1`. That tree was deleted on 2026-07-30,
# correctly: forking shrinks VEP's InputBuffer (maxForkSize = buffer_size /
# (2 * fork)), so `--fork 16` narrows the annotation buffer window ~30x and the
# gene-symbol back-fill drops HGNC_ID for transcripts that fall outside it -- see
# bin/fork_delta.py's module docstring for the whole-genome measurement of that
# effect. It was never a valid golden standard. The path here was not moved when
# the directory was, so every --version 116 run has been silently taking the
# no_gt branch below since 2026-07-30 instead of measuring anything -- a gap that
# went unnoticed because a no_gt summary reads, downstream, as "nothing to
# report" rather than "nothing was measured". Now points at the authoritative,
# unforked ground truth.
GT_DIRS = {
    115: os.path.join(DATA_DIR, "ground_truth_vep"),
    116: os.path.join(DATA_DIR, "ground_truth_vep_116_unforked"),
}

#: The coverage signals printed alongside overall_pct -- a concordance percentage
#: without these is untrustworthy (see module docstring).
_COVERAGE_KEYS = ("join_rate", "only_vepyr", "only_gt", "malformed_vepyr",
                  "malformed_gt", "duplicate_records_vepyr", "duplicate_records_gt")


def summary_path(outdir: str, combo: str, chrom: str | None) -> str:
    suffix = f"_{chrom}" if chrom else ""
    return os.path.join(outdir, f"summary_{combo}{suffix}.json")


def mismatches_path(outdir: str, combo: str, chrom: str | None) -> str:
    suffix = f"_{chrom}" if chrom else ""
    return os.path.join(outdir, "mismatches", f"{combo}{suffix}.tsv")


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="validate.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--combo", required=True, help="combo name, must exist in matrix.tsv")
    ap.add_argument("--version", type=int, required=True, choices=[115, 116])
    ap.add_argument("--vepyr", required=True, help="the vepyr-annotated VCF to check")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chrom", help="restrict the diff to one chromosome (an L1/L2 shard)")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    out_summary = summary_path(args.outdir, args.combo, args.chrom)

    if os.path.exists(out_summary) and os.path.getsize(out_summary) > 0:
        print(f"[validate] {args.combo}: {out_summary} already exists -- resuming (skip)")
        return 0

    matrix = load_matrix(MATRIX_PATH)
    if args.combo not in matrix:
        print(f"[FATAL] combo {args.combo!r} not found in {MATRIX_PATH} "
              f"(known combos: {sorted(matrix)})", file=sys.stderr)
        return 2
    row = matrix[args.combo]
    cache = row[f"cache{args.version}"]
    gt_path = os.path.join(GT_DIRS[args.version], row[f"gt{args.version}"])

    # BEFORE the no_gt short-circuit: an unseeded row is a broken MATRIX, not a
    # missing ground truth, and must not be laundered into a benign-looking
    # status=no_gt summary. It is also what drift.py's flag-aware classification
    # reads -- empty kwargs there would silently reclassify real mismatches.
    try:
        combo_kwargs = kwargs_from_row(row)
    except ValueError as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        return 2

    if not os.path.exists(gt_path):
        summary = {"name": args.combo, "cache": cache, "status": "no_gt"}
        with open(out_summary, "w") as fh:
            json.dump(summary, fh, indent=2)
        print(f"[validate] {args.combo}/{args.version}: no ground truth at {gt_path!r} "
              f"-- wrote {out_summary} (status=no_gt)")
        return 0

    out_mismatches = mismatches_path(args.outdir, args.combo, args.chrom)
    os.makedirs(os.path.dirname(out_mismatches), exist_ok=True)

    summary = diff_files(args.vepyr, gt_path, name=args.combo, cache=cache,
                         combo_kwargs=combo_kwargs, tsv_path=out_mismatches,
                         chrom=args.chrom)

    with open(out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    # A concordance % without the join rate is untrustworthy -- never print
    # overall_pct alone. Every other coverage signal that can silently hollow out
    # the headline number goes on the same line.
    coverage = " ".join(f"{k}={summary[k]}" for k in _COVERAGE_KEYS)
    print(f"[validate] {args.combo}/{args.version} chrom={args.chrom or 'all'}: "
          f"overall_pct={summary['overall_pct']} {coverage}")
    print(f"[validate] wrote {out_summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
