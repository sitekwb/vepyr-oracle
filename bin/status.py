#!/usr/bin/env python3
"""Per-element pipeline state: done / annotated / missing, plus the early-warning
runtime-overrun check.

Usage: status.py --version {115,116}

Reads (from $VEPYR_WORK, unless overridden by --work):
  shards_annotate.tsv, shards_diff.tsv     -- always (both versions run vepyr
                                              annotate + diff; the shard GEOMETRY
                                              is combo/chrom-driven, not
                                              version-driven, so one plan covers
                                              both --version 115 and 116 runs)
  shards_gt.tsv                            -- only for --version 116 (the 115
                                              ground truth already exists as
                                              static files; there is no gt step
                                              for it)

For every shard, prints its state (oracle.status.State: missing / annotated /
done) and, for a `gt`/`annotate` shard whose output already exists, the ACTUAL
elapsed time recorded in `<out>.timing.json` versus the shard's PLANNED
`est_hours` (oracle.shards's safety-adjusted estimate). An element that ran more
than oracle.status.OVERRUN_FACTOR times its planned estimate is flagged -- the
early-warning that the runtime MODEL (sec_per_variant) is wrong, catchable while
the element is still comfortably under the 24h cluster cap, instead of only
being discoverable once a job gets killed at the wall-clock limit.

A DONE element that did NOT overrun prints nothing (there is nothing to act on);
every other state, and every overrun, is listed so the operator's eye goes
straight to what needs attention.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.shards import Step, load_shards
from oracle.status import OVERRUN_FACTOR, State, scan, summarize

#: Overridable so tests can point this CLI at a tmp dir, without touching
#: ~/vepyr/{data,work} at all -- same convention as bin/run_wgs.py / bin/validate.py.
DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))

#: (shards.tsv filename, whether this step applies to a given version) -- the gt
#: step only ever applies to 116 (see module docstring).
_STEP_FILES: tuple[tuple[str, Step, frozenset[int]], ...] = (
    ("shards_annotate.tsv", Step.ANNOTATE, frozenset({115, 116})),
    ("shards_diff.tsv", Step.DIFF, frozenset({115, 116})),
    ("shards_gt.tsv", Step.GT, frozenset({116})),
)


def _resultsdir(step: Step, version: int) -> str:
    """The directory a step's outputs land in. MUST match the directories
    slurm/annotate.sbatch, slurm/diff.sbatch and slurm/gen_gt116.sh actually
    write to -- see those scripts' header comments; this is the canonical
    Python-side definition of that convention."""
    if step is Step.GT:
        return os.path.join(DATA_DIR, "ground_truth_vep_116", "shards")
    return os.path.join(WORK_DIR, f"results_wgs_{version}")


def _load_all_shards(version: int) -> dict[Step, list]:
    by_step: dict[Step, list] = {}
    for filename, step, versions in _STEP_FILES:
        if version not in versions:
            continue
        path = os.path.join(WORK_DIR, filename)
        if not os.path.exists(path):
            print(f"[status] {path} not found -- skipping the {step.value} step "
                  f"(run bin/plan_shards.py first)", file=sys.stderr)
            continue
        by_step[step] = load_shards(path)
    return by_step


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="status.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--version", type=int, required=True, choices=[115, 116])
    args = ap.parse_args(argv)

    by_step = _load_all_shards(args.version)
    if not by_step:
        print(f"[FATAL] no shards.tsv files found for version {args.version} under "
              f"{WORK_DIR} -- run bin/plan_shards.py first", file=sys.stderr)
        return 2

    total_shards = 0
    total_overran = 0
    for step, step_shards in by_step.items():
        resultsdir = _resultsdir(step, args.version)
        statuses = scan(step_shards, resultsdir)
        counts = summarize(statuses)
        total_shards += len(step_shards)
        total_overran += counts.get("overran", 0)

        print(f"\n=== {step.value} (version {args.version}, {len(step_shards)} "
             f"shards, {resultsdir}) ===")
        for st in statuses:
            if st.state == State.DONE and not st.overran:
                continue           # a clean, on-time DONE element needs no attention
            tag = " OVERRAN" if st.overran else ""
            actual = f" actual={st.actual_hours:.1f}h" if st.actual_hours is not None else ""
            print(f"  {st.state:10s}{tag:9s} {st.shard.combo:32s} "
                 f"chrom={st.shard.chrom} chunk={st.shard.chunk} "
                 f"est={st.shard.est_hours:.1f}h{actual}")
        summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()) if k != "overran")
        print(f"  -- {summary}, overran={counts.get('overran', 0)}")

    print(f"\n[status] version {args.version}: {total_shards} total shards across "
         f"{len(by_step)} step(s)")
    if total_overran:
        print(f"[status] WARNING: {total_overran} element(s) ran more than "
             f"{OVERRUN_FACTOR}x their planned est_hours -- the runtime model may be "
             f"wrong. Re-measure sec_per_variant and re-plan before submitting more "
             f"of this shape, rather than trusting the existing shards.tsv.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
