#!/usr/bin/env python3
"""Merge the per-chromosome shard summaries into per-combo WHOLE-GENOME summaries.

`bin/validate.py` diffs one combo x one chromosome per process (the shape that fits
the cluster's 24h wall-clock cap) and leaves one `summary_<combo>_<chrom>.json` behind.
This script folds those shards into the single `summary.json` that `bin/make_report.py`
reads: a LIST of per-combo whole-genome summaries, one entry per combo in
`oracle.matrix.COMBOS` order.

Usage:
  merge_summaries.py <resultsdir> [--chroms 1-22] [-o OUT]

THE COMPLETENESS CONTRACT
-------------------------
`oracle.summary.merge()` is always called WITH `expected_chroms`. Without it, a
chromosome whose SLURM job died without writing its JSON simply vanishes from the
denominator, and a "whole-genome" concordance number gets published from 21
chromosomes with nothing anywhere saying so. `merge()` raises on a missing (or
unexpected, or duplicated) chromosome and this script does NOT catch it: a partial
run must fail loudly rather than round itself up to a whole genome.

`--chroms` exists so a deliberately partial run (a chr22-only smoke check) can be
merged and reported HONESTLY -- the chromosomes actually covered are carried into the
summary's `chroms` key and printed on page 1 of the PDF, so a 2-chromosome report can
never masquerade as a whole-genome one. It is not an escape hatch from the contract:
whatever set you name, the shards must match it EXACTLY.

WHY THIS DISPATCHES ON THE JSON'S `name` AND NOT ON THE FILENAME
---------------------------------------------------------------
The obvious implementation -- `glob(f"summary_{combo}_*.json")` -- is wrong, and
quietly so: combo names are prefixes of one another, so globbing for `hgvs_merged`
also matches `summary_hgvs_merged_pick_allele_21.json`. Every shard of every `pick`
combo would be swallowed into `hgvs_merged`'s merge. (`merge()` would catch it via the
combo-identity check -- but only by killing the run, and only because that guard
happens to exist.) Reading the `name` the shard itself declares is unambiguous by
construction.
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import COMBOS
from oracle.report import coverage
from oracle.summary import merge

#: The autosomes a whole-genome run covers. chrX/chrY/chrM are deliberately NOT in the
#: default: the GIAB HG002 benchmark's confident regions are autosomal, and a run that
#: DOES cover them must say so explicitly with --chroms rather than have them silently
#: appear or silently vanish.
DEFAULT_CHROMS = tuple(str(i) for i in range(1, 23))


def _parse_chroms(spec: str) -> set[str]:
    """"1-22" / "21,22" / "1-4,X" -> the exact chromosome set the shards must cover."""
    chroms: set[str] = set()
    for part in (p.strip() for p in spec.split(",") if p.strip()):
        lo, sep, hi = part.partition("-")
        if sep and lo.isdigit() and hi.isdigit():
            chroms.update(str(i) for i in range(int(lo), int(hi) + 1))
        else:
            chroms.add(part)
    if not chroms:
        raise argparse.ArgumentTypeError(f"no chromosomes in {spec!r}")
    return chroms


def load_shards(resultsdir: str) -> dict[str, list[dict]]:
    """Every `summary_*.json` in the dir, grouped by the combo name it DECLARES.

    (Not by filename -- see the module docstring for the prefix-collision bug that
    makes filename globbing unsafe. `summary.json`, this script's own output, does not
    match the pattern and is never re-read.)
    """
    by_combo: dict[str, list[dict]] = collections.defaultdict(list)
    for path in sorted(glob.glob(os.path.join(resultsdir, "summary_*.json"))):
        with open(path) as fh:
            shard = json.load(fh)
        name = shard.get("name")
        if not name:
            raise ValueError(f"{path}: summary has no 'name' -- cannot tell which combo "
                             f"it belongs to; refusing to guess")
        by_combo[name].append(shard)
    return dict(by_combo)


def merge_combo(name: str, shards: list[dict], expected_chroms: set[str]) -> dict:
    """One combo's shards -> one whole-genome summary (or its `no_gt` stub).

    A combo with NO ground truth has no per-chromosome structure to merge: every shard
    validate.py wrote for it is the same three-key stub. Those are passed through
    unchanged, so the combo still appears in the report (as `no_gt`, contributing no
    evidence) instead of disappearing from it -- but a MIXTURE of stubs and real diffs
    is incoherent (ground truth exists per combo, not per chromosome) and is refused.
    """
    stubs = [s for s in shards if s.get("status") != "ok"]
    if stubs and len(stubs) != len(shards):
        raise ValueError(
            f"combo {name!r}: {len(stubs)} of {len(shards)} shards have status != 'ok' "
            f"(statuses: {sorted({s.get('status') for s in stubs})}). Ground truth exists "
            f"per COMBO, not per chromosome, so a mixture of diffed and un-diffed shards "
            f"means something went wrong mid-run; refusing to merge a partial combo into "
            f"a whole-genome number."
        )
    if stubs:
        stub = dict(stubs[0])
        stub["chroms"] = sorted(filter(None, (s.get("chrom") for s in shards)))
        return stub
    return merge(shards, expected_chroms=expected_chroms)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="merge_summaries.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("resultsdir", help="dir holding summary_<combo>_<chrom>.json")
    ap.add_argument("--chroms", default="1-22", type=_parse_chroms,
                    help="the chromosome set the shards MUST cover exactly "
                         "(default: the 22 autosomes)")
    ap.add_argument("-o", "--out", help="default: <resultsdir>/summary.json")
    args = ap.parse_args(argv)

    out_path = args.out or os.path.join(args.resultsdir, "summary.json")
    by_combo = load_shards(args.resultsdir)

    known = {c.name for c in COMBOS}
    unknown = sorted(set(by_combo) - known)
    if unknown:
        print(f"[FATAL] {args.resultsdir} holds summaries for combo(s) {unknown} that are "
              f"not in oracle.matrix.COMBOS (known: {sorted(known)}). Refusing to guess "
              f"what they are -- a results dir holds exactly one version's runs.",
              file=sys.stderr)
        return 2

    missing = [c.name for c in COMBOS if c.name not in by_combo]
    if missing:
        print(f"[FATAL] no shard summaries at all for combo(s) {missing} in "
              f"{args.resultsdir}. A whole-genome report must not silently omit a combo: "
              f"run bin/validate.py for them, or delete them from oracle.matrix.COMBOS.",
              file=sys.stderr)
        return 2

    merged = []
    for combo in COMBOS:
        summary = merge_combo(combo.name, by_combo[combo.name], args.chroms)
        merged.append(summary)
        cov = coverage(summary)
        trust = "OK" if cov["trustworthy"] else "NOT TRUSTWORTHY"
        # Never print the headline % without the join rate beside it -- the whole
        # point of this pipeline's coverage work. See oracle/report.py.
        print(f"[merge] {combo.name:<28} status={summary.get('status')} "
              f"overall_pct={summary.get('overall_pct')} "
              f"join_rate={summary.get('join_rate')} "
              f"aligned={summary.get('aligned_annotations', 0)} "
              f"chroms={len(summary.get('chroms', []))} -> {trust}")

    with open(out_path, "w") as fh:
        json.dump(merged, fh, indent=2)
    print(f"[merge] wrote {out_path} ({len(merged)} combos)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
