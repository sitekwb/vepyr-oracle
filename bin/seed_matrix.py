#!/usr/bin/env python3
"""Seed work/matrix.tsv FROM the ground truth: recovered VEP flags -> vepyr kwargs.

Run on the cluster. For each combo this reads the `##VEP-command-line=` header that
real VEP wrote into the ground-truth VCF, and DERIVES the vepyr kwargs from it
(oracle.matrix.vep_flags_to_vepyr_kwargs). Both sides land in the same matrix.tsv
row, so the configuration we hand vepyr cannot drift from the one the ground truth
was actually generated with.

It had drifted. The hand-written kwargs ran `hgvs_merged_pick` with vepyr's `pick`
while its ground truth had been generated with `--flag_pick_allele_gene`, and ran all
five pick-family combos with VEP's DEFAULT pick order while their ground truth used
an explicit, different one (`biotype,rank,...` -- not `mane_select,canonical` first).
Different order => different transcript picked => the diff reported "vepyr picked the
wrong transcript" for transcripts we had never told vepyr how to choose.

Prints both sides of every combo side by side, so the agreement is REVIEWABLE rather
than merely asserted.

Exits non-zero if any combo cannot be seeded -- a missing `##VEP-command-line=`, an
unmappable flag, or a cache-flavor disagreement all mean we could NOT reproduce that
ground truth with VEP 116, which would make the whole 116 ladder untrustworthy. Such
a combo's row is written with EMPTY cells rather than guessed ones, so that nothing
downstream can mistake "we don't know" for "no flags needed".
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import (COMBOS, extract_vep_command_line, vep_flags_to_vepyr_kwargs,
                           write_matrix)

#: Overridable so the smoke tests (tests/test_seed_matrix_cli.py) can point this CLI
#: at a tmp dir instead of the real cluster paths, without touching ~/vepyr at all.
DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))

GT_DIR = os.path.join(DATA_DIR, "ground_truth_vep")
OUT = os.path.join(WORK_DIR, "matrix.tsv")


def _format_kwargs(kwargs: dict) -> str:
    """`k=v` pairs, one per line-friendly chunk. The NEEDS_PLUGIN_CACHE sentinel
    prints as `<NEEDS_PLUGIN_CACHE>` (its repr) and NOT as a path -- if it ever shows
    up here looking like a real directory, something has resolved it too early."""
    return "  ".join(f"{k}={v!r}" if isinstance(v, str) else f"{k}={v}"
                     for k, v in kwargs.items())


def main() -> int:
    cmdlines: dict[str, str] = {}
    problems: list[str] = []

    for c in COMBOS:
        gt = os.path.join(GT_DIR, c.gt115)
        if not os.path.exists(gt):
            problems.append(f"{c.name}: ground truth missing: {gt}")
            continue
        cmdline = extract_vep_command_line(gt)
        if not cmdline:
            problems.append(f"{c.name}: no ##VEP-command-line= header in {c.gt115}")
            continue

        # Derive here as well as inside write_matrix() (it is a pure function, so the
        # two cannot disagree) for two reasons: to PRINT the derivation, and to turn a
        # flag we cannot map into a reported problem + an empty row, rather than a
        # traceback that aborts before the other seven combos are even looked at.
        try:
            kwargs = vep_flags_to_vepyr_kwargs(cmdline, combo=c.name,
                                               cache_flavor=c.cache_flavor)
        except ValueError as exc:
            problems.append(f"{c.name}: {exc}")
            continue

        cmdlines[c.name] = cmdline
        print(f"[seed] {c.name}  (cache_flavor={c.cache_flavor})\n"
              f"         VEP  : {cmdline}\n"
              f"         vepyr: {_format_kwargs(kwargs)}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    write_matrix(OUT, cmdlines)
    print(f"[seed] wrote {OUT} ({len(cmdlines)}/{len(COMBOS)} combos seeded from their "
          f"ground-truth headers)")

    if problems:
        print("\n[FATAL] these combos could NOT be seeded from their ground truth "
              "(their matrix.tsv rows are EMPTY -- do not run them):", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
