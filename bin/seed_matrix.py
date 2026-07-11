#!/usr/bin/env python3
"""Seed work/matrix.tsv: vepyr kwargs + the VEP flags recovered from the 115 GT headers.

Run on the cluster. Exits non-zero if any combo's flags cannot be recovered --
a missing ##VEP-command-line= means we could NOT reproduce that GT with VEP 116,
so the whole 116 ladder would be untrustworthy. Fail loudly.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.matrix import COMBOS, extract_vep_command_line, write_matrix

GTDIR = os.path.expanduser("~/vepyr/data/ground_truth_vep")
OUT = os.path.expanduser("~/vepyr/work/matrix.tsv")

flags, problems = {}, []
for c in COMBOS:
    gt = os.path.join(GTDIR, c.gt115)
    if not os.path.exists(gt):
        problems.append(f"{c.name}: ground truth missing: {gt}")
        continue
    cl = extract_vep_command_line(gt)
    if not cl:
        problems.append(f"{c.name}: no ##VEP-command-line= header in {c.gt115}")
        continue
    flags[c.name] = cl
    print(f"[seed] {c.name}\n        {cl}")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
write_matrix(OUT, flags)
print(f"[seed] wrote {OUT} ({len(flags)}/{len(COMBOS)} combos with recovered flags)")

if problems:
    print("\n[FATAL] cannot reproduce these combos with VEP 116:", file=sys.stderr)
    for p in problems:
        print(f"  - {p}", file=sys.stderr)
    sys.exit(1)
