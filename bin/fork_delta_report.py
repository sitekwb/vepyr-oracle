#!/usr/bin/env python3
"""Per-FIELD fork-delta table from bin/merge_summaries.py's merged summary.json.

The question this answers is NOT "how many records differ" -- that number is
compatible with both "only HGNC_ID moved" and "consequence calling moved". It is
"WHICH FIELDS moved, and how many cells of each". So the output is one row per
CSQ field: cells differing / cells compared, the percentage, and the drift
categories behind them.

Reads the LIST of per-combo whole-genome summaries that bin/merge_summaries.py
writes (each entry's `per_field[<field>]` carries `total` / `match` / `pct` /
`by_category`), and prints:

  * one table per combo, fields sorted by cells-differing descending, and
  * a verdict line naming every field that moved anywhere.

Category names are `oracle.drift.Category`'s, which are spelled for the
vepyr-vs-VEP comparison. In a fork-delta run the slots are the ones
bin/fork_delta.py documents:  `vepyr_only` = UNFORKED populated / FORKED empty
(the back-fill hypothesis's predicted signature), `vep_only` = FORKED populated /
UNFORKED empty, `value_diff` = both populated and disagreeing.

Usage:
  fork_delta_report.py SUMMARY_JSON [--markdown] [--all-fields]
"""
from __future__ import annotations

import argparse
import json
import sys

#: Renamed in the header so a reader cannot mistake a fork-delta table for a
#: vepyr-vs-VEP one. See the module docstring for the slot mapping.
_SLOT_NOTE = ("cols: unforked=vepyr slot, forked=gt slot; "
              "vepyr_only => unforked has a value, forked is EMPTY")


def _rows(summary: dict, all_fields: bool) -> list[tuple[str, int, int, float | None, dict]]:
    out = []
    for field, pf in summary.get("per_field", {}).items():
        total, match = pf["total"], pf["match"]
        diff = total - match
        if diff == 0 and not all_fields:
            continue
        out.append((field, diff, total, pf["pct"], pf.get("by_category", {})))
    out.sort(key=lambda r: (-r[1], r[0]))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="fork_delta_report.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("summary_json")
    ap.add_argument("--markdown", action="store_true", help="emit a markdown table")
    ap.add_argument("--all-fields", action="store_true",
                    help="include fields with zero differing cells (default: only movers, "
                         "plus the always-printed totals line)")
    args = ap.parse_args(argv)

    with open(args.summary_json) as fh:
        summaries = json.load(fh)
    if isinstance(summaries, dict):
        summaries = [summaries]

    moved_anywhere: dict[str, int] = {}
    for s in summaries:
        rows = _rows(s, args.all_fields)
        n_fields = len(s.get("per_field", {}))
        cells = sum(pf["total"] for pf in s.get("per_field", {}).values())
        print()
        print(f"### {s['name']}  (cache={s.get('cache')}, "
              f"chroms={len(s.get('chroms', []))}, "
              f"aligned_annotations={s.get('aligned_annotations')}, "
              f"CSQ fields={n_fields}, cells compared={cells})")
        print(f"    join_rate={s.get('join_rate')} only_unforked={s.get('only_vepyr')} "
              f"only_forked={s.get('only_gt')} "
              f"malformed=({s.get('malformed_vepyr')},{s.get('malformed_gt')}) "
              f"dup_records=({s.get('duplicate_records_vepyr')},"
              f"{s.get('duplicate_records_gt')})")
        print(f"    {_SLOT_NOTE}")
        if not rows:
            print("    NO FIELD MOVED -- every compared cell is identical between arms.")
            continue
        if args.markdown:
            print()
            print("| field | cells differing | cells compared | identical % | categories |")
            print("|---|---:|---:|---:|---|")
            for f, d, t, pct, cats in rows:
                print(f"| `{f}` | {d:,} | {t:,} | {pct} | "
                      f"{', '.join(f'{k}={v:,}' for k, v in sorted(cats.items()))} |")
        else:
            print(f"    {'field':<28} {'differ':>14} {'compared':>16} {'identical%':>11}  categories")
            for f, d, t, pct, cats in rows:
                print(f"    {f:<28} {d:>14,} {t:>16,} {pct!s:>11}  "
                      f"{', '.join(f'{k}={v:,}' for k, v in sorted(cats.items()))}")
        for f, d, *_ in rows:
            if d:
                moved_anywhere[f] = moved_anywhere.get(f, 0) + d

    print()
    if not moved_anywhere:
        print("VERDICT INPUT: no field differs in any combo.")
    else:
        print("VERDICT INPUT: fields that moved in at least one combo -- "
              + ", ".join(f"{f} ({n:,} cells)" for f, n in
                          sorted(moved_anywhere.items(), key=lambda kv: -kv[1])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
