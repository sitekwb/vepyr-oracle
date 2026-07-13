#!/usr/bin/env python3
"""Cluster-side triage CLI: slice the 116 mismatch TSVs and probe vepyr's parquet
caches, deciding -- per drifting CSQ field -- whether the drift is a cache
artefact (something missing/extra in a data cache) or a genuine vepyr bug.

Usage:
  triage.py rollup COMBO
  triage.py transcripts COMBO FIELD CATEGORY [--all]
  triage.py in-cache CACHE_DIR ENTITY ID [ID ...]
  triage.py am-cache ID [ID ...]
  triage.py schema CACHE_DIR ENTITY

Paths (overridable via VEPYR_WORK / VEPYR_DATA, same convention as every other
bin/ script -- see bin/status.py / bin/validate.py):
  WORK = ~/vepyr/work
  DATA = ~/vepyr/data
  mismatches at $WORK/results_wgs_116/mismatches/<combo>.tsv -- 116 is the ONLY
    version with real mismatches (VEP 115 is exactly concordant across the whole
    genome; its mismatch TSVs are header-only). There is deliberately no
    --version flag here -- this tool has nothing to do against 115.

`rollup` and `transcripts` reuse oracle.triage (pure stdlib) verbatim to read
and roll up a mismatch TSV -- no pyarrow needed for that half of this script.
`in-cache` / `am-cache` / `schema` instead read vepyr's PARQUET caches, so they
import pyarrow -- LAZILY, inside the functions that need it, so `--help` (and a
plain `rollup`/`transcripts` invocation) still works even if pyarrow were ever
missing. See bin/run_wgs.py's docstring for the same lazy-import convention
around `import vepyr`.

polars 1.42.1 IS installed on the cluster but warns at import time that it is
"Missing required CPU features (avx2, fma, bmi1, bmi2, lzcnt, movbe) ... will
likely result in a crash". This script uses pyarrow (24.0.0, installed)
instead, DELIBERATELY -- do not "fix" this to use polars, and ignore any
cache's READ_WITH_POLARS.md hint.

--- The id-column problem ---
The parquet caches are organised as one directory per entity (e.g.
data/116_GRCh38_merged/transcript/*.parquet), sharded into one file per
chromosome/scaffold -- all files in a directory share one schema. That
schema's id-column NAME is not assumed here. Candidates are tried in order
(see ID_COLUMN_CANDIDATES below); the first one PRESENT in the schema is used.

If NONE of the candidates is present, this PRINTS THE SCHEMA and EXITS
NON-ZERO. It never silently falls back to an empty id-set: an empty id-set
would make `in-cache`/`am-cache` report every id "ABSENT", and that false
verdict is exactly the silent wrongness this whole triage project exists to
eliminate (imagine filing "vepyr has a genuine AlphaMissense gap" off a probe
that was actually just looking at the wrong column and finding nothing).
Run `schema` first when a cache's layout is unfamiliar.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.triage import field_rollup, load_mismatches

#: Overridable so this can be pointed at a scratch dir instead of the real
#: cluster paths -- same convention as bin/status.py / bin/validate.py.
DATA_DIR = os.environ.get("VEPYR_DATA", os.path.expanduser("~/vepyr/data"))
WORK_DIR = os.environ.get("VEPYR_WORK", os.path.expanduser("~/vepyr/work"))
MISMATCHES_DIR = os.path.join(WORK_DIR, "results_wgs_116", "mismatches")

#: Candidate id-column names, tried IN ORDER against a parquet directory's
#: schema; the first one present wins. See the module docstring's "id-column
#: problem" section for why an absent match is a hard failure, never a
#: silent empty-set.
ID_COLUMN_CANDIDATES: tuple[str, ...] = ("stable_id", "transcript_id", "feature", "id", "transcript")

#: data/plugin_cache/plugin/alphamissense/*.parquet -- fixed location, no
#: cache/entity arguments (unlike in-cache, which is generic over both).
AM_CACHE_DIR = os.path.join(DATA_DIR, "plugin_cache", "plugin", "alphamissense")


# ------------------------------------------------------------- mismatches ---

def _mismatches_path(combo: str) -> str:
    return os.path.join(MISMATCHES_DIR, f"{combo}.tsv")


def _require_mismatches(combo: str) -> str:
    path = _mismatches_path(combo)
    if not os.path.exists(path):
        print(f"[FATAL] {path} not found -- is the combo name right, and did "
              f"slurm/diff.sbatch run for it? (VEPYR_WORK={WORK_DIR})", file=sys.stderr)
        raise SystemExit(2)
    return path


def cmd_rollup(args: argparse.Namespace) -> int:
    rows = load_mismatches(_require_mismatches(args.combo))
    rollup = field_rollup(rows)
    ordered = sorted(rollup.items(), key=lambda kv: kv[1]["total"], reverse=True)

    print(f"[rollup] {args.combo}: {len(rows)} mismatch rows across {len(rollup)} field(s)")
    for field, e in ordered:
        by_cat = ", ".join(f"{k}={v}" for k, v in
                           sorted(e["by_category"].items(), key=lambda kv: -kv[1]))
        by_kind = ", ".join(f"{k}={v}" for k, v in
                            sorted(e["by_feature_kind"].items(), key=lambda kv: -kv[1]))
        print(f"  {field:24s} total={e['total']:<8d} by_category=[{by_cat}] "
             f"by_feature_kind=[{by_kind}]")
    return 0


def cmd_transcripts(args: argparse.Namespace) -> int:
    rows = load_mismatches(_require_mismatches(args.combo))
    matched = [r for r in rows if r["field"] == args.field and r["category"] == args.category]
    distinct = sorted({r["feature"] for r in matched})

    print(f"[transcripts] {args.combo} field={args.field} category={args.category}: "
         f"{len(matched)} row(s), {len(distinct)} distinct feature id(s)")
    shown = distinct if args.all else distinct[:20]
    for feat in shown:
        print(f"  {feat if feat else '(empty)'}")
    if not args.all and len(distinct) > 20:
        print(f"  ... and {len(distinct) - 20} more (use --all to print them all)")
    return 0


# --------------------------------------------------------- parquet plumbing ---

def _parquet_files(dir_path: str) -> list[str]:
    if not os.path.isdir(dir_path):
        print(f"[FATAL] {dir_path} is not a directory", file=sys.stderr)
        raise SystemExit(2)
    files = sorted(
        os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(".parquet")
    )
    if not files:
        print(f"[FATAL] no *.parquet files found under {dir_path}", file=sys.stderr)
        raise SystemExit(2)
    return files


def _schema_names(dir_path: str) -> tuple[str, list[str]]:
    """-> (file used, column names). Reads only ONE file's schema: every file in
    an entity/plugin directory is a chrom/scaffold shard of one logical table
    and shares the same schema."""
    import pyarrow.parquet as pq

    files = _parquet_files(dir_path)
    return files[0], pq.read_schema(files[0]).names


def _probe_id_column(names: list[str]) -> str | None:
    return next((c for c in ID_COLUMN_CANDIDATES if c in names), None)


def _require_id_column(dir_path: str) -> tuple[str, list[str]]:
    """-> (id_column, all parquet files in dir_path). Exits loudly -- schema
    printed to stderr, non-zero status -- if no candidate id column is present.
    See module docstring's "id-column problem" section."""
    sample_file, names = _schema_names(dir_path)
    id_col = _probe_id_column(names)
    if id_col is None:
        print(f"[FATAL] none of the candidate id columns {ID_COLUMN_CANDIDATES} "
             f"is present in the schema of {dir_path}", file=sys.stderr)
        print(f"        actual schema (from {sample_file}):", file=sys.stderr)
        for n in names:
            print(f"          {n}", file=sys.stderr)
        raise SystemExit(1)
    return id_col, _parquet_files(dir_path)


def _check_ids_present(dir_path: str, ids: list[str]) -> int:
    import pyarrow.compute as pc
    import pyarrow.dataset as ds

    id_col, files = _require_id_column(dir_path)
    dataset = ds.dataset(files, format="parquet")
    table = dataset.to_table(columns=[id_col], filter=pc.field(id_col).isin(list(ids)))
    found = set(table.column(id_col).to_pylist())

    print(f"[in-cache] {dir_path} (id column: {id_col}, {len(files)} parquet file(s))")
    n_present = sum(i in found for i in ids)
    for i in ids:
        print(f"  {i:32s} {'PRESENT' if i in found else 'ABSENT'}")
    print(f"  -- {n_present}/{len(ids)} present")
    return 0


def cmd_in_cache(args: argparse.Namespace) -> int:
    return _check_ids_present(os.path.join(DATA_DIR, args.cache, args.entity), args.ids)


def cmd_am_cache(args: argparse.Namespace) -> int:
    return _check_ids_present(AM_CACHE_DIR, args.ids)


def cmd_schema(args: argparse.Namespace) -> int:
    dir_path = os.path.join(DATA_DIR, args.cache, args.entity)
    sample_file, names = _schema_names(dir_path)
    id_col = _probe_id_column(names)

    print(f"[schema] {dir_path}")
    print(f"  sample file: {sample_file}")
    for n in names:
        print(f"  {n}" + ("  <- id column" if n == id_col else ""))
    if id_col is None:
        print(f"  (none of the candidate id columns {ID_COLUMN_CANDIDATES} is "
             f"present -- in-cache/am-cache against this directory would fail loudly)")
    return 0


# ------------------------------------------------------------------- main ---

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="triage.py", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("rollup", help="per-field mismatch rollup for one combo")
    p.add_argument("combo")
    p.set_defaults(func=cmd_rollup)

    p = sub.add_parser("transcripts",
                       help="distinct feature ids for one field/category slice")
    p.add_argument("combo")
    p.add_argument("field")
    p.add_argument("category")
    p.add_argument("--all", action="store_true", help="print every id, not just the first 20")
    p.set_defaults(func=cmd_transcripts)

    p = sub.add_parser("in-cache", help="are these ids present in a vepyr parquet cache?")
    p.add_argument("cache", help="e.g. 116_GRCh38_merged")
    p.add_argument("entity", help="e.g. transcript")
    p.add_argument("ids", nargs="+")
    p.set_defaults(func=cmd_in_cache)

    p = sub.add_parser("am-cache",
                       help="do these ids have an AlphaMissense score in the plugin cache?")
    p.add_argument("ids", nargs="+")
    p.set_defaults(func=cmd_am_cache)

    p = sub.add_parser("schema", help="print a parquet cache directory's column names")
    p.add_argument("cache")
    p.add_argument("entity")
    p.set_defaults(func=cmd_schema)

    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
