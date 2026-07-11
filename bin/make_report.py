#!/usr/bin/env python3
"""Render the validation PDF: vepyr vs real Ensembl VEP, whole genome.

Usage: make_report.py <summary.json> <out.pdf> [title-suffix]

`summary.json` is the LIST of per-combo whole-genome summaries written by
`bin/merge_summaries.py`. This script DECIDES NOTHING: every number, ordering and
warning it prints comes from `oracle.report`, which is where the logic is tested.
This file draws boxes.

PAGE ORDER IS A SAFETY PROPERTY, NOT A LAYOUT PREFERENCE
--------------------------------------------------------
The COVERAGE / TRUST page comes BEFORE the pretty heatmaps, because the failure mode
this report exists to prevent is a reader taking `overall_pct` at face value:

    overall_pct = 100.0   join_rate = 0.843      <- real, measured, on chr22

"100.0% concordant, zero mismatches" -- while 15.7% of annotations were never compared
at all. The concordance figure is computed over the annotations the two tools BOTH
emitted; on a `--pick` combo, a variant where vepyr picks the wrong transcript
contributes to neither the numerator nor the denominator, so a `--pick` implementation
that is wrong about nearly everything still scores a perfect 100%.

So: the join rate is rendered in the SAME TABLE ROW as the headline % (never one
without the other), an untrustworthy combo is coloured red and named ON PAGE 1 with the
fraction of annotations that went uncompared, and the full warnings are spelled out on
the page immediately after. A reader who reads ONLY page 1 must not be able to walk
away with a number they think means more than it does.
"""
from __future__ import annotations

import datetime
import json
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages              # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from oracle.report import (CATEGORY_ORDER, coverage, drift_rows,  # noqa: E402
                           headline_risk, rollup_categories, true_pct)

PAGE = (11.0, 8.5)                      # US-letter landscape
RED = "#b40426"
GREEN = "#1a7f37"
GREY = "#666666"
HEADER_BLUE = "#40466e"

#: One-line gloss per drift category -- the report must not assume its reader already
#: knows the taxonomy, least of all that `flag_expected` is the one category that is
#: NOT a bug.
CATEGORY_GLOSS = {
    "vep_only": "VEP populated the field, vepyr left it EMPTY  ->  a real vepyr GAP",
    "value_diff": "both tools populated the field and DISAGREE  ->  the real signal",
    "vepyr_only": "vepyr populated the field, VEP left it empty",
    "order_diff": "same multi-value set, different order  ->  cosmetic",
    "numeric_tol": "both numeric, agree within tolerance  ->  rounding",
    "flag_expected": "the combo's flags never asked VEP for this field  ->  NOT a bug",
}


# --------------------------------------------------------------------------- utils

def _fmt_pct(v: float | None) -> str:
    return "-" if v is None else f"{v:.3f}"


def _fmt_join(v: float | None) -> str:
    return "n/a" if v is None else f"{v * 100:.1f}%"


def _fmt_chroms(chroms: list[str]) -> str:
    """['1','2','3','22'] -> '1-3,22' -- so a 2-chromosome run cannot LOOK whole-genome."""
    if not chroms:
        return "none"
    def key(c: str) -> tuple[int, str]:
        return (int(c), "") if c.isdigit() else (10**6, c)
    runs: list[list[str]] = []
    for c in sorted(chroms, key=key):
        if runs and c.isdigit() and runs[-1][-1].isdigit() and int(c) == int(runs[-1][-1]) + 1:
            runs[-1].append(c)
        else:
            runs.append([c])
    return ",".join(r[0] if len(r) == 1 else f"{r[0]}-{r[-1]}" for r in runs)


#: Vertical advance per point of font size, in figure fractions. The page is 8.5in =
#: 612pt tall, and a line of size S needs ~1.3*S pt of leading, so 1.3/612 ~= 0.00212.
#: (A FIXED advance silently overlaps the larger lines: at 0.019 an 11pt combo heading
#: printed straight through the stats line beneath it.)
_LEADING = 0.00215
_BOTTOM_MARGIN = 0.05


def _text_pages(pdf: PdfPages, title: str, lines: list[tuple[str, str, float]],
                *, subtitle: str = "") -> None:
    """Render (text, colour, size) lines as monospace pages, paginating as needed.

    Pagination is driven by ACCUMULATED HEIGHT, not by a line count: the lines are
    mixed-size, and everything on the trust page must be printed IN FULL. Silently
    truncating a coverage warning to fit a page would be its own kind of dishonesty.
    """
    pages, current, used = [], [], 0.0
    top = 0.90
    budget = top - _BOTTOM_MARGIN
    for line in lines:
        h = line[2] * _LEADING
        if used + h > budget and current:
            pages.append(current)
            current, used = [], 0.0
        current.append(line)
        used += h
    pages.append(current)

    for n, chunk in enumerate(pages):
        fig = plt.figure(figsize=PAGE)
        fig.text(0.05, 0.95, title + (f"  (cont. {n + 1})" if n else ""), size=15,
                 weight="bold")
        y = top
        if subtitle and not n:
            fig.text(0.05, 0.915, subtitle, size=9, color=GREY)
            y = 0.875
        for text, colour, size in chunk:
            fig.text(0.05, y, text, size=size, color=colour, va="top", family="monospace")
            y -= size * _LEADING
        pdf.savefig(fig)
        plt.close(fig)


def _wrap(text: str, width: int = 108, indent: str = "      ") -> list[str]:
    return textwrap.wrap(text, width=width, subsequent_indent=indent) or [""]


def _col_widths(rows: list[list[str]], header: list[str]) -> list[float]:
    """Column widths proportional to the WIDEST cell in each column.

    matplotlib's default is equal widths, which clipped `hgvs_merged_pick_allele_gene`
    -- the row's identity -- straight through its cell border while `cache` sat in a
    column three times wider than it needed. A table whose combo names are unreadable
    cannot be checked by the reader against anything.
    """
    widths = [max(len(str(r[j])) for r in [header, *rows]) for j in range(len(header))]
    total = sum(widths) or 1
    return [w / total for w in widths]


def _table(fig, rows: list[list[str]], header: list[str], *, rect: list[float],
           fontsize: float = 8.5, row_colours: dict[int, str] | None = None,
           cell_colours: dict[tuple[int, int], str] | None = None):
    ax = fig.add_axes(rect)
    ax.axis("off")
    tbl = ax.table(cellText=rows, colLabels=header, loc="center", cellLoc="center",
                   colWidths=_col_widths(rows, header))
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1, 1.55)
    for j in range(len(header)):
        tbl[0, j].set_facecolor(HEADER_BLUE)
        tbl[0, j].set_text_props(color="w", weight="bold")
    for i, colour in (row_colours or {}).items():
        for j in range(len(header)):
            tbl[i + 1, j].set_text_props(color=colour)
    for (i, j), colour in (cell_colours or {}).items():
        tbl[i + 1, j].set_text_props(color=colour, weight="bold")
    return tbl


# --------------------------------------------------------------------------- pages

def page_overview(pdf: PdfPages, summaries: list[dict], covs: list[dict],
                  suffix: str) -> None:
    """Title + the headline table. The join rate sits NEXT TO the headline %, always.

    An untrustworthy combo is red, is stamped `SEE p.2` in its own column, and is named
    again under the table with the fraction of annotations that were never compared --
    a reader who gets no further than this page still cannot mistake a conditional
    number for an unconditional one.
    """
    fig = plt.figure(figsize=PAGE)
    fig.text(0.5, 0.955, "vepyr - Validation vs Ensembl VEP", ha="center", size=20,
             weight="bold")
    fig.text(0.5, 0.917, f"HG002 GIAB benchmark - CSQ per-field concordance {suffix}",
             ha="center", size=12)
    fig.text(0.5, 0.888,
             f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} - "
             f"oracle = real Ensembl VEP ground truth",
             ha="center", size=9, color=GREY)

    header = ["combo", "cache", "overall %", "true % *", "join rate", "aligned annots",
              "status", "trust"]
    rows, row_colours, cell_colours = [], {}, {}
    for i, (s, cov) in enumerate(zip(summaries, covs)):
        ok = s.get("status") == "ok"
        rows.append([
            s.get("name", "?"),
            str(s.get("cache", "")).replace("115_GRCh38_", "").replace("116_GRCh38_", ""),
            _fmt_pct(s.get("overall_pct")) if ok else "-",
            _fmt_pct(true_pct(s)) if ok else "-",
            _fmt_join(cov["join_rate"]) if ok else "-",
            f"{cov['aligned']:,}" if ok else "-",
            s.get("status", "?"),
            "OK" if cov["trustworthy"] else "SEE p.2",
        ])
        if not cov["trustworthy"]:
            row_colours[i] = RED
            cell_colours[(i, 4)] = RED          # the join rate, in red, next to the %
            cell_colours[(i, 7)] = RED
        else:
            cell_colours[(i, 7)] = GREEN
    _table(fig, rows, header, rect=[0.04, 0.50, 0.92, 0.35])

    # The footnotes are ANCHORED to the bottom of the page, not flowed after the
    # warning block: the warning block grows with the number of untrustworthy combos,
    # and a flowed layout pushed "Chromosomes covered" (the line that stops a
    # 2-chromosome run reading as whole-genome) clean off the bottom of the page.
    chroms = sorted({c for s in summaries for c in s.get("chroms", [])},
                    key=lambda c: (int(c), "") if c.isdigit() else (10**6, c))
    fig.text(0.05, 0.075,
             f"Chromosomes covered: {_fmt_chroms(chroms)}  ({len(chroms)} contigs)",
             size=9, color=GREY)
    fig.text(0.05, 0.050,
             "* true % counts flag-semantics expected-diffs as agreement (vepyr computes "
             "HGVSc/HGVSp even when the combo's flags never asked VEP for them).",
             size=8, color=GREY)
    fig.text(0.05, 0.030,
             "  join rate = aligned / (aligned + vepyr-only + VEP-only) annotations. "
             "A concordance % without its join rate is not interpretable.",
             size=8, color=GREY)

    bad = [(s, c) for s, c in zip(summaries, covs) if not c["trustworthy"]]
    if not bad:
        fig.text(0.05, 0.43, "All combos passed the coverage checks (see next page).",
                 size=12, weight="bold", color=GREEN)
        pdf.savefig(fig)
        plt.close(fig)
        return

    y = 0.44
    fig.text(0.05, y,
             f"WARNING: {len(bad)} of {len(summaries)} combos are NOT TRUSTWORTHY.",
             size=13, weight="bold", color=RED)
    y -= 0.030
    fig.text(0.05, y,
             "Their concordance % is CONDITIONAL: it is computed only over the "
             "annotations both tools emitted. Annotations only one tool",
             size=9, color=RED)
    y -= 0.021
    fig.text(0.05, y,
             "emitted (a --pick combo picking a different transcript, a field vepyr "
             "never implemented) are not mismatches - they are INVISIBLE.",
             size=9, color=RED)
    y -= 0.028

    # The risk line is COMPUTED (oracle.report.headline_risk), never phrased here: page 1
    # rendering its own summary of the risk is how "0.0% of annotations NEVER COMPARED"
    # ended up standing in for "two whole CSQ columns were never compared".
    wrapped: list[str] = []
    for s, _cov in bad:
        wrapped += _wrap(f"- {s.get('name')}: overall {_fmt_pct(s.get('overall_pct'))}%, "
                         f"but {headline_risk(s)}", width=118, indent="   ")

    # Every untrustworthy combo is listed HERE, on page 1, whatever it takes: shrink to
    # fit rather than truncate. (8 combos x 2 wrapped lines is the worst case and fits.)
    floor = 0.115                            # keep clear of the anchored footnotes
    step = min(0.019, (y - floor) / max(len(wrapped), 1))
    size = min(8.5, max(5.5, step / 0.0021))
    for chunk in wrapped:
        fig.text(0.07, y, chunk, size=size, color=RED, family="monospace")
        y -= step
    fig.text(0.05, min(y - 0.008, 0.105),
             "Full detail on the COVERAGE / TRUST page that follows.",
             size=9, weight="bold", color=RED)
    pdf.savefig(fig)
    plt.close(fig)


def page_trust(pdf: PdfPages, summaries: list[dict], covs: list[dict]) -> None:
    """The coverage page. Comes BEFORE the heatmaps and prints every warning in full."""
    lines: list[tuple[str, str, float]] = []
    bad = [(s, c) for s, c in zip(summaries, covs) if not c["trustworthy"]]

    if not bad:
        lines += [
            ("Every combo cleared every coverage check:", GREEN, 11),
            ("", "k", 9),
            ("  - join rate >= 99% of annotations compared", GREEN, 9.5),
            ("  - no annotations present on one side only", GREEN, 9.5),
            ("  - no malformed CSQ entries, no duplicate records", GREEN, 9.5),
            ("  - no CSQ column emitted by only one tool", GREEN, 9.5),
            ("  - every combo diffed against real ground truth", GREEN, 9.5),
            ("", "k", 9),
            ("The concordance figures on page 1 can therefore be read at face value.",
             "k", 10),
        ]
    else:
        lines += [
            (f"{len(bad)} of {len(summaries)} combos published a concordance figure that "
             f"is NOT self-standing.", RED, 11),
            ("", "k", 9),
            ("Read every warning below BEFORE quoting any number from page 1.", RED, 10),
            ("", "k", 9),
        ]
        for s, cov in bad:
            lines.append((f"{s.get('name')}  [cache {s.get('cache', '?')}]", RED, 11))
            lines.append(
                (f"   overall {_fmt_pct(s.get('overall_pct'))}%   "
                 f"true {_fmt_pct(true_pct(s))}%   "
                 f"join rate {_fmt_join(cov['join_rate'])}   "
                 f"aligned {cov['aligned']:,}   uncompared {cov['not_compared']:,}",
                 GREY, 9),
            )
            for w in cov["warnings"]:
                for k, chunk in enumerate(_wrap(w)):
                    lines.append((("   !  " if k == 0 else "      ") + chunk, RED, 9))
            lines.append(("", "k", 9))

    ok_combos = [s.get("name") for s, c in zip(summaries, covs) if c["trustworthy"]]
    if bad and ok_combos:
        lines.append((f"Trustworthy combos ({len(ok_combos)}): {', '.join(ok_combos)}",
                      GREEN, 9.5))

    _text_pages(pdf, "COVERAGE / TRUST", lines,
                subtitle="What did we NOT compare? A concordance % speaks only for the "
                         "annotations that were actually joined.")


def page_drift_summary(pdf: PdfPages, summaries: list[dict]) -> None:
    """Drift categories across all combos: real gaps first, non-bugs last."""
    roll = rollup_categories(summaries)
    total = sum(roll.values())

    fig = plt.figure(figsize=PAGE)
    fig.text(0.05, 0.95, "Drift summary - all combos", size=15, weight="bold")
    fig.text(0.05, 0.915,
             f"{total:,} classified mismatches. Ordered by severity: vep_only (real "
             f"vepyr gaps) first, flag_expected (not a bug) last.",
             size=9, color=GREY)

    header = ["category", "mismatches", "% of all", "meaning"]
    rows, cell_colours = [], {}
    for i, cat in enumerate(CATEGORY_ORDER):
        n = roll[str(cat)]
        rows.append([str(cat), f"{n:,}",
                     f"{100 * n / total:.1f}%" if total else "-",
                     CATEGORY_GLOSS[str(cat)]])
        if cat == "vep_only" and n:
            cell_colours[(i, 0)] = RED
            cell_colours[(i, 1)] = RED
    tbl = _table(fig, rows, header, rect=[0.04, 0.55, 0.92, 0.30], fontsize=8.5,
                 cell_colours=cell_colours)
    for i in range(len(rows)):
        tbl[i + 1, 3].set_text_props(ha="left")

    ax = fig.add_axes([0.10, 0.10, 0.80, 0.36])
    cats = [str(c) for c in CATEGORY_ORDER]
    vals = [roll[c] for c in cats]
    colours = [RED if c == "vep_only" else (GREY if c == "flag_expected" else "#4c72b0")
               for c in cats]
    ax.barh(range(len(cats)), vals, color=colours)
    ax.set_yticks(range(len(cats)))
    ax.set_yticklabels(cats, size=9)
    ax.invert_yaxis()
    ax.set_xlabel("mismatching annotations (log scale)", size=9)
    if any(vals):
        ax.set_xscale("symlog")
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:,}", va="center", size=8)
    ax.set_title("vep_only = the vepyr gaps worth fixing", size=10)
    pdf.savefig(fig)
    plt.close(fig)


def page_heatmaps(pdf: PdfPages, ok: list[dict], chunk: int = 40) -> None:
    """Per-field concordance heatmap (the legacy RdYlGn view, kept as-is)."""
    if not ok:
        return
    fields: list[str] = []
    for s in ok:
        for f in s.get("per_field", {}):
            if f not in fields:
                fields.append(f)
    if not fields:
        return
    combos = [s["name"] for s in ok]
    M = np.full((len(combos), len(fields)), np.nan)
    for i, s in enumerate(ok):
        for j, f in enumerate(fields):
            v = s.get("per_field", {}).get(f)
            if v and v["total"]:
                M[i, j] = v["pct"]

    for cs in range(0, len(fields), chunk):
        fl = fields[cs:cs + chunk]
        sub = M[:, cs:cs + chunk]
        fig = plt.figure(figsize=(min(20, 2 + 0.32 * len(fl)), 2 + 0.5 * len(combos)))
        ax = fig.add_subplot(111)
        im = ax.imshow(sub, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100)
        ax.set_xticks(range(len(fl)))
        ax.set_xticklabels(fl, rotation=90, size=6)
        ax.set_yticks(range(len(combos)))
        ax.set_yticklabels(combos, size=8)
        ax.set_title(f"Per-field concordance % (fields {cs + 1}-{cs + len(fl)}) - "
                     f"NB: computed only over ANNOTATIONS THAT JOINED (see p.2)", size=10)
        fig.colorbar(im, ax=ax, fraction=0.02)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def page_combo_drift(pdf: PdfPages, summary: dict, cov: dict, *, top_fields: int = 24,
                     root_causes: int = 6, examples: int = 3) -> None:
    """Per-combo detail: the drift table, then root causes with example variants."""
    rows = drift_rows(summary)
    name = summary.get("name", "?")

    fig = plt.figure(figsize=PAGE)
    fig.text(0.05, 0.95, f"Drift detail - {name}", size=15, weight="bold")
    trust = ("TRUSTWORTHY" if cov["trustworthy"] else
             f"NOT TRUSTWORTHY - {cov['not_compared']:,} annotations never compared")
    fig.text(0.05, 0.915,
             f"overall {_fmt_pct(summary.get('overall_pct'))}%   "
             f"true {_fmt_pct(true_pct(summary))}%   "
             f"join rate {_fmt_join(cov['join_rate'])}   "
             f"aligned {cov['aligned']:,}   [{trust}]",
             size=9, color=GREEN if cov["trustworthy"] else RED, family="monospace")

    if not rows:
        fig.text(0.05, 0.86, "No mismatching fields: every compared annotation agreed.",
                 size=11, color=GREEN)
    else:
        header = ["field", "mismatches", "real *", "concordance %", "dominant cause",
                  "flag-expected?"]
        shown = rows[:top_fields]
        body, cell_colours = [], {}
        for i, r in enumerate(shown):
            body.append([r["field"], f"{r['mismatches']:,}", f"{r['real_mismatches']:,}",
                         _fmt_pct(r["pct"]), str(r["dominant"]),
                         "yes (not a bug)" if r["flag_expected"] else "no"])
            if r["dominant"] == "vep_only":
                cell_colours[(i, 4)] = RED
            if r["flag_expected"]:
                cell_colours[(i, 5)] = GREY
        height = min(0.62, 0.035 * len(shown) + 0.05)
        _table(fig, body, header, rect=[0.04, 0.86 - height, 0.92, height], fontsize=7.5,
               cell_colours=cell_colours)
        foot = 0.84 - height
        fig.text(0.05, foot,
                 "* real = mismatches MINUS flag-semantics expected-diffs. The table is "
                 "sorted by it, so expected-diffs cannot bury a genuine gap.",
                 size=8, color=GREY)
        if len(rows) > top_fields:
            fig.text(0.05, foot - 0.022,
                     f"... and {len(rows) - top_fields} further fields with mismatches "
                     f"(full records in the TSV below).", size=8, color=GREY)
    pdf.savefig(fig)
    plt.close(fig)

    # --- root causes, with example variants ---
    samples = summary.get("samples", {})
    lines: list[tuple[str, str, float]] = []
    for r in rows[:root_causes]:
        field = r["field"]
        tag = "  [flag-expected: NOT a bug]" if r["flag_expected"] else ""
        lines.append((f"{field}  -  {r['mismatches']:,} mismatches "
                      f"({_fmt_pct(r['pct'])}% concordant), dominant cause: "
                      f"{r['dominant']}{tag}",
                      RED if r["dominant"] == "vep_only" and not r["flag_expected"]
                      else "k", 10))
        lines.append((f"      {CATEGORY_GLOSS.get(str(r['dominant']), '')}", GREY, 8.5))
        lines.append((f"      by category: {r['by_category']}", GREY, 8.5))
        for s in samples.get(field, [])[:examples]:
            loc, vepyr_val, vep_val = s[0], s[1], s[2]
            cat = s[3] if len(s) > 3 else "?"
            lines.append((f"      {loc}   [{cat}]", "k", 8.5))
            lines.append((f"         vepyr: {vepyr_val!r}", "k", 8.5))
            lines.append((f"         VEP  : {vep_val!r}", "k", 8.5))
        lines.append(("", "k", 9))

    tsvs = summary.get("mismatches_tsvs") or [summary.get("mismatches_tsv")]
    tsvs = [t for t in tsvs if t]
    lines.append(("Every mismatching annotation (not just these samples) is recorded in:",
                  "k", 9.5))
    lines.append((f"      mismatches/{name}*.tsv"
                  + (f"   ({len(tsvs)} files, one per chromosome)" if len(tsvs) > 1 else ""),
                  "k", 9.5))
    if not cov["trustworthy"]:
        lines.append(("", "k", 9))
        lines.append(("NB: these mismatches are only what the JOIN could see. See the "
                      "COVERAGE / TRUST page:", RED, 9.5))
        for w in cov["warnings"]:
            for k, chunk in enumerate(_wrap(w)):
                lines.append((("   !  " if k == 0 else "      ") + chunk, RED, 8.5))
    _text_pages(pdf, f"Root causes - {name}", lines)


# --------------------------------------------------------------------------- main

def build_report(summaries: list[dict], out_pdf: str, suffix: str = "") -> None:
    covs = [coverage(s) for s in summaries]
    ok = [s for s in summaries if s.get("status") == "ok"]
    with PdfPages(out_pdf) as pdf:
        page_overview(pdf, summaries, covs, suffix)     # 1. headline % + join rate
        page_trust(pdf, summaries, covs)                # 2. BEFORE the pretty pictures
        page_drift_summary(pdf, summaries)              # 3. real gaps first
        page_heatmaps(pdf, ok)                          # 4. per-field heatmaps
        for s, cov in zip(summaries, covs):             # 5. per-combo drift detail
            if s.get("status") == "ok":
                page_combo_drift(pdf, s, cov)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    summary_json, out_pdf = argv[0], argv[1]
    suffix = argv[2] if len(argv) > 2 else ""

    with open(summary_json) as fh:
        summaries = json.load(fh)
    if not isinstance(summaries, list):
        print(f"[FATAL] {summary_json} is not a LIST of per-combo summaries -- run "
              f"bin/merge_summaries.py to produce one.", file=sys.stderr)
        return 2
    if not summaries:
        print(f"[FATAL] {summary_json} contains no combos. There is nothing to report on "
              f"-- refusing to emit an empty PDF that looks like a completed validation.",
              file=sys.stderr)
        return 2

    build_report(summaries, out_pdf, suffix)
    untrusted = [s.get("name") for s in summaries if not coverage(s)["trustworthy"]]
    print(f"[report] wrote {out_pdf} ({len(summaries)} combos)")
    if untrusted:
        print(f"[report] NOT TRUSTWORTHY: {', '.join(untrusted)} -- their concordance % "
              f"is conditional on annotations that were never compared (see page 2)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
