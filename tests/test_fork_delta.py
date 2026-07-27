"""Smoke tests for bin/fork_delta.py -- the UNFORKED-vs-FORKED 116 GT diff CLI.

Three properties, each of which is a way the fork-delta answer could be silently
wrong rather than loudly broken:

1. **Slot orientation.** `oracle.diff.diff_files` is written for vepyr-vs-VEP, so
   its categories are asymmetric: `vepyr_only` and `vep_only` mean opposite
   things. fork_delta.py puts UNFORKED in the vepyr slot and FORKED in the gt
   slot precisely so that the back-fill hypothesis's prediction ("forked LOSES
   the HGNC id") reads as `vepyr_only`. Swap the slots and every conclusion in
   the report inverts while every number stays the same.
2. **No soft failure on a missing arm.** bin/validate.py's `status="no_gt"`
   short-circuit is right for its job and catastrophic for this one: a combo the
   forked arm never had (the matrix grew from 8 to 10 combos between the runs)
   must not come back as "no difference".
3. **Schema compatibility with bin/merge_summaries.py.** The whole point of
   reusing diff_files is that the per-chromosome summaries fold into a
   whole-genome one through the existing, completeness-checked merge. A summary
   carrying an extra or missing key makes `oracle.summary.merge()` refuse.

No cluster, no VEP, no network: `VEPYR_DATA` / `VEPYR_WORK` point at tmp_path.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys

from oracle.matrix import MATRIX_HEADER
from oracle.summary import SHARD_KEYS, merge

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FORK_DELTA = os.path.join(ROOT, "bin", "fork_delta.py")

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr22,length=50818468>\n"
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
    'Ensembl VEP. Format: Allele|Consequence|SYMBOL|HGNC_ID|Feature">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _vcf(*csq_rows: tuple[int, str]) -> str:
    body = "".join(f"chr22\t{pos}\t.\tA\tG\t.\t.\tCSQ={csq}\n" for pos, csq in csq_rows)
    return _HEADER + body


def _setup(tmp_path, *, unforked: str | None, forked: str | None,
           combo: str = "hgvs_merged") -> dict[str, str]:
    data, work = tmp_path / "data", tmp_path / "work"
    work.mkdir(exist_ok=True)
    for arm, content in (("ground_truth_vep_116_unforked", unforked),
                         ("ground_truth_vep_116", forked)):
        shards = data / arm / "shards"
        shards.mkdir(parents=True, exist_ok=True)
        if content is not None:
            (shards / f"{combo}_22.vcf").write_text(content)
    with open(work / "matrix.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MATRIX_HEADER, delimiter="\t")
        w.writeheader()
        w.writerow({
            "name": combo, "cache_flavor": "merged",
            "cache115": "115_GRCh38_merged", "cache116": "116_GRCh38_merged",
            "vepyr_kwargs": json.dumps({"everything": True, "hgvs": True}),
            "vep_flags": "vep --everything --hgvs --merged",
            "gt115": "x.vcf", "gt116": f"{combo}.vcf",
        })
    return dict(os.environ, VEPYR_DATA=str(data), VEPYR_WORK=str(work))


def _run(env, outdir, combo="hgvs_merged", chrom="22"):
    return subprocess.run(
        [sys.executable, FORK_DELTA, "--combo", combo, "--chrom", chrom,
         "--outdir", str(outdir)],
        env=env, capture_output=True, text=True, timeout=60)


def test_forked_losing_a_value_is_reported_as_vepyr_only(tmp_path):
    """UNFORKED has the HGNC id, FORKED lost it -> `vepyr_only`, not `vep_only`.

    This is the orientation the whole report's wording depends on (property 1).
    """
    env = _setup(
        tmp_path,
        unforked=_vcf((100, "G|missense_variant|BRCA2|HGNC:1101|ENST1")),
        forked=_vcf((100, "G|missense_variant|BRCA2||ENST1")),
    )
    out = tmp_path / "out"
    res = _run(env, out)
    assert res.returncode == 0, res.stderr
    summary = json.loads((out / "summary_hgvs_merged_22.json").read_text())
    assert summary["per_field"]["HGNC_ID"] == {
        "total": 1, "match": 0, "pct": 0.0, "by_category": {"vepyr_only": 1},
    }
    # every other field is untouched -- the field-level resolution is the deliverable
    assert summary["per_field"]["Consequence"]["pct"] == 100.0
    assert summary["chrom"] == "22"


def test_a_missing_forked_arm_is_a_hard_failure_not_a_clean_bill(tmp_path):
    """A combo outside the comparable intersection must exit non-zero (property 2)."""
    env = _setup(tmp_path,
                 unforked=_vcf((100, "G|missense_variant|BRCA2|HGNC:1101|ENST1")),
                 forked=None)
    out = tmp_path / "out"
    res = _run(env, out)
    assert res.returncode == 2
    assert "comparable intersection" in res.stderr
    assert not (out / "summary_hgvs_merged_22.json").exists()


def test_summary_merges_through_merge_summaries_schema(tmp_path):
    """The summary is shard-shaped: `oracle.summary.merge()` accepts it (property 3)."""
    env = _setup(
        tmp_path,
        unforked=_vcf((100, "G|missense_variant|BRCA2|HGNC:1101|ENST1"),
                      (200, "G|synonymous_variant|TP53|HGNC:11998|ENST2")),
        forked=_vcf((100, "G|missense_variant|BRCA2||ENST1"),
                    (200, "G|synonymous_variant|TP53|HGNC:11998|ENST2")),
    )
    out = tmp_path / "out"
    assert _run(env, out).returncode == 0
    summary = json.loads((out / "summary_hgvs_merged_22.json").read_text())
    assert set(summary) == SHARD_KEYS
    merged = merge([summary], expected_chroms={"22"})
    assert merged["per_field"]["HGNC_ID"]["total"] == 2
    assert merged["per_field"]["HGNC_ID"]["match"] == 1
    assert merged["chroms"] == ["22"]


def test_resume_is_a_no_op(tmp_path):
    env = _setup(tmp_path,
                 unforked=_vcf((100, "G|missense_variant|BRCA2|HGNC:1101|ENST1")),
                 forked=_vcf((100, "G|missense_variant|BRCA2||ENST1")))
    out = tmp_path / "out"
    assert _run(env, out).returncode == 0
    (out / "summary_hgvs_merged_22.json").write_text('{"sentinel": true}')
    res = _run(env, out)
    assert res.returncode == 0
    assert "resuming (skip)" in res.stdout
    assert json.loads((out / "summary_hgvs_merged_22.json").read_text()) == {"sentinel": True}
