"""bin/seed_matrix.py: the matrix is SEEDED from the ground truth, never hand-written.

seed_matrix reads each ground-truth VCF's `##VEP-command-line=` header and writes
matrix.tsv, deriving the vepyr kwargs from the very flags it recovered. The two sides
of a row therefore cannot disagree -- which is the entire point, because they HAD
disagreed (vepyr got `pick=True` where the GT used `--flag_pick_allele_gene`, and no
pick_order at all where the GT used an explicit non-default one).

It must also PRINT both sides side by side, so a human can eyeball the agreement
rather than take it on faith, and it must exit non-zero on anything it cannot
reproduce -- a combo we cannot reproduce with VEP 116 makes the whole 116 ladder
untrustworthy.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from oracle.matrix import COMBOS, load_matrix
from tests.test_matrix import PICK_ORDER, REAL_CMDLINES

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = os.path.join(ROOT, "bin", "seed_matrix.py")


def _write_gt(gt_dir, combo, cmdline: str | None) -> None:
    """A ground-truth VCF: just enough header for extract_vep_command_line()."""
    header = f"##VEP-command-line='{cmdline}'\n" if cmdline else ""
    (gt_dir / combo.gt115).write_text(
        f"##fileformat=VCFv4.2\n{header}#CHROM\tPOS\tID\tREF\tALT\n")


@pytest.fixture
def cluster(tmp_path):
    """A fake ~/vepyr: data/ground_truth_vep_115_2_unforked/ with all 8 real GT headers."""
    gt_dir = tmp_path / "data" / "ground_truth_vep_115_2_unforked"
    gt_dir.mkdir(parents=True)
    (tmp_path / "work").mkdir()
    for combo in COMBOS:
        _write_gt(gt_dir, combo, REAL_CMDLINES[combo.name])
    return tmp_path


def _run(tmp_path) -> subprocess.CompletedProcess:
    env = {**os.environ,
           "VEPYR_DATA": str(tmp_path / "data"),
           "VEPYR_WORK": str(tmp_path / "work"),
           "PYTHONPATH": ROOT}
    return subprocess.run([sys.executable, SEED], env=env, capture_output=True,
                          text=True, timeout=30)


def test_seeds_a_matrix_whose_two_sides_cannot_disagree(cluster):
    result = _run(cluster)
    assert result.returncode == 0, result.stderr

    rows = load_matrix(str(cluster / "work" / "matrix.tsv"))
    assert len(rows) == len(COMBOS)

    # the headline drift, gone: the GT used --flag_pick_allele_gene, so vepyr gets it
    kwargs = json.loads(rows["hgvs_merged_pick"]["vepyr_kwargs"])
    assert kwargs["flag_pick_allele_gene"] is True
    assert "pick" not in kwargs

    # and every pick-family combo carries the ground truth's explicit pick order
    for name, row in rows.items():
        if "--pick_order" in row["vep_flags"]:
            assert json.loads(row["vepyr_kwargs"])["pick_order"] == PICK_ORDER, name


def test_prints_both_sides_of_every_combo_so_a_human_can_eyeball_them(cluster):
    """A matrix that is merely self-consistent is not enough -- the derived kwargs
    must be REVIEWABLE against the flags they came from, in the same output."""
    result = _run(cluster)
    assert result.returncode == 0, result.stderr

    for combo in COMBOS:
        assert combo.name in result.stdout
    # the VEP side, verbatim
    assert "--flag_pick_allele_gene" in result.stdout
    assert "--pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length" \
        in result.stdout
    # ...and the vepyr side, right next to it (pick_order verbatim -- vepyr's API
    # takes `str | None`, so what is printed is exactly what annotate() receives)
    assert "flag_pick_allele_gene=True" in result.stdout
    assert "pick_order='biotype,rank,mane_select,tsl,canonical,appris,ccds,length'" \
        in result.stdout
    # the plugin sentinel must be visible as a sentinel, never as a plausible path
    assert "plugin_cache_root=<NEEDS_PLUGIN_CACHE>" in result.stdout


def test_missing_ground_truth_exits_nonzero_and_leaves_that_row_empty(cluster):
    os.remove(cluster / "data" / "ground_truth_vep_115_2_unforked" /
              next(c for c in COMBOS if c.name == "hgvs_refseq").gt115)

    result = _run(cluster)
    assert result.returncode != 0
    assert "hgvs_refseq" in result.stderr

    # the row exists but carries NO semantics -- never a plausible-looking `{}`
    rows = load_matrix(str(cluster / "work" / "matrix.tsv"))
    assert rows["hgvs_refseq"]["vepyr_kwargs"] == ""
    assert rows["hgvs_refseq"]["vep_flags"] == ""
    assert rows["hgvs_merged"]["vepyr_kwargs"] != ""   # the others still seed


def test_header_without_a_vep_command_line_exits_nonzero(cluster):
    combo = next(c for c in COMBOS if c.name == "hgvs_merged_am")
    _write_gt(cluster / "data" / "ground_truth_vep_115_2_unforked", combo, None)

    result = _run(cluster)
    assert result.returncode != 0
    assert "hgvs_merged_am" in result.stderr


def test_an_unmappable_flag_exits_nonzero_and_never_writes_a_guessed_row(cluster):
    """The failure mode that matters: a flag we cannot map must NOT quietly become a
    row of kwargs that omits it. Loud, and the cell stays empty."""
    combo = next(c for c in COMBOS if c.name == "hgvs_merged")
    _write_gt(cluster / "data" / "ground_truth_vep_115_2_unforked", combo,
              "vep --everything --hgvs --merged --minimal --cache")

    result = _run(cluster)
    assert result.returncode != 0
    assert "--minimal" in result.stderr and "hgvs_merged" in result.stderr

    rows = load_matrix(str(cluster / "work" / "matrix.tsv"))
    assert rows["hgvs_merged"]["vepyr_kwargs"] == ""
    assert rows["hgvs_merged_pick"]["vepyr_kwargs"] != ""   # unaffected combos survive


def test_a_cache_flavor_disagreement_exits_nonzero(cluster):
    """GT built against the refseq cache, matrix says merged: caught at seed time,
    not discovered as a thousand phantom 'missing transcript' findings later."""
    combo = next(c for c in COMBOS if c.name == "hgvs_merged")
    _write_gt(cluster / "data" / "ground_truth_vep_115_2_unforked", combo,
              "vep --everything --hgvs --refseq --cache")

    result = _run(cluster)
    assert result.returncode != 0
    assert "hgvs_merged" in result.stderr
