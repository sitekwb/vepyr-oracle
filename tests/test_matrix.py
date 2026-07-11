import json
from oracle.matrix import COMBOS, extract_vep_command_line, write_matrix, load_matrix


def test_eight_in_scope_combos_merged_and_refseq_only():
    assert len(COMBOS) == 8
    assert {c.cache_flavor for c in COMBOS} == {"merged", "refseq"}
    assert "hgvs_merged_flag_pick_allele" in {c.name for c in COMBOS}


def test_extract_vep_command_line_from_gt_header(tmp_path):
    p = tmp_path / "gt.vcf"
    p.write_text(
        '##fileformat=VCFv4.2\n'
        "##VEP-command-line='vep --everything --hgvs --pick --cache'\n"
        '#CHROM\tPOS\n')
    assert extract_vep_command_line(str(p)) == "vep --everything --hgvs --pick --cache"


def test_extract_returns_none_when_header_absent(tmp_path):
    p = tmp_path / "gt.vcf"
    p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\n')
    assert extract_vep_command_line(str(p)) is None


def test_matrix_roundtrip(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), {c.name: "vep --everything --hgvs" for c in COMBOS})
    rows = load_matrix(str(p))
    assert len(rows) == 8
    r = rows["hgvs_merged_pick"]
    assert json.loads(r["vepyr_kwargs"])["pick"] is True
    assert r["cache115"] == "115_GRCh38_merged"
    assert r["cache116"] == "116_GRCh38_merged"
    assert rows["hgvs_refseq"]["cache116"] == "116_GRCh38_refseq"
