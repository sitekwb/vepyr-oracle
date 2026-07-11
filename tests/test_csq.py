from oracle.csq import norm_chrom, csq_format_fields, parse_record

HEADER = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence ... Format: Allele|Consequence|Feature|HGVSc">\n'

def test_norm_chrom_strips_chr_prefix():
    assert norm_chrom("chr22") == "22"
    assert norm_chrom("22") == "22"

def test_csq_format_fields_parses_format_list(tmp_path):
    p = tmp_path / "h.vcf"
    p.write_text(HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    assert csq_format_fields(str(p)) == ["Allele", "Consequence", "Feature", "HGVSc"]

def test_parse_record_groups_annotations_by_feature():
    fields = ["Allele", "Consequence", "Feature", "HGVSc"]
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense|ENST1|c.1A>G,G|intron|ENST2|"
    key, feats = parse_record(line, fields, feat_i=2)
    assert key == ("22", 100, "A", "G")
    assert set(feats) == {"ENST1", "ENST2"}
    assert feats["ENST1"]["HGVSc"] == "c.1A>G"
    assert feats["ENST2"]["HGVSc"] == ""
