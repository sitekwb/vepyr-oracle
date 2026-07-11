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
    key, feats, malformed = parse_record(line, fields, feat_i=2)
    assert key == ("22", 100, "A", "G")
    assert set(feats) == {"ENST1", "ENST2"}
    assert feats["ENST1"]["HGVSc"] == "c.1A>G"
    assert feats["ENST2"]["HGVSc"] == ""
    assert malformed == 0

def test_parse_record_drops_and_counts_entries_too_short_for_feature():
    """An entry that cannot carry Feature has no join key -> excluded, counted.

    It must NOT get a synthesised key: the old positional fallback minted the
    same key independently in both files, so unrelated entries cross-paired.
    """
    fields = ["Allele", "Consequence", "Feature", "HGVSc"]
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense|ENST1|c.1A>G,G|upstream"
    key, feats, malformed = parse_record(line, fields, feat_i=2)
    assert set(feats) == {"ENST1"}          # the truncated entry is gone
    assert malformed == 1                   # ...but it is not forgotten
    assert not any(k.isdigit() for k in feats)   # no positional fallback key

def test_parse_record_without_a_feature_column_marks_everything_malformed():
    """No Feature column at all -> nothing is joinable. Count it; never guess."""
    fields = ["Allele", "Consequence"]
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense,G|intron"
    key, feats, malformed = parse_record(line, fields, feat_i=None)
    assert feats == {}
    assert malformed == 2
