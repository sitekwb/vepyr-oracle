from oracle.csq import norm_chrom, csq_format_fields, parse_record

HEADER = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence ... Format: Allele|Consequence|Feature|HGVSc">\n'
FIELDS = ["Allele", "Consequence", "Feature", "HGVSc"]

def test_norm_chrom_strips_chr_prefix():
    assert norm_chrom("chr22") == "22"
    assert norm_chrom("22") == "22"

def test_csq_format_fields_parses_format_list(tmp_path):
    p = tmp_path / "h.vcf"
    p.write_text(HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    assert csq_format_fields(str(p)) == ["Allele", "Consequence", "Feature", "HGVSc"]

def test_parse_record_groups_annotations_by_allele_and_feature():
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense|ENST1|c.1A>G,G|intron|ENST2|"
    key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=0)
    assert key == ("22", 100, "A", "G")
    assert set(feats) == {("G", "ENST1"), ("G", "ENST2")}
    assert feats[("G", "ENST1")]["HGVSc"] == "c.1A>G"
    assert feats[("G", "ENST2")]["HGVSc"] == ""
    assert malformed == 0

def test_parse_record_keeps_one_entry_per_allele_on_a_multiallelic_record():
    """The join key must be (Allele, Feature): one CSQ entry per allele x transcript.

    Keyed by Feature alone, the second allele's entry overwrote the first's --
    silently halving the annotations and leaving the survivor to be diffed against
    the OTHER file's arbitrary survivor, which may be a different allele entirely.
    """
    line = ("chr22\t100\t.\tC\tT,CCGC\t.\t.\t"
            "CSQ=T|missense|ENST_A|c.10C>T,CGC|inframe_insertion|ENST_A|c.10_11insCGC")
    _key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=0)
    assert set(feats) == {("T", "ENST_A"), ("CGC", "ENST_A")}
    assert feats[("T", "ENST_A")]["Consequence"] == "missense"
    assert feats[("CGC", "ENST_A")]["Consequence"] == "inframe_insertion"
    assert malformed == 0

def test_parse_record_keeps_intergenic_entries_of_different_alleles_apart():
    """Intergenic entries have an empty Feature; (Allele, "") keeps them distinct."""
    line = "chr22\t100\t.\tC\tT,CCGC\t.\t.\tCSQ=T|intergenic||,CGC|intergenic||"
    _key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=0)
    assert set(feats) == {("T", ""), ("CGC", "")}
    assert malformed == 0

def test_parse_record_drops_and_counts_entries_too_short_for_feature():
    """An entry that cannot carry Feature has no join key -> excluded, counted.

    It must NOT get a synthesised key: the old positional fallback minted the
    same key independently in both files, so unrelated entries cross-paired.
    """
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense|ENST1|c.1A>G,G|upstream"
    _key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=0)
    assert set(feats) == {("G", "ENST1")}   # the truncated entry is gone
    assert malformed == 1                   # ...but it is not forgotten

def test_parse_record_marks_an_entry_missing_the_allele_column_malformed():
    """Half a key is not a key. Missing EITHER component -> malformed."""
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense|ENST1|c.1A>G,"
    _key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=3)
    assert set(feats) == {("c.1A>G", "ENST1")}   # allele_i=3 -> HGVSc acts as Allele
    assert malformed == 1                        # the empty entry has neither

def test_parse_record_without_a_feature_column_marks_everything_malformed():
    """No Feature column at all -> nothing is joinable. Count it; never guess."""
    line = "chr22\t100\t.\tA\tG\t.\t.\tCSQ=G|missense,G|intron"
    _key, feats, malformed = parse_record(line, ["Allele", "Consequence"],
                                          feat_i=None, allele_i=0)
    assert feats == {}
    assert malformed == 2

def test_parse_record_counts_duplicate_keys_instead_of_overwriting():
    """A repeated (Allele, Feature) is ambiguous: keep the first, count the rest.

    Overwriting would silently discard an annotation; the count makes it visible.
    """
    line = ("chr22\t100\t.\tA\tG\t.\t.\t"
            "CSQ=G|missense|ENST1|c.1A>G,G|synonymous|ENST1|c.3A>G")
    _key, feats, malformed = parse_record(line, FIELDS, feat_i=2, allele_i=0)
    assert set(feats) == {("G", "ENST1")}
    assert feats[("G", "ENST1")]["Consequence"] == "missense"   # first wins
    assert malformed == 1                                       # surplus counted
