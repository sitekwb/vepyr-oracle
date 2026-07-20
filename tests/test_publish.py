"""Testy mapowania nazw i manifestu dla publikacji GT na Drive."""
import pytest

from oracle.publish import drive_filename, parse_md5_manifest, format_md5_manifest


class TestDriveFilename:
    """Nazwa combo (jak na Obelixie) -> nazwa pliku w konwencji folderu 115.2."""

    @pytest.mark.parametrize("combo,expected", [
        ("hgvs_merged", "HG002_annotated_wgs_everything_hgvs_merged.vcf.gz"),
        ("hgvs_merged_am", "HG002_annotated_wgs_everything_hgvs_merged_am.vcf.gz"),
        ("hgvs_merged_pick", "HG002_annotated_wgs_everything_hgvs_merged_pick.vcf.gz"),
        ("hgvs_merged_pick_allele",
         "HG002_annotated_wgs_everything_hgvs_merged_pick_allele.vcf.gz"),
        ("hgvs_merged_pick_allele_gene",
         "HG002_annotated_wgs_everything_hgvs_merged_pick_allele_gene.vcf.gz"),
        ("hgvs_merged_per_gene",
         "HG002_annotated_wgs_everything_hgvs_merged_per_gene.vcf.gz"),
        ("hgvs_merged_flag_pick_allele",
         "HG002_annotated_wgs_everything_hgvs_merged_flag_pick_allele.vcf.gz"),
        ("hgvs_refseq", "HG002_annotated_wgs_everything_hgvs_refseq.vcf.gz"),
    ])
    def test_maps_every_in_scope_combo(self, combo, expected):
        assert drive_filename(combo) == expected

    def test_ensembl_combos_do_not_get_the_hgvs_prefix(self):
        """`everything` i `everything_hgvs` nie sa `hgvs_*` -- ich pliki w folderze
        115.2 nazywaja sie inaczej i zlepienie ich z prefiksem dalo by Markowi plik
        o nazwie obiecujacej cache merged, a zawierajacy ensembl."""
        assert drive_filename("everything") == "HG002_annotated_wgs_everything.vcf.gz"
        assert drive_filename("everything_hgvs") == \
            "HG002_annotated_wgs_everything_hgvs.vcf.gz"

    def test_unknown_combo_raises(self):
        """Cichy fallback nadalby plikowi prawdopodobnie wygladajaca, bledna nazwe."""
        with pytest.raises(ValueError, match="unknown combo"):
            drive_filename("hgvs_merged_typo")

    def test_uncompressed_variant(self):
        assert drive_filename("hgvs_merged", compressed=False) == \
            "HG002_annotated_wgs_everything_hgvs_merged.vcf"


class TestMd5Manifest:
    """Manifest musi byc czytelny przez `md5sum -c` po stronie Marka."""

    def test_format_uses_two_space_separator(self):
        """`md5sum -c` wymaga DWOCH spacji miedzy suma a nazwa (tryb binarny uzywa
        ' *'). Jedna spacja daje 'no properly formatted checksum lines found'."""
        out = format_md5_manifest([("d41d8cd98f00b204e9800998ecf8427e", "a.vcf.gz")])
        assert out == "d41d8cd98f00b204e9800998ecf8427e  a.vcf.gz\n"

    def test_roundtrip(self):
        entries = [("d41d8cd98f00b204e9800998ecf8427e", "a.vcf.gz"),
                   ("0cc175b9c0f1b6a831c399e269772661", "b.vcf.gz")]
        assert parse_md5_manifest(format_md5_manifest(entries)) == entries

    def test_parse_tolerates_binary_star_prefix(self):
        """`md5sum -b` pisze ' *nazwa'. Chcemy to przeczytac, nie wyprodukowac."""
        assert parse_md5_manifest("d41d8cd98f00b204e9800998ecf8427e *a.vcf.gz\n") == \
            [("d41d8cd98f00b204e9800998ecf8427e", "a.vcf.gz")]

    def test_parse_rejects_short_hash(self):
        """Obciety md5 to uszkodzony manifest -- ma krzyknac, nie przejsc dalej."""
        with pytest.raises(ValueError, match="malformed"):
            parse_md5_manifest("d41d8cd9  a.vcf.gz\n")
