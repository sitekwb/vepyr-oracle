"""Combos ensembl-only (`everything`, `everything_hgvs`) w matrix.

Ich `##VEP-command-line=` nie ma ani `--merged`, ani `--refseq` -- cache ensembl
jest domyslny. Cross-check flavoru musi to rozumiec, inaczej oba combos wywalaja sie
przy seedowaniu i komplet 10 plikow nigdy nie powstanie.
"""
import pytest

from oracle.matrix import COMBOS, vep_flags_to_vepyr_kwargs

ENSEMBL_CMDLINE = (
    "vep --cache --offline --dir_cache /data/vep_native_cache_116 "
    "--input_file /data/HG002_normalized.vcf --output_file /out.vcf "
    "--vcf --force_overwrite --everything"
)
ENSEMBL_HGVS_CMDLINE = ENSEMBL_CMDLINE + " --hgvs"


class TestEnsemblCombosRegistered:
    def test_both_ensembl_combos_present(self):
        names = [c.name for c in COMBOS]
        assert "everything" in names
        assert "everything_hgvs" in names

    def test_ten_combos_total(self):
        assert len(COMBOS) == 10

    def test_ensembl_combos_declare_ensembl_flavor(self):
        by_name = {c.name: c for c in COMBOS}
        assert by_name["everything"].cache_flavor == "ensembl"
        assert by_name["everything_hgvs"].cache_flavor == "ensembl"

    def test_ensembl_combos_point_at_the_115_ground_truth(self):
        by_name = {c.name: c for c in COMBOS}
        assert by_name["everything"].gt115 == "HG002_annotated_wgs_everything.vcf"
        assert by_name["everything_hgvs"].gt115 == \
            "HG002_annotated_wgs_everything_hgvs.vcf"

    def test_existing_eight_combos_untouched(self):
        names = [c.name for c in COMBOS]
        for expected in ["hgvs_merged", "hgvs_merged_am", "hgvs_merged_pick",
                         "hgvs_merged_pick_allele", "hgvs_merged_pick_allele_gene",
                         "hgvs_merged_per_gene", "hgvs_merged_flag_pick_allele",
                         "hgvs_refseq"]:
            assert expected in names


class TestEnsemblFlavorCrossCheck:
    def test_no_flavor_flag_means_ensembl(self):
        kwargs = vep_flags_to_vepyr_kwargs(
            ENSEMBL_CMDLINE, combo="everything", cache_flavor="ensembl")
        assert kwargs["everything"] is True

    def test_hgvs_variant_parses(self):
        kwargs = vep_flags_to_vepyr_kwargs(
            ENSEMBL_HGVS_CMDLINE, combo="everything_hgvs", cache_flavor="ensembl")
        assert kwargs["everything"] is True
        assert kwargs["hgvs"] is True

    def test_merged_flag_with_ensembl_flavor_still_fails(self):
        with pytest.raises(ValueError):
            vep_flags_to_vepyr_kwargs(
                ENSEMBL_CMDLINE + " --merged",
                combo="everything", cache_flavor="ensembl")

    def test_missing_flavor_flag_with_merged_declared_still_fails(self):
        with pytest.raises(ValueError):
            vep_flags_to_vepyr_kwargs(
                ENSEMBL_CMDLINE, combo="hgvs_merged", cache_flavor="merged")
