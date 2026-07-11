from oracle.drift import classify, Category

HGVS_OFF = {"everything": True}              # combo WITHOUT hgvs
HGVS_ON  = {"everything": True, "hgvs": True}

def test_flag_expected_wins_when_combo_lacks_the_flag():
    # vepyr computes HGVS even without --hgvs; GT leaves it empty. Not a bug.
    assert classify("HGVSc", HGVS_OFF, "c.1A>G", "") == Category.FLAG_EXPECTED

def test_same_field_is_a_real_drift_when_the_flag_IS_set():
    assert classify("HGVSc", HGVS_ON, "c.1A>G", "") == Category.VEPYR_ONLY

def test_vep_only_is_a_vepyr_gap():
    assert classify("SIFT", HGVS_ON, "", "deleterious(0.01)") == Category.VEP_ONLY

def test_value_diff_when_both_populated_and_differ():
    assert classify("Feature", HGVS_ON, "ENST1", "ENST2") == Category.VALUE_DIFF

def test_numeric_within_tolerance_is_rounding():
    assert classify("gnomADe_AF", HGVS_ON, "0.1234", "0.1235") == Category.NUMERIC_TOL

def test_numeric_beyond_tolerance_is_a_real_value_diff():
    assert classify("gnomADe_AF", HGVS_ON, "0.10", "0.90") == Category.VALUE_DIFF

def test_same_multivalue_set_different_order_is_cosmetic():
    assert classify("Consequence", HGVS_ON, "missense&splice", "splice&missense") == Category.ORDER_DIFF

def test_identical_values_are_not_classified():
    assert classify("Feature", HGVS_ON, "ENST1", "ENST1") is None


# --- IMPORTANT 5: FLAG_EXPECTED absolved REAL bugs. The premise is strictly -------
# --- ONE-DIRECTIONAL (vepyr populates a field the combo's flags never asked VEP ---
# --- for, so VEP leaves it empty). The guard was not: it fired on the FIELD alone, -
# --- before any direction check, so a vepyr REGRESSION and a WRONG value were both -
# --- filed "not a bug" and erased from the headline true %. ------------------------

def test_flag_expected_requires_the_expected_shape_vepyr_populated_vep_empty():
    assert classify("HGVSc", HGVS_OFF, "c.1A>G", "") == Category.FLAG_EXPECTED


def test_vepyr_empty_where_vep_populated_is_a_gap_not_flag_expected():
    """vepyr stopped emitting HGVSc. That is a vepyr GAP, whatever the flags say."""
    assert classify("HGVSc", HGVS_OFF, "", "c.1A>G") == Category.VEP_ONLY


def test_both_populated_and_different_is_a_value_diff_not_flag_expected():
    """vepyr emits a WRONG HGVSc. The flag cannot make a wrong value not-a-bug."""
    assert classify("HGVSc", HGVS_OFF, "c.1A>G", "c.999T>C") == Category.VALUE_DIFF


# --- IMPORTANT 9: the tolerance was ABSOLUTE (1e-4). Rare-variant allele -----------
# --- frequencies live ENTIRELY inside that epsilon, so vepyr reporting AF=0 where --
# --- gnomAD says 8e-5 -- an allele ABSENT vs an allele SEEN -- was filed "rounding". -

def test_zero_versus_a_rare_allele_frequency_is_not_rounding():
    """100% relative error. vepyr says the allele is absent; gnomAD says it is not."""
    assert classify("gnomADe_AF", HGVS_ON, "0", "0.00008") == Category.VALUE_DIFF


def test_rare_afs_an_order_of_magnitude_apart_are_not_rounding():
    """89% relative error, and both values sit inside the old absolute epsilon."""
    assert classify("gnomADe_AF", HGVS_ON, "1e-05", "9e-05") == Category.VALUE_DIFF


def test_last_digit_formatting_noise_on_a_rare_af_is_still_rounding():
    assert classify("gnomADe_AF", HGVS_ON, "0.00008", "0.000080001") == Category.NUMERIC_TOL


def test_the_default_tolerance_is_relative_not_absolute():
    """Exercises the shipped defaults -- nothing used to."""
    assert classify("CADD_PHRED", HGVS_ON, "25.100", "25.1001") == Category.NUMERIC_TOL
    assert classify("CADD_PHRED", HGVS_ON, "25.1", "25.2") == Category.VALUE_DIFF


# --- ...and _as_float was permissive enough to invent numbers out of non-numbers. --

def test_whitespace_padded_value_is_not_treated_as_a_number():
    assert classify("gnomADe_AF", HGVS_ON, " 1 ", "1") == Category.VALUE_DIFF


def test_underscored_digits_are_not_silently_parsed_as_ten():
    """float('1_0') == 10.0, so '1_0' vs '10' used to be 'rounding'."""
    assert classify("gnomADe_AF", HGVS_ON, "1_0", "10") == Category.VALUE_DIFF


def test_infinity_and_nan_are_not_numbers():
    assert classify("CADD_PHRED", HGVS_ON, "inf", "1e400") == Category.VALUE_DIFF
    assert classify("CADD_PHRED", HGVS_ON, "nan", "0") == Category.VALUE_DIFF


# --- MINOR 10: VEP writes "-" for the null case in Amino_acids/Codons/ -----------
# --- Protein_position. If vepyr writes "" there, every such annotation minted a ---
# --- VEP_ONLY -- our HIGHEST-priority category -- burying the real gaps in noise. --

def test_vep_dash_against_vepyr_empty_is_not_a_gap():
    assert classify("Amino_acids", HGVS_ON, "", "-") is None


def test_vepyr_dash_against_vep_empty_is_not_a_drift():
    assert classify("Codons", HGVS_ON, "-", "") is None


def test_dot_is_also_an_empty_sentinel():
    assert classify("Protein_position", HGVS_ON, ".", "") is None
    assert classify("Protein_position", HGVS_ON, "-", ".") is None


def test_a_real_gap_is_still_a_gap_after_sentinel_normalisation():
    assert classify("Amino_acids", HGVS_ON, "-", "M/T") == Category.VEP_ONLY
    assert classify("Amino_acids", HGVS_ON, "M/T", "-") == Category.VEPYR_ONLY
