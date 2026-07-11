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

def test_numeric_within_epsilon_is_rounding():
    assert classify("gnomADe_AF", HGVS_ON, "0.1234", "0.1235", eps=1e-3) == Category.NUMERIC_TOL

def test_numeric_beyond_epsilon_is_a_real_value_diff():
    assert classify("gnomADe_AF", HGVS_ON, "0.10", "0.90", eps=1e-3) == Category.VALUE_DIFF

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
