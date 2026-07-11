"""Drift taxonomy: exactly one category per mismatching (variant, feature, field)."""
from __future__ import annotations
import math
import re
from enum import StrEnum


class Category(StrEnum):
    FLAG_EXPECTED = "flag_expected"   # field not requested by the combo's flags -> not a bug
    VEP_ONLY      = "vep_only"        # VEP populated, vepyr empty -> vepyr GAP (highest priority)
    VEPYR_ONLY    = "vepyr_only"      # vepyr populated, VEP empty
    VALUE_DIFF    = "value_diff"      # both populated, differ -> the real signal
    NUMERIC_TOL   = "numeric_tol"     # both numeric, agree within tolerance -> rounding
    ORDER_DIFF    = "order_diff"      # same multi-value set, different order -> cosmetic


#: field -> the vepyr kwarg that must be set for VEP to populate it.
#: Without the kwarg the GT leaves the field empty while vepyr still computes it.
FIELD_REQUIRES_FLAG: dict[str, str] = {
    "HGVSc": "hgvs",
    "HGVSp": "hgvs",
    "HGVS_OFFSET": "hgvs",
}

#: separators VEP uses inside a single CSQ field for multi-values
_MULTI_SEPS = ("&", ",")

#: Numeric agreement is RELATIVE, with NO absolute floor.
#:
#: `DEFAULT_REL_TOL = 1e-3` (0.1%): VEP prints numeric CSQ values to a handful of
#: significant figures (gnomAD AFs ~4-6 sf, CADD/SIFT/PolyPhen 2-3 sf), so two
#: implementations reading the same source data can differ in the last printed digit
#: and no further. 0.1% sits comfortably above last-digit formatting noise at 4+ sf
#: and far below any difference a scientist would call the same number.
#:
#: `DEFAULT_ABS_TOL = 0.0`: an absolute floor is precisely what made the old
#: `eps = 1e-4` unsafe. ANY floor F silently blesses "vepyr says 0, gnomAD says F/2"
#: -- and rare-variant allele frequencies live ENTIRELY inside 1e-4:
#:     vepyr=0      VEP=0.00008  -> "rounding" (100% relative error)
#:     vepyr=1e-05  VEP=9e-05    -> "rounding" ( 89% relative error)
#: AF=0 (allele absent) vs AF=8e-5 (allele observed) is a categorical difference, not
#: a rounding artefact. With no floor, the only pairs tolerated are those within 0.1%
#: OF EACH OTHER, which can never include zero-vs-nonzero. Values that are exactly
#: equal never reach here (string equality short-circuits), and math.isclose still
#: closes -0.0 against 0.0.
DEFAULT_REL_TOL = 1e-3
DEFAULT_ABS_TOL = 0.0

#: A CSQ value is a number only if it LOOKS like one. `float()` alone accepts
#: `" 1 "` (surrounding whitespace), `"1_0"` (PEP 515 underscores -> 10.0!), `"inf"`
#: and `"nan"`, so `"1_0"` vs `"10"` compared equal and was filed as rounding.
_NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def is_flag_expected(field: str, combo_kwargs: dict, vepyr_val: str, vep_val: str) -> bool:
    """True only for the ONE shape the flag premise actually explains.

    The premise is strictly one-directional: the combo's flags never asked VEP for
    this field, so VEP leaves it EMPTY while vepyr computes it anyway. Anything else
    is a genuine finding, and the guard used to fire on the FIELD ALONE -- before any
    direction or emptiness check -- so on a combo without `hgvs`:

        vepyr=''       VEP='c.1A>G'    -> flag_expected   (a real vepyr GAP)
        vepyr='c.1A>G' VEP='c.999T>C'  -> flag_expected   (a WRONG HGVS)

    were both filed "not a bug" and erased from the headline true %. A vepyr
    regression that stopped emitting HGVSc, or one that emitted the wrong HGVSc,
    was invisible.
    """
    req = FIELD_REQUIRES_FLAG.get(field)
    if req is None or combo_kwargs.get(req):
        return False
    return bool(vepyr_val) and not vep_val


def _as_float(v: str) -> float | None:
    """float(v), but only for strings that actually look like a decimal number."""
    if not isinstance(v, str) or not _NUMERIC_RE.match(v):
        return None
    try:
        return float(v)
    except ValueError:            # pragma: no cover -- the regex already excludes these
        return None


#: All the ways the two tools spell "this field has no value here". VEP writes "-"
#: for the null case in Amino_acids / Codons / Protein_position (and "." is the VCF
#: null); vepyr writes "". Left un-normalised, every such annotation minted a
#: VEP_ONLY -- our HIGHEST-priority category, the one that means "vepyr GAP" -- and
#: buried the real gaps under thousands of spelling differences.
EMPTY_SENTINELS = frozenset({"", "-", "."})


def _norm_empty(v: str) -> str:
    return "" if v in EMPTY_SENTINELS else v


def _same_set_diff_order(a: str, b: str) -> bool:
    for sep in _MULTI_SEPS:
        if sep in a or sep in b:
            pa, pb = a.split(sep), b.split(sep)
            if len(pa) > 1 and sorted(pa) == sorted(pb) and pa != pb:
                return True
    return False


def classify(field: str, combo_kwargs: dict, vepyr_val: str, vep_val: str,
             rel_tol: float = DEFAULT_REL_TOL,
             abs_tol: float = DEFAULT_ABS_TOL) -> Category | None:
    """None when the values agree (not a mismatch).

    Both values are first normalised through EMPTY_SENTINELS, so a field VEP spells
    "-" and vepyr spells "" is a MATCH, not a fabricated gap. (The TSV still records
    the raw strings -- normalisation decides the category, never what the human sees.)
    """
    vepyr_val, vep_val = _norm_empty(vepyr_val), _norm_empty(vep_val)
    if vepyr_val == vep_val:
        return None
    if is_flag_expected(field, combo_kwargs, vepyr_val, vep_val):
        return Category.FLAG_EXPECTED
    if not vepyr_val and vep_val:
        return Category.VEP_ONLY
    if vepyr_val and not vep_val:
        return Category.VEPYR_ONLY
    fa, fb = _as_float(vepyr_val), _as_float(vep_val)
    if fa is not None and fb is not None:
        return (Category.NUMERIC_TOL
                if math.isclose(fa, fb, rel_tol=rel_tol, abs_tol=abs_tol)
                else Category.VALUE_DIFF)
    if _same_set_diff_order(vepyr_val, vep_val):
        return Category.ORDER_DIFF
    return Category.VALUE_DIFF
