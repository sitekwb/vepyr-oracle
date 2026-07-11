"""Drift taxonomy: exactly one category per mismatching (variant, feature, field)."""
from __future__ import annotations
from enum import StrEnum


class Category(StrEnum):
    FLAG_EXPECTED = "flag_expected"   # field not requested by the combo's flags -> not a bug
    VEP_ONLY      = "vep_only"        # VEP populated, vepyr empty -> vepyr GAP (highest priority)
    VEPYR_ONLY    = "vepyr_only"      # vepyr populated, VEP empty
    VALUE_DIFF    = "value_diff"      # both populated, differ -> the real signal
    NUMERIC_TOL   = "numeric_tol"     # both numeric, differ within epsilon -> rounding
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

DEFAULT_EPS = 1e-4


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
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _same_set_diff_order(a: str, b: str) -> bool:
    for sep in _MULTI_SEPS:
        if sep in a or sep in b:
            pa, pb = a.split(sep), b.split(sep)
            if len(pa) > 1 and sorted(pa) == sorted(pb) and pa != pb:
                return True
    return False


def classify(field: str, combo_kwargs: dict, vepyr_val: str, vep_val: str,
             eps: float = DEFAULT_EPS) -> Category | None:
    """None when the values agree (not a mismatch)."""
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
        return Category.NUMERIC_TOL if abs(fa - fb) <= eps else Category.VALUE_DIFF
    if _same_set_diff_order(vepyr_val, vep_val):
        return Category.ORDER_DIFF
    return Category.VALUE_DIFF
