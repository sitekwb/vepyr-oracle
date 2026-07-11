"""The 8 in-scope combos: vepyr kwargs (vepyr side) + VEP CLI flags (116 GT side)."""
from __future__ import annotations
import csv, json
from dataclasses import dataclass
from .csq import open_maybe_gzip

MATRIX_HEADER = ["name", "cache_flavor", "cache115", "cache116",
                 "vepyr_kwargs", "vep_flags", "gt115", "gt116"]

#: The string a NEEDS_PLUGIN_CACHE sentinel becomes once it crosses a JSON boundary
#: (matrix.tsv's vepyr_kwargs column is JSON text). resolve_kwargs() must treat this
#: string identically to the sentinel object itself -- see the NEEDS_PLUGIN_CACHE
#: docstring for why plain string magic is exactly the bug this sentinel replaces.
_PLUGIN_PLACEHOLDER = "@PLUGIN@"


class _NeedsPluginCache:
    """Sentinel: this combo's plugin_cache_root must be resolved before annotate().

    Replaces the old magic string `"@PLUGIN@"`. A previous reviewer flagged that
    string as unsafe: if `resolve_kwargs()` is ever skipped, `vepyr.annotate(
    plugin_cache_root="@PLUGIN@")` runs with a plausible-looking path and,
    depending on how leniently vepyr treats a bad plugin dir, could silently skip
    AlphaMissense annotation entirely -- a combo that LOOKS fine but validates
    nothing. A typed sentinel with no string-like behaviour cannot be mistaken for
    a real value by anything downstream; `resolve_kwargs()` is the only legal way
    to turn it into one, and it raises rather than guesses when it cannot.

    `COMBOS` carries this sentinel as the live Python object. Once it crosses the
    JSON boundary in matrix.tsv (`write_matrix` -> `json.dumps`), it round-trips as
    the plain string `_PLUGIN_PLACEHOLDER` -- `resolve_kwargs()` recognises BOTH
    forms so a combo loaded from disk is exactly as safe as one used in-process.
    """

    def __repr__(self) -> str:
        return "<NEEDS_PLUGIN_CACHE>"

    def __eq__(self, other: object) -> bool:
        # Only ever equal to itself: NEVER to _PLUGIN_PLACEHOLDER. Equality is not
        # how resolve_kwargs() recognises the JSON round-trip form (that check is
        # explicit and separate) -- if it were, `NEEDS_PLUGIN_CACHE == "@PLUGIN@"`
        # would let the raw string sneak in as "the same thing" anywhere this
        # sentinel is compared, quietly reintroducing the magic-string hazard.
        return other is self

    def __hash__(self) -> int:
        return object.__hash__(self)


NEEDS_PLUGIN_CACHE = _NeedsPluginCache()


def _json_default(obj: object) -> str:
    """json.dumps() `default=` hook: the ONLY non-JSON-native value COMBOS ever
    carries is the NEEDS_PLUGIN_CACHE sentinel, which serialises to the string
    placeholder it represents. Anything else is a genuine bug -- re-raise loudly
    rather than let json swallow it as some best-effort str()."""
    if obj is NEEDS_PLUGIN_CACHE:
        return _PLUGIN_PLACEHOLDER
    raise TypeError(f"object of type {type(obj).__name__!r} is not JSON serializable")


@dataclass(frozen=True)
class Combo:
    name: str
    cache_flavor: str          # "merged" | "refseq"
    gt115: str                 # filename inside data/ground_truth_vep/
    kwargs: dict               # passed to vepyr.annotate(**kwargs)


_G = "HG002_annotated_wgs_everything_hgvs"

#: merged+refseq only. The two `everything*` combos are ensembl and are OUT OF SCOPE.
COMBOS: list[Combo] = [
    Combo("hgvs_merged",                  "merged", f"{_G}_merged.vcf",
          dict(everything=True, hgvs=True)),
    Combo("hgvs_merged_am",               "merged", f"{_G}_merged_am.vcf",
          dict(everything=True, hgvs=True, plugin_cache_root=NEEDS_PLUGIN_CACHE)),
    Combo("hgvs_merged_pick",             "merged", f"{_G}_merged_pick.vcf",
          dict(everything=True, hgvs=True, pick=True)),
    Combo("hgvs_merged_pick_allele",      "merged", f"{_G}_merged_pick_allele.vcf",
          dict(everything=True, hgvs=True, pick_allele=True)),
    Combo("hgvs_merged_pick_allele_gene", "merged", f"{_G}_merged_pick_allele_gene.vcf",
          dict(everything=True, hgvs=True, pick_allele_gene=True)),
    Combo("hgvs_merged_per_gene",         "merged", f"{_G}_merged_per_gene.vcf",
          dict(everything=True, hgvs=True, per_gene=True)),
    Combo("hgvs_merged_flag_pick_allele", "merged", f"{_G}_merged_flag_pick_allele.vcf",
          dict(everything=True, hgvs=True, flag_pick_allele=True)),
    Combo("hgvs_refseq",                  "refseq", f"{_G}_refseq.vcf",
          dict(everything=True, hgvs=True)),
]


def extract_vep_command_line(gt_vcf: str) -> str | None:
    """The exact real-VEP invocation recorded in the ground-truth VCF header.

    This is what makes the 116 GT provably the SAME invocation as the 115 GT,
    instead of hand-written flags that could silently drift.
    """
    with open_maybe_gzip(gt_vcf) as f:
        for line in f:
            if line.startswith("##VEP-command-line="):
                return line.split("=", 1)[1].strip().strip("'\"")
            if line.startswith("#CHROM"):
                break
    return None


def write_matrix(path: str, vep_flags_by_combo: dict[str, str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MATRIX_HEADER, delimiter="\t")
        w.writeheader()
        for c in COMBOS:
            w.writerow({
                "name": c.name,
                "cache_flavor": c.cache_flavor,
                "cache115": f"115_GRCh38_{c.cache_flavor}",
                "cache116": f"116_GRCh38_{c.cache_flavor}",
                # default=_json_default: the ONLY non-JSON-native value a combo's
                # kwargs ever carries is NEEDS_PLUGIN_CACHE, which serialises to the
                # "@PLUGIN@" placeholder string -- see that sentinel's docstring.
                "vepyr_kwargs": json.dumps(c.kwargs, default=_json_default),
                "vep_flags": vep_flags_by_combo.get(c.name, ""),
                "gt115": c.gt115,
                "gt116": f"{c.name}.vcf",
            })


def load_matrix(path: str) -> dict[str, dict]:
    with open(path) as fh:
        return {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def resolve_kwargs(combo_kwargs: dict, plugin_cache_root: str | None = None) -> dict:
    """Substitute the plugin-cache sentinel/placeholder with a real path.

    Handles BOTH forms a combo's kwargs can carry `plugin_cache_root` in:
      * the live `NEEDS_PLUGIN_CACHE` object (straight from `COMBOS`), and
      * the string `"@PLUGIN@"` (what that sentinel becomes after a JSON round-trip
        through matrix.tsv -- see `_NeedsPluginCache`'s docstring).

    Raises `ValueError` if either form is present but no `plugin_cache_root` was
    supplied to resolve it with. This must never pass a placeholder through to
    `vepyr.annotate()`: a plausible-looking bogus path can silently skip plugin
    annotation instead of failing, producing a combo that looks fine but validates
    nothing (the exact bug class this sentinel exists to make impossible).

    Combos that do not use the plugin cache at all pass through unchanged -- this
    function only ever touches the `plugin_cache_root` key.
    """
    kwargs = dict(combo_kwargs)
    needs_resolution = (kwargs.get("plugin_cache_root") is NEEDS_PLUGIN_CACHE
                        or kwargs.get("plugin_cache_root") == _PLUGIN_PLACEHOLDER)
    if needs_resolution:
        if not plugin_cache_root:
            raise ValueError(
                "plugin_cache_root is required for this combo but was not supplied "
                "(the combo's kwargs still carry the NEEDS_PLUGIN_CACHE sentinel / "
                f"{_PLUGIN_PLACEHOLDER!r} placeholder). Refusing to pass a bogus "
                "path through to vepyr.annotate(): depending on how leniently vepyr "
                "treats a bad plugin_cache_root, that could silently skip plugin "
                "annotation (e.g. AlphaMissense) instead of failing loudly."
            )
        kwargs["plugin_cache_root"] = plugin_cache_root
    return kwargs
