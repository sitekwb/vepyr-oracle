"""The 8 in-scope combos: ONE source of truth for what each combo means.

Every combo has two sides that MUST describe the same configuration:

  * the VEP side  -- CLI flags, RECOVERED from the `##VEP-command-line=` header
    that real VEP wrote into the ground-truth VCF. Authoritative by construction.
  * the vepyr side -- Python kwargs for `vepyr.annotate()`.

The vepyr side used to be HAND-WRITTEN here, and it had drifted from the ground
truth in ways that produced phantom findings blamed on vepyr:

  * `hgvs_merged_pick` passed `pick=True`, but its GT was generated with
    `--flag_pick_allele_gene` -- an entirely different selection mode (`pick`
    DROPS the non-picked transcripts; `flag_pick_*` keeps them all and merely
    flags one).
  * all five pick-family combos silently dropped the GT's explicit, NON-DEFAULT
    `--pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length`.
    Pick order decides WHICH transcript is picked; VEP's default order leads with
    mane_select/canonical, so vepyr picked a different transcript than the GT and
    the diff reported "vepyr picked the wrong transcript" -- when in fact we had
    never told vepyr the order at all.

So nothing hand-writes a combo's semantics any more: `vep_flags_to_vepyr_kwargs()`
DERIVES the vepyr kwargs from the recovered command line, and `write_matrix()` puts
the derived kwargs and the command line they came from in the same row. The two
sides of matrix.tsv cannot disagree, because one is a pure function of the other.
"""
from __future__ import annotations
import csv, json, shlex
from collections.abc import Callable
from dataclasses import dataclass
from typing import Final
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
    """A combo is an IDENTITY, never a semantics.

    Deliberately has no `kwargs` field: what a combo asks the annotator to DO is
    recovered from its ground truth's `##VEP-command-line=` header and derived by
    `vep_flags_to_vepyr_kwargs()`. Re-adding a `kwargs` field here re-opens the exact
    drift this module was rewritten to eliminate (see the module docstring).
    """
    name: str
    cache_flavor: str          # "merged" | "refseq" -- which CACHE, not which flags
    gt115: str                 # filename inside data/ground_truth_vep/


_G = "HG002_annotated_wgs_everything_hgvs"

#: merged+refseq only. The two `everything*` combos are ensembl and are OUT OF SCOPE.
#: NOTE `hgvs_merged_pick`'s NAME is a historical misnomer: its ground truth was in
#: fact generated with --flag_pick_allele_gene. The header is the truth, not the name.
COMBOS: list[Combo] = [
    Combo("hgvs_merged",                  "merged", f"{_G}_merged.vcf"),
    Combo("hgvs_merged_am",               "merged", f"{_G}_merged_am.vcf"),
    Combo("hgvs_merged_pick",             "merged", f"{_G}_merged_pick.vcf"),
    Combo("hgvs_merged_pick_allele",      "merged", f"{_G}_merged_pick_allele.vcf"),
    Combo("hgvs_merged_pick_allele_gene", "merged", f"{_G}_merged_pick_allele_gene.vcf"),
    Combo("hgvs_merged_per_gene",         "merged", f"{_G}_merged_per_gene.vcf"),
    Combo("hgvs_merged_flag_pick_allele", "merged", f"{_G}_merged_flag_pick_allele.vcf"),
    Combo("hgvs_refseq",                  "refseq", f"{_G}_refseq.vcf"),
]


# --- VEP command line -> vepyr.annotate() kwargs -----------------------------
#
# Three EXPLICIT tables, and nothing else. Every flag a recovered command line can
# carry must appear in exactly one of them, or the parse HARD-FAILS. There is no
# "pass it through", no "looks harmless", no default branch: an unrecognised flag
# could silently change annotation semantics, and silently-wrong is the one failure
# mode this whole pipeline exists to eliminate. A loud ValueError costs a human five
# minutes; a mis-mapped flag costs a false finding that gets filed against vepyr.

#: VEP flag -> vepyr kwarg. Presence of the flag means `kwarg=True`.
_BOOL_FLAG_KWARGS: Final[dict[str, str]] = {
    "everything":            "everything",
    "hgvs":                  "hgvs",
    "pick":                  "pick",
    "pick_allele":           "pick_allele",
    "pick_allele_gene":      "pick_allele_gene",
    "per_gene":              "per_gene",
    "flag_pick":             "flag_pick",
    "flag_pick_allele":      "flag_pick_allele",
    "flag_pick_allele_gene": "flag_pick_allele_gene",
}

#: VEP flags that consume a VALUE and DO change semantics -- parsed, never guessed.
#: (`_VALUE_FLAG_KWARGS` is assembled below, after its parsers are defined.)

#: I/O and environment flags. They say WHERE data comes from and how the process is
#: run -- not what the annotation MEANS -- so vepyr gets them from the harness
#: (cache_dir, reference_fasta, output_vcf, workers) rather than from this table.
#: `--merged`/`--refseq` are NOT here: they select the cache, which the matrix already
#: tracks as `cache_flavor`, so they are cross-checked against it instead (below).
_IGNORED_BOOL_FLAGS: Final[frozenset[str]] = frozenset({
    "cache", "offline", "vcf", "force_overwrite", "no_stats",
})

#: ...same, but these consume a value token, which must be swallowed with them.
_IGNORED_VALUE_FLAGS: Final[frozenset[str]] = frozenset({
    "database", "dir_cache", "dir_plugins", "fasta", "input_file", "output_file",
    "assembly", "fork",
})

#: --merged / --refseq: the CACHE the ground truth was built against.
_CACHE_FLAVOR_FLAGS: Final[frozenset[str]] = frozenset({"merged", "refseq"})

#: the only plugin in scope. Any other plugin changes what gets annotated, and vepyr
#: would have to be told about it -- so it must fail loudly rather than be dropped.
_KNOWN_PLUGIN: Final[str] = "AlphaMissense"


def _next_token_as_value(tokens: list[str], i: int, flag: str,
                         combo: str) -> tuple[str, int]:
    """`tokens[i]` as the value of `--flag`, plus the index to resume parsing at.

    A value that looks like a flag (`--pick_order --cache`) means the command line is
    not what we think it is -- raise rather than swallow `--cache` as a pick order.
    """
    if i >= len(tokens) or tokens[i].startswith("-"):
        raise ValueError(
            f"combo {combo!r}: --{flag} expects a value but none follows it in the "
            f"recovered VEP command line")
    return tokens[i], i + 1


def _parse_pick_order(value: str, combo: str) -> str:
    """`--pick_order biotype,rank,...` -> `"biotype,rank,..."` -- the RAW comma string.

    DO NOT split this into a list. Verified against the real API on the cluster:

        inspect.signature(vepyr.annotate).parameters["pick_order"]
          annotation : str | None
          default    : None

    vepyr takes the comma string exactly as VEP's `--pick_order` spells it, so the
    value is passed through verbatim and only VALIDATED here (never reshaped).

    Splitting it is not a harmless nicety, it is the bug: `pick_order` decides WHICH
    transcript gets picked, and if vepyr were to ignore a value it cannot use rather
    than reject it, we would silently fall back to VEP's DEFAULT order (mane_select /
    canonical first) while the ground truth used the explicit order (biotype, rank
    first). Different order => different transcript => a diff that reports "vepyr
    picked the wrong transcript" for an order WE never passed it -- indistinguishable
    from a real engine bug. That phantom finding is what this whole module fixes; see
    the module docstring.
    """
    order = value.strip()
    if not [t for t in order.split(",") if t.strip()]:
        raise ValueError(f"combo {combo!r}: --pick_order has no terms: {value!r}")
    return order


def _parse_plugin(value: str, combo: str) -> _NeedsPluginCache:
    """`--plugin AlphaMissense,file=/vep/data/AlphaMissense_hg38.tsv.gz` -> sentinel.

    The plugin's `file=` path is a VEP-side path (inside the VEP container/cluster
    image); vepyr instead wants its own `plugin_cache_root`, which only the caller
    knows. So the flag maps to NEEDS_PLUGIN_CACHE, and `resolve_kwargs()` is forced
    to supply a real path before `annotate()` ever sees it.
    """
    name = value.split(",", 1)[0].strip()
    if name != _KNOWN_PLUGIN:
        raise ValueError(
            f"combo {combo!r}: unknown VEP plugin {name!r} in --plugin {value!r}. "
            f"Only {_KNOWN_PLUGIN} is mapped to a vepyr kwarg (plugin_cache_root). "
            f"Refusing to drop it: a plugin the ground truth ran and vepyr did not "
            f"produces annotations on one side only, which the diff would blame on "
            f"vepyr. Add an explicit mapping for it (after checking vepyr supports "
            f"it) if this plugin is now in scope.")
    return NEEDS_PLUGIN_CACHE


#: VEP flag -> (vepyr kwarg, parser). These consume the next token (or an `=value`).
_VALUE_FLAG_KWARGS: Final[dict[str, tuple[str, Callable[[str, str], object]]]] = {
    "pick_order": ("pick_order", _parse_pick_order),
    "plugin":     ("plugin_cache_root", _parse_plugin),
}


def vep_flags_to_vepyr_kwargs(vep_command_line: str, *, combo: str = "<unknown>",
                              cache_flavor: str | None = None) -> dict:
    """Derive `vepyr.annotate()` kwargs from a ground truth's recovered VEP flags.

    This is the ONLY way a combo's vepyr-side semantics come into existence. Hand
    writing them is what let `hgvs_merged_pick` run vepyr with `pick=True` against a
    ground truth generated with `--flag_pick_allele_gene`, and let all five
    pick-family combos run with VEP's DEFAULT pick order against a ground truth
    generated with an explicit, different one.

    Args:
        vep_command_line: the value of the GT VCF's `##VEP-command-line=` header,
            with or without its leading program token (`vep`, `./vep`, ...).
        combo: the combo's name -- used only to make failures nameable.
        cache_flavor: the combo's declared cache flavor ("merged" | "refseq"). When
            given, the command line's own `--merged`/`--refseq` MUST agree with it,
            otherwise the ground truth was built against a different cache than the
            one we annotate with, and every "missing transcript" downstream is an
            artefact of that mismatch. `None` skips the cross-check (parsing a
            command line that belongs to no combo).

    Returns:
        kwargs suitable for `vepyr.annotate(**kwargs)`, possibly carrying the
        NEEDS_PLUGIN_CACHE sentinel -- pass them through `resolve_kwargs()` first.

    Raises:
        ValueError: on ANY flag not in one of the explicit tables above, on a
            cache-flavor disagreement, or on a malformed value. Never guesses.
    """
    tokens = shlex.split(vep_command_line)
    if tokens and not tokens[0].startswith("-"):
        tokens = tokens[1:]                      # the program token: `vep`, `./vep`, ...

    kwargs: dict = {}
    declared_flavor: str | None = None
    i = 0
    while i < len(tokens):
        token = tokens[i]
        i += 1
        if not token.startswith("-"):
            raise ValueError(
                f"combo {combo!r}: stray positional argument {token!r} in the recovered "
                f"VEP command line -- it may be the value of a flag this parser does "
                f"not know consumes one, which would mean the flag itself was "
                f"misread. Refusing to guess.")

        # `--flag VALUE` and `--flag=VALUE` are both legal VEP; a flag that takes a
        # value consumes the next token, so the value can never be mistaken for a flag.
        flag, has_inline, inline_value = token.lstrip("-").partition("=")
        value = ""
        if flag in _VALUE_FLAG_KWARGS or flag in _IGNORED_VALUE_FLAGS:
            if has_inline:
                value = inline_value
            else:
                value, i = _next_token_as_value(tokens, i, flag, combo)

        if flag in _BOOL_FLAG_KWARGS:
            kwargs[_BOOL_FLAG_KWARGS[flag]] = True
        elif flag in _VALUE_FLAG_KWARGS:
            kwarg, parse = _VALUE_FLAG_KWARGS[flag]
            kwargs[kwarg] = parse(value, combo)
        elif flag in _CACHE_FLAVOR_FLAGS:
            if declared_flavor is not None and declared_flavor != flag:
                raise ValueError(
                    f"combo {combo!r}: the recovered VEP command line declares both "
                    f"--merged and --refseq; they select different caches and cannot "
                    f"both be true")
            declared_flavor = flag
        elif flag in _IGNORED_VALUE_FLAGS or flag in _IGNORED_BOOL_FLAGS:
            pass                                  # I/O or environment: not semantics
            # (its value, if any, was already consumed above)
        else:
            raise ValueError(
                f"combo {combo!r}: unrecognised VEP flag --{flag} in the recovered "
                f"command line. It is NOT ignored: an unmapped flag could silently "
                f"change annotation semantics (which transcripts, which fields), and "
                f"the resulting diff would blame vepyr for a configuration WE never "
                f"passed it. Add --{flag} to oracle.matrix's explicit tables: to "
                f"_BOOL_FLAG_KWARGS/_VALUE_FLAG_KWARGS if vepyr.annotate() has a "
                f"matching kwarg, or to _IGNORED_*_FLAGS only if it is purely about "
                f"I/O or the environment.")

    if cache_flavor is not None:
        if declared_flavor is None:
            raise ValueError(
                f"combo {combo!r}: the recovered VEP command line specifies neither "
                f"--merged nor --refseq, so its ground truth was built against the "
                f"plain Ensembl cache -- but the combo declares cache_flavor="
                f"{cache_flavor!r}. Annotating against a different cache than the "
                f"ground truth used makes every transcript difference meaningless.")
        if declared_flavor != cache_flavor:
            raise ValueError(
                f"combo {combo!r}: the recovered VEP command line says --{declared_flavor} "
                f"but the combo declares cache_flavor={cache_flavor!r}. The ground truth "
                f"was built against the {declared_flavor} cache; the matrix would have us "
                f"annotate against the {cache_flavor} one.")
    return kwargs


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


def write_matrix(path: str, vep_cmdlines_by_combo: dict[str, str]) -> None:
    """Write matrix.tsv, DERIVING each combo's vepyr kwargs from its VEP command line.

    The `vepyr_kwargs` column is a pure function of the `vep_flags` column sitting
    beside it in the same row -- that is the entire point, and it is why nothing here
    accepts hand-written kwargs.

    A combo whose command line could not be recovered gets EMPTY cells in both
    columns, never a plausible-looking `{}`: empty kwargs would silently annotate with
    vepyr's defaults and look fine. Empty makes every consumer choke (bin/run_wgs.py's
    json.loads, slurm/gen_gt116.sh's empty-flags check), which is the intended outcome
    -- bin/seed_matrix.py exits non-zero in that case anyway.

    Raises:
        ValueError: if a recovered command line cannot be mapped (see
            `vep_flags_to_vepyr_kwargs`). Better no matrix than a wrong one.
    """
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MATRIX_HEADER, delimiter="\t")
        w.writeheader()
        for c in COMBOS:
            cmdline = (vep_cmdlines_by_combo.get(c.name) or "").strip()
            kwargs_json = ""
            if cmdline:
                kwargs = vep_flags_to_vepyr_kwargs(cmdline, combo=c.name,
                                                   cache_flavor=c.cache_flavor)
                # default=_json_default: the ONLY non-JSON-native value derived kwargs
                # ever carry is NEEDS_PLUGIN_CACHE, which serialises to the "@PLUGIN@"
                # placeholder string -- see that sentinel's docstring.
                kwargs_json = json.dumps(kwargs, default=_json_default)
            w.writerow({
                "name": c.name,
                "cache_flavor": c.cache_flavor,
                "cache115": f"115_GRCh38_{c.cache_flavor}",
                "cache116": f"116_GRCh38_{c.cache_flavor}",
                "vepyr_kwargs": kwargs_json,
                "vep_flags": cmdline,
                "gt115": c.gt115,
                "gt116": f"{c.name}.vcf",
            })


def load_matrix(path: str) -> dict[str, dict]:
    with open(path) as fh:
        return {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def kwargs_from_row(row: dict) -> dict:
    """The combo's vepyr kwargs, from a matrix.tsv row loaded by `load_matrix()`.

    An EMPTY `vepyr_kwargs` cell means `write_matrix()` could not derive that combo's
    semantics from its ground truth (missing `##VEP-command-line=` header, unmappable
    flag, cache-flavor disagreement -- bin/seed_matrix.py reports which, and exits
    non-zero). Refuse it. Falling back to `{}` would annotate with vepyr's DEFAULTS
    against a ground truth generated with `--everything --hgvs [--pick_order ...]`,
    and every resulting difference would be filed as a vepyr bug.

    Raises:
        ValueError: if the row was never seeded, or its kwargs cell is not JSON.
    """
    cell = (row.get("vepyr_kwargs") or "").strip()
    name = row.get("name", "<unknown>")
    if not cell:
        raise ValueError(
            f"combo {name!r} has an EMPTY vepyr_kwargs cell in matrix.tsv: its ground "
            f"truth's VEP command line could not be recovered or mapped, so we do not "
            f"know what this combo MEANS. Refusing to annotate with default kwargs "
            f"against a ground truth that was generated with flags -- re-run "
            f"bin/seed_matrix.py and fix what it reports.")
    try:
        return json.loads(cell)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"combo {name!r} has a malformed vepyr_kwargs cell in matrix.tsv "
            f"({cell!r}): {exc}. Re-run bin/seed_matrix.py.") from exc


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
