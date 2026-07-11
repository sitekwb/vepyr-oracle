import json
import pytest
from oracle.matrix import (COMBOS, NEEDS_PLUGIN_CACHE, extract_vep_command_line,
                           load_matrix, resolve_kwargs, write_matrix)


def test_eight_in_scope_combos_merged_and_refseq_only():
    assert len(COMBOS) == 8
    assert {c.cache_flavor for c in COMBOS} == {"merged", "refseq"}
    assert "hgvs_merged_flag_pick_allele" in {c.name for c in COMBOS}


def test_extract_vep_command_line_from_gt_header(tmp_path):
    p = tmp_path / "gt.vcf"
    p.write_text(
        '##fileformat=VCFv4.2\n'
        "##VEP-command-line='vep --everything --hgvs --pick --cache'\n"
        '#CHROM\tPOS\n')
    assert extract_vep_command_line(str(p)) == "vep --everything --hgvs --pick --cache"


def test_extract_returns_none_when_header_absent(tmp_path):
    p = tmp_path / "gt.vcf"
    p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\n')
    assert extract_vep_command_line(str(p)) is None


def test_matrix_roundtrip(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), {c.name: "vep --everything --hgvs" for c in COMBOS})
    rows = load_matrix(str(p))
    assert len(rows) == 8
    r = rows["hgvs_merged_pick"]
    assert json.loads(r["vepyr_kwargs"])["pick"] is True
    assert r["cache115"] == "115_GRCh38_merged"
    assert r["cache116"] == "116_GRCh38_merged"
    assert rows["hgvs_refseq"]["cache116"] == "116_GRCh38_refseq"


# --- resolve_kwargs / NEEDS_PLUGIN_CACHE sentinel ----------------------------
#
# COMBOS used to encode hgvs_merged_am's plugin path as the magic string
# "@PLUGIN@", substituted at runtime. Forgotten substitution meant
# vepyr.annotate(plugin_cache_root="@PLUGIN@") ran with a plausible-looking string
# that could silently skip AlphaMissense annotation -- a combo that looks fine but
# validates nothing. A typed sentinel cannot be mistaken for a real path: it must
# be resolved (or must raise), never silently passed through.

def test_hgvs_merged_am_uses_the_typed_sentinel_not_a_magic_string():
    am = next(c for c in COMBOS if c.name == "hgvs_merged_am")
    assert am.kwargs["plugin_cache_root"] is NEEDS_PLUGIN_CACHE
    assert am.kwargs["plugin_cache_root"] != "@PLUGIN@"  # a real string would defeat the point


def test_resolve_kwargs_raises_when_the_sentinel_object_is_unresolved():
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs({"plugin_cache_root": NEEDS_PLUGIN_CACHE})


def test_resolve_kwargs_raises_when_the_round_tripped_string_is_unresolved():
    # after write_matrix -> load_matrix -> json.loads, the sentinel object is gone;
    # only the string "@PLUGIN@" survives the JSON round-trip. Must raise all the same.
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs({"plugin_cache_root": "@PLUGIN@"})


def test_resolve_kwargs_raises_on_empty_string_plugin_cache_root():
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs({"plugin_cache_root": NEEDS_PLUGIN_CACHE}, plugin_cache_root="")


def test_resolve_kwargs_substitutes_the_sentinel_object():
    out = resolve_kwargs({"everything": True, "plugin_cache_root": NEEDS_PLUGIN_CACHE},
                         plugin_cache_root="/data/plugin_cache")
    assert out == {"everything": True, "plugin_cache_root": "/data/plugin_cache"}


def test_resolve_kwargs_substitutes_the_round_tripped_string():
    out = resolve_kwargs({"everything": True, "plugin_cache_root": "@PLUGIN@"},
                         plugin_cache_root="/data/plugin_cache")
    assert out == {"everything": True, "plugin_cache_root": "/data/plugin_cache"}


def test_resolve_kwargs_leaves_non_plugin_combos_untouched():
    kw = {"everything": True, "hgvs": True, "pick": True}
    out = resolve_kwargs(kw)
    assert out == kw
    assert out is not kw  # must not hand back (or mutate) the caller's dict


def test_resolve_kwargs_does_not_mutate_the_input():
    kw = {"plugin_cache_root": NEEDS_PLUGIN_CACHE}
    resolve_kwargs(kw, plugin_cache_root="/data/plugin_cache")
    assert kw["plugin_cache_root"] is NEEDS_PLUGIN_CACHE


def test_write_matrix_round_trips_the_sentinel_as_the_at_plugin_string(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), {c.name: "vep --everything --hgvs" for c in COMBOS})
    rows = load_matrix(str(p))
    kwargs = json.loads(rows["hgvs_merged_am"]["vepyr_kwargs"])
    assert kwargs["plugin_cache_root"] == "@PLUGIN@"
    # and that round-tripped string still resolves/raises exactly like the sentinel
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs(kwargs)
    assert resolve_kwargs(kwargs, plugin_cache_root="/x")["plugin_cache_root"] == "/x"
