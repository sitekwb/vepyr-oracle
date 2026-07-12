"""matrix.py: ONE source of truth for what a combo means.

The bug this file exists to prevent, in full:

The pipeline compares vepyr against real VEP. Two sides must describe the SAME
configuration -- the VEP side (CLI flags, RECOVERED from the `##VEP-command-line=`
header of the ground-truth VCF) and the vepyr side (Python kwargs). The vepyr side
used to be HAND-WRITTEN in COMBOS, and it had drifted from the ground truth:

  * `hgvs_merged_pick` passed `pick=True` -- the GT was generated with
    `--flag_pick_allele_gene`. A COMPLETELY DIFFERENT selection mode.
  * all five pick-family combos silently dropped the GT's explicit non-default
    `--pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length`.
    Pick order decides WHICH transcript is picked; VEP's default order leads with
    mane_select/canonical, so vepyr picked a different transcript than the GT and
    the diff blamed vepyr for a transcript WE never told it how to choose. A chr22
    run produced exactly such a phantom finding.

The fix: the vepyr kwargs are DERIVED from the recovered VEP command line. Nothing
hand-writes a combo's semantics any more, so the two sides cannot disagree.
"""
import json

import pytest

from oracle.matrix import (COMBOS, NEEDS_PLUGIN_CACHE, extract_vep_command_line,
                           load_matrix, resolve_kwargs, vep_flags_to_vepyr_kwargs,
                           write_matrix)

# --- fixtures: the REAL command lines, recovered from the real ground truth ---
#
# Verbatim from the `##VEP-command-line=` headers of the 115 GT VCFs (recovered by
# running bin/seed_matrix.py against them on the cluster). These are the contract:
# whatever these say, the derived vepyr kwargs must mean the same thing.

_IO = ("--cache --offline --database 0 "
       "--dir_cache /data/vep_cache --dir_plugins /data/vep_plugins "
       "--fasta /data/Homo_sapiens.GRCh38.dna.primary_assembly.fa "
       "--input_file /data/HG002.vcf.gz --output_file /data/out.vcf "
       "--vcf --force_overwrite --no_stats --assembly GRCh38 --fork 8")

_PICK_ORDER = "--pick_order biotype,rank,mane_select,tsl,canonical,appris,ccds,length"

#: what every pick-family combo's `--pick_order` must become on the vepyr side
PICK_ORDER = ["biotype", "rank", "mane_select", "tsl", "canonical", "appris",
              "ccds", "length"]

REAL_CMDLINES: dict[str, str] = {
    "hgvs_merged":
        f"vep --everything --hgvs --merged {_IO}",
    "hgvs_merged_am":
        f"vep --everything --hgvs --merged "
        f"--plugin AlphaMissense,file=/data/plugin_cache/AlphaMissense_hg38.tsv.gz {_IO}",
    # NOTE: the GT for the combo NAMED "..._pick" was generated with
    # --flag_pick_allele_gene, NOT --pick. The name lies; the header does not.
    "hgvs_merged_pick":
        f"vep --everything --hgvs --merged --flag_pick_allele_gene {_PICK_ORDER} {_IO}",
    "hgvs_merged_pick_allele":
        f"vep --everything --hgvs --merged --pick_allele {_PICK_ORDER} {_IO}",
    "hgvs_merged_pick_allele_gene":
        f"vep --everything --hgvs --merged --pick_allele_gene {_PICK_ORDER} {_IO}",
    "hgvs_merged_per_gene":
        f"vep --everything --hgvs --merged --per_gene {_PICK_ORDER} {_IO}",
    "hgvs_merged_flag_pick_allele":
        f"vep --everything --hgvs --merged --flag_pick_allele {_PICK_ORDER} {_IO}",
    "hgvs_refseq":
        f"vep --everything --hgvs --refseq {_IO}",
}

#: the kwargs each real command line MUST derive to, spelled out in full
EXPECTED_KWARGS: dict[str, dict] = {
    "hgvs_merged":
        dict(everything=True, hgvs=True),
    "hgvs_merged_am":
        dict(everything=True, hgvs=True, plugin_cache_root=NEEDS_PLUGIN_CACHE),
    "hgvs_merged_pick":
        dict(everything=True, hgvs=True, flag_pick_allele_gene=True, pick_order=PICK_ORDER),
    "hgvs_merged_pick_allele":
        dict(everything=True, hgvs=True, pick_allele=True, pick_order=PICK_ORDER),
    "hgvs_merged_pick_allele_gene":
        dict(everything=True, hgvs=True, pick_allele_gene=True, pick_order=PICK_ORDER),
    "hgvs_merged_per_gene":
        dict(everything=True, hgvs=True, per_gene=True, pick_order=PICK_ORDER),
    "hgvs_merged_flag_pick_allele":
        dict(everything=True, hgvs=True, flag_pick_allele=True, pick_order=PICK_ORDER),
    "hgvs_refseq":
        dict(everything=True, hgvs=True),
}


def _flavor(combo_name: str) -> str:
    return next(c.cache_flavor for c in COMBOS if c.name == combo_name)


def _derive(combo_name: str) -> dict:
    return vep_flags_to_vepyr_kwargs(REAL_CMDLINES[combo_name], combo=combo_name,
                                     cache_flavor=_flavor(combo_name))


# --- COMBOS: identity only, never semantics ---------------------------------

def test_eight_in_scope_combos_merged_and_refseq_only():
    assert len(COMBOS) == 8
    assert {c.cache_flavor for c in COMBOS} == {"merged", "refseq"}
    assert "hgvs_merged_flag_pick_allele" in {c.name for c in COMBOS}


def test_combos_carry_no_hand_written_kwargs():
    """The whole point: a combo is a NAME + a cache + a GT file. Its SEMANTICS come
    from the recovered VEP command line and from nowhere else. If a `kwargs` field
    ever reappears on Combo, hand-written drift becomes possible again."""
    for c in COMBOS:
        assert not hasattr(c, "kwargs"), (
            f"{c.name} carries hand-written kwargs -- that is the bug this module "
            f"was rewritten to make impossible; derive them from vep_flags instead")


def test_every_combo_in_the_fixture_table_and_vice_versa():
    """A new combo with no recovered command line would otherwise sail through."""
    assert {c.name for c in COMBOS} == set(REAL_CMDLINES) == set(EXPECTED_KWARGS)


# --- vep_flags_to_vepyr_kwargs: the real command lines ----------------------

@pytest.mark.parametrize("combo", sorted(REAL_CMDLINES))
def test_real_command_lines_derive_the_expected_kwargs(combo):
    assert _derive(combo) == EXPECTED_KWARGS[combo]


def test_hgvs_merged_pick_derives_flag_pick_allele_gene_NOT_pick():
    """THE headline bug: the combo is named "_pick" but its GT used
    --flag_pick_allele_gene. Passing vepyr `pick=True` runs a different selection
    mode entirely -- one transcript kept and the rest DROPPED, versus all transcripts
    kept with one flagged."""
    kwargs = _derive("hgvs_merged_pick")
    assert kwargs["flag_pick_allele_gene"] is True
    assert "pick" not in kwargs


@pytest.mark.parametrize("combo", [c for c, cl in REAL_CMDLINES.items()
                                   if "--pick_order" in cl])
def test_every_pick_family_combo_carries_the_ground_truths_pick_order(combo):
    """Regression: NO combo may silently lack pick_order when its command line
    specified one. A dropped pick_order means vepyr picks with VEP's DEFAULT order
    (mane_select/canonical first) while the GT picked with the explicit order
    (biotype,rank first) -- a different transcript, reported as a vepyr bug."""
    assert _derive(combo)["pick_order"] == PICK_ORDER


def test_no_recovered_flag_is_silently_dropped_from_any_combo():
    """The general form of the pick_order regression: every semantic flag present in
    a command line must show up in the derived kwargs (under its mapped name)."""
    aliases = {"flag_pick_allele_gene": "flag_pick_allele_gene", "plugin": "plugin_cache_root"}
    for combo, cmdline in REAL_CMDLINES.items():
        derived = _derive(combo)
        for semantic in ("everything", "hgvs", "pick", "pick_allele", "pick_allele_gene",
                         "per_gene", "flag_pick", "flag_pick_allele",
                         "flag_pick_allele_gene", "pick_order", "plugin"):
            if f"--{semantic} " in f"{cmdline} ":
                key = aliases.get(semantic, semantic)
                assert key in derived, (
                    f"{combo}: the GT command line specifies --{semantic} but the "
                    f"derived vepyr kwargs do not carry {key!r}: {derived}")


def test_pick_order_is_a_list_of_terms_not_a_comma_string():
    """UNVERIFIED against vepyr's real API -- a list is the natural Python shape and
    must be confirmed on the first real run; if vepyr wants the raw comma string,
    THIS is the single place to change it."""
    order = _derive("hgvs_merged_per_gene")["pick_order"]
    assert isinstance(order, list)
    assert order[0] == "biotype" and order[1] == "rank"   # NOT VEP's default order
    assert order != "biotype,rank,mane_select,tsl,canonical,appris,ccds,length"


def test_plugin_becomes_the_typed_sentinel_not_a_parsed_path():
    """The plugin's own `file=` path is a VEP-side path (inside the VEP container);
    vepyr wants its own plugin_cache_root, which only the caller knows -- so it must
    arrive as the sentinel that resolve_kwargs() is forced to fill in."""
    kwargs = _derive("hgvs_merged_am")
    assert kwargs["plugin_cache_root"] is NEEDS_PLUGIN_CACHE


def test_io_and_environment_flags_are_ignored_not_mapped():
    kwargs = _derive("hgvs_merged")
    assert kwargs == {"everything": True, "hgvs": True}


def test_merged_and_refseq_are_not_kwargs_they_are_the_cache():
    for name in ("hgvs_merged", "hgvs_refseq"):
        derived = _derive(name)
        assert "merged" not in derived and "refseq" not in derived


def test_flags_may_use_the_equals_form():
    kwargs = vep_flags_to_vepyr_kwargs(
        "vep --everything --hgvs --merged --pick_order=biotype,rank --cache",
        combo="x", cache_flavor="merged")
    assert kwargs == {"everything": True, "hgvs": True, "pick_order": ["biotype", "rank"]}


def test_command_line_without_the_leading_program_token_still_parses():
    assert vep_flags_to_vepyr_kwargs("--everything --hgvs --merged", combo="x",
                                     cache_flavor="merged") == {"everything": True,
                                                                "hgvs": True}


def test_cache_flavor_may_be_omitted_when_there_is_no_combo_to_check_against():
    assert vep_flags_to_vepyr_kwargs("vep --everything --merged") == {"everything": True}


# --- hard failures: silently-wrong is the failure mode we are eliminating ----

def test_unknown_flag_raises_naming_the_flag_and_the_combo():
    """An unrecognised flag could silently change annotation semantics (e.g.
    --minimal, --allele_number, --nearest). Refusing to guess is the whole point."""
    with pytest.raises(ValueError) as exc:
        vep_flags_to_vepyr_kwargs("vep --everything --merged --minimal",
                                  combo="hgvs_merged", cache_flavor="merged")
    assert "--minimal" in str(exc.value)
    assert "hgvs_merged" in str(exc.value)


def test_unknown_flag_raises_even_in_the_equals_form():
    with pytest.raises(ValueError, match="--custom"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged --custom=/data/x.bed",
                                  combo="c", cache_flavor="merged")


def test_refseq_command_line_with_a_merged_combo_raises():
    """--merged/--refseq select the CACHE, which the matrix tracks as cache_flavor.
    If the two disagree, the GT was made against a different cache than we annotate
    with -- every 'missing transcript' after that is an artefact."""
    with pytest.raises(ValueError) as exc:
        vep_flags_to_vepyr_kwargs(f"vep --everything --hgvs --refseq {_IO}",
                                  combo="hgvs_merged", cache_flavor="merged")
    assert "refseq" in str(exc.value) and "merged" in str(exc.value)


def test_merged_command_line_with_a_refseq_combo_raises():
    with pytest.raises(ValueError) as exc:
        vep_flags_to_vepyr_kwargs(f"vep --everything --hgvs --merged {_IO}",
                                  combo="hgvs_refseq", cache_flavor="refseq")
    assert "hgvs_refseq" in str(exc.value)


def test_command_line_declaring_both_merged_and_refseq_raises():
    with pytest.raises(ValueError, match="both"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged --refseq",
                                  combo="c", cache_flavor="merged")


def test_command_line_declaring_no_cache_flavor_raises_for_a_flavored_combo():
    """No --merged and no --refseq means the GT used the plain ENSEMBL cache, while
    the matrix says merged. Silent, and fatal to every transcript comparison."""
    with pytest.raises(ValueError, match="neither --merged nor --refseq"):
        vep_flags_to_vepyr_kwargs("vep --everything --hgvs --cache",
                                  combo="hgvs_merged", cache_flavor="merged")


def test_unknown_plugin_raises_rather_than_being_mapped_to_the_am_sentinel():
    with pytest.raises(ValueError, match="SpliceAI"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged --plugin SpliceAI,snv=/x.gz",
                                  combo="c", cache_flavor="merged")


def test_value_flag_with_a_missing_value_raises():
    with pytest.raises(ValueError, match="pick_order"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged --pick_order",
                                  combo="c", cache_flavor="merged")


def test_stray_positional_argument_raises():
    with pytest.raises(ValueError, match="stray|positional"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged /data/leftover.vcf",
                                  combo="c", cache_flavor="merged")


def test_empty_pick_order_raises():
    with pytest.raises(ValueError, match="pick_order"):
        vep_flags_to_vepyr_kwargs("vep --everything --merged --pick_order=",
                                  combo="c", cache_flavor="merged")


# --- extract_vep_command_line -----------------------------------------------

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


# --- write_matrix: the two sides cannot disagree, BY CONSTRUCTION ------------

def test_matrix_roundtrip(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), REAL_CMDLINES)
    rows = load_matrix(str(p))
    assert len(rows) == 8
    r = rows["hgvs_merged_pick"]
    assert r["cache115"] == "115_GRCh38_merged"
    assert r["cache116"] == "116_GRCh38_merged"
    assert rows["hgvs_refseq"]["cache116"] == "116_GRCh38_refseq"
    assert r["vep_flags"] == REAL_CMDLINES["hgvs_merged_pick"]


def test_matrix_vepyr_kwargs_column_is_derived_from_the_vep_flags_column(tmp_path):
    """The single property the whole rewrite buys: for EVERY row, the vepyr side is a
    pure function of the VEP side sitting right next to it in the same row."""
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), REAL_CMDLINES)
    for name, row in load_matrix(str(p)).items():
        derived = vep_flags_to_vepyr_kwargs(row["vep_flags"], combo=name,
                                            cache_flavor=row["cache_flavor"])
        # the sentinel is "@PLUGIN@" once it has crossed the JSON boundary
        expected = {k: ("@PLUGIN@" if v is NEEDS_PLUGIN_CACHE else v)
                    for k, v in derived.items()}
        assert json.loads(row["vepyr_kwargs"]) == expected


def test_written_matrix_never_hands_vepyr_pick_for_the_flag_pick_combo(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), REAL_CMDLINES)
    kwargs = json.loads(load_matrix(str(p))["hgvs_merged_pick"]["vepyr_kwargs"])
    assert kwargs["flag_pick_allele_gene"] is True and "pick" not in kwargs


def test_written_matrix_carries_pick_order_for_every_pick_family_row(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), REAL_CMDLINES)
    rows = load_matrix(str(p))
    checked = 0
    for name, row in rows.items():
        if "--pick_order" not in row["vep_flags"]:
            continue
        assert json.loads(row["vepyr_kwargs"])["pick_order"] == PICK_ORDER, name
        checked += 1
    assert checked == 5, "all five pick-family combos must have been checked"


def test_write_matrix_leaves_kwargs_empty_when_a_command_line_is_unrecoverable(tmp_path):
    """A combo whose GT header is missing has NO semantics we can trust. Writing an
    empty-but-valid `{}` would look fine and annotate with vepyr's defaults; leave the
    cell empty instead so any consumer chokes on it (and seed_matrix.py exits 1)."""
    p = tmp_path / "matrix.tsv"
    partial = {k: v for k, v in REAL_CMDLINES.items() if k != "hgvs_refseq"}
    write_matrix(str(p), partial)
    rows = load_matrix(str(p))
    assert rows["hgvs_refseq"]["vepyr_kwargs"] == ""
    assert rows["hgvs_refseq"]["vep_flags"] == ""
    assert rows["hgvs_merged"]["vepyr_kwargs"] != ""


def test_write_matrix_propagates_a_bad_command_line_instead_of_writing_it(tmp_path):
    p = tmp_path / "matrix.tsv"
    with pytest.raises(ValueError, match="hgvs_merged"):
        write_matrix(str(p), {**REAL_CMDLINES,
                              "hgvs_merged": "vep --everything --merged --minimal"})


# --- resolve_kwargs / NEEDS_PLUGIN_CACHE sentinel ----------------------------
#
# The AM combo's plugin path used to be the magic string "@PLUGIN@", substituted at
# runtime. Forgotten substitution meant vepyr.annotate(plugin_cache_root="@PLUGIN@")
# ran with a plausible-looking string that could silently skip AlphaMissense -- a
# combo that looks fine but validates nothing. A typed sentinel cannot be mistaken
# for a real path: it must be resolved (or must raise), never silently passed on.

def test_derived_am_kwargs_use_the_typed_sentinel_not_a_magic_string():
    kwargs = _derive("hgvs_merged_am")
    assert kwargs["plugin_cache_root"] is NEEDS_PLUGIN_CACHE
    assert kwargs["plugin_cache_root"] != "@PLUGIN@"  # a real string defeats the point


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


def test_resolve_kwargs_preserves_pick_order(tmp_path):
    """pick_order must survive the sentinel substitution untouched -- it is the value
    the whole drift bug turned on."""
    kw = _derive("hgvs_merged_pick")
    assert resolve_kwargs(kw)["pick_order"] == PICK_ORDER


def test_write_matrix_round_trips_the_sentinel_as_the_at_plugin_string(tmp_path):
    p = tmp_path / "matrix.tsv"
    write_matrix(str(p), REAL_CMDLINES)
    rows = load_matrix(str(p))
    kwargs = json.loads(rows["hgvs_merged_am"]["vepyr_kwargs"])
    assert kwargs["plugin_cache_root"] == "@PLUGIN@"
    # and that round-tripped string still resolves/raises exactly like the sentinel
    with pytest.raises(ValueError, match="plugin_cache_root is required"):
        resolve_kwargs(kwargs)
    assert resolve_kwargs(kwargs, plugin_cache_root="/x")["plugin_cache_root"] == "/x"
