"""oracle/status.py: per-element pipeline state + the runtime-overrun early warning."""
from __future__ import annotations

import json

from oracle.shards import Level, Shard, Step
from oracle.status import (OVERRUN_FACTOR, State, annotate_out_path, diff_summary_path,
                           element_status, scan, summarize)


def _shard(step: Step, level: Level, chrom: str, chunk: int | None = None,
          n_variants: int = 1000, est_hours: float = 10.0) -> Shard:
    return Shard(combo="hgvs_merged", step=step, level=level, chrom=chrom, chunk=chunk,
                n_variants=n_variants, est_hours=est_hours, est_hours_raw=est_hours / 1.4,
                safety_factor=1.4)


# --- path helpers ----------------------------------------------------------------

def test_annotate_out_path_l0_has_no_chrom_suffix():
    assert annotate_out_path("/r", "hgvs_merged", "all", None) == "/r/vepyr_hgvs_merged.vcf"


def test_annotate_out_path_l1_has_chrom_suffix():
    assert annotate_out_path("/r", "hgvs_merged", "22", None) == "/r/vepyr_hgvs_merged_22.vcf"


def test_annotate_out_path_l2_has_chrom_and_chunk_suffix():
    assert (annotate_out_path("/r", "hgvs_merged", "1", 3)
           == "/r/vepyr_hgvs_merged_1_c3.vcf")


def test_diff_summary_path_matches_validate_py_convention():
    assert diff_summary_path("/r", "hgvs_merged", "all") == "/r/summary_hgvs_merged.json"
    assert diff_summary_path("/r", "hgvs_merged", "22") == "/r/summary_hgvs_merged_22.json"


# --- annotate/gt shard state -------------------------------------------------------

def test_annotate_shard_missing_when_no_output(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22")
    st = element_status(s, str(tmp_path))
    assert st.state == State.MISSING
    assert st.actual_hours is None
    assert not st.overran


def test_annotate_shard_done_when_output_exists_and_nonempty(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22")
    out = tmp_path / "vepyr_hgvs_merged_22.vcf"
    out.write_text("not empty\n")
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE


def test_annotate_shard_missing_when_output_exists_but_is_empty(tmp_path):
    """Mirrors bin/run_wgs.py's own resume check (exists AND non-empty) -- an
    empty file left by a killed job must not read as done."""
    s = _shard(Step.ANNOTATE, Level.L1, "22")
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("")
    st = element_status(s, str(tmp_path))
    assert st.state == State.MISSING


def test_l2_annotate_shard_out_path_includes_chunk(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L2, "1", chunk=4)
    out = tmp_path / "vepyr_hgvs_merged_1_c4.vcf"
    out.write_text("x")
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert st.out_path == str(out)


def test_gt_shard_out_path_carries_no_vepyr_prefix(tmp_path):
    """slurm/gen_gt116.sh writes real-VEP shards, not vepyr's -- a "vepyr_"
    prefix on them would be actively misleading."""
    s = _shard(Step.GT, Level.L1, "22")
    out = tmp_path / "hgvs_merged_22.vcf"
    out.write_text("x")
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert st.out_path == str(out)
    assert not (tmp_path / "vepyr_hgvs_merged_22.vcf").exists()


# --- diff shard state ---------------------------------------------------------------

def test_diff_shard_missing_when_neither_upstream_nor_summary_exist(tmp_path):
    s = _shard(Step.DIFF, Level.L1, "22")
    st = element_status(s, str(tmp_path))
    assert st.state == State.MISSING


def test_diff_shard_annotated_when_upstream_vcf_exists_but_no_summary_yet(tmp_path):
    s = _shard(Step.DIFF, Level.L1, "22")
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    st = element_status(s, str(tmp_path))
    assert st.state == State.ANNOTATED


def test_diff_shard_done_when_summary_exists(tmp_path):
    s = _shard(Step.DIFF, Level.L1, "22")
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    (tmp_path / "summary_hgvs_merged_22.json").write_text('{"status": "ok"}')
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE


def test_diff_shard_l0_uses_no_chrom_suffix(tmp_path):
    s = _shard(Step.DIFF, Level.L0, "all")
    (tmp_path / "summary_hgvs_merged.json").write_text('{"status": "ok"}')
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert st.out_path == str(tmp_path / "summary_hgvs_merged.json")


# --- overrun early warning -----------------------------------------------------------

def test_no_overrun_when_actual_is_within_the_factor(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=10.0)
    out = tmp_path / "vepyr_hgvs_merged_22.vcf"
    out.write_text("x")
    (tmp_path / "vepyr_hgvs_merged_22.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 10.0 * 3600 * 1.1}))     # 1.1x -- under OVERRUN_FACTOR
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert not st.overran
    assert st.actual_hours == 11.0


def test_overrun_flagged_when_actual_exceeds_the_factor(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=10.0)
    out = tmp_path / "vepyr_hgvs_merged_22.vcf"
    out.write_text("x")
    (tmp_path / "vepyr_hgvs_merged_22.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 10.0 * 3600 * (OVERRUN_FACTOR + 0.01)}))
    st = element_status(s, str(tmp_path))
    assert st.overran


def test_overrun_boundary_is_strictly_greater_than_not_equal(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=10.0)
    out = tmp_path / "vepyr_hgvs_merged_22.vcf"
    out.write_text("x")
    (tmp_path / "vepyr_hgvs_merged_22.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 10.0 * 3600 * OVERRUN_FACTOR}))     # exactly the factor
    st = element_status(s, str(tmp_path))
    assert not st.overran


def test_missing_timing_file_is_not_an_overrun(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=10.0)
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert st.actual_hours is None
    assert not st.overran


def test_malformed_timing_json_does_not_crash_the_scan(tmp_path):
    s = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=10.0)
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    (tmp_path / "vepyr_hgvs_merged_22.vcf.timing.json").write_text("{not valid json")
    st = element_status(s, str(tmp_path))
    assert st.state == State.DONE
    assert st.actual_hours is None
    assert not st.overran


def test_diff_shards_never_carry_a_timing_comparison(tmp_path):
    """Diff shards have no timing.json of their own (run_wgs.py, not validate.py,
    writes timing) -- overrun detection only ever applies to gt/annotate shards."""
    s = _shard(Step.DIFF, Level.L1, "22", est_hours=1.0)
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    (tmp_path / "summary_hgvs_merged_22.json").write_text('{"status": "ok"}')
    st = element_status(s, str(tmp_path))
    assert st.actual_hours is None
    assert not st.overran


# --- scan() / summarize() -------------------------------------------------------------

def test_scan_and_summarize_count_states_and_overruns(tmp_path):
    done_ok = _shard(Step.ANNOTATE, Level.L1, "21", est_hours=10.0)
    (tmp_path / "vepyr_hgvs_merged_21.vcf").write_text("x")
    (tmp_path / "vepyr_hgvs_merged_21.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 10.0 * 3600}))

    done_overran = _shard(Step.ANNOTATE, Level.L1, "22", est_hours=1.0)
    (tmp_path / "vepyr_hgvs_merged_22.vcf").write_text("x")
    (tmp_path / "vepyr_hgvs_merged_22.vcf.timing.json").write_text(
        json.dumps({"elapsed_s": 1.0 * 3600 * (OVERRUN_FACTOR + 1)}))

    missing = _shard(Step.ANNOTATE, Level.L1, "20", est_hours=5.0)

    statuses = scan([done_ok, done_overran, missing], str(tmp_path))
    counts = summarize(statuses)
    assert counts[State.DONE] == 2
    assert counts[State.MISSING] == 1
    assert counts["overran"] == 1
