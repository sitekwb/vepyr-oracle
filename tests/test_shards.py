import pytest
from oracle.shards import (BAND_HI_H, DEFAULT_SAFETY_FACTOR, Level, Shard, Step, WalltimeError,
                           load_shards, plan_all, plan_combo, write_shards)

COUNTS = {str(i): 100_000 for i in range(1, 11)}   # 10 chroms x 100k = 1M variants

def test_l0_when_the_whole_wgs_fits_in_the_band():
    # 1M * 0.01s = 2.8h -> fits
    shards = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.01)
    assert len(shards) == 1
    assert shards[0].level is Level.L0
    assert shards[0].est_hours < 16

def test_l1_when_wgs_overflows_but_each_chrom_fits():
    # whole = 55.5h (too big); per chrom = 5.5h (fits)
    shards = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.2)
    assert len(shards) == 10
    assert all(s.level is Level.L1 for s in shards)
    assert all(s.est_hours <= 16 for s in shards)

def test_l2_when_a_single_chrom_overflows_the_band():
    # per chrom = 27.8h -> must chunk
    shards = plan_combo("c", Step.GT, COUNTS, sec_per_variant=1.0)
    assert all(s.level is Level.L2 for s in shards)
    assert all(s.est_hours <= 16 for s in shards)
    assert len(shards) > 10

def test_diff_step_never_emits_l2_and_raises_instead():
    # a chrom that would need L2 -> for the diff step this is a hard error,
    # because summary.merge() requires chromosome-atomic shards.
    with pytest.raises(WalltimeError, match="chromosome-atomic"):
        plan_combo("c", Step.DIFF, COUNTS, sec_per_variant=1.0)

def test_diff_step_l1_is_fine():
    shards = plan_combo("c", Step.DIFF, COUNTS, sec_per_variant=0.2)
    assert all(s.level is Level.L1 for s in shards)

def test_hard_fail_when_a_chunk_cannot_fit_under_the_24h_cap():
    with pytest.raises(WalltimeError):
        plan_combo("c", Step.GT, {"1": 100_000}, sec_per_variant=1.0, max_chunks_per_chrom=1)

def test_roundtrip(tmp_path):
    p = tmp_path / "shards.tsv"
    shards = plan_all(["a", "b"], Step.GT, COUNTS, sec_per_variant=0.2)
    write_shards(str(p), shards)
    back = load_shards(str(p))
    assert back == shards


# --- safety margin -----------------------------------------------------------
#
# The downside is ASYMMETRIC: under-provisioning costs a job KILLED at the 24h cap
# after burning up to 24h of compute (repeatedly, on a preemption-prone cluster);
# over-provisioning costs a few more array elements, which run concurrently anyway.
# So the planner plans against `sec_per_variant * safety_factor`, not the raw rate.

def test_safety_factor_splits_finer_than_the_raw_model():
    lo = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.5, safety_factor=1.0)
    hi = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.5, safety_factor=1.4)
    assert len(hi) > len(lo)
    assert all(s.est_hours <= BAND_HI_H for s in lo + hi)


def test_a_rate_that_would_have_fit_at_16h_unadjusted_now_splits():
    # 100k * 0.55s = 15.3h -> fits the band on the raw model, so L1 with no margin...
    raw = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.55, safety_factor=1.0)
    assert all(s.level is Level.L1 for s in raw)
    # ...but 15.3h * 1.4 = 21.4h, which would overrun the band and leaves almost no
    # headroom under the 24h cap -- so with the margin it must chunk to L2 instead.
    adj = plan_combo("c", Step.GT, COUNTS, sec_per_variant=0.55, safety_factor=1.4)
    assert all(s.level is Level.L2 for s in adj)
    assert all(s.est_hours <= BAND_HI_H for s in adj)


def test_est_hours_is_adjusted_and_est_hours_raw_is_the_unadjusted_model():
    (s,) = plan_combo("c", Step.GT, {"1": 100_000}, sec_per_variant=0.1, safety_factor=1.4)
    assert s.est_hours_raw == pytest.approx(100_000 * 0.1 / 3600)      # 2.78h
    assert s.est_hours == pytest.approx(s.est_hours_raw * 1.4)         # 3.89h
    assert s.safety_factor == 1.4


def test_safety_factor_defaults_to_the_conservative_bias():
    assert DEFAULT_SAFETY_FACTOR == 1.4
    (s,) = plan_combo("c", Step.GT, {"1": 100_000}, sec_per_variant=0.1)
    assert s.safety_factor == DEFAULT_SAFETY_FACTOR


def test_the_plan_records_the_margin_it_assumed(tmp_path):
    p = tmp_path / "shards.tsv"
    shards = plan_all(["a"], Step.GT, COUNTS, sec_per_variant=0.2, safety_factor=1.25)
    write_shards(str(p), shards)
    back = load_shards(str(p))
    assert back == shards
    assert all(s.safety_factor == 1.25 for s in back)   # self-documenting on disk
    assert all(s.est_hours_raw < s.est_hours for s in back)


# --- per-combo measured rates ------------------------------------------------
#
# One global sec_per_variant does not transfer across combos: hgvs_merged_am runs the
# AlphaMissense plugin, hgvs_refseq does not. A MISSING measurement must raise, never
# fall back to a global average -- that fallback is exactly how a slow combo blows the
# cap unnoticed.

def test_plan_all_accepts_a_single_rate_for_every_combo():
    shards = plan_all(["a", "b"], Step.GT, COUNTS, sec_per_variant=0.2)
    assert {s.combo for s in shards} == {"a", "b"}
    a = [s for s in shards if s.combo == "a"]
    b = [s for s in shards if s.combo == "b"]
    assert len(a) == len(b)


def test_plan_all_accepts_per_combo_rates():
    # "slow" is measured 5x more expensive per variant -> it must be sharded finer
    shards = plan_all(["fast", "slow"], Step.GT, COUNTS,
                      sec_per_variant={"fast": 0.2, "slow": 1.0})
    fast = [s for s in shards if s.combo == "fast"]
    slow = [s for s in shards if s.combo == "slow"]
    assert all(s.level is Level.L1 for s in fast)
    assert all(s.level is Level.L2 for s in slow)
    assert len(slow) > len(fast)
    assert all(s.est_hours <= BAND_HI_H for s in shards)


def test_plan_all_raises_when_a_combo_has_no_measured_rate():
    with pytest.raises(ValueError, match="no measured sec_per_variant"):
        plan_all(["fast", "slow"], Step.GT, COUNTS, sec_per_variant={"fast": 0.2})
