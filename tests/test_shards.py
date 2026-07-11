import pytest
from oracle.shards import plan_combo, plan_all, Level, Step, WalltimeError, write_shards, load_shards

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
