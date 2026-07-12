"""Guard: no ACTIVE pipeline code path may point the ANNOTATION input at the
raw, un-normalized HG002 benchmark VCF.

CONTEXT (the CRITICAL bug this guards against): every ground-truth VCF was
built with `--input_file .../HG002_normalized.vcf`, i.e. the benchmark AFTER
`bcftools norm -m -any` split every multi-allelic site into one row per ALT
(see slurm/prep_input.sh's header comment for the record-count proof:
4,048,342 raw + 47,781 multi-allelic splits = 4,096,123 ground-truth records).
Feeding vepyr (or a fresh real-VEP 116 run) the RAW, un-split benchmark
instead would still run to completion -- it would just silently break the
(chrom,pos,ref,alt) join on all 47,781 multi-allelic sites: a joint row
`C -> T,CCGC` never matches its split counterparts `C -> T` / `C -> CCGC`.
Those variants would then silently drop out of every downstream comparison --
exactly the ones multi-ALT handling is weakest on, and exactly the ones this
whole validation pipeline most needs to measure.

The raw benchmark's filename is legitimately allowed to appear in exactly ONE
place: slurm/prep_input.sh, where it is the *source* `bcftools norm` reads
FROM to produce HG002_normalized.vcf.gz in the first place. Anywhere else it
appears in the active pipeline is either a fresh copy of the same bug, or a
deliberate change that must update this allowlist explicitly -- never
silently.

`legacy/` is excluded from this scan on purpose: it is frozen, dead code
(explicitly excluded from deploy.sh's rsync to the cluster, and not exercised
by any other test in this suite) -- not part of the pipeline this guard
protects, and rewriting it is out of scope for this fix.
"""
from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RAW_BENCHMARK = "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
NORMALIZED_INPUT = "HG002_normalized.vcf.gz"

#: Files allowed to mention the raw benchmark's filename -- see module
#: docstring. Anything else that mentions it is exactly the bug this guard
#: exists to catch.
_ALLOWLIST = {os.path.join(ROOT, "slurm", "prep_input.sh")}

#: Directories that make up the ACTIVE pipeline (deployed to the cluster by
#: deploy.sh). legacy/ and tests/ are deliberately excluded -- see module
#: docstring.
_SCAN_DIRS = ["bin", "slurm", "oracle"]
_SCAN_EXTS = (".py", ".sh", ".sbatch")


def _pipeline_files() -> list[str]:
    out = []
    for d in _SCAN_DIRS:
        base = os.path.join(ROOT, d)
        for dirpath, _dirnames, filenames in os.walk(base):
            for fn in filenames:
                if fn.endswith(_SCAN_EXTS):
                    out.append(os.path.join(dirpath, fn))
    return sorted(out)


def test_pipeline_scan_finds_the_files_this_guard_depends_on():
    """Pins the scan itself: if directory layout ever changes such that this
    finds nothing, the guard below would trivially (and silently) pass."""
    files = _pipeline_files()
    assert os.path.join(ROOT, "bin", "run_wgs.py") in files
    assert os.path.join(ROOT, "slurm", "prep_input.sh") in files
    assert os.path.join(ROOT, "slurm", "gen_gt116.sh") in files


def test_raw_benchmark_filename_appears_only_in_the_allowlisted_normalizer():
    offenders = []
    for path in _pipeline_files():
        if path in _ALLOWLIST:
            continue
        with open(path) as fh:
            if RAW_BENCHMARK in fh.read():
                offenders.append(os.path.relpath(path, ROOT))
    assert offenders == [], (
        f"these ACTIVE pipeline files reference the RAW benchmark "
        f"({RAW_BENCHMARK!r}) outside slurm/prep_input.sh -- that silently "
        f"breaks the (chrom,pos,ref,alt) join on 47,781 multi-allelic sites "
        f"(see this test module's docstring): {offenders}"
    )


def test_run_wgs_whole_wgs_input_points_at_the_normalized_file():
    path = os.path.join(ROOT, "bin", "run_wgs.py")
    src = open(path).read()
    assert f'"{NORMALIZED_INPUT}"' in src, (
        "bin/run_wgs.py's WHOLE_WGS_INPUT must point at HG002_normalized.vcf.gz")
    assert RAW_BENCHMARK not in src


def test_gen_gt116_whole_wgs_input_points_at_the_normalized_file():
    path = os.path.join(ROOT, "slurm", "gen_gt116.sh")
    src = open(path).read()
    assert NORMALIZED_INPUT in src, (
        "slurm/gen_gt116.sh's WHOLE_WGS_INPUT must point at HG002_normalized.vcf.gz "
        "-- the 116 GT must be built from the SAME input as the 115 GT")
    assert RAW_BENCHMARK not in src


def test_prep_input_is_the_one_file_allowed_to_read_the_raw_benchmark():
    path = os.path.join(ROOT, "slurm", "prep_input.sh")
    src = open(path).read()
    assert RAW_BENCHMARK in src, (
        "slurm/prep_input.sh must read the raw benchmark -- it is the "
        "normalization SOURCE; if this stops being true the allowlist above "
        "is stale and should be revisited, not silently left in place")
