"""Per-element pipeline state: done / annotated / missing, plus the early-warning
overrun check.

THE PROBLEM THIS SOLVES: every shard in `shards.tsv` carries a SAFETY-ADJUSTED
`est_hours` (see oracle/shards.py) -- but that estimate is a MODEL, and the
whole point of the 12-16h band + 1.4x margin is to survive the model being
somewhat wrong. It is NOT designed to survive the model being wrong and nobody
noticing until a job hits the 24h cap and gets killed after burning a day of
compute. `scan()` below is the early-warning: for every element whose output
already exists and carries a `<out>.timing.json` (written by `bin/run_wgs.py`),
it compares the ACTUAL elapsed time against the PLANNED `est_hours` and flags
anything that ran materially over -- catchable while the run is still healthy,
not after the fact.

No pysam import at module level: this only ever reads `shards.tsv` rows and
stats/reads small JSON files on disk, so it stays laptop-testable like the rest
of `oracle/`.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

from .shards import Shard, Step

#: How far ACTUAL elapsed time has to exceed the PLANNED (safety-adjusted)
#: est_hours before an element is flagged as an overrun.
#:
#: 1.25 is deliberately a check on the ADJUSTED estimate -- the number the
#: planner actually decided against (raw rate * safety_factor, default 1.4) --
#: not the raw model. The raw model is EXPECTED to be beaten some of the time;
#: that slack is what the margin exists to absorb. Flagging at 1.25x the
#: ADJUSTED estimate means the RAW model was overrun by ~1.75x (1.25 * 1.4) --
#: a systematic divergence, not noise -- which is exactly the situation that
#: would eventually walk an element into the 24h cap if it keeps recurring
#: across a shard set and nobody re-measures sec_per_variant.
OVERRUN_FACTOR = 1.25


class State:
    MISSING = "missing"       # no output yet for this shard's step
    ANNOTATED = "annotated"   # (diff shards only) the upstream vepyr/GT VCF exists,
                              # but this shard's own diff has not run yet
    DONE = "done"             # this shard's terminal output exists and is non-empty


@dataclass(frozen=True)
class ElementStatus:
    shard: Shard
    state: str
    out_path: str
    timing_path: str | None
    actual_hours: float | None
    overran: bool


def annotate_out_path(resultsdir: str, combo: str, chrom: str, chunk: int | None,
                      *, prefix: str = "vepyr_") -> str:
    """The vepyr/real-VEP output VCF path for one (combo, chrom[, chunk]) shard.

    MUST match what `slurm/annotate.sbatch` and `slurm/gen_gt116.sh` actually
    write -- see those scripts' header comments. This is the single canonical
    Python-side definition of that naming convention; the sbatch scripts are
    bash and cannot import it, so their own comments point back here and must
    be kept in sync by hand (a real, acknowledged limitation of the bash/Python
    split -- see the task report).

    `prefix` defaults to `"vepyr_"` (the `annotate` step's convention, set by
    `slurm/annotate.sbatch`). The `gt` step's real-VEP shards (written by
    `slurm/gen_gt116.sh`) carry NO tool prefix -- `"vepyr_"` on a real-VEP
    output would be actively misleading, and the directory itself
    (`ground_truth_vep_116/shards/`) already supplies that context. Callers
    dispatch on `shard.step`; see `bin/status.py`.
    """
    if chrom == "all":
        return os.path.join(resultsdir, f"{prefix}{combo}.vcf")
    sfx = f"_{chrom}" + (f"_c{chunk}" if chunk is not None else "")
    return os.path.join(resultsdir, f"{prefix}{combo}{sfx}.vcf")


def diff_summary_path(resultsdir: str, combo: str, chrom: str) -> str:
    """MUST match `bin/validate.py::summary_path()` exactly (diff shards are
    never L2, so there is never a chunk suffix to account for here)."""
    sfx = "" if chrom == "all" else f"_{chrom}"
    return os.path.join(resultsdir, f"summary_{combo}{sfx}.json")


def _timing_path(vcf_path: str) -> str:
    return f"{vcf_path}.timing.json"


def _exists_nonempty(path: str) -> bool:
    return os.path.exists(path) and os.path.getsize(path) > 0


def element_status(shard: Shard, resultsdir: str) -> ElementStatus:
    """One shard's state (+ overrun flag), given the directory its step's
    outputs land in.

    `gt` / `annotate` shards: DONE means the output VCF exists and is
    non-empty (mirrors `bin/run_wgs.py`'s own resume check exactly, so this
    reports the same thing a re-run of that CLI would decide). There is no
    intermediate ANNOTATED state for these two steps.

    `diff` shards: MISSING until `bin/validate.py` has written
    `summary_<combo>[_<chrom>].json` (-> DONE). ANNOTATED means the upstream
    vepyr output this diff needs already exists but the diff itself has not
    run yet -- the state a diff array element sits in right after its
    `aftercorr` dependency on the annotate array clears, before it is
    scheduled. Diff shards are never L2 (oracle.shards enforces this), so
    their upstream lookup never needs a chunk suffix.
    """
    if shard.step is Step.DIFF:
        summary = diff_summary_path(resultsdir, shard.combo, shard.chrom)
        upstream = annotate_out_path(resultsdir, shard.combo, shard.chrom, None)
        if _exists_nonempty(summary):
            state = State.DONE
        elif _exists_nonempty(upstream):
            state = State.ANNOTATED
        else:
            state = State.MISSING
        return ElementStatus(shard=shard, state=state, out_path=summary,
                             timing_path=None, actual_hours=None, overran=False)

    # shard.chunk is already None unless shard.level is L2 (oracle.shards's own
    # Shard contract), so no extra level check is needed here. The gt step's
    # real-VEP shards carry no "vepyr_" prefix (see annotate_out_path's docstring).
    prefix = "" if shard.step is Step.GT else "vepyr_"
    out = annotate_out_path(resultsdir, shard.combo, shard.chrom, shard.chunk, prefix=prefix)
    timing = _timing_path(out)
    exists = _exists_nonempty(out)
    state = State.DONE if exists else State.MISSING

    actual_hours: float | None = None
    overran = False
    if exists and os.path.exists(timing):
        try:
            with open(timing) as fh:
                data = json.load(fh)
        except (json.JSONDecodeError, OSError):
            data = None            # a malformed/unreadable timing file must not crash the scan
        if data is not None:
            elapsed_s = data.get("elapsed_s")
            if elapsed_s is not None:
                actual_hours = elapsed_s / 3600.0
                overran = actual_hours > shard.est_hours * OVERRUN_FACTOR

    return ElementStatus(shard=shard, state=state, out_path=out,
                         timing_path=timing if exists else None,
                         actual_hours=actual_hours, overran=overran)


def scan(shards: list[Shard], resultsdir: str) -> list[ElementStatus]:
    return [element_status(s, resultsdir) for s in shards]


def summarize(statuses: list[ElementStatus]) -> dict[str, int]:
    """State counts, plus `"overran"` -- counted separately because it is
    ORTHOGONAL to state (an overrun is only ever observed on an element that is
    already DONE; a MISSING or ANNOTATED element has no completed timing to
    compare)."""
    counts: dict[str, int] = {}
    for st in statuses:
        counts[st.state] = counts.get(st.state, 0) + 1
    counts["overran"] = sum(1 for st in statuses if st.overran)
    return counts
