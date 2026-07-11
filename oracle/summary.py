"""Per-combo summary: additive counters so per-chrom shards merge cleanly."""
from __future__ import annotations
import collections
from typing import Any

MAX_SAMPLES_PER_FIELD_CATEGORY = 5


def new_accumulator(fields: list[str]) -> dict:
    return {
        "fields": list(fields),
        "total": collections.Counter(),
        "mismatch": collections.Counter(),
        "by_category": collections.defaultdict(collections.Counter),  # field -> cat -> n
        "samples": collections.defaultdict(list),                      # field -> [[loc, v, g]]
        "sample_seen": collections.Counter(),                          # (field, cat) -> n
        "aligned": 0,
    }


def record_match(acc: dict, field: str) -> None:
    acc["total"][field] += 1


def record_mismatch(acc: dict, field: str, category: str, sample: list[str]) -> None:
    acc["total"][field] += 1
    acc["mismatch"][field] += 1
    acc["by_category"][field][str(category)] += 1
    key = (field, str(category))
    if acc["sample_seen"][key] < MAX_SAMPLES_PER_FIELD_CATEGORY:
        acc["sample_seen"][key] += 1
        acc["samples"][field].append(list(sample) + [str(category)])


def _pct(match: int, total: int) -> float | None:
    return round(100 * match / total, 3) if total else None


def finalize(acc: dict, name: str, cache: str, **extra: Any) -> dict:
    per_field = {}
    for f in acc["fields"]:
        tot = acc["total"][f]
        mm = acc["mismatch"][f]
        per_field[f] = {
            "total": tot,
            "match": tot - mm,
            "pct": _pct(tot - mm, tot),
            "by_category": dict(acc["by_category"][f]),
        }
    T = sum(acc["total"].values())
    M = sum(acc["mismatch"].values())
    return {
        "name": name, "cache": cache, "status": "ok",
        "aligned_annotations": acc["aligned"],
        "overall_pct": _pct(T - M, T),
        "per_field": per_field,
        "samples": {f: acc["samples"][f] for f in acc["fields"] if acc["samples"][f]},
        **extra,
    }


#: How merge() combines each key of a shard summary. Every key finalize() emits
#: MUST be listed in exactly one class below -- merge() refuses a shard carrying an
#: unclassified key, and refuses one MISSING a classified key.
#:
#: This table is the fix for the bug that made every coverage counter evaporate at
#: merge: merge() used to rebuild the dict from scratch and copy 7 keys by hand, so
#: `only_vepyr`, `only_gt`, `join_rate`, `malformed_*`, `shared_fields`, the
#: field-asymmetry lists, `chrom` and `mismatches_tsv` were all silently deleted at
#: exactly the step that produces the published whole-genome number.
#:
#: (Why a dict + this table and not a dataclass: a summary crosses a process
#: boundary as JSON -- one file per shard, written by a cluster job, globbed by the
#: report step -- so a dataclass buys a (de)serialisation layer and still would not
#: stop a field being forgotten inside merge()'s own combining logic. An exhaustive
#: partition of the key space DOES stop it, and additionally rejects a stale
#: half-schema JSON left over from an earlier run.)
IDENTITY_KEYS = ("name", "cache")          # must be identical across shards
STATUS_KEY = "status"
SUM_KEYS = ("aligned_annotations", "only_vepyr", "only_gt",
            "malformed_vepyr", "malformed_gt")
EQUAL_KEYS = ("shared_fields", "vepyr_only_fields", "gt_only_fields")
COLLECT_KEYS = {"chrom": "chroms", "mismatches_tsv": "mismatches_tsvs"}
DERIVED_KEYS = ("overall_pct", "per_field", "samples", "join_rate")

SHARD_KEYS = (set(IDENTITY_KEYS) | {STATUS_KEY} | set(SUM_KEYS) | set(EQUAL_KEYS)
              | set(COLLECT_KEYS) | set(DERIVED_KEYS))


def _check_schema(summaries: list[dict]) -> None:
    """Every shard must carry EXACTLY the classified key set -- no more, no less."""
    for i, s in enumerate(summaries):
        unknown = sorted(set(s) - SHARD_KEYS)
        if unknown:
            raise ValueError(
                f"merge(): shard {i} ({s.get('chrom')!r}) carries key(s) {unknown} "
                f"that merge() does not know how to combine. Classify them in "
                f"summary.py (SUM_KEYS / EQUAL_KEYS / COLLECT_KEYS / DERIVED_KEYS) -- "
                f"refusing to merge rather than silently drop them from the "
                f"whole-genome summary the report is built from."
            )
        missing = sorted(SHARD_KEYS - set(s))
        if missing:
            raise ValueError(
                f"merge(): shard {i} ({s.get('chrom')!r}) is missing key(s) {missing}; "
                f"refusing to merge a half-schema summary (a stale JSON from an "
                f"earlier run would otherwise contribute counters we cannot account for)"
            )


def _check_shard_set(summaries: list[dict], expected_chroms: set[str] | None) -> None:
    """The shard SET itself must be exactly right: no duplicates, no failures, no holes.

    Every one of these is a live failure mode of the sharded WGS run, not a
    hypothetical:

    * DUPLICATES -- the cluster's 24h wall-clock cap actively manufactures them: a
      job times out, is resubmitted, and the glob picks up BOTH outputs. Summed
      twice, chr22's 900 matches + 100 mismatches turned 1500 comparisons into 2500
      and 93.333% into 92.0% with no exception raised.
    * FAILED SHARDS -- merge() used to hardcode ``"status": "ok"`` and never look at
      its inputs, laundering a shard that failed into a clean whole-genome summary.
    * HOLES -- if chr7's job dies and its JSON is never written, merge() would
      cheerfully publish a "whole-genome" number computed from 21 chromosomes.
    """
    seen: dict[str, int] = {}
    for i, s in enumerate(summaries):
        chrom = s["chrom"]
        if chrom in seen:
            raise ValueError(
                f"merge(): chromosome {chrom!r} appears in shard {seen[chrom]} AND "
                f"shard {i}. Two ways to get here, and merge() cannot tell them apart "
                f"-- both are wrong:\n"
                f"  (a) a resubmitted job's output was globbed alongside the original "
                f"(the 24h wall-clock cap manufactures this). Merging both "
                f"DOUBLE-COUNTS every one of its annotations and silently moves the "
                f"published concordance. De-duplicate the shard JSONs -- keep exactly "
                f"one per chromosome -- and re-run.\n"
                f"  (b) someone produced SUB-CHROMOSOMAL (region) diff shards. "
                f"merge() requires chromosome-atomic shards: sub-chromosomal splits "
                f"must be concatenated at the VCF level (bcftools concat) BEFORE "
                f"diffing, never summed as separate summaries afterwards."
            )
        seen[chrom] = i
        if s["status"] != "ok":
            raise ValueError(
                f"merge(): shard {i} (chrom {chrom!r}) has status {s['status']!r}, "
                f"not 'ok'. Refusing to merge a shard that did not complete -- its "
                f"counters are partial, and folding them in would understate coverage "
                f"while reading as a clean whole-genome result."
            )

    if expected_chroms is not None:
        got = set(seen)
        missing = sorted(expected_chroms - got)
        unexpected = sorted(got - expected_chroms)
        if missing or unexpected:
            raise ValueError(
                f"merge(): shard set does not cover the expected chromosomes -- "
                f"missing={missing} unexpected={unexpected}. A whole-genome number "
                f"must not be published from a partial shard set (a job that died "
                f"without writing its JSON would otherwise just vanish from the "
                f"denominator)."
            )


def _check_mergeable(summaries: list[dict]) -> None:
    """Shards must describe the SAME combo over the SAME fields, else the summed
    counters are silently wrong. Fail loudly rather than narrow to shard 0."""
    base = summaries[0]
    base_fields = set(base["per_field"])
    for i, s in enumerate(summaries[1:], start=1):
        if s["name"] != base["name"]:
            raise ValueError(
                f"merge(): shard {i} is a different combo "
                f"({s['name']!r} != {base['name']!r} from shard 0); "
                f"refusing to merge shards from different combo runs"
            )
        if s["cache"] != base["cache"]:
            raise ValueError(
                f"merge(): shard {i} ({s['name']!r}) is from a different cache "
                f"({s['cache']!r} != {base['cache']!r} from shard 0); "
                f"refusing to merge -- the result would sum counters across caches "
                f"and be stamped with only shard 0's cache"
            )
        fields = set(s["per_field"])
        if fields != base_fields:
            missing = sorted(base_fields - fields)
            extra = sorted(fields - base_fields)
            raise ValueError(
                f"merge(): shard {i} ({s['name']!r}) has a different field set "
                f"than shard 0 -- missing={missing} extra={extra}; "
                f"refusing to merge (summed counters would be silently wrong)"
            )
        for k in EQUAL_KEYS:
            # These describe the CSQ schema of the run, not a per-shard tally. If two
            # shards disagree they were produced from different VCF headers and the
            # merged summary cannot honestly claim either value.
            if s[k] != base[k]:
                raise ValueError(
                    f"merge(): shard {i} ({s['chrom']!r}) disagrees with shard 0 on "
                    f"{k!r} ({s[k]!r} != {base[k]!r}); the shards were produced from "
                    f"different CSQ headers -- refusing to merge"
                )


def merge(summaries: list[dict], *, expected_chroms: set[str] | None = None) -> dict:
    """Merge per-shard summaries of the SAME combo into one whole-genome summary.

    Everything a shard carries survives: counters are SUMMED, `join_rate` is
    RECOMPUTED from the summed counters (averaging per-shard rates is not
    associative -- a shard that joined 1/1 and a shard that joined 1/101 average to
    50%, while the truth is 2%), schema facts are asserted equal and carried, and
    the per-shard `chrom` / `mismatches_tsv` are collected into `chroms` /
    `mismatches_tsvs`.

    `expected_chroms` is the completeness contract: pass the chromosomes the run was
    supposed to cover (the report step passes chr1..chr22) and merge() refuses a
    shard set that does not match it exactly. Without it, a chromosome whose job
    died simply disappears from the denominator.

    THE SHARD CONTRACT (enforced, not assumed -- an unwritten assumption is how
    every bug in this module got here):

    * Shards are **chromosome-atomic**. One summary JSON per chromosome, or one for
      the whole genome. `chrom` is the shard's identity, and it is what the
      duplicate check keys on.
    * Sub-chromosomal splitting is a **VCF-level** concern that happens strictly
      BEFORE diffing. The pipeline's other sharded steps may go sub-chromosomal --
      real-VEP ground-truth generation must, since a whole-genome VEP run blows the
      24h cap -- but their region shards are stitched back together with `bcftools
      concat` and never reach this function. `diff_files()` (the only producer of
      summary JSONs) runs whole-genome or per-chromosome.
    * Region shards are therefore NOT supported here, deliberately: nothing produces
      them, and summing region summaries would need a shard identity of
      (chrom, start, end) plus an overlap check to stay honest. If one ever appears,
      the duplicate-chrom guard rejects it with an error that explains this.
    """
    if not summaries:
        raise ValueError("merge() needs at least one summary")
    _check_schema(summaries)
    _check_shard_set(summaries, expected_chroms)
    _check_mergeable(summaries)
    base = summaries[0]
    fields = list(base["per_field"])
    per_field = {
        f: {"total": 0, "match": 0, "pct": None, "by_category": collections.Counter()}
        for f in fields
    }
    samples: dict[str, list] = {f: [] for f in fields}
    summed = dict.fromkeys(SUM_KEYS, 0)
    collected: dict[str, list] = {out: [] for out in COLLECT_KEYS.values()}

    for s in summaries:
        for k in SUM_KEYS:
            summed[k] += s[k]
        for src_key, out_key in COLLECT_KEYS.items():
            collected[out_key].append(s[src_key])
        for f in fields:
            src = s["per_field"].get(f)
            if not src:
                continue
            per_field[f]["total"] += src["total"]
            per_field[f]["match"] += src["match"]
            per_field[f]["by_category"].update(src.get("by_category", {}))
            samples[f].extend(s.get("samples", {}).get(f, []))

    T = M = 0
    for f in fields:
        pf = per_field[f]
        pf["by_category"] = dict(pf["by_category"])
        pf["pct"] = _pct(pf["match"], pf["total"])
        T += pf["total"]
        M += pf["total"] - pf["match"]

    aligned = summed["aligned_annotations"]
    seen = aligned + summed["only_vepyr"] + summed["only_gt"]

    return {
        "name": base["name"], "cache": base["cache"],
        # DERIVED, never hardcoded: _check_shard_set() has already refused anything
        # that is not "ok", so this can only ever restate what the shards said.
        "status": base["status"],
        **summed,
        "join_rate": round(aligned / seen, 4) if seen else None,
        "overall_pct": _pct(T - M, T),
        "per_field": per_field,
        "samples": {f: v[:MAX_SAMPLES_PER_FIELD_CATEGORY * 3] for f, v in samples.items() if v},
        **{k: base[k] for k in EQUAL_KEYS},
        "chroms": sorted(collected["chroms"]),
        "mismatches_tsvs": collected["mismatches_tsvs"],
    }
