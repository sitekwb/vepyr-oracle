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


def merge(summaries: list[dict]) -> dict:
    """Merge per-shard summaries of the SAME combo into one."""
    if not summaries:
        raise ValueError("merge() needs at least one summary")
    _check_mergeable(summaries)
    base = summaries[0]
    fields = list(base["per_field"])
    per_field = {
        f: {"total": 0, "match": 0, "pct": None, "by_category": collections.Counter()}
        for f in fields
    }
    samples: dict[str, list] = {f: [] for f in fields}
    aligned = 0
    for s in summaries:
        aligned += s.get("aligned_annotations", 0)
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

    return {
        "name": base["name"], "cache": base["cache"], "status": "ok",
        "aligned_annotations": aligned,
        "overall_pct": _pct(T - M, T),
        "per_field": per_field,
        "samples": {f: v[:MAX_SAMPLES_PER_FIELD_CATEGORY * 3] for f, v in samples.items() if v},
    }
