"""The 8 in-scope combos: vepyr kwargs (vepyr side) + VEP CLI flags (116 GT side)."""
from __future__ import annotations
import csv, json
from dataclasses import dataclass
from .csq import open_maybe_gzip

MATRIX_HEADER = ["name", "cache_flavor", "cache115", "cache116",
                 "vepyr_kwargs", "vep_flags", "gt115", "gt116"]


@dataclass(frozen=True)
class Combo:
    name: str
    cache_flavor: str          # "merged" | "refseq"
    gt115: str                 # filename inside data/ground_truth_vep/
    kwargs: dict               # passed to vepyr.annotate(**kwargs)


_G = "HG002_annotated_wgs_everything_hgvs"

#: merged+refseq only. The two `everything*` combos are ensembl and are OUT OF SCOPE.
COMBOS: list[Combo] = [
    Combo("hgvs_merged",                  "merged", f"{_G}_merged.vcf",
          dict(everything=True, hgvs=True)),
    Combo("hgvs_merged_am",               "merged", f"{_G}_merged_am.vcf",
          dict(everything=True, hgvs=True, plugin_cache_root="@PLUGIN@")),
    Combo("hgvs_merged_pick",             "merged", f"{_G}_merged_pick.vcf",
          dict(everything=True, hgvs=True, pick=True)),
    Combo("hgvs_merged_pick_allele",      "merged", f"{_G}_merged_pick_allele.vcf",
          dict(everything=True, hgvs=True, pick_allele=True)),
    Combo("hgvs_merged_pick_allele_gene", "merged", f"{_G}_merged_pick_allele_gene.vcf",
          dict(everything=True, hgvs=True, pick_allele_gene=True)),
    Combo("hgvs_merged_per_gene",         "merged", f"{_G}_merged_per_gene.vcf",
          dict(everything=True, hgvs=True, per_gene=True)),
    Combo("hgvs_merged_flag_pick_allele", "merged", f"{_G}_merged_flag_pick_allele.vcf",
          dict(everything=True, hgvs=True, flag_pick_allele=True)),
    Combo("hgvs_refseq",                  "refseq", f"{_G}_refseq.vcf",
          dict(everything=True, hgvs=True)),
]


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


def write_matrix(path: str, vep_flags_by_combo: dict[str, str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MATRIX_HEADER, delimiter="\t")
        w.writeheader()
        for c in COMBOS:
            w.writerow({
                "name": c.name,
                "cache_flavor": c.cache_flavor,
                "cache115": f"115_GRCh38_{c.cache_flavor}",
                "cache116": f"116_GRCh38_{c.cache_flavor}",
                "vepyr_kwargs": json.dumps(c.kwargs),
                "vep_flags": vep_flags_by_combo.get(c.name, ""),
                "gt115": c.gt115,
                "gt116": f"{c.name}.vcf",
            })


def load_matrix(path: str) -> dict[str, dict]:
    with open(path) as fh:
        return {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
