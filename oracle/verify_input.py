"""Key-set identity check: does OUR normalized input carry the exact same
(chrom, pos, ref, alt) variant keys as the ground truth real VEP was run
against?

THE BUG THIS GUARDS AGAINST: the pipeline joins vepyr output against real-VEP
ground truth on (chrom, pos, ref, alt). Every ground-truth VCF was built with
`--input_file .../HG002_normalized.vcf` -- the raw HG002 benchmark AFTER
`bcftools norm -m -any` split every multi-allelic site into one row per ALT
(see slurm/prep_input.sh's header for the record-count proof). Feeding vepyr
the RAW, un-split benchmark instead still "works" (nothing crashes) but
silently breaks the join key on every one of the 47,781 multi-allelic sites:
a joint record `REF=C ALT=T,CCGC` and its split counterparts `REF=C ALT=T` /
`REF=C ALT=CCGC` are three DIFFERENT (chrom,pos,ref,alt) keys, none of which
match each other -- so those variants silently drop out of the comparison,
and multi-ALT handling is exactly the weak spot we most need to measure.

This module answers ONE narrow, decisive question: for a given chromosome, is
the SET of (chrom,pos,ref,alt) keys in our candidate normalized file
IDENTICAL to the set real VEP actually saw? A record-count match (GATE 1 in
slurm/prep_input.sh) only proves the two files have the same number of rows
-- it says nothing about whether those rows carry the SAME variants. Two
files can agree on count and disagree entirely on content (e.g. if
left-alignment shifted every indel's position by one base, or normalized its
REF/ALT padding differently). Key-set identity (GATE 2) is the only check
that actually proves the downstream join will work.

DELIBERATELY LITERAL KEYS -- ALT IS NEVER SPLIT ON COMMA HERE: if one side
still carries an un-split multi-allelic row (`ALT="T,CCGC"`), that string is
kept AS THE KEY VERBATIM, never exploded into two keys. Exploding on comma
would make this check blind to precisely the bug it exists to catch: a
joint-vs-split representation mismatch would then produce the SAME exploded
key set on both sides, defeating the whole point. A real representation
mismatch (joint vs. split, or an un-left-aligned indel) must show up as a
literal set difference, not be normalized away before the comparison runs.

MEMORY CHARACTERISTIC -- READ THIS BEFORE POINTING THIS AT MORE THAN ONE
CHROMOSOME AT A TIME: `read_keys()` streams its input file line-by-line
(never loads the whole file into memory) but ACCUMULATES every
matching-chromosome key into a Python set. For chr22 (~50k variants) that is
trivial. For the whole genome (~4.1M variants) it would hold roughly that
many tuples in memory PER SIDE -- which is why this module's contract is
deliberately per-chromosome (`chrom` is a required, explicit argument,
never "all"), and why slurm/prep_input.sh's GATE 2 runs it against chr22
alone rather than genome-wide. A genome-wide check, if ever needed, should
loop this one chromosome at a time (comparing and discarding each
chromosome's sets before moving to the next) rather than building one
4M-entry set per side.

No pysam import at module level: this is a plain text/gzip scan (via
`oracle.csq.open_maybe_gzip`), so it stays laptop-testable like the rest of
`oracle/`.
"""
from __future__ import annotations

from dataclasses import dataclass

from .csq import norm_chrom, open_maybe_gzip

#: (chrom, pos, ref, alt) -- chrom has already been passed through
#: norm_chrom(), so a "chr22"-prefixed file and a bare-"22" file compare equal.
VariantKey = tuple[str, int, str, str]


def read_keys(vcf_path: str, chrom: str) -> set[VariantKey]:
    """Every (chrom, pos, ref, alt) key on `chrom` in `vcf_path`, as a set.

    Streams the file once; only keys matching `chrom` (after norm_chrom()
    normalization on BOTH the file's own contig column and the requested
    `chrom`, so "chr22" and "22" name the same chromosome) are ever held in
    memory -- see the module docstring's memory-characteristic note for why
    this is a per-chromosome contract, not a whole-genome one.

    ALT is kept EXACTLY as written in the file -- a still-joint multi-allelic
    ALT (`"T,CCGC"`) is one literal key component, never split into two keys.
    See the module docstring for why that is the point, not an oversight.
    """
    target = norm_chrom(chrom)
    keys: set[VariantKey] = set()
    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            this_chrom = norm_chrom(cols[0])
            if this_chrom != target:
                continue
            keys.add((this_chrom, int(cols[1]), cols[3], cols[4]))
    return keys


@dataclass(frozen=True)
class KeySetComparison:
    """The result of comparing two files' variant key sets for one chromosome."""

    chrom: str
    ours_total: int
    gt_total: int
    only_ours: frozenset[VariantKey]
    only_gt: frozenset[VariantKey]

    @property
    def overlap(self) -> int:
        """Keys present on BOTH sides -- the count the downstream join can
        actually use."""
        return self.ours_total - len(self.only_ours)

    @property
    def identical(self) -> bool:
        """True iff the two files carry EXACTLY the same variant keys on this
        chromosome -- the only condition under which the (chrom,pos,ref,alt)
        join used downstream is safe."""
        return not self.only_ours and not self.only_gt


def compare_key_sets(ours: set[VariantKey], gt: set[VariantKey], *,
                     chrom: str) -> KeySetComparison:
    """Compare two key sets for one chromosome.

    Pure set arithmetic, deliberately separated from file I/O (`read_keys`)
    so the comparison itself is trivially unit-testable against small,
    hand-built sets -- see tests/test_verify_input.py.
    """
    return KeySetComparison(
        chrom=chrom,
        ours_total=len(ours),
        gt_total=len(gt),
        only_ours=frozenset(ours - gt),
        only_gt=frozenset(gt - ours),
    )


def format_examples(keys: frozenset[VariantKey], limit: int = 10) -> list[str]:
    """Up to `limit` keys, sorted for deterministic, diffable output."""
    return [f"{c}:{p} {r}>{a}" for c, p, r, a in sorted(keys)[:limit]]
