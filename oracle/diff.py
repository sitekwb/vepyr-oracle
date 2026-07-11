"""Streaming merge-join diff of two CSQ-annotated VCFs, with drift classification."""
from __future__ import annotations
import csv
from .csq import open_maybe_gzip, csq_format_fields, chrom_ranks, norm_chrom, parse_record
from .drift import classify
from .summary import new_accumulator, record_match, record_mismatch, finalize

#: `allele` is load-bearing, not decoration: on a multi-allelic record a
#: (chrom,pos,ref,alt,feature) tuple no longer identifies a single annotation.
TSV_HEADER = ["chrom", "pos", "ref", "alt", "allele", "feature", "field",
              "vepyr_val", "vep_val", "category"]


def _reader(path, fields, feat_i, allele_i, chrom=None, stats=None):
    """Yield (key, feats); tally keyless (malformed) CSQ entries into `stats`."""
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rec = parse_record(line, fields, feat_i, allele_i)
            if rec is None:
                continue
            key, feats, malformed = rec
            if chrom is not None and key[0] != chrom:
                continue
            if stats is not None:
                stats["malformed"] += malformed
            yield key, feats


def _checked(gen, stream: str, rank):
    """Fail loudly if a stream is not strictly ascending in (contig_rank, pos).

    The merge-join advances two pointers in lockstep and can only do that over
    globally sorted streams. On an unsorted stream it does NOT crash -- it walks
    straight past records that would have matched, dumping them into
    only_vepyr/only_gt and starving the per-field counters. The report then looks
    like annotation drift when it is really a plumbing bug. Concatenated per-chrom
    shards and karyotypic-vs-lexicographic contig collation both produce exactly
    this, so guard the invariant instead of trusting it.
    """
    prev = None
    for pos_key, buf in gen:
        if prev is not None and rank(pos_key) <= rank(prev):
            raise ValueError(
                f"{stream} not sorted by (contig_rank,pos): "
                f"{pos_key[0]}:{pos_key[1]} came after {prev[0]}:{prev[1]}. "
                f"The merge-join requires both VCFs strictly ascending in "
                f"(contig_rank,pos), where contig rank comes from the ground-truth "
                f"VCF's ##contig header order -- so a contig-collation mismatch "
                f"(karyotypic vs lexicographic), a contig missing from that header, "
                f"or concatenated per-chrom shards will trip this. Refusing to "
                f"continue: diffing an unsorted stream silently undercounts matches."
            )
        prev = pos_key
        yield pos_key, buf


def _record_key(ref: str, alt: str) -> tuple[str, tuple[str, ...]]:
    """Outer record key: (REF, sorted ALT alleles).

    ALT is a comma-separated LIST whose order carries no meaning -- `T,CCGC` and
    `CCGC,T` are the same variant. Keyed by the raw string, a tool that happened to
    list the alleles in the other order failed to join at all and the whole record
    fell into only_vepyr/only_gt, silently hollowing out the sample. Multi-ALT is
    exactly where the two tools are most likely to disagree on representation, so
    this must not be left to chance. Sorting makes listing order un-observable.

    Only the OUTER key is normalised: the per-entry CSQ `Allele` still disambiguates
    annotations *within* the record, and the raw ALT string is kept for output.
    """
    return ref, tuple(sorted(alt.split(",")))


def _grouped(gen):
    """Collapse consecutive records at the same (chrom,pos).

    -> (chrom, pos), {record_key: (raw_alt, feats)}
    """
    cur_key = None
    buf: dict = {}
    for key, feats in gen:
        chrom, pos, ref, alt = key
        pos_key = (chrom, pos)
        if cur_key is not None and pos_key != cur_key:
            yield cur_key, buf
            buf = {}
        cur_key = pos_key
        buf[_record_key(ref, alt)] = (alt, feats)
    if cur_key is not None:
        yield cur_key, buf


def _annotations(buf: dict) -> int:
    """Total CSQ annotations across every record in a (chrom,pos) group."""
    return sum(len(feats) for _raw_alt, feats in buf.values())


def diff_files(vepyr_vcf: str, gt_vcf: str, *, name: str, cache: str,
               combo_kwargs: dict, tsv_path: str, chrom: str | None = None,
               rel_tol: float | None = None, abs_tol: float | None = None) -> dict:
    vf = csq_format_fields(vepyr_vcf)
    gf = csq_format_fields(gt_vcf)
    shared = [f for f in vf if f in gf]
    # A field only one side emits cannot be value-compared, so it is dropped from
    # the diff -- but the drop must not be silent: a field vepyr never implements
    # would otherwise leave a report that is 100% green on the fields it does emit.
    vepyr_only_fields = sorted(set(vf) - set(gf))
    gt_only_fields = sorted(set(gf) - set(vf))
    v_idx = {f: i for i, f in enumerate(vf)}
    g_idx = {f: i for i, f in enumerate(gf)}
    v_feat, v_allele = v_idx.get("Feature"), v_idx.get("Allele")
    g_feat, g_allele = g_idx.get("Feature"), g_idx.get("Allele")

    ranks = chrom_ranks(gt_vcf)
    rank_src = gt_vcf
    if not ranks:
        ranks, rank_src = chrom_ranks(vepyr_vcf), vepyr_vcf

    def rank(pk: tuple[str, int]) -> tuple[int, str, int]:
        """Merge-join order key: (contig_rank, contig_name, pos).

        An unranked contig gets NO default. The old `ranks.get(c, 9999)` gave every
        contig missing from the ##contig header the SAME rank, and the merge-join's
        equal-branch fires on rank equality without ever comparing contig NAMES --
        so two DIFFERENT unranked contigs (a vepyr-only chrUn_A and a GT-only
        chrUn_B) became "the same contig", were joined, and minted fabricated
        value_diffs at a join_rate of 1.0. The quieter variant: norm_chrom("chrM")
        is "M" but a header saying `##contig=<ID=MT>` ranks "MT", so "M" sorted last
        behind the sentinel, tripped no guard, and the ENTIRE MITOCHONDRION vanished
        into only_vepyr/only_gt while the report read 100%.

        The contig NAME is part of the key so that even a rank collision can never
        make two distinct contigs compare equal.
        """
        c, pos = pk
        if c not in ranks:
            raise ValueError(
                f"contig {c!r} has no rank: it is absent from the ##contig header of "
                f"{rank_src}. The merge-join orders BOTH streams by that header, so "
                f"an unranked contig has no defined position in the stream and cannot "
                f"be joined or ordered. Refusing to guess: a default rank makes every "
                f"unranked contig compare EQUAL to every other (joining chrUn_A "
                f"against chrUn_B and fabricating value_diffs), and hides a naming "
                f"mismatch such as chrM-vs-MT by silently sweeping the whole contig "
                f"into only_vepyr/only_gt. Add {c!r} to the ground-truth VCF's "
                f"##contig header, or restrict the run with chrom=."
            )
        return ranks[c], c, pos

    cf = norm_chrom(chrom) if chrom else None
    acc = new_accumulator(shared)
    only_v = only_g = 0
    # keyless CSQ entries: excluded from the join (they cannot be identified),
    # but counted so their existence is never invisible in the report.
    v_stats, g_stats = {"malformed": 0}, {"malformed": 0}

    kw = {k: v for k, v in (("rel_tol", rel_tol), ("abs_tol", abs_tol)) if v is not None}

    with open(tsv_path, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(TSV_HEADER)

        vg = _checked(_grouped(_reader(vepyr_vcf, vf, v_feat, v_allele, cf, v_stats)),
                      "vepyr_vcf", rank)
        gg = _checked(_grouped(_reader(gt_vcf, gf, g_feat, g_allele, cf, g_stats)),
                      "gt_vcf", rank)
        V, G = next(vg, None), next(gg, None)
        while V and G:
            if rank(V[0]) < rank(G[0]):
                only_v += _annotations(V[1]); V = next(vg, None)
            elif rank(G[0]) < rank(V[0]):
                only_g += _annotations(G[1]); G = next(gg, None)
            else:
                (chrom_, pos_) = V[0]
                for rk_ in set(V[1]) & set(G[1]):
                    raw_alt, vfeats = V[1][rk_]
                    _gt_raw_alt, gfeats = G[1][rk_]
                    ref = rk_[0]
                    # keys are (allele, feature): one annotation per allele x transcript
                    for key in set(vfeats) & set(gfeats):
                        allele, ft = key
                        acc["aligned"] += 1
                        for fld in shared:
                            a = vfeats[key].get(fld, "")
                            b = gfeats[key].get(fld, "")
                            cat = classify(fld, combo_kwargs, a, b, **kw)
                            if cat is None:
                                record_match(acc, fld)
                            else:
                                loc = f"{chrom_}:{pos_} {ref}>{raw_alt} {allele}|{ft}"
                                record_mismatch(acc, fld, cat, [loc, a, b])
                                writer.writerow([chrom_, pos_, ref, raw_alt, allele, ft,
                                                 fld, a, b, str(cat)])
                    # annotations present on one side only, INSIDE a joined record:
                    # a transcript vepyr emits and VEP does not (or vice versa) was
                    # previously uncounted entirely -- it hid inside a perfect score.
                    only_v += len(set(vfeats) - set(gfeats))
                    only_g += len(set(gfeats) - set(vfeats))
                for rk_ in set(V[1]) - set(G[1]):
                    only_v += len(V[1][rk_][1])
                for rk_ in set(G[1]) - set(V[1]):
                    only_g += len(G[1][rk_][1])
                V, G = next(vg, None), next(gg, None)
        while V:
            only_v += _annotations(V[1]); V = next(vg, None)
        while G:
            only_g += _annotations(G[1]); G = next(gg, None)

    # What fraction of all annotations seen did we actually get to compare?
    # overall_pct only speaks for the annotations that JOINED; on its own it cannot
    # distinguish "vepyr agrees with VEP" from "we compared almost nothing and the
    # scraps agreed". Never raise on a low rate -- some non-overlap is legitimate --
    # but the report must be able to say "we compared X% of annotations".
    seen = acc["aligned"] + only_v + only_g
    join_rate = round(acc["aligned"] / seen, 4) if seen else None

    return finalize(acc, name=name, cache=cache,
                    only_vepyr=only_v, only_gt=only_g, join_rate=join_rate,
                    shared_fields=len(shared),
                    vepyr_only_fields=vepyr_only_fields,
                    gt_only_fields=gt_only_fields,
                    malformed_vepyr=v_stats["malformed"],
                    malformed_gt=g_stats["malformed"],
                    mismatches_tsv=tsv_path,
                    chrom=chrom or "all")
