"""Streaming merge-join diff of two CSQ-annotated VCFs, with drift classification."""
from __future__ import annotations
import csv
from .csq import open_maybe_gzip, csq_format_fields, chrom_ranks, norm_chrom, parse_record
from .drift import classify
from .summary import new_accumulator, record_match, record_mismatch, finalize

TSV_HEADER = ["chrom", "pos", "ref", "alt", "feature", "field",
              "vepyr_val", "vep_val", "category"]


def _reader(path, fields, feat_i, chrom=None):
    with open_maybe_gzip(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            rec = parse_record(line, fields, feat_i)
            if rec and (chrom is None or rec[0][0] == chrom):
                yield rec


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


def _grouped(gen):
    """Collapse consecutive records at the same (chrom,pos) into {(ref,alt): feats}."""
    cur_key = None
    buf: dict = {}
    for key, feats in gen:
        pos_key = (key[0], key[1])
        if cur_key is None:
            cur_key, buf = pos_key, {(key[2], key[3]): feats}
        elif pos_key == cur_key:
            buf[(key[2], key[3])] = feats
        else:
            yield cur_key, buf
            cur_key, buf = pos_key, {(key[2], key[3]): feats}
    if cur_key is not None:
        yield cur_key, buf


def diff_files(vepyr_vcf: str, gt_vcf: str, *, name: str, cache: str,
               combo_kwargs: dict, tsv_path: str, chrom: str | None = None,
               eps: float | None = None) -> dict:
    vf = csq_format_fields(vepyr_vcf)
    gf = csq_format_fields(gt_vcf)
    shared = [f for f in vf if f in gf]
    v_feat = {f: i for i, f in enumerate(vf)}.get("Feature")
    g_feat = {f: i for i, f in enumerate(gf)}.get("Feature")

    ranks = chrom_ranks(gt_vcf) or chrom_ranks(vepyr_vcf)
    rank = lambda pk: (ranks.get(pk[0], 9999), pk[1])

    cf = norm_chrom(chrom) if chrom else None
    acc = new_accumulator(shared)
    only_v = only_g = 0

    kw = {} if eps is None else {"eps": eps}

    with open(tsv_path, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(TSV_HEADER)

        vg = _checked(_grouped(_reader(vepyr_vcf, vf, v_feat, cf)), "vepyr_vcf", rank)
        gg = _checked(_grouped(_reader(gt_vcf, gf, g_feat, cf)), "gt_vcf", rank)
        V, G = next(vg, None), next(gg, None)
        while V and G:
            if rank(V[0]) < rank(G[0]):
                only_v += len(V[1]); V = next(vg, None)
            elif rank(G[0]) < rank(V[0]):
                only_g += len(G[1]); G = next(gg, None)
            else:
                (chrom_, pos_) = V[0]
                for ra in set(V[1]) & set(G[1]):
                    vfeats, gfeats = V[1][ra], G[1][ra]
                    for ft in set(vfeats) & set(gfeats):
                        acc["aligned"] += 1
                        for fld in shared:
                            a = vfeats[ft].get(fld, "")
                            b = gfeats[ft].get(fld, "")
                            cat = classify(fld, combo_kwargs, a, b, **kw)
                            if cat is None:
                                record_match(acc, fld)
                            else:
                                record_mismatch(acc, fld, cat,
                                                [f"{chrom_}:{pos_} {ra[0]}>{ra[1]} {ft}", a, b])
                                writer.writerow([chrom_, pos_, ra[0], ra[1], ft, fld, a, b, str(cat)])
                only_v += len(set(V[1]) - set(G[1]))
                only_g += len(set(G[1]) - set(V[1]))
                V, G = next(vg, None), next(gg, None)
        while V:
            only_v += len(V[1]); V = next(vg, None)
        while G:
            only_g += len(G[1]); G = next(gg, None)

    return finalize(acc, name=name, cache=cache,
                    only_vepyr=only_v, only_gt=only_g,
                    shared_fields=len(shared), mismatches_tsv=tsv_path,
                    chrom=chrom or "all")
