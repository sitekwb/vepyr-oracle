#!/usr/bin/env python3
"""vepyr 0.2.0 validation vs real-VEP ground truth. Usage: validate.py <chrom|all> <outdir>"""
import sys, os, json, time, gzip, re, collections
import vepyr
DATA=os.path.expanduser("~/vepyr/data"); FASTA=f"{DATA}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INPUT=f"{DATA}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"; GTDIR=f"{DATA}/ground_truth_vep"; PLUGIN=f"{DATA}/plugin_cache"
MATRIX=[
 ("everything","HG002_annotated_wgs_everything.vcf","115_GRCh38_ensembl",dict(everything=True)),
 ("everything_hgvs","HG002_annotated_wgs_everything_hgvs.vcf","115_GRCh38_ensembl",dict(everything=True,hgvs=True)),
 ("hgvs_merged","HG002_annotated_wgs_everything_hgvs_merged.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True)),
 ("hgvs_merged_am","HG002_annotated_wgs_everything_hgvs_merged_am.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,plugin_cache_root=PLUGIN)),
 ("hgvs_merged_pick","HG002_annotated_wgs_everything_hgvs_merged_pick.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,pick=True)),
 ("hgvs_merged_pick_allele","HG002_annotated_wgs_everything_hgvs_merged_pick_allele.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,pick_allele=True)),
 ("hgvs_merged_pick_allele_gene","HG002_annotated_wgs_everything_hgvs_merged_pick_allele_gene.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,pick_allele_gene=True)),
 ("hgvs_merged_per_gene","HG002_annotated_wgs_everything_hgvs_merged_per_gene.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,per_gene=True)),
 ("hgvs_merged_flag_pick_allele","HG002_annotated_wgs_everything_hgvs_merged_flag_pick_allele.vcf","115_GRCh38_merged",dict(everything=True,hgvs=True,flag_pick_allele=True)),
 ("hgvs_refseq","HG002_annotated_wgs_everything_hgvs_refseq.vcf","115_GRCh38_refseq",dict(everything=True,hgvs=True)),
]
def op(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
def norm(c): return c[3:] if c.startswith('chr') else c
def csq_fields(path):
    with op(path) as f:
        for l in f:
            if l.startswith('##INFO=<ID=CSQ'):
                m=re.search(r'Format:\s*([^">]+)',l); return m.group(1).split('|') if m else []
            if l.startswith('#CHROM'): break
    return []
def chrom_ranks(path):
    r={}; i=0
    with op(path) as f:
        for l in f:
            if l.startswith('##contig'):
                m=re.search(r'ID=([^,>]+)',l)
                if m: r[norm(m.group(1))]=i; i+=1
            elif l.startswith('#CHROM'): break
    return r
def parse(line, fields, feat_i):
    c=line.rstrip('\n').split('\t')
    if len(c)<8: return None
    m=re.search(r'(?:^|;)CSQ=([^;\t]+)',c[7])
    if not m: return None
    feats={}
    for e in m.group(1).split(','):
        p=e.split('|'); ft=p[feat_i] if feat_i is not None and feat_i<len(p) else str(len(feats))
        feats[ft]={fields[i]:(p[i] if i<len(p) else '') for i in range(len(fields))}
    return (norm(c[0]), int(c[1]), c[3], c[4]), feats
def reader(path, fields, feat_i, cf=None):
    with op(path) as f:
        for l in f:
            if l.startswith('#'): continue
            r=parse(l, fields, feat_i)
            if r and (cf is None or r[0][0]==cf): yield r
def grouped(gen):
    cp=None; buf=None
    for key,feats in gen:
        k=(key[0],key[1])
        if cp is None: cp=k; buf={(key[2],key[3]):feats}
        elif k==cp: buf[(key[2],key[3])]=feats
        else: yield cp,buf; cp=k; buf={(key[2],key[3]):feats}
    if cp is not None: yield cp,buf
def result(shared, tot, mm, samp, aligned, only_v, only_g, vkeys, gkeys):
    pf={f:{"match":tot[f]-mm[f],"total":tot[f],"pct":round(100*(tot[f]-mm[f])/tot[f],3) if tot[f] else None} for f in shared}
    T=sum(tot.values()); M=sum(mm.values())
    return dict(vepyr_records=vkeys, gt_records=gkeys, only_vepyr=only_v, only_gt=only_g,
        aligned_annotations=aligned, shared_fields=len(shared),
        overall_pct=round(100*(T-M)/T,3) if T else None, per_field=pf,
        samples={f:samp[f][:3] for f in shared if mm[f]})
def diff_stream(vfile, gfile, cf=None):
    vf=csq_fields(vfile); gf=csq_fields(gfile); shared=[x for x in vf if x in gf]
    vfi={f:i for i,f in enumerate(vf)}.get('Feature'); gfi={f:i for i,f in enumerate(gf)}.get('Feature')
    ranks=chrom_ranks(gfile) or chrom_ranks(vfile)
    rk=lambda cp:(ranks.get(cp[0],9999),cp[1])
    vg=grouped(reader(vfile,vf,vfi,cf)); gg=grouped(reader(gfile,gf,gfi,cf))
    tot=collections.Counter(); mm=collections.Counter(); samp=collections.defaultdict(list)
    aligned=only_v=only_g=vk=gk=0
    V=next(vg,None); G=next(gg,None)
    while V and G:
        if rk(V[0])<rk(G[0]): only_v+=len(V[1]); vk+=len(V[1]); V=next(vg,None)
        elif rk(G[0])<rk(V[0]): only_g+=len(G[1]); gk+=len(G[1]); G=next(gg,None)
        else:
            vk+=len(V[1]); gk+=len(G[1])
            for ra in set(V[1])&set(G[1]):
                vf2=V[1][ra]; gf2=G[1][ra]
                for ft in set(vf2)&set(gf2):
                    aligned+=1
                    for fld in shared:
                        a=vf2[ft].get(fld,''); b=gf2[ft].get(fld,''); tot[fld]+=1
                        if a!=b:
                            mm[fld]+=1
                            if len(samp[fld])<5: samp[fld].append([f"{V[0][0]}:{V[0][1]} {ra[0]}>{ra[1]} {ft}",a,b])
            only_v+=len(set(V[1])-set(G[1])); only_g+=len(set(G[1])-set(V[1]))
            V=next(vg,None); G=next(gg,None)
    while V: only_v+=len(V[1]); vk+=len(V[1]); V=next(vg,None)
    while G: only_g+=len(G[1]); gk+=len(G[1]); G=next(gg,None)
    return result(shared,tot,mm,samp,aligned,only_v,only_g,vk,gk)
def subset(chrom, out):
    import pysam
    vin=pysam.VariantFile(INPUT); vout=pysam.VariantFile(out,'wz',header=vin.header); n=0
    for r in vin.fetch(chrom): vout.write(r); n+=1
    vout.close(); vin.close(); pysam.tabix_index(out,preset='vcf',force=True); return n
def main():
    chrom=sys.argv[1]; outdir=sys.argv[2]; os.makedirs(outdir,exist_ok=True)
    cf=None if chrom=='all' else norm(chrom)
    if chrom=='all': inp=INPUT
    else: inp=f"{outdir}/input_{chrom}.vcf.gz"; print(f"[subset] {chrom}: {subset(chrom,inp)} variants",flush=True)
    results=[]
    for name,gtf,cache,kw in MATRIX:
        gt=f"{GTDIR}/{gtf}"; cdir=f"{DATA}/{cache}"; out=f"{outdir}/vepyr_{name}.vcf"
        r={"name":name,"cache":cache,"gt":gtf}
        if not os.path.exists(gt): r["status"]="no_gt"; results.append(r); print(f"[skip] {name}",flush=True); continue
        kwargs=dict(reference_fasta=FASTA,output_vcf=out,show_progress=False,workers=8,**kw); t0=time.time()
        try:
            res=vepyr.annotate(inp,cdir,**kwargs); r["annotate_s"]=round(time.time()-t0,1)
            print(f"[annotate] {name}: {res} ({r['annotate_s']}s)",flush=True)
        except Exception as e:
            r["status"]="annotate_fail"; r["error"]=f"{type(e).__name__}: {e}"[:300]; results.append(r)
            print(f"[ANNOTATE_FAIL] {name}: {r['error']}",flush=True)
            json.dump(results,open(f"{outdir}/summary.json","w"),indent=2); continue
        try:
            t1=time.time(); r.update(diff_stream(out,gt,cf)); r["diff_s"]=round(time.time()-t1,1); r["status"]="ok"
            print(f"[diff] {name}: overall={r.get('overall_pct')}% aligned={r.get('aligned_annotations')} ({r['diff_s']}s)",flush=True)
        except Exception as e:
            import traceback; r["status"]="diff_fail"; r["error"]=f"{type(e).__name__}: {e}"[:300]; traceback.print_exc()
            print(f"[DIFF_FAIL] {name}: {r['error']}",flush=True)
        results.append(r); json.dump(results,open(f"{outdir}/summary.json","w"),indent=2)
    print("[DONE]",f"{outdir}/summary.json",flush=True)
if __name__=='__main__': main()
