import sys, gzip, re, collections
def op(p): return gzip.open(p,'rt') if p.endswith('.gz') else open(p)
def csq_fields(path):
    with op(path) as f:
        for l in f:
            if l.startswith('##INFO=<ID=CSQ'):
                m=re.search(r'Format:\s*([^">]+)', l); return m.group(1).split('|') if m else []
            if l.startswith('#CHROM'): break
    return []
def norm_chrom(c): return c[3:] if c.startswith('chr') else c
def load(path):
    fields=csq_fields(path); idx={f:i for i,f in enumerate(fields)}
    feat_i=idx.get('Feature'); rec=collections.defaultdict(dict)
    with op(path) as f:
        for l in f:
            if l.startswith('#'): continue
            c=l.rstrip('\n').split('\t')
            if len(c)<8: continue
            key=(norm_chrom(c[0]),c[1],c[3],c[4])
            m=re.search(r'(?:^|;)CSQ=([^;\t]+)',c[7])
            if not m: continue
            for e in m.group(1).split(','):
                parts=e.split('|'); feat=parts[feat_i] if feat_i is not None and feat_i<len(parts) else str(len(rec[key]))
                rec[key][feat]={fields[i]:(parts[i] if i<len(parts) else '') for i in range(len(fields))}
    return fields, rec
def main(vepyr_vcf, gt_vcf, sample=8):
    vf,vr=load(vepyr_vcf); gf,gr=load(gt_vcf)
    shared=[f for f in vf if f in gf]
    keys=set(vr)&set(gr)
    print(f"vepyr fields={len(vf)} gt fields={len(gf)} shared={len(shared)}")
    print(f"records: vepyr={len(vr)} gt={len(gr)} common_keys={len(keys)} only_vepyr={len(set(vr)-set(gr))} only_gt={len(set(gr)-set(vr))}")
    tot=collections.Counter(); mm=collections.Counter(); samples=collections.defaultdict(list)
    aligned=feat_only=0
    for k in keys:
        for feat in set(vr[k])&set(gr[k]):
            aligned+=1
            for fld in shared:
                a=vr[k][feat].get(fld,''); b=gr[k][feat].get(fld,'')
                tot[fld]+=1
                if a!=b:
                    mm[fld]+=1
                    if len(samples[fld])<sample: samples[fld].append((k,feat,a,b))
    print(f"aligned (key+Feature) annotations = {aligned}")
    print("field concordance (match% over aligned):")
    for fld in shared:
        t=tot[fld] or 1; ok=t-mm[fld]
        print(f"  {fld:24s} {100*ok/t:6.2f}%  ({ok}/{t})")
    print("\nsample mismatches (field: key feat vepyr|gt):")
    for fld in shared:
        if mm[fld]:
            print(f"  [{fld}]")
            for k,feat,a,b in samples[fld][:3]:
                print(f"     {k} {feat}: {a!r} | {b!r}")
if __name__=='__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv)>3 else 8)
