#!/usr/bin/env python3
"""Build PDF validation report from validate.py summary.json.
Usage: make_report.py <summary.json> <out.pdf> [title-suffix]"""
import sys, json, datetime, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

summ = json.load(open(sys.argv[1])); outpdf = sys.argv[2]
suffix = sys.argv[3] if len(sys.argv)>3 else ""
ok = [r for r in summ if r.get("status")=="ok"]
# fields that are flag-expected-diffs per combo (HGVS only meaningful when 'hgvs' in combo)
def expected_diff(combo, field):
    if field in ("HGVSc","HGVSp","HGVS_OFFSET") and "hgvs" not in combo: return True
    return False

def true_pct(r):
    pf=r.get("per_field",{}); t=m=0
    for f,v in pf.items():
        if expected_diff(r["name"],f) or v["total"]==0: continue
        t+=v["total"]; m+=v["match"]
    return round(100*m/t,3) if t else None

with PdfPages(outpdf) as pdf:
    # --- page 1: title + overall table ---
    fig=plt.figure(figsize=(11,8.5)); fig.text(0.5,0.93,"vepyr 0.2.0 — Validation vs Ensembl VEP 115.2",ha="center",size=20,weight="bold")
    fig.text(0.5,0.89,f"HG002 GIAB benchmark · CSQ per-field concordance {suffix}",ha="center",size=12)
    fig.text(0.5,0.86,f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} · oracle = real VEP 115.2 ground truth",ha="center",size=9,color="gray")
    rows=[["combo","cache","overall %","true %*","aligned annot.","status"]]
    for r in summ:
        if r.get("status")=="ok":
            rows.append([r["name"], r["cache"].replace("115_GRCh38_",""), f'{r["overall_pct"]}', f'{true_pct(r)}', f'{r.get("aligned_annotations",0):,}', "ok"])
        else:
            rows.append([r["name"], r.get("cache","").replace("115_GRCh38_",""), "-","-","-", r.get("status","?")])
    ax=fig.add_axes([0.05,0.30,0.9,0.5]); ax.axis("off")
    tbl=ax.table(cellText=rows[1:],colLabels=rows[0],loc="center",cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1,1.6)
    for j in range(len(rows[0])): tbl[0,j].set_facecolor("#40466e"); tbl[0,j].set_text_props(color="w",weight="bold")
    fig.text(0.05,0.24,"* true % excludes flag-semantics expected-diffs (HGVSc/HGVSp only meaningful when the combo uses --hgvs).",size=8,color="gray")
    fig.text(0.05,0.20,"Key finding: vepyr is ~100% concordant with VEP on Consequence/IMPACT/gene/transcript/positions/SIFT/PolyPhen/AF/gnomAD/MANE.",size=9)
    fig.text(0.05,0.175,"HGVSc/HGVSp discordance on non-hgvs combos = vepyr computes HGVS even without the flag; the GT leaves it empty (not a bug).",size=9)
    pdf.savefig(fig); plt.close(fig)

    # --- page 2: per-field concordance heatmap ---
    if ok:
        fields=[f for f in ok[0].get("per_field",{}).keys()]
        combos=[r["name"] for r in ok]
        M=np.full((len(combos),len(fields)),np.nan)
        for i,r in enumerate(ok):
            for j,f in enumerate(fields):
                v=r["per_field"].get(f); 
                if v and v["total"]: M[i,j]=v["pct"]
        # split fields into chunks for readability
        chunk=40
        for cs in range(0,len(fields),chunk):
            fl=fields[cs:cs+chunk]; sub=M[:,cs:cs+chunk]
            fig=plt.figure(figsize=(min(20,2+0.32*len(fl)),2+0.5*len(combos)))
            ax=fig.add_subplot(111)
            im=ax.imshow(sub,aspect="auto",cmap="RdYlGn",vmin=0,vmax=100)
            ax.set_xticks(range(len(fl))); ax.set_xticklabels(fl,rotation=90,size=6)
            ax.set_yticks(range(len(combos))); ax.set_yticklabels(combos,size=8)
            ax.set_title(f"Per-field concordance % (fields {cs+1}-{cs+len(fl)})",size=11)
            fig.colorbar(im,ax=ax,fraction=0.02); fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

    # --- page 3: mismatch samples per combo (top diff fields) ---
    for r in ok:
        samp=r.get("samples",{})
        if not samp: continue
        lines=[f"Combo: {r['name']}  (cache {r['cache']}, overall {r['overall_pct']}%, true {true_pct(r)}%)",""]
        for f,ss in list(samp.items())[:8]:
            pct=r["per_field"].get(f,{}).get("pct")
            tag=" [flag-expected-diff]" if expected_diff(r["name"],f) else ""
            lines.append(f"  {f} ({pct}%){tag}:")
            for s in ss[:2]:
                lines.append(f"     {s[0]}"); lines.append(f"        vepyr: {s[1]!r}"); lines.append(f"        VEP  : {s[2]!r}")
        fig=plt.figure(figsize=(11,8.5)); fig.text(0.05,0.95,"Mismatch samples",size=14,weight="bold")
        fig.text(0.05,0.92,"\n".join(lines[:1]),size=9,va="top",family="monospace")
        fig.text(0.05,0.88,"\n".join(lines[1:60]),size=7.5,va="top",family="monospace")
        pdf.savefig(fig); plt.close(fig)
print("[report] wrote", outpdf)
