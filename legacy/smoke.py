import vepyr, os, time, pysam
DATA=os.path.expanduser("~/vepyr/data"); W=os.path.expanduser("~/vepyr/work/smoke"); os.makedirs(W,exist_ok=True)
vin=pysam.VariantFile(f"{DATA}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
outp=f"{W}/hg002_chr22.vcf.gz"; vout=pysam.VariantFile(outp,"wz",header=vin.header); n=0
for r in vin.fetch("chr22"): vout.write(r); n+=1
vout.close(); pysam.tabix_index(outp, preset="vcf", force=True)
print(f"[smoke] chr22 input variants: {n} (bgzipped+tabixed)", flush=True)
t0=time.time()
try:
    res=vepyr.annotate(outp, f"{DATA}/115_GRCh38_ensembl", everything=True,
        reference_fasta=f"{DATA}/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        output_vcf=f"{W}/vepyr_chr22.vcf", show_progress=False, workers=8)
    print(f"[smoke] annotate OK: {res}  elapsed={time.time()-t0:.0f}s", flush=True)
except Exception as e:
    import traceback; print(f"[smoke] ANNOTATE_FAIL: {type(e).__name__}: {e}", flush=True); traceback.print_exc(); raise
os.system(f"echo '=== CSQ header ==='; grep -m1 'ID=CSQ' {W}/vepyr_chr22.vcf | cut -c1-600")
os.system(f"echo '=== 2 records (trunc) ==='; grep -v '^#' {W}/vepyr_chr22.vcf | head -2 | cut -c1-260")
os.system(f"echo -n '=== annotated CSQ-record count: '; grep -vc '^#' {W}/vepyr_chr22.vcf")
print("[smoke] DONE", flush=True)
