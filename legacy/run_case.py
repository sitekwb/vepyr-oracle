import sys, json, time, vepyr
vcf, cache_dir, kw_json = sys.argv[1], sys.argv[2], sys.argv[3]
kw = json.loads(kw_json)
t0 = time.time()
res = vepyr.annotate(vcf, cache_dir, **kw)
print(f"RESULT_PATH {res}  elapsed={time.time()-t0:.1f}s", flush=True)
