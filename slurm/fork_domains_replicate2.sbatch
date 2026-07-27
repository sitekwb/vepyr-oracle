#!/usr/bin/env bash
#SBATCH --job-name=fdomrep2
#SBATCH --partition=cpu
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=%x-%j.out
#
# Round 2 of the DOMAINS adjudication (see fork_domains_replicate.sbatch).
#
# Round 1 ran a 190 kb slice and found unforked == forked == a THIRD ordering,
# i.e. --fork did not reproduce the production difference -- but a 190 kb slice
# fills VEP's InputBuffer very differently from a whole chromosome, and the
# hypothesised mechanism IS buffer-scoped, so that alone does not settle it.
#
# This round makes the input production-LIKE (chr2:1-30Mb, ~40k variants, enough
# to fill and recycle the buffer many times) and adds the mechanism probe the
# hypothesis actually names: an UNFORKED run whose --buffer_size is set to what
# `--fork 16` leaves each worker (default 5000 / (2*16) = 156). If DOMAINS
# ordering is buffer-scoped, THAT is what flips it -- with no fork involved.
set -euo pipefail
source "$HOME/vepyr/work/slurm/lib_common.sh"
load_apptainer
DATA="$HOME/vepyr/data"; WORK="$HOME/vepyr/work"
VEP_SIF="$HOME/bvp/vep116.sif"; BCFTOOLS_SIF="$HOME/aidiva-onb/bcftools.sif"
FASTA="$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CACHE="$DATA/vep_native_cache_116"
OUT="$WORK/fork_delta/domains_replicate2"; mkdir -p "$OUT"
SLICE="$OUT/input_chr2_30Mb.vcf.gz"
if [ ! -s "$SLICE" ]; then
    apptainer exec -B "$WORK:$WORK" "$BCFTOOLS_SIF" \
        bcftools view -r chr2:1-30000000 -O z -o "$SLICE" "$WORK/input/input_2.vcf.gz"
fi
echo "[rep2] slice records: $(apptainer exec -B "$WORK:$WORK" "$BCFTOOLS_SIF" bcftools view -H "$SLICE" | wc -l)"

run_vep() {
    local tag="$1"; shift
    local o="$OUT/rep_${tag}.vcf"
    if [ -s "$o" ]; then echo "[rep2] $tag exists -- skip"; return 0; fi
    local t0=$SECONDS
    apptainer exec -B "$DATA:$DATA" -B "$WORK:$WORK" "$VEP_SIF" \
        vep --cache --offline --dir_cache "$CACHE" --fasta "$FASTA" \
            --input_file "$SLICE" --output_file "$o" --vcf --force_overwrite \
            --assembly GRCh38 --everything --hgvs --merged --no_stats "$@"
    echo "[rep2] $tag done in $((SECONDS - t0))s"
}

run_vep forked16 --fork 16
run_vep unforked_buf156 --buffer_size 156
run_vep unforked

echo "=========== DOMAINS for ENST00000308624 at the two production-divergent loci ==========="
for tag in unforked unforked_buf156 forked16; do
    f="$OUT/rep_${tag}.vcf"
    idx=$(grep -m1 '^##INFO=<ID=CSQ' "$f" | sed 's/.*Format: //;s/">.*//' \
          | tr '|' '\n' | grep -n '^DOMAINS$' | cut -d: -f1)
    fidx=$(grep -m1 '^##INFO=<ID=CSQ' "$f" | sed 's/.*Format: //;s/">.*//' \
           | tr '|' '\n' | grep -n '^Feature$' | cut -d: -f1)
    for pos in 1083634 1267458; do
        d=$(awk -F'\t' -v p="$pos" '$2==p' "$f" | tr ';' '\n' | grep '^CSQ=' \
            | sed 's/^CSQ=//' | tr ',' '\n' \
            | awk -F'|' -v t=ENST00000308624 -v i="$idx" -v fi="$fidx" '$fi==t {print $i; exit}')
        echo "$tag  chr2:$pos  DOMAINS=$d"
    done
done

echo "=========== whole-slice pairwise identity of the three runs ==========="
# Body-only comparison (headers legitimately differ: they record the command line).
for pair in "unforked unforked_buf156" "unforked forked16" "unforked_buf156 forked16"; do
    set -- $pair
    a="$OUT/rep_$1.vcf"; b="$OUT/rep_$2.vcf"
    n=$(diff <(grep -v '^#' "$a") <(grep -v '^#' "$b") | grep -c '^<' || true)
    echo "$1 vs $2 : $n differing records (of $(grep -vc '^#' "$a"))"
done
