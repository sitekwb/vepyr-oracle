#!/usr/bin/env bash
#SBATCH --job-name=fdomrep
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=%x-%j.out
#
# ADJUDICATION EXPERIMENT for the fork-delta finding on DOMAINS.
#
# The whole-genome fork delta found exactly TWO (variant, transcript) loci where
# the DOMAINS field carries the SAME set of domains in a DIFFERENT order between
# the unforked and forked GT (oracle.drift.Category.ORDER_DIFF). Two runs -- one
# per arm -- cannot tell apart:
#   (a) `--fork` reorders DOMAINS, vs
#   (b) DOMAINS ordering is nondeterministic PER PROCESS (Perl hash ordering is
#       seeded per interpreter), so ANY two runs differ, forked or not.
# Those have opposite consequences for trusting the forked GT, so guessing is not
# an option. This runs the SAME tiny region FOUR times -- three UNFORKED, one
# FORKED -- and prints the DOMAINS value each time. If the three unforked runs
# disagree among THEMSELVES, the answer is (b) and fork is exonerated for DOMAINS.
set -euo pipefail
source "$HOME/vepyr/work/slurm/lib_common.sh"
load_apptainer
DATA="$HOME/vepyr/data"; WORK="$HOME/vepyr/work"
VEP_SIF="$HOME/bvp/vep116.sif"; BCFTOOLS_SIF="$HOME/aidiva-onb/bcftools.sif"
FASTA="$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CACHE="$DATA/vep_native_cache_116"
OUT="$WORK/fork_delta/domains_replicate"; mkdir -p "$OUT"
REGION="chr2:1080000-1270000"
SLICE="$OUT/input_slice.vcf.gz"
if [ ! -s "$SLICE" ]; then
  apptainer exec -B "$WORK:$WORK" "$BCFTOOLS_SIF" \
    bcftools view -r "$REGION" -O z -o "$SLICE" "$WORK/input/input_2.vcf.gz"
fi
echo "[rep] slice records: $(apptainer exec -B "$WORK:$WORK" "$BCFTOOLS_SIF" bcftools view -H "$SLICE" | wc -l)"
run_vep() {  # $1=tag  $2...=extra flags
  local tag="$1"; shift
  local o="$OUT/rep_${tag}.vcf"
  [ -s "$o" ] && { echo "[rep] $tag exists"; return 0; }
  apptainer exec -B "$DATA:$DATA" -B "$WORK:$WORK" "$VEP_SIF" \
    vep --cache --offline --dir_cache "$CACHE" --fasta "$FASTA" \
        --input_file "$SLICE" --output_file "$o" --vcf --force_overwrite \
        --assembly GRCh38 --everything --hgvs --merged --no_stats "$@"
}
run_vep unforked_a
run_vep unforked_b
run_vep unforked_c
run_vep forked_16 --fork 16
echo "=========== DOMAINS for ENST00000308624 at the two divergent loci ==========="
for tag in unforked_a unforked_b unforked_c forked_16; do
  for pos in 1083634 1267458; do
    v=$(awk -v p="$pos" -F"\t" '$2==p' "$OUT/rep_${tag}.vcf" \
        | tr ";" "\n" | grep "^CSQ=" | sed "s/^CSQ=//" | tr "," "\n" \
        | awk -F"|" -v t="ENST00000308624" '$7==t {print $0}' \
        | awk -F"|" '{print $(NF)}' | head -1)
    # DOMAINS is not the last field; pull it by header index instead
    idx=$(grep -m1 "^##INFO=<ID=CSQ" "$OUT/rep_${tag}.vcf" | sed "s/.*Format: //;s/\">.*//" \
          | tr "|" "\n" | grep -n "^DOMAINS$" | cut -d: -f1)
    d=$(awk -v p="$pos" -F"\t" '$2==p' "$OUT/rep_${tag}.vcf" \
        | tr ";" "\n" | grep "^CSQ=" | sed "s/^CSQ=//" | tr "," "\n" \
        | awk -F"|" -v t="ENST00000308624" -v i="$idx" '$7==t {print $i}' | head -1)
    echo "$tag  chr2:$pos  DOMAINS=$d"
  done
done
