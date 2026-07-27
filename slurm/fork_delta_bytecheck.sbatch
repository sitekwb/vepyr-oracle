#!/usr/bin/env bash
#SBATCH --job-name=fdbyte
#SBATCH --partition=cpu
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out
#
# INDEPENDENT cross-check of the fork-delta result, deliberately NOT using
# oracle.diff at all.
#
# The per-field table is computed by parsing the CSQ string. That is the right
# resolution for the question, but it is blind by construction to anything
# OUTSIDE CSQ (QUAL, FILTER, the FORMAT/sample columns, any other INFO key) and
# it can only report on fields it knows how to split out. If the engine were
# silently dropping a field, the table would look clean. So: compare the two
# arms' merged VCFs as raw TEXT, header excluded (the header legitimately differs
# -- it records the command line, which is the one thing that IS different).
#
# CONSTANT MEMORY, DELIBERATELY. The first version of this script used `diff`,
# which loads both inputs; on the 29 GB hgvs_merged pair it was OOM-killed, and
# because the pipeline ended in `grep -c ... || true` it reported "0 differing
# lines" -- a FALSE clean bill produced by the check dying, which is exactly the
# silent-wrongness this whole pipeline exists to eliminate. Both arms are
# coordinate-sorted from the same input and have identical record counts
# (verified: 4 096 123 in every combo of both arms), so line i corresponds to
# line i and a POSITIONAL, streaming comparison is both valid and O(1) memory.
# `paste` with \001 as the delimiter cannot collide with VCF text.
#
# Expectations, if the CSQ-level finding is complete and correct:
#   hgvs_refseq  -> 0 differing lines (the diff found zero differing cells)
#   hgvs_merged  -> differing lines > 0, but ZERO once HGNC: tokens and the two
#                   known DOMAINS orderings are normalised away.
set -euo pipefail
FORKED="$HOME/vepyr/data/ground_truth_vep_116"
UNFORKED="$HOME/vepyr/data/ground_truth_vep_116_unforked"

for combo in hgvs_refseq hgvs_merged; do
    a="$UNFORKED/$combo.vcf"; b="$FORKED/$combo.vcf"
    echo "=== $combo ==="
    na=$(grep -vc '^#' "$a"); nb=$(grep -vc '^#' "$b")
    echo "  records: unforked=$na forked=$nb"
    if [ "$na" != "$nb" ]; then
        echo "  [FATAL] record counts differ -- positional comparison is invalid" >&2
        exit 1
    fi

    raw=$(paste -d $'\001' <(grep -v '^#' "$a") <(grep -v '^#' "$b") \
          | awk -F'\001' '$1!=$2 {n++} END {print n+0}')
    echo "  differing lines (raw):                       $raw of $na"

    # Strip every HGNC: token from both sides. What survives is everything the
    # HGNC_ID finding does NOT explain.
    noh=$(paste -d $'\001' <(grep -v '^#' "$a" | sed 's/HGNC:[0-9]*//g') \
                           <(grep -v '^#' "$b" | sed 's/HGNC:[0-9]*//g') \
          | awk -F'\001' '$1!=$2 {n++} END {print n+0}')
    echo "  differing lines after stripping HGNC tokens: $noh"

    if [ "$noh" -gt 0 ]; then
        echo "  --- residual differing lines (first 3, truncated to 200 chars/side):"
        paste -d $'\001' <(grep -v '^#' "$a" | sed 's/HGNC:[0-9]*//g') \
                         <(grep -v '^#' "$b" | sed 's/HGNC:[0-9]*//g') \
          | awk -F'\001' '$1!=$2 {print "  UNFORKED " substr($1,1,200);
                                  print "  FORKED   " substr($2,1,200); c++}
                          c>=3 {exit}'
    fi
done
