#!/bin/bash
# Reference-guided scaffolding against Kita1 haplotype 1 — the primary reference
# designated in Kudoh et al. 2026 PLoS Genetics 22(4):e1012123
# (172 contigs, 414.4 Mb, contig N50 33.9 Mb, six T2T chromosomes).
# Reference file used here is the chromosome-only subset (chr1-chr10).
# Usage: 23_ragtag_scaffolding.sh <Kita1_hap1_chr.fasta> <genome.fa> [out_dir] [threads]
# Requires ragtag.py and minimap2 on PATH (for example, an activated conda env).
set -e
REF=${1:?Usage: 23_ragtag_scaffolding.sh <Kita1_hap1_chr.fasta> <genome.fa> [out_dir] [threads]}
QRY=${2:?Usage: 23_ragtag_scaffolding.sh <Kita1_hap1_chr.fasta> <genome.fa> [out_dir] [threads]}
OUT=${3:-ragtag_out_Kita1}
THREADS=${4:-16}

echo "=== RagTag vs Kita1 hap1 start: $(date) ==="
ragtag.py scaffold "$REF" "$QRY" -o "$OUT" -t $THREADS --mm2-params "-x asm20 -t $THREADS"
echo "=== done: $(date) ==="

echo "=== stats ==="
cat "$OUT"/ragtag.scaffold.stats

awk '/^>/ { if (s) print len; s=1; len=0; next } { len += length($0) } END { if (s) print len }' \
  "$OUT"/ragtag.scaffold.fasta | sort -nr > "$OUT"/scaffold_lengths.txt
awk '
{ L[NR]=$1; total+=$1 }
END { n=NR; half=total/2; run=0
  for (i=1;i<=n;i++){ run+=L[i]; if(!n50 && run>=half){n50=L[i]; l50=i} }
  printf "sequences\t%d\ntotal_bp\t%d\nlongest\t%d\nN50\t%d\nL50\t%d\n", n,total,L[1],n50,l50 }' "$OUT"/scaffold_lengths.txt

echo "=== anchored bp per chromosome ==="
awk '$5=="W" && $1 ~ /_RagTag/ {a[$1]+=($3-$2+1); n[$1]++} END{for(c in a) printf "%s\t%d\t%d\n", c,a[c],n[c]}' \
  "$OUT"/ragtag.scaffold.agp | sort -k1,1V
