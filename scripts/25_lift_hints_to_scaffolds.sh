#!/bin/bash
# Lift the repaired hints from original contig coordinates onto the
# Kita1-scaffolded pseudochromosomes, using the RagTag AGP.
# Usage: 25_lift_hints_to_scaffolds.sh [work_dir]
#   work_dir holds hints_fixed_*.gff and ragtag_out_Kita1/ (default: current directory).
# Requires ragtag.py and minimap2 on PATH (for example, an activated conda env).
set -e
cd "${1:-.}"

AGP=ragtag_out_Kita1/ragtag.scaffold.agp

for f in hints_fixed_introns hints_fixed_exonparts; do
  echo "=== lifting $f ==="
  ragtag.py updategff "$f".gff "$AGP" > "$f".scaffolded.gff
  echo "  input : $(wc -l < $f.gff) lines"
  echo "  output: $(wc -l < $f.scaffolded.gff) lines"
done

echo
echo "=== sanity: contig names before/after ==="
echo "--- before ---"
cut -f1 hints_fixed_introns.gff | sort -u | head -3
echo "--- after ---"
cut -f1 hints_fixed_introns.scaffolded.gff | sort -u | head -12
echo
echo "=== how many hints landed on pseudochromosomes vs unplaced ==="
awk '{ if ($1 ~ /_RagTag$/) chr++; else un++ } END { printf "on pseudochromosomes: %d\nunplaced contigs:     %d\n", chr, un }' \
  hints_fixed_introns.scaffolded.gff
echo
echo "=== sample of lifted intron hints ==="
head -3 hints_fixed_introns.scaffolded.gff
