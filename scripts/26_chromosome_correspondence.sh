#!/bin/bash
# Determine whether the Kita1 and Waka1 D. tokoro assemblies use the same
# chromosome numbering, by aligning the two references to each other.
# Usage: 26_chromosome_correspondence.sh [work_dir]
#   work_dir holds ref_Dtokoro_Kita1_hap1_chr.fasta and ref_Dtokoro_Waka1_chr.fasta.
# Requires ragtag.py and minimap2 on PATH (for example, an activated conda env).
set -e
cd "${1:-.}"

KITA=ref_Dtokoro_Kita1_hap1_chr.fasta
WAKA=ref_Dtokoro_Waka1_chr.fasta

if [ ! -s kita_vs_waka.paf ]; then
  minimap2 -x asm10 -t 8 "$WAKA" "$KITA" > kita_vs_waka.paf 2> kita_vs_waka.log
fi

echo "=== best Waka1 match for each Kita1 chromosome (by total aligned bp) ==="
# col1 query(Kita1), col6 target(Waka1), col10 residue matches
awk '{ key=$1"\t"$6; m[key]+=$10 }
     END { for (k in m) print k"\t"m[k] }' kita_vs_waka.paf \
| sort -k1,1V -k3,3nr \
| awk '!seen[$1]++ { printf "Kita1 %-8s -> Waka1 %-8s  (%d bp aligned)\n", $1, $2, $3 }'

echo
echo "=== full matrix (top 3 per Kita1 chromosome) ==="
awk '{ key=$1"\t"$6; m[key]+=$10 }
     END { for (k in m) print k"\t"m[k] }' kita_vs_waka.paf \
| sort -k1,1V -k3,3nr \
| awk '{ c[$1]++; if (c[$1]<=3) printf "  %-8s %-8s %d\n", $1, $2, $3 }'
