#!/bin/bash
set -euo pipefail

# ============================================================
# 16_cp_align_to_reference.sh
# Align assembled contigs to a reference chloroplast genome
# using minimap2 and generate a sorted/indexed BAM.
#
# Usage:
#   bash 16_cp_align_to_reference.sh <contigs.fa> <cp_reference.fa> <output_prefix> [threads]
#
# Example:
#   bash 16_cp_align_to_reference.sh whole_genome.h2.fa D_tokoro_cp.fa dtenuipes_cp 8
# ============================================================

CONTIGS="${1:-whole_genome.h2.fa}"
REFERENCE="${2:-D_tokoro_cp.fa}"
PREFIX="${3:-dtenuipes_cp}"
THREADS="${4:-8}"

SAM_OUTPUT="${PREFIX}.sam"
SORTED_BAM="${PREFIX}.sorted.bam"

command -v minimap2 >/dev/null 2>&1 || { echo "Error: minimap2 not found."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found."; exit 1; }

[[ -f "${CONTIGS}" ]] || { echo "Error: contigs file not found: ${CONTIGS}"; exit 1; }
[[ -f "${REFERENCE}" ]] || { echo "Error: reference file not found: ${REFERENCE}"; exit 1; }

echo "Aligning contigs to chloroplast reference..."
minimap2 -t "${THREADS}" -ax asm5 "${REFERENCE}" "${CONTIGS}" > "${SAM_OUTPUT}"

echo "Converting, sorting, and indexing BAM..."
samtools view -@ "${THREADS}" -bS "${SAM_OUTPUT}" | samtools sort -@ "${THREADS}" -o "${SORTED_BAM}"
samtools index "${SORTED_BAM}"

rm -f "${SAM_OUTPUT}"

echo "Done."
echo "Sorted BAM: ${SORTED_BAM}"
echo "BAM index: ${SORTED_BAM}.bai"
