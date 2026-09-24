#!/bin/bash
set -euo pipefail

# ============================================================
# 14_stringtie_transcripts.sh
# Assemble transcripts from RNA-seq BAM using StringTie
#
# Usage:
#   bash 14_stringtie_transcripts.sh <aligned.bam> <output.gtf> [threads]
#
# Example:
#   bash 14_stringtie_transcripts.sh annotation_output/rna_aligned.sorted.bam annotation_output/transcripts.gtf 16
#
# Notes:
# - Optional but useful for transcript structure support
# - Can be kept as an auxiliary annotation product
# ============================================================

INPUT_BAM="${1:-annotation_output/rna_aligned.sorted.bam}"
OUTPUT_GTF="${2:-annotation_output/transcripts.gtf}"
THREADS="${3:-16}"

mkdir -p "$(dirname "${OUTPUT_GTF}")"

stringtie "${INPUT_BAM}" \
    -p "${THREADS}" \
    -o "${OUTPUT_GTF}"

echo "StringTie transcript assembly completed."
echo "Transcript GTF: ${OUTPUT_GTF}"
