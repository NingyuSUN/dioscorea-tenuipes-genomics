#!/bin/bash
set -euo pipefail

# ============================================================
# 10_te_repeatmasker.sh
# Soft-mask transposable elements with RepeatMasker
#
# Usage:
#   bash 10_te_repeatmasker.sh <genome.fa> <output_dir> [species] [threads]
#
# Example:
#   bash 10_te_repeatmasker.sh D_tenuipes.fa repeatmasker_out arabidopsis 12
#
# Notes:
# - Output softmasked genome is typically *.masked
# - Replace species/library as needed
# - For BRAKER with --softmasking, use the masked genome
# ============================================================

GENOME="${1:-D_tenuipes.fa}"
OUTPUT_DIR="${2:-repeatmasker_out}"
SPECIES="${3:-arabidopsis}"
THREADS="${4:-12}"

mkdir -p "${OUTPUT_DIR}"

RepeatMasker \
    -species "${SPECIES}" \
    -pa "${THREADS}" \
    -s \
    -nolow \
    -norna \
    -xsmall \
    -dir "${OUTPUT_DIR}" \
    "${GENOME}"

echo "RepeatMasker completed."
echo "Softmasked genome should be in: ${OUTPUT_DIR}"
