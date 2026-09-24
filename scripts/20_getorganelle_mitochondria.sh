#!/bin/bash
set -euo pipefail

# ============================================================
# 20_getorganelle_mitochondria.sh
# Assemble plant mitochondrial genome using GetOrganelle.
#
# Usage:
#   bash 20_getorganelle_mitochondria.sh <R1.fq.gz> <R2.fq.gz> <output_dir> [threads]
#
# Example:
#   bash 20_getorganelle_mitochondria.sh Dt_Izu_1.fq Dt_Izu_2.fq mito_getorganelle 16
#
# Notes:
# - Recommended for public repo because GetOrganelle is widely
#   used for organelle genome assembly.
# - For plant mitochondria, use -F embplant_mt
# ============================================================

R1="${1:-reads_R1.fq.gz}"
R2="${2:-reads_R2.fq.gz}"
OUTPUT_DIR="${3:-mito_getorganelle}"
THREADS="${4:-16}"

command -v get_organelle_from_reads.py >/dev/null 2>&1 || { echo "Error: GetOrganelle not found."; exit 1; }

[[ -f "${R1}" ]] || { echo "Error: forward reads not found: ${R1}"; exit 1; }
[[ -f "${R2}" ]] || { echo "Error: reverse reads not found: ${R2}"; exit 1; }

get_organelle_from_reads.py     -1 "${R1}"     -2 "${R2}"     -o "${OUTPUT_DIR}"     -F embplant_mt     -t "${THREADS}"

echo "GetOrganelle mitochondrial assembly completed."
echo "Output directory: ${OUTPUT_DIR}"
