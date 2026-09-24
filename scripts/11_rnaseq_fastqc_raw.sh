#!/bin/bash
set -euo pipefail

# ============================================================
# 11_rnaseq_fastqc_raw.sh
# Run FastQC on raw paired-end RNA-seq reads
#
# Usage:
#   bash 11_rnaseq_fastqc_raw.sh <input_dir> <output_dir> [threads]
#
# Example:
#   bash 11_rnaseq_fastqc_raw.sh RNA_Seq RNA_Seq/fastqc_raw 8
#
# Assumes input files are named:
#   sample_1.fq.gz
#   sample_2.fq.gz
# ============================================================

INPUT_DIR="${1:-RNA_Seq}"
OUTPUT_DIR="${2:-RNA_Seq/fastqc_raw}"
THREADS="${3:-8}"

mkdir -p "${OUTPUT_DIR}"

for R1 in "${INPUT_DIR}"/*_1.fq.gz; do
    [[ -e "$R1" ]] || continue
    R2="${R1/_1.fq.gz/_2.fq.gz}"

    if [[ -f "$R1" && -f "$R2" ]]; then
        fastqc -t "${THREADS}" -o "${OUTPUT_DIR}" "$R1" "$R2"
        echo "Processed: $(basename "$R1") and $(basename "$R2")"
    else
        echo "Warning: missing pair for $R1, skipping."
    fi
done

echo "Raw RNA-seq FastQC completed."
