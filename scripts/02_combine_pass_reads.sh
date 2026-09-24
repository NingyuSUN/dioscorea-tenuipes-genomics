#!/bin/bash
set -euo pipefail

GUPPY_OUT_DIR="${1:-guppy_out}"
OUTPUT_FASTQ="${2:-combined_pass.fastq.gz}"

PASS_DIR="${GUPPY_OUT_DIR}/pass"

mapfile -t PASS_FILES < <(find "${PASS_DIR}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

cat "${PASS_FILES[@]}" > "${OUTPUT_FASTQ}"

echo "Combined ${#PASS_FILES[@]} files into: ${OUTPUT_FASTQ}"
