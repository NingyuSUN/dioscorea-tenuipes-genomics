#!/bin/bash
set -euo pipefail

DRAFT_ASSEMBLY="${1:-flye_out/assembly.fasta}"
READS="${2:-combined_pass.fastq.gz}"
OUTPUT_DIR="${3:-medaka_out}"
THREADS="${4:-30}"
MODEL="${5:-}"

mkdir -p "${OUTPUT_DIR}"

if [[ -n "${MODEL}" ]]; then
    medaka_consensus \
        -i "${READS}" \
        -d "${DRAFT_ASSEMBLY}" \
        -o "${OUTPUT_DIR}" \
        -t "${THREADS}" \
        -m "${MODEL}"
else
    medaka_consensus \
        -i "${READS}" \
        -d "${DRAFT_ASSEMBLY}" \
        -o "${OUTPUT_DIR}" \
        -t "${THREADS}"
fi
