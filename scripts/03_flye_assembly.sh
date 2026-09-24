#!/bin/bash
set -euo pipefail

READS="${1:-combined_pass.fastq.gz}"
OUTPUT_DIR="${2:-flye_out}"
GENOME_SIZE="${3:-400m}"
THREADS="${4:-30}"

flye \
    --nano-raw "${READS}" \
    -g "${GENOME_SIZE}" \
    -t "${THREADS}" \
    --scaffold \
    -o "${OUTPUT_DIR}" \
    --keep-haplotypes \
    --asm-coverage 50
