#!/bin/bash
set -euo pipefail

RAW_INPUT_DIR="${1:-raw_fast5}"
OUTPUT_DIR="${2:-guppy_out}"
THREADS="${3:-30}"

GUPPY_BIN="${GUPPY_BIN:-guppy_basecaller}"
GUPPY_CONFIG="${GUPPY_CONFIG:-dna_r10.4.1_e8.2_400bps_sup.cfg}"
MIN_QSCORE="${MIN_QSCORE:-8}"
DEVICE="${DEVICE:-auto}"

mkdir -p "${OUTPUT_DIR}"

"${GUPPY_BIN}" \
    --config "${GUPPY_CONFIG}" \
    -i "${RAW_INPUT_DIR}" \
    -s "${OUTPUT_DIR}" \
    --min_qscore "${MIN_QSCORE}" \
    -r \
    -x "${DEVICE}" \
    --num_callers "${THREADS}"

find "${OUTPUT_DIR}/pass" -type f -name "*.fastq" -exec gzip -f {} \; 2>/dev/null || true
find "${OUTPUT_DIR}/fail" -type f -name "*.fastq" -exec gzip -f {} \; 2>/dev/null || true

echo "Guppy basecalling completed."
