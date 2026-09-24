#!/bin/bash
set -euo pipefail

ASSEMBLY1="${1:-flye_out/assembly.fasta}"
ASSEMBLY2="${2:-medaka_out/consensus.fasta}"
REFERENCE_FA="${3:-reference.fa}"
REFERENCE_GFF="${4:-reference.gff3}"
OUTPUT_DIR="${5:-quast_out}"
THREADS="${6:-30}"

mkdir -p "${OUTPUT_DIR}"

quast.py \
    "${ASSEMBLY1}" \
    "${ASSEMBLY2}" \
    -t "${THREADS}" \
    -r "${REFERENCE_FA}" \
    -g "${REFERENCE_GFF}" \
    -o "${OUTPUT_DIR}" \
    -b
