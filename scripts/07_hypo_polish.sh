#!/bin/bash
set -euo pipefail

DRAFT_ASSEMBLY="${1:-medaka_out/consensus.fasta}"
ILLUMINA_LIST="${2:-mapping/il_names.txt}"
SHORT_BAM="${3:-mapping/mapped_sr.sorted.bam}"
LONG_BAM="${4:-mapping/mapped_lg.sorted.bam}"
OUTPUT_ASSEMBLY="${5:-final_assembly.fa}"
GENOME_SIZE="${6:-400m}"
COVERAGE="${7:-50}"
THREADS="${8:-24}"

hypo \
    -d "${DRAFT_ASSEMBLY}" \
    -r @"${ILLUMINA_LIST}" \
    -s "${GENOME_SIZE}" \
    -c "${COVERAGE}" \
    -b "${SHORT_BAM}" \
    -B "${LONG_BAM}" \
    -p 48 \
    -t "${THREADS}" \
    -o "${OUTPUT_ASSEMBLY}"
