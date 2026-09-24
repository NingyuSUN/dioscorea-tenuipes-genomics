#!/bin/bash
set -euo pipefail

REF="${1:-medaka_out/consensus.fasta}"
NANOPORE_READS="${2:-combined_pass.fastq.gz}"
R1="${3:-illumina_qc/trimmed/trimmed_R1.fastq.gz}"
R2="${4:-illumina_qc/trimmed/trimmed_R2.fastq.gz}"
OUTPUT_DIR="${5:-mapping}"
THREADS="${6:-24}"

mkdir -p "${OUTPUT_DIR}"

LONG_BAM="${OUTPUT_DIR}/mapped_lg.sorted.bam"
SHORT_BAM="${OUTPUT_DIR}/mapped_sr.sorted.bam"

minimap2 -ax map-ont "${REF}" "${NANOPORE_READS}" | \
    samtools sort -@ "${THREADS}" -o "${LONG_BAM}"
samtools index "${LONG_BAM}"

bwa index "${REF}"

bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
    samtools sort -@ "${THREADS}" -o "${SHORT_BAM}"
samtools index "${SHORT_BAM}"

printf "%s\n%s\n" "${R1}" "${R2}" > "${OUTPUT_DIR}/il_names.txt"
