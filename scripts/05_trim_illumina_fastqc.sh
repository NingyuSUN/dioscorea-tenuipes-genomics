#!/bin/bash
set -euo pipefail

R1="${1:-reads_R1.fastq.gz}"
R2="${2:-reads_R2.fastq.gz}"
OUTPUT_DIR="${3:-illumina_qc}"
ADAPTERS="${4:-TruSeq3-PE-2.fa}"
THREADS="${5:-16}"

RAW_FASTQC_DIR="${OUTPUT_DIR}/raw_fastqc"
TRIMMED_DIR="${OUTPUT_DIR}/trimmed"
TRIMMED_FASTQC_DIR="${OUTPUT_DIR}/trimmed_fastqc"

mkdir -p "${RAW_FASTQC_DIR}" "${TRIMMED_DIR}" "${TRIMMED_FASTQC_DIR}"

TRIMMED_R1="${TRIMMED_DIR}/trimmed_R1.fastq.gz"
TRIMMED_R2="${TRIMMED_DIR}/trimmed_R2.fastq.gz"
UNPAIRED_R1="${TRIMMED_DIR}/trimmed_R1.unpaired.fastq.gz"
UNPAIRED_R2="${TRIMMED_DIR}/trimmed_R2.unpaired.fastq.gz"

fastqc -t 4 -o "${RAW_FASTQC_DIR}" "${R1}" "${R2}"

trimmomatic PE -threads "${THREADS}" -phred33 \
    "${R1}" "${R2}" \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" \
    "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10:8:true \
    HEADCROP:10 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    MAXINFO:40:0.5

fastqc -t 4 -o "${TRIMMED_FASTQC_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"
