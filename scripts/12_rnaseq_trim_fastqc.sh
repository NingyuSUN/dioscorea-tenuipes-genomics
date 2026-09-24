#!/bin/bash
set -euo pipefail

# ============================================================
# 12_rnaseq_trim_fastqc.sh
# Trim paired-end RNA-seq reads with Trimmomatic and run FastQC
#
# Usage:
#   bash 12_rnaseq_trim_fastqc.sh <input_dir> <trimmed_dir> <fastqc_dir> <adapters.fa> [threads] [max_jobs]
#
# Example:
#   bash 12_rnaseq_trim_fastqc.sh RNA_Seq RNA_Seq/trimmomatic RNA_Seq/fastqc_trimmomatic TruSeq3-PE-2.fa 16 4
#
# Assumes input files are named:
#   sample_1.fq.gz
#   sample_2.fq.gz
# ============================================================

INPUT_DIR="${1:-RNA_Seq}"
TRIMMED_DIR="${2:-RNA_Seq/trimmomatic}"
FASTQC_DIR="${3:-RNA_Seq/fastqc_trimmomatic}"
ADAPTERS_FILE="${4:-TruSeq3-PE-2.fa}"
THREADS="${5:-16}"
MAX_JOBS="${6:-4}"

MINLEN="${MINLEN:-36}"
SLIDINGWINDOW="${SLIDINGWINDOW:-4:20}"
HEADCROP="${HEADCROP:-10}"
MAXINFO="${MAXINFO:-40:0.5}"
LOG_FILE="${TRIMMED_DIR}/trimmomatic.log"

mkdir -p "${TRIMMED_DIR}" "${FASTQC_DIR}"

log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" | tee -a "${LOG_FILE}"
}

run_trimmomatic() {
    local sample="$1"
    local r1_file="$2"
    local r2_file="$3"

    local sample_dir="${TRIMMED_DIR}/${sample}"
    local trimmed_r1="${sample_dir}/${sample}_R1.trimmomatic.fastq.gz"
    local trimmed_r2="${sample_dir}/${sample}_R2.trimmomatic.fastq.gz"
    local unpaired_r1="${sample_dir}/${sample}_R1.unpaired.fq.gz"
    local unpaired_r2="${sample_dir}/${sample}_R2.unpaired.fq.gz"
    local trimlog="${sample_dir}/trim.log"
    local summary="${sample_dir}/trim.stat"

    mkdir -p "${sample_dir}"

    log "Trimming sample: ${sample}"
    trimmomatic PE -threads "${THREADS}" -phred33 \
        -trimlog "${trimlog}" \
        -summary "${summary}" \
        "${r1_file}" "${r2_file}" \
        "${trimmed_r1}" "${unpaired_r1}" \
        "${trimmed_r2}" "${unpaired_r2}" \
        ILLUMINACLIP:"${ADAPTERS_FILE}":2:30:10:8:true \
        HEADCROP:"${HEADCROP}" \
        SLIDINGWINDOW:"${SLIDINGWINDOW}" \
        MINLEN:"${MINLEN}" \
        MAXINFO:"${MAXINFO}"

    log "Running FastQC on trimmed reads for ${sample}"
    fastqc -t 4 -o "${FASTQC_DIR}" \
        "${trimmed_r1}" \
        "${trimmed_r2}"
}

JOB_COUNT=0
log "Starting RNA-seq trimming pipeline"

for R1 in "${INPUT_DIR}"/*_1.fq.gz; do
    [[ -e "$R1" ]] || continue
    sample=$(basename "$R1" _1.fq.gz)
    R2="${INPUT_DIR}/${sample}_2.fq.gz"

    if [[ -f "$R1" && -f "$R2" ]]; then
        run_trimmomatic "$sample" "$R1" "$R2" &
        ((JOB_COUNT+=1))
        if [[ "${JOB_COUNT}" -ge "${MAX_JOBS}" ]]; then
            wait
            JOB_COUNT=0
        fi
    else
        log "Skipping ${sample}: paired files missing"
    fi
done

wait
log "RNA-seq trimming pipeline completed."
