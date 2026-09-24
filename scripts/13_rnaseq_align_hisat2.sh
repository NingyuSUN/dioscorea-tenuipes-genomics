#!/bin/bash
set -euo pipefail

# ============================================================
# 13_rnaseq_align_hisat2.sh
# Build HISAT2 index and align trimmed RNA-seq reads to genome
#
# Usage:
#   bash 13_rnaseq_align_hisat2.sh <genome.fa> <trimmed_dir> <output_dir> [threads]
#
# Example:
#   bash 13_rnaseq_align_hisat2.sh D_tenuipes.softmasked.fa RNA_Seq/trimmomatic annotation_output 16
#
# Notes:
# - This script merges all sample BAMs into one coordinate-sorted BAM
# - Input trimmed reads are expected under:
#     trimmed_dir/sample/sample_R1.trimmomatic.fastq.gz
#     trimmed_dir/sample/sample_R2.trimmomatic.fastq.gz
# ============================================================

GENOME="${1:-D_tenuipes.softmasked.fa}"
TRIMMED_DIR="${2:-RNA_Seq/trimmomatic}"
OUTPUT_DIR="${3:-annotation_output}"
THREADS="${4:-16}"

INDEX_PREFIX="${OUTPUT_DIR}/hisat2_index/D_tenuipes"
BAM_DIR="${OUTPUT_DIR}/hisat2_bam"
MERGED_BAM="${OUTPUT_DIR}/rna_aligned.sorted.bam"

mkdir -p "${OUTPUT_DIR}" "$(dirname "${INDEX_PREFIX}")" "${BAM_DIR}"

echo "Building HISAT2 index..."
hisat2-build "${GENOME}" "${INDEX_PREFIX}"

echo "Aligning RNA-seq samples..."
BAM_LIST=()

for sample_dir in "${TRIMMED_DIR}"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")

    R1="${sample_dir}/${sample}_R1.trimmomatic.fastq.gz"
    R2="${sample_dir}/${sample}_R2.trimmomatic.fastq.gz"
    SAM_OUT="${BAM_DIR}/${sample}.sam"
    BAM_OUT="${BAM_DIR}/${sample}.sorted.bam"

    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "Aligning ${sample} ..."
        hisat2 -p "${THREADS}" -x "${INDEX_PREFIX}" -1 "${R1}" -2 "${R2}" -S "${SAM_OUT}"
        samtools view -@ "${THREADS}" -Sb "${SAM_OUT}" | samtools sort -@ "${THREADS}" -o "${BAM_OUT}"
        samtools index "${BAM_OUT}"
        rm -f "${SAM_OUT}"
        BAM_LIST+=("${BAM_OUT}")
    else
        echo "Skipping ${sample}: trimmed read pair not found"
    fi
done

if [[ "${#BAM_LIST[@]}" -eq 0 ]]; then
    echo "Error: no BAM files generated."
    exit 1
fi

echo "Merging BAM files..."
samtools merge -@ "${THREADS}" -o "${MERGED_BAM}" "${BAM_LIST[@]}"
samtools index "${MERGED_BAM}"

echo "RNA-seq alignment completed."
echo "Merged BAM: ${MERGED_BAM}"
