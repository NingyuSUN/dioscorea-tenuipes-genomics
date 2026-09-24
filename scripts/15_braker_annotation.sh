#!/bin/bash
set -euo pipefail

# ============================================================
# 15_braker_annotation.sh
# Run BRAKER gene annotation using RNA-seq BAM + protein evidence
#
# Usage:
#   bash 15_braker_annotation.sh <softmasked_genome.fa> <rna_bam> <protein1.fa,protein2.fa,...> [species_name] [threads]
#
# Example:
#   bash 15_braker_annotation.sh \
#       repeatmasker_out/D_tenuipes.fa.masked \
#       annotation_output/rna_aligned.sorted.bam \
#       Dalata.faa,Dtokoro.fasta,Drotundata.faa \
#       dioscorea_tenuipes \
#       32
#
# Notes:
# - Protein files are comma-separated and will be merged automatically
# - Genome should be softmasked for --softmasking
# - D. tokoro can be included as one of the protein evidence sources
# - BRAKER writes outputs to the current working directory unless configured otherwise
# ============================================================

GENOME="${1:-D_tenuipes.softmasked.fa}"
RNA_BAM="${2:-annotation_output/rna_aligned.sorted.bam}"
PROTEIN_LIST="${3:-Dalata.faa,Dtokoro.fasta,Drotundata.faa}"
SPECIES_NAME="${4:-dioscorea_tenuipes}"
THREADS="${5:-32}"

WORK_DIR="${WORK_DIR:-braker_work}"
MERGED_PROTEIN="${WORK_DIR}/protein_evidence.fa"

mkdir -p "${WORK_DIR}"

IFS=',' read -r -a PROT_FILES <<< "${PROTEIN_LIST}"

: > "${MERGED_PROTEIN}"
for prot in "${PROT_FILES[@]}"; do
    if [[ -f "${prot}" ]]; then
        cat "${prot}" >> "${MERGED_PROTEIN}"
    else
        echo "Warning: protein file not found, skipping: ${prot}"
    fi
done

if [[ ! -s "${MERGED_PROTEIN}" ]]; then
    echo "Error: no valid protein evidence files found."
    exit 1
fi

cd "${WORK_DIR}"

braker.pl \
    --species="${SPECIES_NAME}" \
    --genome="${GENOME}" \
    --bam="${RNA_BAM}" \
    --prot_seq="${MERGED_PROTEIN}" \
    --softmasking \
    --threads "${THREADS}"

echo "BRAKER annotation completed."
echo "BRAKER working directory: ${WORK_DIR}"
