#!/bin/bash
set -euo pipefail

# ============================================================
# 17_cp_consensus_from_bam.sh
# Build a chloroplast consensus sequence from a reference-
# guided BAM using bcftools and report low-coverage sites.
#
# Usage:
#   bash 17_cp_consensus_from_bam.sh <cp_reference.fa> <sorted.bam> <output_prefix> [min_depth]
#
# Example:
#   bash 17_cp_consensus_from_bam.sh D_tokoro_cp.fa dtenuipes_cp.sorted.bam tenuipes_cp 10
# ============================================================

REFERENCE="${1:-D_tokoro_cp.fa}"
BAM="${2:-dtenuipes_cp.sorted.bam}"
PREFIX="${3:-tenuipes_cp}"
MIN_DEPTH="${4:-10}"

BCF="${PREFIX}.variants.bcf"
VCF_GZ="${PREFIX}.variants.vcf.gz"
CONSENSUS="${PREFIX}.consensus.fasta"
COVERAGE="${PREFIX}.coverage.txt"
LOW_COV="${PREFIX}.low_coverage_regions.txt"

command -v bcftools >/dev/null 2>&1 || { echo "Error: bcftools not found."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "Error: tabix not found."; exit 1; }

[[ -f "${REFERENCE}" ]] || { echo "Error: reference not found: ${REFERENCE}"; exit 1; }
[[ -f "${BAM}" ]] || { echo "Error: BAM not found: ${BAM}"; exit 1; }

echo "Generating pileup and calling variants..."
bcftools mpileup -f "${REFERENCE}" "${BAM}" -Ou |     bcftools call -mv -Oz -o "${VCF_GZ}"

tabix -p vcf "${VCF_GZ}"

echo "Generating consensus sequence..."
bcftools consensus -f "${REFERENCE}" -o "${CONSENSUS}" "${VCF_GZ}"

echo "Checking low-coverage sites..."
samtools depth -a "${BAM}" > "${COVERAGE}"
awk -v min="${MIN_DEPTH}" '$3 < min {print $1 "\t" $2 "\t" $3}' "${COVERAGE}" > "${LOW_COV}"

echo "Done."
echo "Consensus chloroplast genome: ${CONSENSUS}"
echo "Low coverage report: ${LOW_COV}"
