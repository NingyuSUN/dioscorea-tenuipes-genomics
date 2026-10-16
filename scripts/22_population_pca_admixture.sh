\
#!/usr/bin/env bash
# HISTORICAL: chromosome-label replacement is not a coordinate transformation. Use scripts 32–34 for the current contig-coordinate workflow.
set -euo pipefail

# =========================================================
# Script: 22_population_pca_admixture.sh
# Purpose:
#   Convert filtered SNP VCF to PLINK format, project contigs
#   to pseudo-chromosomes, run PCA and ADMIXTURE.
#
# Notes:
#   - Set UPDIR, REF_GENOME and CONTIGS in the environment before running
#   - Requires a filtered VCF from the GATK pipeline
# =========================================================

UPDIR="${UPDIR:?Set UPDIR to the directory holding the filtered VCF}"
WORKDIR="${UPDIR}/contigs_to_chromosome"
REF_GENOME="${REF_GENOME:?Set REF_GENOME to the D. tokoro Kita1 FASTA}"
CONTIGS="${CONTIGS:?Set CONTIGS to the D. tenuipes contig FASTA}"
VCF_FILE="${VCF_FILE:-${UPDIR}/final_filtered_biallelic_renamed.vcf.gz}"

PLINK_DIR="${WORKDIR}/plink"
ADMIXTURE_DIR="${WORKDIR}/admixture"

mkdir -p "${ADMIXTURE_DIR}" "${WORKDIR}" "${PLINK_DIR}"

log() {
    echo "[$(date +'%F %T')] $*"
}

check_tools() {
    local tools=(minimap2 plink admixture awk sort)
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "ERROR: Required tool '$tool' not found."
            exit 1
        fi
    done
}

check_inputs() {
    for f in "$REF_GENOME" "$CONTIGS" "$VCF_FILE"; do
        if [[ ! -f "$f" ]]; then
            log "ERROR: Missing input file: $f"
            exit 1
        fi
    done
}

align_contigs() {
    log "Step 1: Aligning contigs to reference genome..."
    cd "$WORKDIR"

    if [[ ! -f "sister_species.mmi" ]]; then
        minimap2 -d sister_species.mmi "$REF_GENOME"
    fi

    minimap2 -ax asm5 sister_species.mmi "$CONTIGS" > "${WORKDIR}/contigs_to_ref.sam"

    if [[ ! -s "${WORKDIR}/contigs_to_ref.sam" ]]; then
        log "ERROR: Alignment failed or SAM file is empty."
        exit 1
    fi
}

extract_contig_mapping() {
    log "Step 2: Extracting contig-to-chromosome mapping..."
    awk 'BEGIN{OFS="\t"} !/^@/ && $3 != "*" {print $1, $3, $4}' "${WORKDIR}/contigs_to_ref.sam" \
        | sort -k1,1 -k3,3n \
        | awk '!seen[$1]++' \
        > "${WORKDIR}/contigs_to_chromosome.txt"

    if [[ ! -s "${WORKDIR}/contigs_to_chromosome.txt" ]]; then
        log "ERROR: No contig-to-chromosome mappings were extracted."
        exit 1
    fi
}

make_renaming_table() {
    log "Step 3: Creating pseudo-chromosome renaming table..."
    awk 'BEGIN{OFS="\t"} {print $1, "Pseudo_"$2}' "${WORKDIR}/contigs_to_chromosome.txt" \
        > "${WORKDIR}/contig_renaming.txt"
}

vcf_to_plink() {
    log "Step 4: Converting VCF to PLINK..."
    plink --vcf "$VCF_FILE" --make-bed --out "${PLINK_DIR}/population" --allow-extra-chr
}

update_bim() {
    log "Step 5: Updating chromosome labels in BIM file..."
    awk '
        BEGIN{OFS="\t"}
        NR==FNR {map[$1]=$2; next}
        {
            if ($1 in map) $1=map[$1]
            else if ($1 ~ /^Pseudo_/) $1="11"
            print
        }
    ' "${WORKDIR}/contig_renaming.txt" "${PLINK_DIR}/population.bim" > "${PLINK_DIR}/population_fixed.bim"

    mv "${PLINK_DIR}/population_fixed.bim" "${PLINK_DIR}/population.bim"
}

verify_labels() {
    log "Step 6: Verifying chromosome labels..."
    awk '{print $1}' "${PLINK_DIR}/population.bim" | sort | uniq -c
}

ld_prune_and_pca() {
    log "Step 7: LD pruning..."
    plink --bfile "${PLINK_DIR}/population" \
        --indep-pairwise 50 5 0.2 \
        --out "${PLINK_DIR}/population_pruned" \
        --allow-extra-chr

    log "Step 8: Running PCA..."
    plink --bfile "${PLINK_DIR}/population" \
        --extract "${PLINK_DIR}/population_pruned.prune.in" \
        --pca \
        --out "${PLINK_DIR}/population_pruned_pca" \
        --allow-extra-chr
}

run_admixture() {
    log "Step 9: Running ADMIXTURE for K=2 to 5..."
    cd "${ADMIXTURE_DIR}"

    for K in 2 3 4 5; do
        mkdir -p "${ADMIXTURE_DIR}/K${K}"
        admixture --cv "${PLINK_DIR}/population.bed" "$K" 2>&1 | tee "${ADMIXTURE_DIR}/K${K}/log${K}.out"

        if compgen -G "${PLINK_DIR}/population.${K}.*" > /dev/null; then
            mv "${PLINK_DIR}/population.${K}."* "${ADMIXTURE_DIR}/K${K}/"
        else
            log "Warning: No ADMIXTURE output files found for K=${K}"
        fi
    done

    log "Step 10: Collecting cross-validation errors..."
    grep -h "CV error" "${ADMIXTURE_DIR}"/K*/log*.out > "${ADMIXTURE_DIR}/cross_validation_errors.txt" || true
}

main() {
    check_tools
    check_inputs
    align_contigs
    extract_contig_mapping
    make_renaming_table
    vcf_to_plink
    update_bim
    verify_labels
    ld_prune_and_pca
    run_admixture
    log "Pipeline completed successfully."
}

main "$@"
