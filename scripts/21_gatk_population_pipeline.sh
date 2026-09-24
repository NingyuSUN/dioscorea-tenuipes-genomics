\
#!/usr/bin/env bash
set -euo pipefail

# =========================================================
# Script: 21_gatk_population_pipeline.sh
# Purpose:
#   Population variant calling and downstream SNP processing
#   for Dioscorea tenuipes using BWA + SAMtools + GATK.
#
# Main steps:
#   1. Map trimmed reads to reference
#   2. Sort/index BAM
#   3. Mark duplicates
#   4. Per-sample GVCF calling
#   5. GenomicsDBImport
#   6. Joint genotyping in interval batches
#   7. Merge batch VCFs
#   8. Variant filtration
#   9. Select filtered biallelic SNPs
#   10. STRUCTURE / phylogeny export
#
# Notes:
#   - Designed for HPC/SLURM environments
#   - Set WORK_DIR and REFERENCE in the environment before running
# =========================================================

# -----------------------------
# Configuration
# -----------------------------
WORK_DIR="${WORK_DIR:?Set WORK_DIR to the directory holding trimmomatic/<sample>/ reads}"
REFERENCE="${REFERENCE:?Set REFERENCE to the D. tenuipes reference FASTA}"
TRIMMED_DIR="$WORK_DIR/trimmomatic"

SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-8}"
THREADS="${THREADS:-$((SLURM_CPUS_PER_TASK * 2))}"
JAVA_MEM="${JAVA_MEM:-100G}"
MAX_JOBS="${MAX_JOBS:-8}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
GATK_OUTPUT_DIR="$WORK_DIR/gatk_run_${TIMESTAMP}"
LOG_FILE="$GATK_OUTPUT_DIR/pipeline_progress.log"

mkdir -p \
  "$GATK_OUTPUT_DIR"/{mapped_sorted,metrics,gvcf,tmp,genomics_db,joint_calls,filtered_calls,structure,phylogeny}

log() {
    echo "[$(date +'%F %T')] $*" | tee -a "$LOG_FILE"
}

check_tools() {
    local tools=(bwa samtools gatk vcftools plink)
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "ERROR: Required tool '$tool' not found."
            exit 1
        fi
    done

    if ! command -v vcf2phylip.py >/dev/null 2>&1 && ! command -v vcf2phylip >/dev/null 2>&1; then
        log "ERROR: vcf2phylip.py (or vcf2phylip) not found."
        exit 1
    fi

    if ! command -v raxmlHPC-PTHREADS-AVX2 >/dev/null 2>&1 && ! command -v raxmlHPC >/dev/null 2>&1; then
        log "ERROR: RAxML executable not found."
        exit 1
    fi
}

prepare_reference_indexes() {
    log "Preparing reference indexes if needed..."

    if [[ ! -f "${REFERENCE}.bwt" ]]; then
        log "Building BWA index..."
        bwa index "$REFERENCE"
    fi

    if [[ ! -f "${REFERENCE}.fai" ]]; then
        log "Building FASTA index..."
        samtools faidx "$REFERENCE"
    fi

    local dict="${REFERENCE%.*}.dict"
    if [[ ! -f "$dict" ]]; then
        log "Building sequence dictionary..."
        gatk CreateSequenceDictionary -R "$REFERENCE" -O "$dict"
    fi
}

create_interval_list() {
    log "Creating genome interval list from reference FASTA index..."
    awk '{print $1 ":1-" $2}' "${REFERENCE}.fai" > "$GATK_OUTPUT_DIR/genome_intervals_fixed.list"

    if [[ ! -s "$GATK_OUTPUT_DIR/genome_intervals_fixed.list" ]]; then
        log "ERROR: Interval list is empty."
        exit 1
    fi
}

process_gatk() {
    local sample="$1"
    local fq1="$TRIMMED_DIR/$sample/${sample}_R1.trimmomatic.fastq.gz"
    local fq2="$TRIMMED_DIR/$sample/${sample}_R2.trimmomatic.fastq.gz"
    local sorted_bam="$GATK_OUTPUT_DIR/mapped_sorted/${sample}.sorted.bam"
    local dedup_bam="$GATK_OUTPUT_DIR/mapped_sorted/${sample}.dedup.bam"
    local gvcf_out="$GATK_OUTPUT_DIR/gvcf/${sample}.g.vcf.gz"

    log "[START] Processing sample: $sample"

    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        log "ERROR: Missing FASTQ files for sample $sample. Skipping."
        return 1
    fi

    log "Aligning reads for $sample"
    bwa mem -t "$THREADS" -Y -R "@RG\tID:${sample}\tSM:${sample}\tPL:Illumina" \
        "$REFERENCE" "$fq1" "$fq2" \
        | samtools sort -@ "$THREADS" -m 2G -o "$sorted_bam" -

    samtools index "$sorted_bam"

    log "Marking duplicates for $sample"
    gatk --java-options "-Xmx${JAVA_MEM}" MarkDuplicates \
        -I "$sorted_bam" \
        -O "$dedup_bam" \
        -M "$GATK_OUTPUT_DIR/metrics/${sample}.metrics.txt" \
        --TMP_DIR "$GATK_OUTPUT_DIR/tmp"

    samtools index "$dedup_bam"

    log "Running HaplotypeCaller for $sample"
    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$dedup_bam" \
        -O "$gvcf_out" \
        -ERC GVCF \
        --native-pair-hmm-threads "$THREADS"

    log "[DONE] Sample $sample completed"
}

run_per_sample_gvcf() {
    log "Starting per-sample GVCF generation..."

    local current_jobs=0
    local sample

    while IFS= read -r sample; do
        process_gatk "$sample" &
        current_jobs=$((current_jobs + 1))

        if [[ "$current_jobs" -ge "$MAX_JOBS" ]]; then
            wait -n
            current_jobs=$((current_jobs - 1))
        fi
    done < <(find "$TRIMMED_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)

    wait

    if [[ $(find "$GATK_OUTPUT_DIR/gvcf" -name "*.g.vcf.gz" | wc -l) -eq 0 ]]; then
        log "ERROR: No GVCFs were generated."
        exit 1
    fi

    log "All GVCFs are complete."
}

create_sample_map() {
    log "Creating sample map..."
    find "$GATK_OUTPUT_DIR/gvcf" -name "*.g.vcf.gz" \
        | sort \
        | awk -F'/' '{sample=$NF; sub(/\.g\.vcf\.gz$/, "", sample); print sample "\t" $0}' \
        > "$GATK_OUTPUT_DIR/sample_map.txt"

    if [[ ! -s "$GATK_OUTPUT_DIR/sample_map.txt" ]]; then
        log "ERROR: sample_map.txt is empty."
        exit 1
    fi
}

run_genomicsdbimport() {
    log "Step 1: GenomicsDBImport - STARTED"
    gatk --java-options "-Xmx${JAVA_MEM}" GenomicsDBImport \
        --genomicsdb-workspace-path "$GATK_OUTPUT_DIR/genomics_db" \
        --sample-name-map "$GATK_OUTPUT_DIR/sample_map.txt" \
        --batch-size 5 \
        -L "$GATK_OUTPUT_DIR/genome_intervals_fixed.list" \
        --tmp-dir "$GATK_OUTPUT_DIR/tmp"
    log "Step 1: GenomicsDBImport - COMPLETED"
}

run_joint_genotyping_batches() {
    log "Step 2: Joint Genotyping - STARTED"

    split -l 100 "$GATK_OUTPUT_DIR/genome_intervals_fixed.list" "$GATK_OUTPUT_DIR/interval_batch_"

    local batch
    for batch in "$GATK_OUTPUT_DIR"/interval_batch_*; do
        [[ -f "$batch" ]] || continue
        local base
        base="$(basename "$batch")"
        log "Processing batch: $base"

        gatk --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs \
            -R "$REFERENCE" \
            -V "gendb://$GATK_OUTPUT_DIR/genomics_db" \
            -L "$batch" \
            -O "$GATK_OUTPUT_DIR/joint_calls/${base}.vcf.gz"
    done

    log "Step 2: Joint Genotyping - COMPLETED"
}

merge_joint_calls() {
    log "Step 3: Merging Joint Calls - STARTED"

    mapfile -t vcfs < <(find "$GATK_OUTPUT_DIR/joint_calls" -name "*.vcf.gz" | sort)
    if [[ ${#vcfs[@]} -eq 0 ]]; then
        log "ERROR: No batch VCFs found for merging."
        exit 1
    fi

    local inputs=()
    local v
    for v in "${vcfs[@]}"; do
        inputs+=("-I" "$v")
    done

    gatk --java-options "-Xmx${JAVA_MEM}" MergeVcfs \
        "${inputs[@]}" \
        -O "$GATK_OUTPUT_DIR/joint_calls_merged.vcf.gz"

    log "Step 3: Merging Joint Calls - COMPLETED"
}

filter_variants() {
    log "Step 4: Variant Filtration - STARTED"
    gatk --java-options "-Xmx${JAVA_MEM}" VariantFiltration \
        -R "$REFERENCE" \
        -V "$GATK_OUTPUT_DIR/joint_calls_merged.vcf.gz" \
        -O "$GATK_OUTPUT_DIR/filtered_calls/filtered_calls.vcf.gz" \
        --filter-name "QD2" --filter-expression "QD < 2.0" \
        --filter-name "FS60" --filter-expression "FS > 60.0" \
        --filter-name "MQ30" --filter-expression "MQ < 30.0" \
        --filter-name "SOR4" --filter-expression "SOR > 4.0" \
        --filter-name "ReadPosRankSum-6" --filter-expression "ReadPosRankSum < -6.0" \
        --filter-name "MQRankSum-10" --filter-expression "MQRankSum < -10.0"
    log "Step 4: Variant Filtration - COMPLETED"

    log "Step 5: Selecting High-Quality Biallelic SNPs - STARTED"
    gatk --java-options "-Xmx${JAVA_MEM}" SelectVariants \
        -R "$REFERENCE" \
        -V "$GATK_OUTPUT_DIR/filtered_calls/filtered_calls.vcf.gz" \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --exclude-filtered true \
        -O "$GATK_OUTPUT_DIR/final_filtered_biallelic.vcf.gz"
    log "Step 5: Selecting High-Quality Biallelic SNPs - COMPLETED"
}

export_structure_and_phylogeny_inputs() {
    log "Step 6: Converting VCF to STRUCTURE format - STARTED"
    vcftools --gzvcf "$GATK_OUTPUT_DIR/final_filtered_biallelic.vcf.gz" \
        --plink \
        --out "$GATK_OUTPUT_DIR/structure/structure_data"

    plink --file "$GATK_OUTPUT_DIR/structure/structure_data" \
        --recode structure \
        --out "$GATK_OUTPUT_DIR/structure/structure_input"
    log "Step 6: Converting VCF to STRUCTURE format - COMPLETED"

    log "Step 7: LD pruning before phylogeny - STARTED"
    plink --vcf "$GATK_OUTPUT_DIR/final_filtered_biallelic.vcf.gz" \
        --indep-pairwise 50 10 0.2 \
        --out "$GATK_OUTPUT_DIR/phylogeny/pruned"

    plink --vcf "$GATK_OUTPUT_DIR/final_filtered_biallelic.vcf.gz" \
        --extract "$GATK_OUTPUT_DIR/phylogeny/pruned.prune.in" \
        --recode vcf \
        --out "$GATK_OUTPUT_DIR/phylogeny/pruned_filtered"
    log "Step 7: LD pruning - COMPLETED"

    log "Step 8: Converting VCF to PHYLIP - STARTED"
    if command -v vcf2phylip.py >/dev/null 2>&1; then
        vcf2phylip.py -i "$GATK_OUTPUT_DIR/phylogeny/pruned_filtered.vcf" \
            -o "$GATK_OUTPUT_DIR/phylogeny/phylip_matrix.phy"
    else
        vcf2phylip -i "$GATK_OUTPUT_DIR/phylogeny/pruned_filtered.vcf" \
            -o "$GATK_OUTPUT_DIR/phylogeny/phylip_matrix.phy"
    fi
    log "Step 8: Converting VCF to PHYLIP - COMPLETED"

    log "Step 9: Running RAxML - STARTED"
    if command -v raxmlHPC-PTHREADS-AVX2 >/dev/null 2>&1; then
        raxmlHPC-PTHREADS-AVX2 \
            -T 32 \
            -s "$GATK_OUTPUT_DIR/phylogeny/phylip_matrix.phy" \
            -n tenuipes_tree \
            -m GTRGAMMA \
            -p 12345 \
            -# 100 \
            -w "$GATK_OUTPUT_DIR/phylogeny"
    else
        raxmlHPC \
            -s "$GATK_OUTPUT_DIR/phylogeny/phylip_matrix.phy" \
            -n tenuipes_tree \
            -m GTRGAMMA \
            -p 12345 \
            -# 100 \
            -w "$GATK_OUTPUT_DIR/phylogeny"
    fi
    log "Step 9: RAxML - COMPLETED"
}

main() {
    check_tools
    prepare_reference_indexes
    create_interval_list
    run_per_sample_gvcf
    create_sample_map
    run_genomicsdbimport
    run_joint_genotyping_batches
    merge_joint_calls
    filter_variants
    export_structure_and_phylogeny_inputs

    log "Pipeline completed successfully."
    log "Final SNP dataset: $GATK_OUTPUT_DIR/final_filtered_biallelic.vcf.gz"
}

main "$@"
