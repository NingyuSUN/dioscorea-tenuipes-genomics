# Scripts 01-15: historical nuclear assembly, polishing, QC, repeat masking and
# BRAKER annotation. TEMPLATE: the original HPC runs used hand-edited paths and
# an unrecorded environment (see environment.lock.template.yml). This module
# chains the scripts with their documented positional arguments; it has been
# dry-run checked, not executed. Enabled only with config/assembly.yaml.
#
# Not wrapped: organelle scripts 16-20 and 35-43 (current plastome comes from
# the audited GetOrganelle recovery, docs/PLASTOME_ANNOTATION_INDEL_2026-09-14.md),
# server-specific 21-31, and one-off probes 48-52.

A = os.path.abspath(config["assembly"]["outdir"])
AC = config["assembly"]
AT = AC["threads"]
AS = f"{REPO}/scripts"

ASSEMBLY_TARGETS = [
    f"{A}/busco_out/busco_run",
    f"{A}/braker_work/braker/braker.gtf",
    f"{A}/annotation_output/transcripts.gtf",
]
if AC.get("quast_reference_fa"):
    ASSEMBLY_TARGETS.append(f"{A}/quast_out/report.tsv")


rule assembly_all:
    input:
        ASSEMBLY_TARGETS,


rule asm_basecall:
    input:
        AC["fast5_dir"],
    output:
        directory(f"{A}/guppy_out"),
    log:
        f"{A}/logs/01_guppy_basecalling.log",
    params:
        scripts=AS,
    threads: AT["guppy"]
    shell:
        "bash {params.scripts}/01_guppy_basecalling.sh {input} {output} {threads} > {log} 2>&1"


rule asm_combine_pass:
    input:
        rules.asm_basecall.output,
    output:
        f"{A}/combined_pass.fastq.gz",
    log:
        f"{A}/logs/02_combine_pass_reads.log",
    params:
        scripts=AS,
    shell:
        "bash {params.scripts}/02_combine_pass_reads.sh {input} {output} > {log} 2>&1"


rule asm_flye:
    input:
        rules.asm_combine_pass.output,
    output:
        f"{A}/flye_out/assembly.fasta",
    log:
        f"{A}/logs/03_flye_assembly.log",
    params:
        scripts=AS,
        out=f"{A}/flye_out",
        size=AC["genome_size"],
    threads: AT["flye"]
    shell:
        "bash {params.scripts}/03_flye_assembly.sh {input} {params.out} {params.size} {threads} > {log} 2>&1"


rule asm_medaka:
    input:
        draft=rules.asm_flye.output,
        reads=rules.asm_combine_pass.output,
    output:
        f"{A}/medaka_out/consensus.fasta",
    log:
        f"{A}/logs/04_medaka_polish.log",
    params:
        scripts=AS,
        out=f"{A}/medaka_out",
        model=AC.get("medaka_model", ""),
    threads: AT["medaka"]
    shell:
        "bash {params.scripts}/04_medaka_polish.sh {input.draft} {input.reads} {params.out} {threads} {params.model} > {log} 2>&1"


rule asm_trim_illumina:
    input:
        r1=AC["illumina_r1"],
        r2=AC["illumina_r2"],
        adapters=AC["adapters"],
    output:
        r1=f"{A}/illumina_qc/trimmed/trimmed_R1.fastq.gz",
        r2=f"{A}/illumina_qc/trimmed/trimmed_R2.fastq.gz",
    log:
        f"{A}/logs/05_trim_illumina_fastqc.log",
    params:
        scripts=AS,
        out=f"{A}/illumina_qc",
    threads: AT["trimmomatic"]
    shell:
        "bash {params.scripts}/05_trim_illumina_fastqc.sh {input.r1} {input.r2} {params.out} {input.adapters} {threads} > {log} 2>&1"


rule asm_map_for_hypo:
    input:
        ref=rules.asm_medaka.output,
        ont=rules.asm_combine_pass.output,
        r1=rules.asm_trim_illumina.output.r1,
        r2=rules.asm_trim_illumina.output.r2,
    output:
        long=f"{A}/mapping/mapped_lg.sorted.bam",
        short=f"{A}/mapping/mapped_sr.sorted.bam",
        names=f"{A}/mapping/il_names.txt",
    log:
        f"{A}/logs/06_map_reads_for_hypo.log",
    params:
        scripts=AS,
        out=f"{A}/mapping",
    threads: AT["mapping"]
    shell:
        "bash {params.scripts}/06_map_reads_for_hypo.sh {input.ref} {input.ont} {input.r1} {input.r2} {params.out} {threads} > {log} 2>&1"


rule asm_hypo:
    input:
        draft=rules.asm_medaka.output,
        names=rules.asm_map_for_hypo.output.names,
        short=rules.asm_map_for_hypo.output.short,
        long=rules.asm_map_for_hypo.output.long,
    output:
        f"{A}/final_assembly.fa",
    log:
        f"{A}/logs/07_hypo_polish.log",
    params:
        scripts=AS,
        outdir=A,
        size=AC["genome_size"],
        coverage=AC["hypo_coverage"],
    threads: AT["hypo"]
    # Hypo writes its aux/ directory into the working directory.
    shell:
        "cd {params.outdir} && bash {params.scripts}/07_hypo_polish.sh {input.draft} {input.names} {input.short} {input.long} "
        "{output} {params.size} {params.coverage} {threads} > {log} 2>&1"


rule asm_busco:
    input:
        rules.asm_hypo.output,
    output:
        directory(f"{A}/busco_out/busco_run"),
    log:
        f"{A}/logs/08_busco_assessment.log",
    params:
        scripts=AS,
        lineage=AC["busco_lineage"],
        out=f"{A}/busco_out",
    # BUSCO genome mode is memory-hungry; keep threads low on ~32 GB machines.
    threads: AT["busco"]
    shell:
        "bash {params.scripts}/08_busco_assessment.sh {input} {params.lineage} {params.out} {threads} > {log} 2>&1"


rule asm_quast:
    input:
        flye=rules.asm_flye.output,
        medaka=rules.asm_medaka.output,
        ref=AC.get("quast_reference_fa") or [],
        gff=AC.get("quast_reference_gff") or [],
    output:
        f"{A}/quast_out/report.tsv",
    log:
        f"{A}/logs/09_quast_assessment.log",
    params:
        scripts=AS,
        out=f"{A}/quast_out",
    threads: AT["quast"]
    shell:
        "bash {params.scripts}/09_quast_assessment.sh {input.flye} {input.medaka} {input.ref} {input.gff} {params.out} {threads} > {log} 2>&1"


rule asm_repeatmasker:
    input:
        rules.asm_hypo.output,
    output:
        f"{A}/repeatmasker_out/final_assembly.fa.masked",
    log:
        f"{A}/logs/10_te_repeatmasker.log",
    params:
        scripts=AS,
        out=f"{A}/repeatmasker_out",
        species=AC["repeatmasker_species"],
    threads: AT["repeatmasker"]
    shell:
        "bash {params.scripts}/10_te_repeatmasker.sh {input} {params.out} {params.species} {threads} > {log} 2>&1"


rule asm_rnaseq_fastqc_raw:
    input:
        AC["rnaseq_dir"],
    output:
        directory(f"{A}/RNA_Seq/fastqc_raw"),
    log:
        f"{A}/logs/11_rnaseq_fastqc_raw.log",
    params:
        scripts=AS,
    threads: AT["fastqc"]
    shell:
        "bash {params.scripts}/11_rnaseq_fastqc_raw.sh {input} {output} {threads} > {log} 2>&1"


rule asm_rnaseq_trim:
    input:
        reads=AC["rnaseq_dir"],
        adapters=AC["adapters"],
        raw_qc=rules.asm_rnaseq_fastqc_raw.output,
    output:
        trimmed=directory(f"{A}/RNA_Seq/trimmomatic"),
        qc=directory(f"{A}/RNA_Seq/fastqc_trimmomatic"),
    log:
        f"{A}/logs/12_rnaseq_trim_fastqc.log",
    params:
        scripts=AS,
        jobs=AC["rnaseq_parallel_samples"],
        per_job=max(1, AT["trimmomatic"] // AC["rnaseq_parallel_samples"]),
    threads: AT["trimmomatic"]
    shell:
        "bash {params.scripts}/12_rnaseq_trim_fastqc.sh {input.reads} {output.trimmed} {output.qc} {input.adapters} "
        "{params.per_job} {params.jobs} > {log} 2>&1"


rule asm_hisat2:
    input:
        genome=rules.asm_repeatmasker.output,
        trimmed=rules.asm_rnaseq_trim.output.trimmed,
    output:
        f"{A}/annotation_output/rna_aligned.sorted.bam",
    log:
        f"{A}/logs/13_rnaseq_align_hisat2.log",
    params:
        scripts=AS,
        out=f"{A}/annotation_output",
    threads: AT["hisat2"]
    shell:
        "bash {params.scripts}/13_rnaseq_align_hisat2.sh {input.genome} {input.trimmed} {params.out} {threads} > {log} 2>&1"


rule asm_stringtie:
    input:
        rules.asm_hisat2.output,
    output:
        f"{A}/annotation_output/transcripts.gtf",
    log:
        f"{A}/logs/14_stringtie_transcripts.log",
    params:
        scripts=AS,
    threads: AT["stringtie"]
    shell:
        "bash {params.scripts}/14_stringtie_transcripts.sh {input} {output} {threads} > {log} 2>&1"


rule asm_braker:
    input:
        genome=rules.asm_repeatmasker.output,
        bam=rules.asm_hisat2.output,
        proteins=AC["proteins"],
    output:
        f"{A}/braker_work/braker/braker.gtf",
    log:
        f"{A}/logs/15_braker_annotation.log",
    params:
        scripts=AS,
        proteins=",".join(os.path.abspath(p) for p in AC["proteins"]),
        species=AC["braker_species"],
        work=f"{A}/braker_work",
    threads: AT["braker"]
    shell:
        "WORK_DIR={params.work} bash {params.scripts}/15_braker_annotation.sh {input.genome} {input.bam} "
        "{params.proteins} {params.species} {threads} > {log} 2>&1"
