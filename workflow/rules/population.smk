# Scripts 32 -> 33 -> 34: masked panels, PCA / ADMIXTURE / WC84 FST, report.
# Filters and limitations: docs/RUN_POPULATION_ANALYSIS.md.

PA = f"{R}/population_analysis"
ST = f"{PA}/structure"
BCF_PANELS = ["main14", "strict14", "complete14", "supplement17"]
SAMPLE_SETS = ["main14", "supplement17", "AIC", "HCO", "IZU"]
PRUNED = ["bed", "bim", "fam"]

POPULATION_TARGETS = [f"{PA}/report/SUMMARY_COMPLETE"]


rule pop_filter_variants:
    input:
        vcf=config["inputs"]["population_vcf"],
        manifest=config["inputs"]["sample_manifest"],
    output:
        done=protected(f"{PA}/FILTER_COMPLETE"),
        masked=protected(f"{PA}/masked17.bcf"),
        panels=protected(expand(f"{PA}/{{p}}.bcf", p=BCF_PANELS)),
        samples=protected(expand(f"{PA}/{{s}}.samples.txt", s=SAMPLE_SETS)),
    log:
        f"{LOGS}/32_filter_population_variants.log",
    params:
        py=PY,
        out=PA,
        script=script("32_filter_population_variants.py"),
    threads: 3
    resources:
        mem_mb=4000,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(PA) + "{params.py} -u {params.script} --vcf {input.vcf} "
        "--manifest {input.manifest} --out {params.out} > {log} 2>&1"


rule pop_structure:
    input:
        done=rules.pop_filter_variants.output.done,
    output:
        done=protected(f"{ST}/STRUCTURE_COMPLETE"),
        pca=protected(f"{ST}/PCA_COMPLETE"),
        pruned=protected(expand(f"{ST}/{{p}}/pruned.{{e}}", p=BCF_PANELS, e=PRUNED)),
        main14_pca=protected(expand(f"{ST}/main14/pca.{{e}}", e=["eigenvec", "eigenval", "log"])),
    log:
        f"{LOGS}/33_run_population_structure.log",
    params:
        py=PY,
        data=PA,
        out=ST,
        script=script("33_run_population_structure.py"),
    # Up to six concurrent jobs: ADMIXTURE -j2 or PLINK --threads 2 at 5000 MB.
    threads: 12
    resources:
        mem_mb=16000,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(ST) + "{params.py} -u {params.script} --data {params.data} "
        "--out {params.out} --tools " + POPGEN_BIN + " > {log} 2>&1"


rule pop_report:
    input:
        done=rules.pop_structure.output.done,
        audit=config["inputs"]["audit_package"],
    output:
        protected(f"{PA}/report/SUMMARY_COMPLETE"),
    log:
        f"{LOGS}/34_summarize_population_results.log",
    params:
        py=PY,
        data=PA,
        analysis=ST,
        out=f"{PA}/report",
        script=script("34_summarize_population_results.py"),
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(f"{PA}/report") + "{params.py} {params.script} --data {params.data} "
        "--analysis {params.analysis} --audit {input.audit} --out {params.out} > {log} 2>&1"
