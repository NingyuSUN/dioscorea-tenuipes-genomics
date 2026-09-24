# Scripts 53 -> 55 -> 54: 38 PCA panels and 114 FST estimates for sample and
# candidate-relative influence.

IN = f"{R}/influence"

INFLUENCE_TARGETS = [f"{R}/influence_summary/SUMMARY_COMPLETE"]


rule infl_pca_panels:
    input:
        rules.pop_filter_variants.output.done,
        rules.pop_structure.output.done,
    output:
        done=protected(f"{IN}/PCA_COMPLETE"),
        panels=protected(f"{IN}/panels.json"),
        bed=protected(multiext(f"{IN}/source_own14", ".bed", ".bim", ".fam")),
    log:
        f"{LOGS}/53_run_influence_panels.log",
    params:
        py=PY,
        data=PA,
        out=IN,
        script=script("53_run_influence_panels.py"),
    # Four concurrent jobs, each PLINK --threads 2 --memory 3500.
    threads: 8
    resources:
        mem_mb=14500,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(IN) + "{params.py} -u {params.script} --data {params.data} "
        "--out {params.out} --tools " + POPGEN_BIN + " > {log} 2>&1"


rule infl_fst:
    input:
        bed=rules.infl_pca_panels.output.bed,
        panels=rules.infl_pca_panels.output.panels,
    output:
        done=protected(f"{IN}/fst/FST_COMPLETE"),
        summary=protected(f"{IN}/fst/fst_summary.tsv"),
    log:
        f"{LOGS}/55_compute_influence_fst.log",
    params:
        py=PY,
        prefix=f"{IN}/source_own14",
        out=f"{IN}/fst",
        script=script("55_compute_influence_fst.py"),
    threads: 1
    resources:
        mem_mb=6000,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(f"{IN}/fst") + "{params.py} -u {params.script} --bed-prefix {params.prefix} "
        "--panels {input.panels} --out {params.out} > {log} 2>&1"


rule infl_stage_fst:
    """Script 54 reads fst_summary.tsv from its --root, next to panels.json."""
    input:
        rules.infl_fst.output.summary,
    output:
        f"{IN}/fst_summary.tsv",
    shell:
        "ln -s fst/fst_summary.tsv {output}"


rule infl_summarize:
    input:
        rules.infl_pca_panels.output.done,
        rules.infl_stage_fst.output,
    output:
        protected(f"{R}/influence_summary/SUMMARY_COMPLETE"),
    log:
        f"{LOGS}/54_summarize_influence.log",
    params:
        apy=APY,
        root=IN,
        out=f"{R}/influence_summary",
        script=script("54_summarize_influence.py"),
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/relationships_py.yaml"
    shell:
        ENV + fresh(f"{R}/influence_summary") + "{params.apy} {params.script} --root {params.root} "
        "--out {params.out} > {log} 2>&1"
