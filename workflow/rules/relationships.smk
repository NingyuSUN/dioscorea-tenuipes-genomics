# Scripts 44 -> 46 and staged 45 -> 47/47b -> 49: 17-individual networks, NJ
# trees, 147 stability panels and relatedness checks.
#
# Layout. RP is script 44's output (and 46 writes RP/relatedness). RS is the
# staging root that scripts 45 and 49 expect: data/<panel>/ holding symlinks to
# the 22 PLINK panels, relatedness -> RP/relatedness, then results/,
# independent_R/ and network_cross_validation.json.

RP = f"{R}/relationship_panels"
RS = f"{R}/relationships"
NEW_PANELS = ["strict17", "complete17"] + [f"loo_{i:02d}" for i in range(17)]
OLD_PANELS = ["supplement17", "main14", "strict14", "complete14"]
PLINK_BASE = "--allow-extra-chr --keep-allele-order --threads 2 --memory 4000"

RELATIONSHIP_TARGETS = [f"{RS}/VALIDATION.json"]


rule rel_prepare_panels:
    input:
        done=rules.pop_filter_variants.output.done,
        masked=rules.pop_filter_variants.output.masked,
    output:
        done=protected(f"{RP}/COMPLETE"),
        manifest=protected(f"{RP}/panels.json"),
        pruned=protected(expand(f"{RP}/{{p}}/pruned.{{e}}", p=NEW_PANELS, e=PRUNED)),
    log:
        f"{LOGS}/44_prepare_relationship_panels.log",
    params:
        py=PY,
        data=PA,
        out=RP,
        script=script("44_prepare_relationship_panels.py"),
    # Three concurrent panels, each PLINK --threads 2 --memory 4000.
    threads: 6
    resources:
        mem_mb=12500,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + fresh(RP) + "{params.py} -u {params.script} --data {params.data} "
        "--out {params.out} --plink " + POPGEN_BIN + "/plink > {log} 2>&1"


rule rel_relatedness:
    input:
        rules.rel_prepare_panels.output.done,
        expand(f"{ST}/supplement17/pruned.{{e}}", e=PRUNED),
    output:
        protected(f"{RP}/relatedness/COMPLETE"),
    log:
        f"{LOGS}/46_relatedness.log",
    params:
        rscript=RSCRIPT,
        rp=RP,
        pa=PA,
        script=script("46_relatedness.R"),
    threads: 2
    resources:
        mem_mb=8000,
    conda:
        "../envs/relationships_r.yaml"
    shell:
        ENV + "{params.rscript} {params.script} {params.rp} {params.pa} > {log} 2>&1"


rule rel_stage:
    """Assemble RS/data from both panel sources and link relatedness."""
    input:
        rules.rel_prepare_panels.output.done,
        rules.rel_relatedness.output,
        rules.pop_structure.output.pruned,
    output:
        touch(f"{RS}/data/STAGED"),
    run:
        os.makedirs(f"{RS}/data", exist_ok=True)
        links = [(f"{RP}/panels.json", f"{RS}/data/panels.json"), (f"{RP}/relatedness", f"{RS}/relatedness")]
        for p in OLD_PANELS:
            links += [(f"{ST}/{p}/pruned.{e}", f"{RS}/data/{p}/pruned.{e}") for e in PRUNED]
        for p in NEW_PANELS:
            links += [(f"{RP}/{p}/{f}", f"{RS}/data/{p}/{f}") for f in ["pruned.bed", "pruned.bim", "pruned.fam", "panel.json", "flat.mdist"]]
        for src, dst in links:
            if not os.path.exists(src):
                raise WorkflowError(f"Missing staging source {src}")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if os.path.lexists(dst):
                if os.path.realpath(dst) != os.path.realpath(src):
                    raise WorkflowError(f"{dst} exists and points elsewhere")
                continue
            os.symlink(src, dst)


rule rel_primary17_flat_distance:
    """Independent PLINK distance for supplement17, compared by script 49."""
    input:
        expand(f"{ST}/supplement17/pruned.{{e}}", e=PRUNED),
        rules.rel_stage.output,
    output:
        f"{RS}/data/primary17_flat.mdist",
    log:
        f"{LOGS}/primary17_flat_distance.log",
    params:
        st=ST,
        prefix=f"{RS}/data/primary17_flat",
    threads: 2
    resources:
        mem_mb=4000,
    conda:
        "../envs/popgen.yaml"
    shell:
        ENV + POPGEN_BIN + "/plink " + PLINK_BASE + " --bfile {params.st}/supplement17/pruned "
        "--distance square 1-ibs flat-missing --out {params.prefix} > {log} 2>&1"


rule rel_analyze:
    input:
        rules.rel_stage.output,
    output:
        stability=protected(f"{RS}/results/panel_stability.tsv"),
        half=protected(f"{RS}/results/half_contig_stability.tsv"),
        subsampling=protected(f"{RS}/results/baseline_split_subsampling.tsv"),
        pairwise=protected(f"{RS}/results/pairwise_qc.tsv"),
        distance=protected(f"{RS}/results/supplement17/distance.tsv"),
        nexus=protected(f"{RS}/results/supplement17/network.nex"),
        splits=protected(f"{RS}/results/supplement17/network_splits.tsv"),
    log:
        f"{LOGS}/45_analyze_relationships.log",
    params:
        apy=APY,
        rs=RS,
        script=script("45_analyze_relationships.py"),
        replicates=config["relationships"]["replicates"],
    threads: 1
    resources:
        mem_mb=8000,
    conda:
        "../envs/relationships_py.yaml"
    shell:
        ENV + "PYTHONHASHSEED=0 {params.apy} {params.script} --root {params.rs} "
        "--replicates {params.replicates} > {log} 2>&1"


rule rel_r_nj_network:
    """ape NJ cross-check for all 147 panels; automatic-order phangorn network."""
    input:
        rules.rel_analyze.output.stability,
    output:
        cross=f"{RS}/results/cross_implementation.tsv",
        splits=f"{RS}/results/supplement17/phangorn_splits.tsv",
    log:
        f"{LOGS}/47_validate_and_plot.log",
    params:
        rscript=RSCRIPT,
        rs=RS,
        script=script("47_validate_and_plot.R"),
    threads: 1
    resources:
        mem_mb=4000,
    conda:
        "../envs/relationships_r.yaml"
    shell:
        ENV + "{params.rscript} {params.script} {params.rs}/results > {log} 2>&1"


rule rel_network_cycle:
    input:
        rules.rel_analyze.output.nexus,
    output:
        f"{RS}/results/supplement17/cycle.txt",
    params:
        apy=APY,
        helper=helper("network_crosscheck.py"),
    conda:
        "../envs/relationships_py.yaml"
    shell:
        ENV + "{params.apy} {params.helper} cycle --nexus {input} --out {output}"


rule rel_r_fixed_order:
    input:
        cycle=rules.rel_network_cycle.output,
        distance=rules.rel_analyze.output.distance,
    output:
        f"{RS}/results/supplement17/phangorn_same_order_splits.tsv",
    log:
        f"{LOGS}/47b_network_fixed_order.log",
    params:
        rscript=RSCRIPT,
        rs=RS,
        script=script("47b_network_fixed_order.R"),
    conda:
        "../envs/relationships_r.yaml"
    shell:
        ENV + "{params.rscript} {params.script} {params.rs}/results/supplement17 > {log} 2>&1"


rule rel_collect_independent_r:
    """Copy the R outputs to independent_R/, where script 49 reads them."""
    input:
        cross=rules.rel_r_nj_network.output.cross,
        auto=rules.rel_r_nj_network.output.splits,
        same=rules.rel_r_fixed_order.output,
    output:
        cross=f"{RS}/independent_R/cross_implementation.tsv",
        auto=f"{RS}/independent_R/supplement17/phangorn_splits.tsv",
        same=f"{RS}/independent_R/supplement17/phangorn_same_order_splits.tsv",
    params:
        rs=RS,
    shell:
        "set -euo pipefail; cp {input.cross} {output.cross}; cp {input.auto} {output.auto}; "
        "cp {input.same} {output.same}; "
        "for f in nj_unrooted.pdf nj_unrooted.png network_phangorn.pdf phangorn_network.rds phangorn_structure.txt; do "
        '[ -e {params.rs}/results/supplement17/$f ] && cp {params.rs}/results/supplement17/$f {params.rs}/independent_R/supplement17/; done; true'


rule rel_network_compare:
    input:
        distance=rules.rel_analyze.output.distance,
        primary=rules.rel_analyze.output.splits,
        same=rules.rel_collect_independent_r.output.same,
        auto=rules.rel_collect_independent_r.output.auto,
    output:
        f"{RS}/network_cross_validation.json",
    params:
        apy=APY,
        helper=helper("network_crosscheck.py"),
    conda:
        "../envs/relationships_py.yaml"
    shell:
        ENV + "{params.apy} {params.helper} compare --distance {input.distance} --primary {input.primary} "
        "--same-order {input.same} --automatic {input.auto} --out {output}"


rule rel_validate:
    input:
        rules.rel_stage.output,
        rules.rel_primary17_flat_distance.output,
        rules.rel_collect_independent_r.output.cross,
        rules.rel_network_compare.output,
        rules.rel_analyze.output.stability,
    output:
        protected(f"{RS}/VALIDATION.json"),
    log:
        f"{LOGS}/49_summarize_validate.log",
    params:
        apy=APY,
        rs=RS,
        script=script("49_summarize_validate.py"),
    threads: 1
    resources:
        mem_mb=8000,
    conda:
        "../envs/relationships_py.yaml"
    shell:
        ENV + "{params.apy} {params.script} --root {params.rs} > {log} 2>&1"
