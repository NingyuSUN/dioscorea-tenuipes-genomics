# Prepared-input wrappers for selected later-analysis entry points.
# Historical review scripts remain independent: workflow/research/README.md.
import json
import sys
sys.path.insert(0, os.path.join(REPO, "workflow", "scripts"))
from research_contracts import nuclear_inputs, cp_ils_inputs, cp_tree_inputs

RC = config["research"]
RR = os.path.realpath(RC["outdir"])
RPY = RC.get("python") or "python"
RSCRIPTS = os.path.join(REPO, "scripts", "research")
RN = f"{RR}/nuclear"
RI = f"{RR}/chloroplast_ils"
PHASE = RC.get("chloroplast_ils_phase", "full")
if PHASE not in ("pilot", "full"):
    raise WorkflowError("chloroplast_ils_phase must be pilot or full")


def plan_files(plan, key):
    # Read declared plan metadata only. No genotypes are read during DAG creation.
    try:
        with open(plan) as handle:
            prefix = json.load(handle)[key]
    except (OSError, KeyError, ValueError) as exc:
        raise WorkflowError(f"Cannot read prepared population plan {plan}: {exc}")
    return [prefix + "." + ext for ext in ("bed", "bim", "fam")]


RESEARCH_TARGETS = {
    "nuclear": [f"{RN}/concordance_independent/SUMMARY.json"],
    "population": [f"{RR}/spatial_fst/RUN_RECEIPT.json", f"{RR}/balanced_fst/RUN_RECEIPT.json"],
    "chloroplast_ils": [f"{RI}/results_{PHASE}/VALIDATION.json"],
    "chloroplast_trees": [f"{RR}/chloroplast_trees/WORKFLOW_COMPLETE.json"],
}


rule research_all:
    input:
        [p for module in RC.get("modules", []) for p in RESEARCH_TARGETS[module]],


rule research_nuclear_all:
    input: RESEARCH_TARGETS["nuclear"]


rule research_population_all:
    input: RESEARCH_TARGETS["population"]


rule research_chloroplast_ils_all:
    input: RESEARCH_TARGETS["chloroplast_ils"]


rule research_chloroplast_trees_all:
    input: RESEARCH_TARGETS["chloroplast_trees"]


rule research_stage_nuclear:
    input:
        files=lambda wc: nuclear_inputs(RC["nuclear_input_root"]),
        stage=helper("stage_research_inputs.py"),
        contract=helper("research_contracts.py"),
    output: protected(f"{RN}/STAGED.json")
    params:
        python=RPY, stage=helper("stage_research_inputs.py"),
        source=RC["nuclear_input_root"], dest=RN,
    shell:
        "{params.python:q} {params.stage:q} --source {params.source:q} --destination {params.dest:q} --kind nuclear"


rule research_locus_trees:
    input:
        staged=rules.research_stage_nuclear.output,
        script=f"{RSCRIPTS}/nuclear/run_loci.py",
        guard=helper("research_run_guard.py"),
    output: protected(f"{RN}/WORKFLOW_LOCUS_COMPLETE.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/nuclear/run_loci.py", root=RN,
        iqtree=RC.get("iqtree") or "iqtree2",
        native=f"{RN}/LOCUS_ALL_STATUS.json",
    threads: RC.get("locus_workers", 3)
    resources: mem_mb=6000
    log: f"{RR}/logs/locus_trees.log"
    shell:
        "{params.python:q} {params.script:q} --root {params.root:q} --iqtree {params.iqtree:q} --workers {threads} --mode all > {log:q} 2>&1 && "
        "{params.python:q} {input.guard:q} complete --native {params.native:q} --completion {output:q} --kind loci"


rule research_gene_concordance:
    input:
        completion=rules.research_locus_trees.output,
        script=f"{RSCRIPTS}/nuclear/gene_concordance.py",
        reference=RC["reference_tree"],
    output: protected(f"{RN}/concordance_independent/SUMMARY.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/nuclear/gene_concordance.py", root=RN,
        fresh=fresh(f"{RN}/concordance_independent"),
    log: f"{RR}/logs/gene_concordance.log"
    shell:
        "{params.fresh} {params.python:q} {params.script:q} --root {params.root:q} --reference-tree {input.reference:q} > {log:q} 2>&1"


rule research_spatial_fst:
    input:
        plan=RC["spatial_plan"],
        data=lambda wc: plan_files(RC["spatial_plan"], "input_prefix"),
        script=f"{RSCRIPTS}/population/spatial_fst.py",
        core=f"{RSCRIPTS}/population/influence_fst_core.py",
    output: protected(f"{RR}/spatial_fst/RUN_RECEIPT.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/population/spatial_fst.py", out=f"{RR}/spatial_fst",
        fresh=fresh(f"{RR}/spatial_fst"), chunk=RC.get("chunk_sites", 100000),
    log: f"{RR}/logs/spatial_fst.log"
    shell:
        "{params.fresh} {params.python:q} {params.script:q} --plan {input.plan:q} --out {params.out:q} --chunk-sites {params.chunk} > {log:q} 2>&1"


rule research_balanced_fst:
    input:
        plan=RC["balanced_plan"],
        data=lambda wc: plan_files(RC["balanced_plan"], "source_prefix"),
        script=f"{RSCRIPTS}/population/balanced_fst.py",
        core=f"{RSCRIPTS}/population/influence_fst_core.py",
    output: protected(f"{RR}/balanced_fst/RUN_RECEIPT.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/population/balanced_fst.py", out=f"{RR}/balanced_fst",
        fresh=fresh(f"{RR}/balanced_fst"), chunk=RC.get("chunk_sites", 100000),
    log: f"{RR}/logs/balanced_fst.log"
    shell:
        "{params.fresh} {params.python:q} {params.script:q} --plan {input.plan:q} --out {params.out:q} --chunk-sites {params.chunk} > {log:q} 2>&1"


rule research_stage_cp_ils:
    input:
        files=lambda wc: cp_ils_inputs(RC["chloroplast_ils_input_root"]),
        stage=helper("stage_research_inputs.py"),
        contract=helper("research_contracts.py"),
    output: protected(f"{RI}/STAGED.json")
    params:
        python=RPY, stage=helper("stage_research_inputs.py"),
        source=RC["chloroplast_ils_input_root"], dest=RI,
    shell:
        "{params.python:q} {params.stage:q} --source {params.source:q} --destination {params.dest:q} --kind chloroplast_ils"


rule research_cp_ils_simulate:
    input:
        staged=rules.research_stage_cp_ils.output,
        script=f"{RSCRIPTS}/chloroplast_ils/simulate.py",
    output: protected(f"{RI}/results_{PHASE}/RUN_RECEIPT.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/chloroplast_ils/simulate.py", root=RI, phase=PHASE,
        fresh=fresh(f"{RI}/results_{PHASE}"),
    log: f"{RR}/logs/chloroplast_ils_simulate.log"
    shell:
        "{params.fresh} {params.python:q} {params.script:q} --root {params.root:q} --phase {params.phase:q} > {log:q} 2>&1"


rule research_cp_ils_validate:
    input:
        receipt=rules.research_cp_ils_simulate.output,
        script=f"{RSCRIPTS}/chloroplast_ils/validate.py",
    output: protected(f"{RI}/results_{PHASE}/VALIDATION.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/chloroplast_ils/validate.py", root=RI, phase=PHASE,
    log: f"{RR}/logs/chloroplast_ils_validate.log"
    shell:
        "{params.python:q} {params.script:q} --root {params.root:q} --phase {params.phase:q} > {log:q} 2>&1"


rule research_cp_trees:
    input:
        files=lambda wc: cp_tree_inputs(RC["cp_inputs"]),
        script=f"{RSCRIPTS}/chloroplast/run_cp_trees.py",
        guard=helper("research_run_guard.py"),
        contract=helper("research_contracts.py"),
    output: protected(f"{RR}/chloroplast_trees/WORKFLOW_COMPLETE.json")
    params:
        python=RPY, script=f"{RSCRIPTS}/chloroplast/run_cp_trees.py", out=f"{RR}/chloroplast_trees",
        fresh=fresh(f"{RR}/chloroplast_trees"), inputs=RC["cp_inputs"],
        native=f"{RR}/chloroplast_trees/RUN_STATUS.json",
    threads: 4
    resources: mem_mb=6144
    log: f"{RR}/logs/chloroplast_trees.log"
    shell:
        "set -euo pipefail; {params.python:q} {input.guard:q} cp-resources --inputs {params.inputs:q} --allocated-threads {threads}; "
        "{params.fresh} {params.python:q} {params.script:q} --inputs {params.inputs:q} --out {params.out:q} --mode full --execute > {log:q} 2>&1 && "
        "{params.python:q} {input.guard:q} complete --native {params.native:q} --completion {output:q} --kind cp_trees"
