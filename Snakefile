# Snakemake entry point for the Dioscorea tenuipes repository.
#
# Wraps the existing numbered scripts in scripts/ without changing them.
# Usage, safety rules and the mapping to historical runs: workflow/README.md.
#
#   snakemake -n                                      # dry run, population modules
#   snakemake --profile workflow/profiles/default     # run
#   snakemake --configfile config/config.yaml config/assembly.yaml -n assembly_all

from snakemake.utils import min_version

min_version("8.0")


configfile: "config/config.yaml"


include: "workflow/rules/common.smk"
include: "workflow/rules/population.smk"
include: "workflow/rules/relationships.smk"
include: "workflow/rules/influence.smk"


if "assembly" in config:

    include: "workflow/rules/assembly.smk"


if "research" in config:
    include: "workflow/rules/research.smk"


MODULE_TARGETS = {
    "population": POPULATION_TARGETS,
    "relationships": RELATIONSHIP_TARGETS,
    "influence": INFLUENCE_TARGETS,
}


rule all:
    input:
        [t for m in config["targets"] for t in MODULE_TARGETS[m]],
    default_target: True


rule population_all:
    input:
        POPULATION_TARGETS,


rule relationships_all:
    input:
        RELATIONSHIP_TARGETS,


rule influence_all:
    input:
        INFLUENCE_TARGETS,
