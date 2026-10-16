# Dioscorea tenuipes · genome and comparative analysis

Code for nuclear genome assembly and annotation, population genomics, organelle validation, and subsequent nuclear–chloroplast comparative analyses.

**[Every script and its purpose](docs/SCRIPT_CATALOG.md)** · **[Software & package versions](software_versions.md)** · **[Snakemake workflows](workflow/README.md)** · **[Project index](PROJECT_INDEX.md)**

![Overview of assembly, population, organelle and comparative workflows](docs/figures/pipeline_overview.svg)

## Scripts by purpose

| Analysis | Scripts | What they do |
|---|---|---|
| Assembly and polishing | [01–09](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis) | Basecall reads, assemble contigs, polish the sequence, and measure assembly quality. |
| Annotation and scaffold review | [10–15, 23–29](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis) | Process RNA-seq, annotate genes/repeats, and audit reference-guided scaffolds and transcript evidence. |
| Population variants and structure | [21–22, 30–34](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis) | Call/audit variants, mask genotypes, construct panels, and estimate PCA, ADMIXTURE and FST. |
| Individual relationships | [44–49](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis) | Estimate relatedness, unrooted trees and networks; cross-check the implementations. |
| Sample and filter sensitivity | [50–55](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis), [research/population](docs/SCRIPT_CATALOG.md#later-population-sensitivity) | Check public-read processing, sample influence, allele balance, genomic regions, balanced sampling and genotype likelihoods. |
| Organelle sequence evidence | [16–20, 35–43](docs/SCRIPT_CATALOG.md#assembly-annotation-organelles-and-population-analysis) | Recover candidate sequences, annotate the plastome, and check sequence differences against reads. |
| Nuclear markers and phylogeny | [research/nuclear](docs/SCRIPT_CATALOG.md#nuclear-markers-and-gene-concordance) | Audit marker identity, build amino-acid matrices, estimate gene trees and count concordant/alternative relationships. |
| Chloroplast comparisons | [research/chloroplast](docs/SCRIPT_CATALOG.md#chloroplast-comparison-and-topology-controls) | Prepare expanded plastome panels and compare topology, source choice and alignment/masking controls. |
| Matched *D. panthaica* comparisons | [research/panthaica](docs/SCRIPT_CATALOG.md#d-panthaica-matched-marker-comparison) | Compare matched marker sets with and without the additional taxon and test locus/site sensitivity. |
| Nuclear lineage-sorting models | [research/nuclear_ils](docs/SCRIPT_CATALOG.md#nuclear-multispecies-coalescent-simulations) | Generate conditional nuclear gene histories and re-estimate their quartet relationships. |
| Chloroplast lineage-sorting models | [research/chloroplast_ils](docs/SCRIPT_CATALOG.md#conditional-chloroplast-genealogies-and-distance-scenarios) | Simulate linked plastid genealogies and sequence-distance scenarios, then check them analytically. |

The catalogue gives a separate purpose for **each script**, including helpers and historical reference code. [The later-analysis source map](workflow/research/SOURCE_MAP.csv) records source stages, hashes and portability changes.

## Population workflow

![Population filtering, structure, relationships and sensitivity workflow](docs/figures/pipeline_population.svg)

The root [Snakefile](Snakefile) connects population filtering/structure (32–34), relationships (44–49) and sample influence (53–55). Assembly/annotation (01–15) is an optional template. [The rule map](workflow/README.md) lists the precise coverage.

```bash
# Plan the existing population workflow; no jobs run with -n.
snakemake all -n --cores 4 --configfile config/config.yaml

# Plan later analyses from separately prepared inputs.
snakemake research_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml
```

Paths in the configuration files are examples. Later targets cover locus trees/concordance, regional and balanced-sample FST, and conditional chloroplast simulation/validation. The bounded chloroplast tree driver is a separate optional target. See [later workflow contracts](workflow/research/README.md) for commands and boundaries.

## Comparative analyses

![Nuclear, chloroplast, matched-marker and conditional genealogy workflows](docs/figures/pipeline_comparative.svg)

Nuclear and chloroplast relationships are evaluated separately. Conditional simulations describe outcomes under stated model assumptions; they do not give a calibrated probability of hybridization or prove a single evolutionary mechanism.

## Versions and verification

[software_versions.md](software_versions.md) records the software and package versions actually found in run logs and environment receipts, including stage-specific differences. [The evidence table](docs/provenance/software_versions.tsv) supplies version excerpts and source hashes. Environment templates and this release's code-check environment are identified separately.

The [2026-09-24 verification record](docs/WORKFLOW_VERIFICATION.md) describes the earlier population replay. [The current release checks](docs/PUBLIC_CODE_RELEASE_2026-10-16.md) cover code, small artificial tests, documentation and dry-run planning. No research analysis was rerun for this update.

A few older scripts are retained for inspection: **22** relabels chromosomes without a full coordinate transformation; **24** reconstructs incomplete hints and is unsuitable for annotation; **27** is a repeat-analysis template whose resource settings need review. Current contig-coordinate population analyses use **32–55**; transcript-structure evidence uses **29**. Reference-guided scaffolds remain candidate structures. Mitochondrial assembly attempts do not constitute a complete mitochondrial reference.

## Repository layout

```text
scripts/                 Numbered scripts and shared helpers
scripts/research/        Later scientific scripts, grouped by purpose
Snakefile                Snakemake entry point
config/                  Population, assembly and later-analysis configuration
workflow/                Rules, environments, contracts and dry-run checks
docs/SCRIPT_CATALOG.md   Purpose of every script
docs/figures/            Three SVG diagrams and their editable generator
software_versions.md    Stage-specific software and package versions
PROJECT_INDEX.md         Current release state and history
```

Large sequencing data and manuscript drafts are distributed separately. Code is released under the [MIT License](LICENSE).
