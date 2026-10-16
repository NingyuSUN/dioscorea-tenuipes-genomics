# Software versions

Versions below come from archived run logs, runtime receipts, environment exports, and repository configuration. Different analysis stages used different environments. A recorded version is specific to its stage; it is not a claim that every script ran with that version.

[The version evidence index](docs/provenance/software_versions.tsv) records each tool, version, stage, evidence type, source SHA-256, and a short exact version excerpt. Paths beginning `docs/evidence/` and archive-member paths are logical identifiers for source files retained in the private archive, not links to files distributed here. For an `archive.tar.gz!member` identifier, the SHA-256 covers the uncompressed member bytes. Other hashes cover the named file bytes. Only version excerpts are included publicly.

## Recorded analysis environments

| Analysis stage | Recorded software | Evidence and scope |
|---|---|---|
| Initial plastid phylogeny | MAFFT **7.525**; IQ-TREE **2.0.7** | Actual alignment and inference log headers, 2026-09-15. |
| Expanded marker phylogeny | IQ-TREE **2.4.0** | P16 inference receipt records completed executions. |
| Whole-plastome comparison | MAFFT **7.526**; IQ-TREE **2.4.0**; BLAST+ `blastn` **2.17.0+**; Python **3.13.15**; Biopython **1.88** | P19 tool-version capture. These versions do not replace the earlier phylogeny versions. |
| Nuclear marker binding | miniprot **0.18-r281** | Mapping log version header. |
| Protein completeness assessment | BUSCO **5.8.3**; HMMER `hmmsearch` **3.4**; Python **3.13.15** | BUSCO result summary. The lineage was `embryophyta_odb10`, dated 2024-01-08, with 1,614 markers; this is a dataset version, not a software version. |
| Original gene annotation | BRAKER **3.0.8** | Original 2025 annotation log retained in the 2026 inventory. |
| Genome-to-organelle alignment | minimap2 **2.31-r1302** | Alignment log version header. |
| Plastome annotation / indel checks | Python **3.10.16**; Biopython **1.79** | Execution-environment receipt linked to the annotation scripts. |
| Public-read validation | Python **3.12.14**; NumPy **2.5.3**; BWA **0.7.17-r1188**; samtools **1.13**; HTSlib **1.13+ds**; fastp **0.24.0** | Archived version-command outputs. BWA prints its version through its help output. |
| Population filtering / production environment | bcftools **1.13**; HTSlib **1.13+ds**; PLINK **1.90b6.21**; ADMIXTURE **1.3.0**; VCFtools **0.1.16** | bcftools/HTSlib version capture; PLINK/ADMIXTURE/VCFtools exact package export. An installed-package export alone does not prove execution of every package. |
| Population relationships / influence plots | Python **3.12.14**; NumPy **2.5.3**; SciPy **1.18.1**; pandas **3.0.5**; matplotlib **3.11.2**; Biopython **1.88**; splitspy **0.0.10** | Stage version receipts. splitspy supplies circular ordering only; its native weight solver is not used. |
| Network / tree cross-checks | R **4.5.3**; ape **5.8-1**; phangorn **2.12.1** | R session information. |
| Population relatedness | R **4.5.3**; SNPRelate **1.44.0**; GENESIS **2.40.0**; GWASTools **1.56.0**; BiocParallel **1.44.0**; gdsfmt **1.46.0** | R session information. The evidence index also records loaded R dependency namespaces; loading a namespace does not establish direct analytical use. |
| PCAngsd build / runtime audit (P4) | Python **3.12.14**; PCAngsd **1.36.4**; NumPy **2.5.3**; SciPy **1.18.1**; Cython **3.2.4** | Historical runtime manifest. This row documents the audited tool environment, not completion of a production genotype-likelihood analysis. |
| Nuclear multispecies-coalescent (MSC) panel generation | msprime **1.4.4** | Completed generation state in the ILS evidence archive. The separate environment check records Python **3.11.16**, NumPy **2.4.6**, and tskit **1.0.3**, but is explicitly an import/signature check. |
| Chloroplast incomplete-lineage-sorting (ILS) simulations | Python **3.11.16**; NumPy **2.4.6**; msprime **1.4.4** | Full simulation run receipt. Its status records simulation completion with validation pending at that checkpoint; the version record is not a scientific validation claim. |

## Configured environments

| File | Meaning |
|---|---|
| `workflow/envs/relationships_py.yaml` | Analysis-package versions match the archived relationships receipts. Python is constrained to 3.12, not an exact patch version. |
| `workflow/envs/relationships_r.yaml` | Main R package versions match archived session information. R is constrained to 4.5; the archived runtime was 4.5.3. |
| `workflow/envs/popgen.yaml` | PLINK, ADMIXTURE, and VCFtools are pinned. `bcftools>=1.13` is a lower bound, not an exact historical version. |
| `workflow/envs/research-requirements.txt` | Exact Python package versions used for the current research-code engineering checks; not a historical analysis lock. |
| `workflow/envs/ci.yaml` | Code-check environment constraints. They do not describe all historical analyses. |
| `environment.yml` | General assembly/annotation dependency specification, mostly without exact version pins. |
| `environment.lock.template.yml` | Best-effort reconstruction template, not a verified historical lockfile. For example, its minimap2 2.28, samtools 1.21, and BWA 0.7.18 differ from recorded later-run versions above. |

The evidence index preserves configured version constraints separately from observed versions. These files do not establish one environment containing every tool listed here.

## Versions not established by the indexed records

Exact historical runtime versions were not established for Guppy, Flye, Medaka, FastQC, Trimmomatic, Hypo, QUAST, RepeatMasker, HISAT2, StringTie, AUGUSTUS, GetOrganelle, GATK, trimAl, Seq-Gen, or ete3. Configuration-template versions must not be substituted for missing runtime evidence. Scripts 18 and 19 name `NOVOPlasty4.3.5.pl`; **4.3.5 is a filename indication**, not a captured runtime version.

The scikit-allel **1.3.3** reference identifies documentation inspected for the Hudson FST implementation. The archived methods explicitly distinguish that reference from an installed scikit-allel run; it is therefore excluded from the installed analysis environments above. Likewise, a basecalling model/configuration filename is not the Guppy executable version.

## Current code-check environment

The 2026-10-16 public-code engineering environment recorded Python **3.12.14**, Snakemake **9.27.0**, NumPy **2.5.4**, SciPy **1.18.1**, pandas **3.0.6**, matplotlib **3.11.2**, Biopython **1.88**, pysam **0.24.1**, splitspy **0.0.10**, PyYAML **6.0.3**, msprime **1.4.4**, and tskit **1.0.3**. Temporary bcftools **1.22** and HTSlib **1.22.1** were used for the existing code tests. The [sanitized environment receipt](docs/provenance/current_validation_environment.json) records these versions. They describe code checks in this update, not historical research runs or a new scientific analysis.

## Historical-script portability update

The 22 configured historical stages use the recorded versions above; no new
scientific runtime was measured. [Their dependency table](workflow/portable/DEPENDENCIES.md)
now identifies Python packages, external tools, lineage/data bundles and Linux
resource controls by stage. The same recorded engineering environment is used
for artificial-input checks. PCAngsd build hashes remain explicit; installing a
package with the same version does not establish an identical compiled build.
