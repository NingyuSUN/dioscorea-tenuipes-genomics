# Dependencies of the configured historical scripts

Install the Python packages in [research-requirements.txt](../envs/research-requirements.txt)
for the shared engineering environment. External executable paths are configured
per stage; they are not assumed to live in a particular Conda environment or
server directory. The [recorded version table](../../software_versions.md) is the
source for historical version claims. A package specification is not proof that
a complete historical environment has been reproduced.

| Stages | Python packages | Executables, data dependencies and platform requirements |
|---|---|---|
| Population full panel and comparison | NumPy, SciPy; PCAngsd for the full-panel fit | PLINK for the hard-call arm; the explicit three-genotype reference is bundled. The PCAngsd Python and compiled modules and the PLINK binary are hash-bound. |
| Nuclear recovery and concatenated tree | Biopython | MEGAHIT, BUSCO with its own runtime and lineage dataset, MAFFT, IQ-TREE; the configured BUSCO interpreter must contain its BUSCO dependencies. |
| Nuclear locus preparation / CDS translation | Biopython | Prepared matrices, partition summaries, reference tree, genomic sequence and marker evidence; no tree inference during preparation. |
| Chloroplast matrix/source controls | Biopython where tree parsing is used | IQ-TREE binary identity is recorded by preparation; later execution also needs the exact prepared models and matrices. |
| Read coverage / de novo plastome validation | pysam | minimap2, samtools and paired reads; mapping tools are explicit configuration inputs. |
| CP98/97 tree drivers | Biopython in bundled helpers | IQ-TREE; Linux cgroup v2, systemd user scopes, taskset, timeout and configured scope helper; frozen admission/evidence receipts. Memory, swap and affinity checks are retained. |
| Panthaica RNA / marker / alignment stages | Biopython in the configured BUSCO runtime | Cutadapt, rnaSPAdes/SPAdes runtime, BUSCO, HMMER and MAFFT as required by the stage; original read identity, lineage and matched-marker evidence remain required. |
| Nuclear conditional history fitting / audit | NumPy, SciPy for fitting; tskit for the panel audit | Frozen panel manifests and inference evidence; Linux cgroup v2 and a matching scope receipt for bounded execution. The audit consumes existing inference outputs. |
| Chloroplast conditional preparation / figures | NumPy, SciPy, Matplotlib, Biopython; msprime/tskit for the separate simulator | Hash-bound guide/evidence files and prepared simulation summaries. Plotting does not rerun simulation. |

PCAngsd **1.36.4** is recorded in the historical runtime audit. Its compiled
module hashes depend on the actual build, so the configuration requires all four
reviewed module digests together when replacing the original build. Never copy
placeholder hashes into a production plan or disable these identity checks.

The optional workflow engine is Snakemake (minimum 8; this release's checks use
9.27.0). Tests use fake executables and artificial inputs where tools are absent.
The release therefore does not claim that every external tool is installed or
that real-data execution has been validated on a new machine.
