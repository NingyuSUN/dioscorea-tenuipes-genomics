# Script purposes and dependencies

The complete [script catalogue](docs/SCRIPT_CATALOG.md) describes every numbered script, later research entry point, and supporting module. [Software versions](software_versions.md) separates observed versions from configured constraints.

| Scripts / module | Main software and packages |
|---|---|
| 01–09: assembly and polishing | Guppy, Flye, Medaka, FastQC, Trimmomatic/Java, minimap2, BWA, samtools, Hypo, BUSCO, QUAST |
| 10–15, 23–29: annotation and scaffolding | RepeatMasker, RepeatModeler, HISAT2, StringTie, BRAKER, RagTag; Python |
| 16–20, 35–43: organelles | minimap2, samtools, bcftools, NOVOPlasty/Perl, GetOrganelle; Python/Biopython |
| 21–22: historical population processing | BWA, GATK/Java, samtools, minimap2, PLINK, ADMIXTURE |
| 30–34: population auditing and structure | Python, pysam, NumPy, pandas, matplotlib; HTSlib/C compiler, bcftools, PLINK, ADMIXTURE, VCFtools |
| 44–49: relationships | NumPy, SciPy, pandas, matplotlib, Biopython, splitspy; R/SNPRelate/GENESIS/GWASTools/BiocParallel/ape/phangorn |
| 48, 50–52: public-data probes | Supplied Entrez adapter for 48, ENA/NCBI endpoints; fastp, BWA, samtools; Python |
| 53–55 and research/population | NumPy, SciPy, pandas, matplotlib; PLINK, VCFtools and bcftools where invoked; explicit genotype-likelihood scripts also use their stage configurations |
| research/nuclear | Python, Biopython, NumPy; BUSCO/HMMER, miniprot, MAFFT, IQ-TREE |
| research/chloroplast | Python, Biopython, NumPy/SciPy where imported; MAFFT, IQ-TREE; recovery drivers additionally use samtools/bcftools |
| research/panthaica | Python, Biopython, NumPy; MAFFT, IQ-TREE; historical recovery uses stage-specific RNA-seq/mapping tools |
| research/nuclear_ils | Python, NumPy, msprime/tskit, Biopython; IQ-TREE for re-estimation |
| research/chloroplast_ils | Python, NumPy, msprime/tskit; SciPy and Biopython for independent checks |
| Snakemake wrappers | Snakemake >=8, PyYAML, a Unix shell; scientific dependencies installed separately |

The later source files preserve assertions about input hashes, sample ordering and marker definitions. Those assertions are part of the historical methods, not general-purpose defaults to remove when a new dataset fails them. Resource-bounded drivers require Linux cgroup v2 in addition to the software packages.
