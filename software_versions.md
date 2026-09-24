# Software versions and dependency notes

This file is derived from the uploaded shell scripts (`01`–`20`). It records the tools that are explicitly invoked by the code and notes whether the version is pinned by the scripts.

## Important limitation

The scripts **do not record exact software versions** for most tools. That means this file is **code-matched**, but not a guaranteed reconstruction of the exact historical environment used for the original runs.

## Tools explicitly invoked by the scripts

| Tool | Where used | Version pinned in script? | Notes |
|---|---|---:|---|
| Guppy / `guppy_basecaller` | `01_guppy_basecalling.sh` | No | Script sets `GUPPY_CONFIG=dna_r10.4.1_e8.2_400bps_sup.cfg`, `MIN_QSCORE=8`, `DEVICE=auto`. Guppy itself is external/manual. |
| `gzip` | `01_guppy_basecalling.sh` | No | Used to compress Guppy outputs. |
| `find`, `cat` | `02_combine_pass_reads.sh` | No | Shell/coreutils only. |
| Flye | `03_flye_assembly.sh` | No | Uses `--nano-raw`, `--scaffold`, `--keep-haplotypes`, `--asm-coverage 50`. |
| Medaka / `medaka_consensus` | `04_medaka_polish.sh` | No | Optional model can be supplied as argument 5. |
| FastQC | `05_trim_illumina_fastqc.sh`, `11_rnaseq_fastqc_raw.sh`, `12_rnaseq_trim_fastqc.sh` | No | Uses `-t` thread argument. |
| Trimmomatic | `05_trim_illumina_fastqc.sh`, `12_rnaseq_trim_fastqc.sh` | No | Requires Java runtime. Adapter file is user-supplied. |
| minimap2 | `06_map_reads_for_hypo.sh`, `16_cp_align_to_reference.sh` | No | ONT mapping in script 06 and assembly/reference alignment in script 16. |
| samtools | `06_map_reads_for_hypo.sh`, `13_rnaseq_align_hisat2.sh`, `16_cp_align_to_reference.sh`, `17_cp_consensus_from_bam.sh` | No | Sorting, indexing, view, merge, depth. |
| BWA | `06_map_reads_for_hypo.sh` | No | `bwa index` and `bwa mem`. |
| Hypo | `07_hypo_polish.sh` | No | Uses both short-read and long-read BAM inputs. |
| BUSCO | `08_busco_assessment.sh` | No | Requires lineage dataset such as `embryophyta_odb10`. |
| QUAST / `quast.py` | `09_quast_assessment.sh` | No | Uses reference FASTA and GFF. |
| RepeatMasker | `10_te_repeatmasker.sh` | No | Called with `-species`, `-pa`, `-xsmall`, `-nolow`, `-norna`. |
| HISAT2 / `hisat2-build`, `hisat2` | `13_rnaseq_align_hisat2.sh` | No | Builds genome index and aligns trimmed RNA-seq reads. |
| StringTie | `14_stringtie_transcripts.sh` | No | Transcript assembly from merged BAM. |
| BRAKER / `braker.pl` | `15_braker_annotation.sh` | No | Runs with RNA BAM, protein evidence, and `--softmasking`. |
| `bcftools` | `17_cp_consensus_from_bam.sh` | No | `mpileup`, `call`, `consensus`. |
| `tabix` | `17_cp_consensus_from_bam.sh` | No | Indexes compressed VCF. Usually from htslib. |
| `awk` | `17_cp_consensus_from_bam.sh` | No | Generates low-coverage table. |
| NOVOPlasty | `18_novoplasty_chloroplast.sh`, `19_novoplasty_mitochondria.sh` | No | Invoked as `perl NOVOPlasty4.3.5.pl -c config.txt`; external/manual. |
| Perl | `18_novoplasty_chloroplast.sh`, `19_novoplasty_mitochondria.sh` | No | Required to execute NOVOPlasty. |
| GetOrganelle / `get_organelle_from_reads.py` | `20_getorganelle_mitochondria.sh` | No | Plant mitochondrial mode `-F embplant_mt`. |

## Implicit or likely secondary dependencies

These are not always called directly in the scripts, but are commonly required because of the software chosen:

- **Java** for Trimmomatic.
- **RMBlast** and **TRF** for RepeatMasker in many standard installations.
- **AUGUSTUS** and related BRAKER backend requirements for BRAKER runs.
- **Python** runtime for several tools distributed through conda.

## External/manual tools not safely managed inside a generic conda environment

### Guppy
The basecalling script uses:

- executable: `guppy_basecaller`
- config: `dna_r10.4.1_e8.2_400bps_sup.cfg`

Because Guppy is typically distributed by Oxford Nanopore with platform-specific installation constraints, it is listed as **external/manual** rather than embedded in `environment.yml`.

### NOVOPlasty
Scripts 18 and 19 call a local Perl script:

- `NOVOPlasty4.3.5.pl`

That indicates a manually downloaded program rather than a package reliably pinned by the uploaded scripts.

## Suggested reproducibility practice

For stronger reproducibility, record the following from the real runtime environment you used on the cluster:

```bash
conda list
samtools --version
minimap2 --version
flye --version
medaka --version
hisat2 --version
stringtie --version
busco --version
quast.py --version
bcftools --version
RepeatMasker -version
braker.pl --version
get_organelle_from_reads.py --version
```

Then replace the placeholder or unpinned entries in `environment.lock.template.yml` with the actual versions you used.
