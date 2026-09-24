# Per-script dependency table

This table is derived directly from the uploaded shell scripts.

| Script | Main purpose | Tools explicitly invoked | Key default outputs |
|---|---|---|---|
| `01_guppy_basecalling.sh` | ONT basecalling | `guppy_basecaller`, `find`, `gzip` | `guppy_out/`, gzipped pass/fail FASTQ |
| `02_combine_pass_reads.sh` | Combine passed ONT reads | `find`, `cat` | `combined_pass.fastq.gz` |
| `03_flye_assembly.sh` | Long-read assembly | `flye` | `flye_out/` |
| `04_medaka_polish.sh` | Medaka polishing | `medaka_consensus` | `medaka_out/` |
| `05_trim_illumina_fastqc.sh` | Illumina QC and trimming | `fastqc`, `trimmomatic` | `illumina_qc/raw_fastqc/`, `illumina_qc/trimmed/`, `illumina_qc/trimmed_fastqc/` |
| `06_map_reads_for_hypo.sh` | Map ONT + Illumina reads to polished assembly | `minimap2`, `samtools`, `bwa` | `mapping/mapped_lg.sorted.bam`, `mapping/mapped_sr.sorted.bam`, `mapping/il_names.txt` |
| `07_hypo_polish.sh` | Hybrid polishing with Hypo | `hypo` | `final_assembly.fa` |
| `08_busco_assessment.sh` | Assembly completeness assessment | `busco` | `busco_out/` |
| `09_quast_assessment.sh` | Assembly quality comparison | `quast.py` | `quast_out/` |
| `10_te_repeatmasker.sh` | TE masking | `RepeatMasker` | `repeatmasker_out/` and masked genome |
| `11_rnaseq_fastqc_raw.sh` | Raw RNA-seq QC | `fastqc` | `RNA_Seq/fastqc_raw/` |
| `12_rnaseq_trim_fastqc.sh` | RNA-seq trimming + QC | `trimmomatic`, `fastqc`, `tee`, `date` | `RNA_Seq/trimmomatic/`, `RNA_Seq/fastqc_trimmomatic/` |
| `13_rnaseq_align_hisat2.sh` | RNA-seq genome alignment | `hisat2-build`, `hisat2`, `samtools` | `annotation_output/hisat2_index/`, `annotation_output/hisat2_bam/`, `annotation_output/rna_aligned.sorted.bam` |
| `14_stringtie_transcripts.sh` | Transcript assembly | `stringtie` | `annotation_output/transcripts.gtf` |
| `15_braker_annotation.sh` | Gene annotation | `braker.pl`, `cat`, shell built-ins | `braker_work/` |
| `16_cp_align_to_reference.sh` | Align contigs to cp reference | `minimap2`, `samtools` | `<prefix>.sorted.bam`, `<prefix>.sorted.bam.bai` |
| `17_cp_consensus_from_bam.sh` | Reference-guided cp consensus | `bcftools`, `samtools`, `tabix`, `awk` | `<prefix>.variants.vcf.gz`, `<prefix>.consensus.fasta`, coverage reports |
| `18_novoplasty_chloroplast.sh` | cp de novo assembly | `perl`, `NOVOPlasty4.3.5.pl` | NOVOPlasty output from config |
| `19_novoplasty_mitochondria.sh` | mt de novo assembly | `perl`, `NOVOPlasty4.3.5.pl` | NOVOPlasty output from config |
| `20_getorganelle_mitochondria.sh` | mt assembly with GetOrganelle | `get_organelle_from_reads.py` | `mito_getorganelle/` |

## Minimal tool set by workflow block

### Nuclear assembly and polishing
- Guppy
- Flye
- Medaka
- minimap2
- samtools
- bwa
- Hypo

### Assembly evaluation
- BUSCO
- QUAST

### Repeat and gene annotation
- RepeatMasker
- HISAT2
- StringTie
- BRAKER

### Organelle analyses
- minimap2
- samtools
- bcftools
- tabix
- NOVOPlasty
- GetOrganelle

### General runtime utilities
- Perl
- Java (for Trimmomatic)
- common POSIX tools (`find`, `cat`, `awk`, `tee`, `date`, `gzip`)
