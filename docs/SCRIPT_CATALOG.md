# Script catalogue

Every analysis script below has a specific purpose. `historical_reference` means the script retains a stage-specific directory/configuration contract and is kept for method review. `entrypoint` means a command-line entry point is present; it does not imply a newly reproduced scientific result.

## Assembly, annotation, organelles and population analysis

| Script | Purpose |
|---|---|
| [01_guppy_basecalling.sh](../scripts/01_guppy_basecalling.sh) | Basecall Nanopore signal files and separate pass/fail reads. |
| [02_combine_pass_reads.sh](../scripts/02_combine_pass_reads.sh) | Combine passed Nanopore reads into one compressed FASTQ. |
| [03_flye_assembly.sh](../scripts/03_flye_assembly.sh) | Assemble nuclear contigs from Nanopore long reads with Flye. |
| [04_medaka_polish.sh](../scripts/04_medaka_polish.sh) | Polish the draft assembly using Nanopore reads with Medaka. |
| [05_trim_illumina_fastqc.sh](../scripts/05_trim_illumina_fastqc.sh) | Assess and trim paired Illumina reads before hybrid polishing. |
| [06_map_reads_for_hypo.sh](../scripts/06_map_reads_for_hypo.sh) | Map long and short reads to the draft for hybrid polishing. |
| [07_hypo_polish.sh](../scripts/07_hypo_polish.sh) | Polish the assembly with both long- and short-read alignments. |
| [08_busco_assessment.sh](../scripts/08_busco_assessment.sh) | Measure assembly completeness against a BUSCO lineage dataset. |
| [09_quast_assessment.sh](../scripts/09_quast_assessment.sh) | Compare assembly continuity and reference-alignment statistics with QUAST. |
| [10_te_repeatmasker.sh](../scripts/10_te_repeatmasker.sh) | Soft-mask repeats with a selected RepeatMasker species library. |
| [11_rnaseq_fastqc_raw.sh](../scripts/11_rnaseq_fastqc_raw.sh) | Assess raw RNA-seq read quality. |
| [12_rnaseq_trim_fastqc.sh](../scripts/12_rnaseq_trim_fastqc.sh) | Trim RNA-seq reads and assess post-trimming quality. |
| [13_rnaseq_align_hisat2.sh](../scripts/13_rnaseq_align_hisat2.sh) | Align RNA-seq libraries to the genome with HISAT2 and merge alignments. |
| [14_stringtie_transcripts.sh](../scripts/14_stringtie_transcripts.sh) | Assemble transcript structures from RNA-seq alignments with StringTie. |
| [15_braker_annotation.sh](../scripts/15_braker_annotation.sh) | Run BRAKER gene annotation using RNA-seq and protein evidence. |
| [16_cp_align_to_reference.sh](../scripts/16_cp_align_to_reference.sh) | Align assembly sequences to a chloroplast reference. |
| [17_cp_consensus_from_bam.sh](../scripts/17_cp_consensus_from_bam.sh) | Call a reference-guided chloroplast consensus and report depth. |
| [18_novoplasty_chloroplast.sh](../scripts/18_novoplasty_chloroplast.sh) | Attempt de novo chloroplast assembly with NOVOPlasty. |
| [19_novoplasty_mitochondria.sh](../scripts/19_novoplasty_mitochondria.sh) | Attempt de novo mitochondrial assembly with NOVOPlasty; not a complete mitochondrial release. |
| [20_getorganelle_mitochondria.sh](../scripts/20_getorganelle_mitochondria.sh) | Attempt plant mitochondrial assembly with GetOrganelle. |
| [21_gatk_population_pipeline.sh](../scripts/21_gatk_population_pipeline.sh) | Map population reads, call per-sample variants, and perform joint GATK genotyping. |
| [22_population_pca_admixture.sh](../scripts/22_population_pca_admixture.sh) | Historical PCA/ADMIXTURE route with chromosome relabeling; labels are not transformed coordinates. Superseded by 32–34. |
| [23_ragtag_scaffolding.sh](../scripts/23_ragtag_scaffolding.sh) | Order and orient contigs against a reference with RagTag; produces candidate scaffolds. |
| [24_repair_stringtie_hints.py](../scripts/24_repair_stringtie_hints.py) | Historical reconstruction from incomplete transcript hints; unsuitable for annotation. Use 29. |
| [25_lift_hints_to_scaffolds.sh](../scripts/25_lift_hints_to_scaffolds.sh) | Lift transcript hints through scaffold coordinates for comparison. |
| [26_chromosome_correspondence.sh](../scripts/26_chromosome_correspondence.sh) | Compare contig-to-reference chromosome assignments. |
| [27_repeatmodeler_softmask.sh](../scripts/27_repeatmodeler_softmask.sh) | Template for building a repeat library and soft-masking; review CPU/batch settings before use. |
| [28_audit_annotation_and_scaffolding.py](../scripts/28_audit_annotation_and_scaffolding.py) | Audit retained annotation, BUSCO, scaffold assignments, and lifted hints. |
| [29_prepare_rna_evidence.py](../scripts/29_prepare_rna_evidence.py) | Derive traceable intron evidence from complete StringTie transcript structures. |
| [30_audit_population_vcf.py](../scripts/30_audit_population_vcf.py) | Compare exact VCF alleles/genotypes and audit sample-level missingness and depth. |
| [31_diagnose_population_samples.py](../scripts/31_diagnose_population_samples.py) | Diagnose population samples using deterministic genomic windows. |
| [32_filter_population_variants.py](../scripts/32_filter_population_variants.py) | Mask low-quality genotypes and construct explicitly defined population panels. |
| [33_run_population_structure.py](../scripts/33_run_population_structure.py) | Run PCA, repeated ADMIXTURE, and pairwise Weir–Cockerham FST. |
| [34_summarize_population_results.py](../scripts/34_summarize_population_results.py) | Summarize population filtering, structure, and differentiation in tables and figures. |
| [35_recover_organelle_sequences.py](../scripts/35_recover_organelle_sequences.py) | Recover graph-supported organelle sequences without filling gaps from a reference. |
| [36_validate_organelle_evidence.sh](../scripts/36_validate_organelle_evidence.sh) | Map sequencing evidence to recovered organelle candidates. |
| [37_summarize_organelle_evidence.py](../scripts/37_summarize_organelle_evidence.py) | Summarize plastome read support and unjoined mitochondrial homology candidates. |
| [38_compare_published_plastome.py](../scripts/38_compare_published_plastome.py) | Retrieve a versioned published plastome and prepare sequence comparison. |
| [39_summarize_plastome_comparison.py](../scripts/39_summarize_plastome_comparison.py) | Validate alignment columns and summarize plastome sequence differences. |
| [40_annotate_plastome.py](../scripts/40_annotate_plastome.py) | Project published plastome annotations through a verified alignment. |
| [41_validate_plastome_indels.py](../scripts/41_validate_plastome_indels.py) | Check plastome indels using diagnostic read-junction k-mers. |
| [42_refine_plastome_annotation.py](../scripts/42_refine_plastome_annotation.py) | Refine annotation boundaries and repeat-indel evidence. |
| [43_verify_plastome_annotation.py](../scripts/43_verify_plastome_annotation.py) | Cross-check delivered FASTA, GenBank, GFF, CDS, and proteins. |
| [44_prepare_relationship_panels.py](../scripts/44_prepare_relationship_panels.py) | Build strict, complete-call, and leave-one-out relationship panels. |
| [45_analyze_relationships.py](../scripts/45_analyze_relationships.py) | Estimate pairwise distances, unrooted NJ trees, and descriptive networks. |
| [46_relatedness.R](../scripts/46_relatedness.R) | Estimate relatedness with KING, PC-AiR, and PC-Relate. |
| [47_validate_and_plot.R](../scripts/47_validate_and_plot.R) | Independently check and plot trees/networks in R. |
| [47b_network_fixed_order.R](../scripts/47b_network_fixed_order.R) | Fit network split weights in R using a fixed circular ordering. |
| [48_discover_public_dt.py](../scripts/48_discover_public_dt.py) | Discover public accessions using a supplied Entrez API adapter and ENA. |
| [49_summarize_validate.py](../scripts/49_summarize_validate.py) | Validate relationship outputs and assemble summary figures/tables. |
| [50_probe_public_wgs.py](../scripts/50_probe_public_wgs.py) | Run a bounded public-read QC and mapping pilot. |
| [51_summarize_wgs_probe.py](../scripts/51_summarize_wgs_probe.py) | Summarize public-read pilot QC and mapping with explicit denominators. |
| [52_probe_trim_sensitivity.py](../scripts/52_probe_trim_sensitivity.py) | Compare raw and trimmed alignments on matched read identifiers. |
| [53_run_influence_panels.py](../scripts/53_run_influence_panels.py) | Construct leave-one-out and candidate-relative PCA sensitivity panels. |
| [54_summarize_influence.py](../scripts/54_summarize_influence.py) | Summarize how sample exclusions affect PCA and FST. |
| [55_compute_influence_fst.py](../scripts/55_compute_influence_fst.py) | Compute chunked Weir–Cockerham FST sensitivity from masked PLINK genotypes. |

## Shared libraries and workflow helpers

| Script | Purpose |
|---|---|
| [influence_fst_core.py](../scripts/influence_fst_core.py) | Compute diploid biallelic Weir–Cockerham FST components. |
| [relationship_core.py](../scripts/relationship_core.py) | Read genotype panels, compute IBS distances, and check unrooted tree/network summaries. |
| [wgs_probe_core.py](../scripts/wgs_probe_core.py) | Bound downloads and validate FASTQ structure for public-read pilots. |
| [compare_population_genotypes.c](../scripts/lib/compare_population_genotypes.c) | Stream exact-allele genotype comparisons through HTSlib. |
| [network_crosscheck.py](../workflow/scripts/network_crosscheck.py) | Compare network split weights under shared and independent circular orderings. |
| [stage_research_inputs.py](../workflow/scripts/stage_research_inputs.py) | Copy declared prepared inputs into a separate research output root without running analysis. |
| [build_pipeline_figures.py](../docs/figures/build_pipeline_figures.py) | Regenerate the three pipeline SVG diagrams deterministically. |


| [plot_organelle_summary.py](../scripts/lib/plot_organelle_summary.py) | Plot plastome structure, differences and retained-read support from comparison tables. |
| [check_plastome_differences.py](../scripts/lib/check_plastome_differences.py) | Count quality-filtered read bases at plastome alignment differences. |
| [research_contracts.py](../workflow/scripts/research_contracts.py) | Enumerate contained, explicitly declared research input files and their binary dependencies. |
| [research_run_guard.py](../workflow/scripts/research_run_guard.py) | Check allocated resources and native driver status before recording workflow completion. |

## Later population sensitivity

| Script | Purpose | Role |
|---|---|---|
| [balanced_fst.py](../scripts/research/population/balanced_fst.py) | Compare Weir–Cockerham and Hudson FST across all balanced IZU sample subsets. | entrypoint |
| [explicit3.py](../scripts/research/population/explicit3.py) | Compute the explicit three-genotype likelihood reference model for covariance analysis. | library |
| [extract_original_panel.py](../scripts/research/population/extract_original_panel.py) | Extract original-panel genotype likelihoods while preserving sample and site identity. | entrypoint |
| [historical/compare_full_panel.py](../scripts/research/population/historical/compare_full_panel.py) | Summarize covariance and subspace differences across fixed-panel model arms. | historical_reference |
| [historical/run_full_panel.py](../scripts/research/population/historical/run_full_panel.py) | Fit and compare fixed-dimensional genotype-likelihood covariance models on the full SNP panel. | historical_reference |
| [historical/run_plink_arm.py](../scripts/research/population/historical/run_plink_arm.py) | Compute the hard-call PLINK comparison arm for the same population panel. | historical_reference |
| [influence_fst_core.py](../scripts/research/population/influence_fst_core.py) | Compute population allele-count masks and Weir–Cockerham FST components. | library |
| [panel_identity.py](../scripts/research/population/panel_identity.py) | Verify sample and variant identity when binding the original SNP panel to genotype likelihoods. | library |
| [pl_diagnostic.py](../scripts/research/population/pl_diagnostic.py) | Audit three-genotype likelihoods and their alignment to the selected SNP panel. | entrypoint |
| [prepare_ab_beds.py](../scripts/research/population/prepare_ab_beds.py) | Prepare genotype panels under prespecified heterozygote allele-balance filters. | entrypoint |
| [run_ab_sensitivity.py](../scripts/research/population/run_ab_sensitivity.py) | Compare PCA and FST across allele-balance filters on matched population panels. | entrypoint |
| [spatial_fst.py](../scripts/research/population/spatial_fst.py) | Quantify regional and leave-region-out sensitivity of pairwise population FST. | entrypoint |
| [test_reference.py](../scripts/research/population/test_reference.py) | Test the explicit three-genotype likelihood reference on synthetic cases. | test |

## Nuclear markers and gene concordance

| Script | Purpose | Role |
|---|---|---|
| [audit_iqtree_gcf.py](../scripts/research/nuclear/audit_iqtree_gcf.py) | Cross-check independent branch concordance counts against IQ-TREE gCF output. | entrypoint |
| [audit_miniprot_stops.py](../scripts/research/nuclear/audit_miniprot_stops.py) | Check miniprot stop-codon flags against the nuclear marker translations. | library |
| [audit_paired_fastq.py](../scripts/research/nuclear/audit_paired_fastq.py) | Validate paired FASTQ structure, sequence alphabets and expected read counts. | entrypoint |
| [audit_quartet_projection.py](../scripts/research/nuclear/audit_quartet_projection.py) | Independently verify focal-quartet projection inputs and bootstrap counts. | entrypoint |
| [build_busco_matrices.py](../scripts/research/nuclear/build_busco_matrices.py) | Build single-copy BUSCO amino-acid matrices with fixed occupancy and marker exclusions. | entrypoint |
| [gene_concordance.py](../scripts/research/nuclear/gene_concordance.py) | Count concordant, alternative and unresolved gene-tree splits for each reference branch. | entrypoint |
| [historical/genomic_cds_translation.py](../scripts/research/nuclear/historical/genomic_cds_translation.py) | Verify genomic CDS translation and marker binding for the hime2 reference. | historical_reference |
| [historical/prepare_loci.py](../scripts/research/nuclear/historical/prepare_loci.py) | Split the frozen concatenated matrix into loci while omitting all-missing taxon rows. | historical_reference |
| [historical/recover_nuclear.py](../scripts/research/nuclear/historical/recover_nuclear.py) | Recover BUSCO nuclear markers from the selected genome assembly. | historical_reference |
| [historical/run_nuclear_main.py](../scripts/research/nuclear/historical/run_nuclear_main.py) | Coordinate the historical expanded nuclear panel and concatenated IQ-TREE analysis. | historical_reference |
| [marker_structure.py](../scripts/research/nuclear/marker_structure.py) | Compare candidate historical marker segments with current gene-model structure. | entrypoint |
| [project_locus_quartets.py](../scripts/research/nuclear/project_locus_quartets.py) | Project existing nuclear point and bootstrap trees onto the fixed focal quartet. | entrypoint |
| [run_loci.py](../scripts/research/nuclear/run_loci.py) | Infer one supported IQ-TREE phylogeny per retained nuclear locus. | entrypoint |
| [test_gene_concordance.py](../scripts/research/nuclear/test_gene_concordance.py) | Test split canonicalization and decisive-locus concordance on synthetic trees. | test |
| [warning_sensitivity.py](../scripts/research/nuclear/warning_sensitivity.py) | Recalculate branch concordance after applying recorded gene-tree warning exclusions. | entrypoint |

## Chloroplast comparison and topology controls

| Script | Purpose | Role |
|---|---|---|
| [historical/audit_hachi_nuclear_coverage.py](../scripts/research/chloroplast/historical/audit_hachi_nuclear_coverage.py) | Measure nuclear-marker read coverage for the Hachi sample. | historical_reference |
| [historical/compare_published_edo.py](../scripts/research/chloroplast/historical/compare_published_edo.py) | Compare the published Edo plastome with read-supported chloroplast states. | historical_reference |
| [historical/prepare_inputs.py](../scripts/research/chloroplast/historical/prepare_inputs.py) | Prepare the historical paired chloroplast matrix and partition-model controls. | historical_reference |
| [historical/prepare_source_inputs.py](../scripts/research/chloroplast/historical/prepare_source_inputs.py) | Prepare fixed-column chloroplast source variants for paired AU comparisons. | historical_reference |
| [historical/published_edo_controls.py](../scripts/research/chloroplast/historical/published_edo_controls.py) | Evaluate source and callable-site controls for the published Edo plastome. | historical_reference |
| [historical/recover_hachi_cp_callable.py](../scripts/research/chloroplast/historical/recover_hachi_cp_callable.py) | Recover callable chloroplast bases from competitive-reference read alignments. | library |
| [historical/run_cp97.py](../scripts/research/chloroplast/historical/run_cp97.py) | Evaluate the historical chloroplast topology after removing the Hachi tip. | historical_reference |
| [historical/run_expanded_trees.py](../scripts/research/chloroplast/historical/run_expanded_trees.py) | Execute the historical expanded 98-tip chloroplast comparison on admitted fixed-column matrices. | historical_reference |
| [historical/source_au_common.py](../scripts/research/chloroplast/historical/source_au_common.py) | Validate chloroplast alignments and canonicalize unrooted source-control topologies. | library |
| [historical/test_preparation_only.py](../scripts/research/chloroplast/historical/test_preparation_only.py) | Test the historical chloroplast matrix preparation without tree inference. | historical_test |
| [historical/validate_denovo_cp_reads.py](../scripts/research/chloroplast/historical/validate_denovo_cp_reads.py) | Validate de novo chloroplast contigs against paired read evidence. | historical_reference |
| [prepare_expanded_cp.py](../scripts/research/chloroplast/prepare_expanded_cp.py) | Admit validated additional plastomes into the fixed-column chloroplast matrix. | entrypoint |
| [recover_hachi_cp_callable.py](../scripts/research/chloroplast/recover_hachi_cp_callable.py) | Recover callable chloroplast bases from competitive-reference read alignments. | entrypoint |
| [run_cp_trees.py](../scripts/research/chloroplast/run_cp_trees.py) | Infer paired chloroplast trees on identical retained alignment columns. | linux_guarded_entrypoint |
| [run_source_au.py](../scripts/research/chloroplast/run_source_au.py) | Run matched chloroplast-source searches and approximately unbiased topology comparisons. | linux_guarded_entrypoint |
| [source_au_common.py](../scripts/research/chloroplast/source_au_common.py) | Validate chloroplast alignments and canonicalize unrooted source-control topologies. | library |
| [test_source_au.py](../scripts/research/chloroplast/test_source_au.py) | Validate archived source-control matrices, candidate unions and AU numeric checks. | historical_test |

## D. panthaica matched-marker comparison

| Script | Purpose | Role |
|---|---|---|
| [historical/exec_verified_scope.py](../scripts/research/panthaica/historical/exec_verified_scope.py) | Verify local resource limits for historical panthaica recovery steps. | runtime_guard |
| [historical/prepare_alignment_controls.py](../scripts/research/panthaica/historical/prepare_alignment_controls.py) | Add panthaica to frozen nuclear profiles and construct matched-column controls. | historical_reference |
| [historical/run_marker_recovery_draft.py](../scripts/research/panthaica/historical/run_marker_recovery_draft.py) | Recover BUSCO markers from the assembled panthaica transcriptome. | historical_reference |
| [historical/run_rna_pipeline_draft.py](../scripts/research/panthaica/historical/run_rna_pipeline_draft.py) | Perform the historical RNA read-QC and transcriptome-assembly stages. | historical_reference |
| [historical/test_alignment_controls.py](../scripts/research/panthaica/historical/test_alignment_controls.py) | Test preservation of original taxa and columns when adding panthaica. | historical_test |
| [historical/test_marker_recovery_draft.py](../scripts/research/panthaica/historical/test_marker_recovery_draft.py) | Test panthaica marker admission and sequence-translation checks. | historical_test |
| [historical/test_rna_pipeline_draft.py](../scripts/research/panthaica/historical/test_rna_pipeline_draft.py) | Test RNA-stage receipt and paired-read validation. | historical_test |
| [matched27/exec_verified_scope.py](../scripts/research/panthaica/matched27/exec_verified_scope.py) | Verify local resource limits before the matched-locus comparison. | runtime_guard |
| [matched27/prepare_matched27.py](../scripts/research/panthaica/matched27/prepare_matched27.py) | Prepare paired matched-locus, callable-site and leave-one-locus-out panthaica inputs. | entrypoint |
| [matched27/run_matched27.py](../scripts/research/panthaica/matched27/run_matched27.py) | Infer paired with/without-panthaica trees on the same 27 nuclear loci. | linux_guarded_entrypoint |
| [matched27/sensitivity/exec_verified_scope.py](../scripts/research/panthaica/matched27/sensitivity/exec_verified_scope.py) | Verify local resource limits before a panthaica sensitivity run. | runtime_guard |
| [matched27/sensitivity/run_sensitivity.py](../scripts/research/panthaica/matched27/sensitivity/run_sensitivity.py) | Evaluate paired callable-site masking and leave-one-locus-out panthaica sensitivity. | linux_guarded_entrypoint |

## Nuclear multispecies-coalescent simulations

| Script | Purpose | Role |
|---|---|---|
| [historical/audit_panel.py](../scripts/research/nuclear_ils/historical/audit_panel.py) | Audit conditional simulation manifests, generated alignments and inference jobs. | historical_reference |
| [historical/fit_history.py](../scripts/research/nuclear_ils/historical/fit_history.py) | Fit a fixed conditional species-history scenario to the retained nuclear evidence. | historical_reference |
| [inference/reestimate_panel.py](../scripts/research/nuclear_ils/inference/reestimate_panel.py) | Re-estimate gene trees from simulated alignments and summarize quartet outcomes. | linux_guarded_entrypoint |
| [inference/test_artificial.py](../scripts/research/nuclear_ils/inference/test_artificial.py) | Test simulated-panel inference admission and parsing with synthetic inputs. | test |
| [inference/vendor/engineering_parser.py](../scripts/research/nuclear_ils/inference/vendor/engineering_parser.py) | Parse IQ-TREE inference outputs and verify simulated-panel scientific invariants. | library |
| [inference/vendor/exec_verified_scope.py](../scripts/research/nuclear_ils/inference/vendor/exec_verified_scope.py) | Verify inherited local CPU and memory limits before an inference command. | runtime_guard |
| [inference/vendor/original_quartet_projector.py](../scripts/research/nuclear_ils/inference/vendor/original_quartet_projector.py) | Project inferred gene trees to the same fixed focal quartet as the observations. | library |
| [simulator/exec_verified_scope.py](../scripts/research/nuclear_ils/simulator/exec_verified_scope.py) | Verify inherited local CPU and memory limits before the simulation command. | runtime_guard |
| [simulator/simulate_panel.py](../scripts/research/nuclear_ils/simulator/simulate_panel.py) | Generate conditional nuclear genealogies and amino-acid alignments with preserved missingness. | linux_guarded_entrypoint |

## Conditional chloroplast genealogies and distance scenarios

| Script | Purpose | Role |
|---|---|---|
| [historical/make_report_figures.py](../scripts/research/chloroplast_ils/historical/make_report_figures.py) | Plot conditional plastid topology and distance-scenario summaries. | historical_reference |
| [historical/prepare_inputs.py](../scripts/research/chloroplast_ils/historical/prepare_inputs.py) | Prepare the frozen nuclear guide and observed chloroplast quartet/distance inputs. | historical_reference |
| [simulate.py](../scripts/research/chloroplast_ils/simulate.py) | Simulate conditional plastid genealogies and JC69 distance scenarios. | entrypoint |
| [validate.py](../scripts/research/chloroplast_ils/validate.py) | Validate simulated plastid topology frequencies and analytic distance probabilities. | entrypoint |

The [source map](../workflow/research/SOURCE_MAP.csv) records the origin, input/output hashes, and portability edits for each later-analysis source file. Tests in `tests/` cover small artificial inputs and engineering invariants; they are not research datasets. Software evidence is in [software_versions.md](../software_versions.md).

## Engineering tests

| Script | Purpose |
|---|---|
| [test_influence_fst.py](../tests/test_influence_fst.py) | Check influence fst with small artificial inputs and failure cases. |
| [test_influence_pca.py](../tests/test_influence_pca.py) | Check influence pca with small artificial inputs and failure cases. |
| [test_plastome_annotation.py](../tests/test_plastome_annotation.py) | Check plastome annotation with small artificial inputs and failure cases. |
| [test_population_audit.py](../tests/test_population_audit.py) | Check population audit with small artificial inputs and failure cases. |
| [test_population_diagnostics.py](../tests/test_population_diagnostics.py) | Check population diagnostics with small artificial inputs and failure cases. |
| [test_population_filters.py](../tests/test_population_filters.py) | Check population filters with small artificial inputs and failure cases. |
| [test_population_summary.py](../tests/test_population_summary.py) | Check population summary with small artificial inputs and failure cases. |
| [test_relationships.py](../tests/test_relationships.py) | Check relationships with small artificial inputs and failure cases. |
| [test_research_staging.py](../tests/test_research_staging.py) | Check research staging with small artificial inputs and failure cases. |
| [test_rna_evidence.py](../tests/test_rna_evidence.py) | Check rna evidence with small artificial inputs and failure cases. |
| [test_wgs_probe.py](../tests/test_wgs_probe.py) | Check wgs probe with small artificial inputs and failure cases. |
| [dryrun.sh](../workflow/tests/dryrun.sh) | Plan population/assembly jobs and check frozen-directory refusal without executing analysis. |
| [research_dryrun.py](../workflow/tests/research_dryrun.py) | Plan later-analysis jobs and check input/output contracts without executing analysis. |
