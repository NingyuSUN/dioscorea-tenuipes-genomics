# Population analysis: filters, methods and limits (scripts 32–34)

Two samples, AIC-M6-S11 and AIC-F20-S13, are excluded because their original reads/BAM could not be recovered for inspection. Main analysis: AIC 4, HCO 4, IZU 6. The public individuals Mie1, Yuga and Utsu are supplementary individuals, not a pooled fourth population.

## Inputs and coordinates

- Population VCF: 20,625,692 biallelic SNPs, 19 samples. SHA256 of the file used: `43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9`.
- All positions stay on the original 9,423-contig reference. Script 22's chromosome relabeling is not reused.
- Sample manifest: TSV with `sample` and `group` columns, in VCF sample order.

## Tools

PLINK 1.90b6.21, ADMIXTURE 1.3.0, VCFtools 0.1.16, bcftools 1.13 with the setGT/fill-tags plugins, and Python 3.12 with numpy, pandas, matplotlib and scipy. The exact conda packages of the environment used are listed in `workflow/envs/popgen_production.explicit.txt`.

## Run

The Snakemake workflow runs scripts 32–34 in order (`snakemake population_all`; see `workflow/README.md`). To run them by hand:

```bash
python -u scripts/32_filter_population_variants.py --vcf population.vcf.gz --manifest sample_manifest.tsv --out out/population_analysis
python -u scripts/33_run_population_structure.py --data out/population_analysis --out out/population_analysis/structure --tools <dir with plink, admixture, vcftools>
python scripts/34_summarize_population_results.py --data out/population_analysis --analysis out/population_analysis/structure --audit <dir with main.sample_qc.bcftools.tsv> --out out/population_analysis/report
```

Each step refuses to overwrite its existing output directory.

## Filters fixed before seeing population results

- **Primary genotype filter**:
  - diploid complete GT, DP 5–100 inclusive, GQ ≥20;
  - missing DP/GQ and partial GT become `./.`; no genotype is filled as reference;
  - the high-depth cap is a QC heuristic and does not certify the absence of paralogs.
- **Site filters**: PASS biallelic SNP; ≤10% missing in the selected sample set; minor allele count ≥1 after recalculating AC/AN from the retained genotypes. With 14 or 17 samples, ≤10% allows at most one missing individual per site.
- **Stricter sensitivity**: DP 8–60 inclusive, GQ ≥30, same site missingness.
- **Complete-case sensitivity**: primary genotype thresholds with zero missingness in all 14 study individuals.
- **Supplementary panel**: 17 individuals, primary thresholds, its own site filtering. Differences in the 17-individual plot therefore reflect both sample addition and panel selection; this is not a formal projection experiment.

**PCA/ADMIXTURE only**:

- MAC ≥2.
- Space SNPs at least 1 kb apart within each original contig, then PLINK `--indep-pairwise 50kb 1 0.2`. This reduces local LD; residual LD across contig boundaries cannot be assessed from a contig assembly.
- One deterministic hash-selected SNP per represented contig gives an additional PCA sensitivity check. It does not establish physical independence of contigs.

**PCA reporting**:

- Raw eigenvalues are preserved.
- Figures omit explained-variance percentages, because relationship matrices built with missing data have small negative eigenvalues. Positive-spectrum fractions in the JSON output are labelled as such.
- Missing calls stay missing; PLINK applies its internal handling.
- Missingness and mean depth are checked for correlation with PC1/PC2. These descriptive checks do not prove freedom from confounding.
- PCA is exploratory with 4/4/6 individuals and unresolved relatedness and library metadata.

**ADMIXTURE**:

- K = 1–4, three fixed seeds (20260914–20260916), five-fold CV, panels primary14, strict14 and supplement17.
- Mean CV is compared across seeds only within a panel.
- For each K, the plot shows the run with the highest full-data likelihood. Component permutations are aligned across seeds and panels using the shared study individuals.
- Components are statistical clusters, not proof of ancestral populations, migration or hybridization.
- ADMIXTURE's input parser needs numeric chromosome codes, so after pruning a parser-only copy of the BIM with numeric contig codes is made. Genotypes and positions are unchanged and the mapping is saved. This copy must never be used for coordinate-aware analyses.

**FST**:

- Pairwise Weir–Cockerham estimator in VCFtools, panels primary14 and strict14.
- Uses the full filtered MAC ≥1 panel (no PCA thinning). Reports weighted overall estimates and 100-kb windows on the original contigs. Negative estimates are retained.
- These are descriptive estimates conditional on the retained SNPs, not selection tests, and no confidence interval is derived from window averages.
- SNP discovery, missingness, small sample size and unresolved sample dependence can all affect the estimates.
- Public samples are not pooled for FST.

## Scope boundaries

No all-sites/callable-region denominator is available. Variant-only window π, Tajima's D, demographic models and runs of homozygosity are not valid results from this workflow. It provides no mapping-level QC, reference-bias test, ploidy confirmation or raw-read diagnosis. Sex inferred from F/M labels is unconfirmed and must not become a biological covariate without metadata verification.

## Primary documentation

- [bcftools setGT](https://samtools.github.io/bcftools/howtos/plugin.setGT.html)
- [PLINK LD pruning](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK PCA](https://www.cog-genomics.org/plink/1.9/strat)
- [ADMIXTURE manual](https://dalexander.github.io/admixture/admixture-manual.pdf). The hosted manual is newer than 1.3.0, so the installed CLI behavior was tested directly.
- [VCFtools FST](https://vcftools.github.io/man_latest.html)
