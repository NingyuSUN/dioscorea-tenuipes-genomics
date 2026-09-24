# Workflow verification (2026-09-24)

Before this workflow existed, scripts 32–55 had been run one at a time. The Snakemake workflow was then run end to end from the same population VCF, into a new directory, and every output was compared with those original runs.

## Run

- 19/19 jobs completed in about 40 minutes on a 36-core, 31 GB RAM Linux server, with the default profile.
- Same software versions as the original runs:
  - PLINK 1.90b6.21, ADMIXTURE 1.3.0, VCFtools 0.1.16, bcftools 1.13;
  - R 4.5.3 with SNPRelate 1.44.0 and GENESIS 2.40.0;
  - Python 3.12 with numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.2, biopython 1.88 and splitspy 0.0.10.

## Results of the comparison

**Byte-identical to the original runs:**

- the genotype records of all 24 BCF files;
- all PLINK panels;
- PCA eigenvectors and eigenvalues;
- all ADMIXTURE Q matrices and CV errors;
- all FST estimates, including the 114 sample-influence estimates;
- all population-report tables.

**Different only in run metadata:**

- logs (paths, timestamps, PLINK's printed seed, ADMIXTURE iteration times);
- BCF headers (command line and date) and their indices;
- JSON files that record input paths or header-dependent hashes;
- SVG metadata.

**Equal to floating-point precision:** relationship and influence summaries, with maximum differences of 3–7 × 10⁻¹⁶. The validation checks produced identical results.

**One substantive difference:** the Neighbor-Net split weights of one leave-one-out panel.

## Neighbor-Net circular ordering and near-ties

For that panel the distance matrix and NJ tree were byte-identical, but splitspy's Neighbor-Net returned a different circular ordering. The NNLS split weights therefore differed by up to 7 × 10⁻⁴.

- **Uninitialised memory was ruled out.** Filling splitspy's internal work matrix with NaN, 0 or 10⁶ changes none of the 147 orderings, and repeated runs on a fixed matrix are deterministic.
- **The cause is near-ties in the agglomeration.** Adding symmetric relative noise of about 10⁻¹⁵ to the distances changes the ordering in 31–38 of the 147 panels, mostly half-contig subsets. The main 17-individual network is not affected. The last floating-point bits, which can differ between platforms, decide the join order.
- **Effect on conclusions.** Weights were recomputed under every alternative ordering found.
  - Each locality's own split keeps a positive weight: minimum 0.0287, against a largest shift of 0.0013.
  - Network fit changes by at most 0.0011 percentage points.

  The result that locality splits are retained in every panel does not depend on how the ties are broken. Network weights of the sensitivity panels should not be reported to more than about three decimals.

See also `NETWORK_SOLVER_AUDIT.md`: splitspy's native weight solver is not used.

## What this shows

The computation is reproducible from the population VCF with the recorded software. It is not new biological evidence.
