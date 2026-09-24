# Population input audit (scripts 30–31)

Checks that the population VCF matches the joint-call VCF it was derived from, and diagnoses sample quality, before any population analysis. All source VCFs are read-only; use a new output directory for every run. Positions are on the original contigs, not scaffold coordinates.

## Inputs and versions

- Main VCF `Dtenuipes_population.vcf.gz`: 4,560,509,276 bytes, SHA256 `43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9`.
- Joint VCF `joint_calls_merged.vcf.gz`: 6,822,793,877 bytes, SHA256 `156db9911bf5790dfb0dc5a33b2f728f6d77884ee529148872428a358a4746a0`.
- Reference `genome.fa`: SHA256 `36fecc3e8b3054728d9779d22d61957062fd914bf24cd04f09e27d363d97c87b`.
- bcftools 1.13 / HTSlib 1.13. Python comparison: pysam 0.22.1; the diagnostics also need numpy, and the figures matplotlib.

## Full comparison

```bash
bcftools stats --threads 2 -s - Dtenuipes_population.vcf.gz > out/main.stats.txt
bcftools stats --threads 2 -s - joint_calls_merged.vcf.gz > out/joint.stats.txt
bcftools norm -f genome.fa -c e -Ou Dtenuipes_population.vcf.gz > /dev/null 2> out/reference_check.log
./compare_population_genotypes joint_calls_merged.vcf.gz Dtenuipes_population.vcf.gz \
  out/full_comparison > out/full_comparison.json 2> out/progress.log
```

How the C helper compares records:

- It uses HTSlib to parse both files completely.
- It merges coordinate groups, then matches REF and the complete ordered ALT list within each group.
- It verifies old/new sample aliases and contig order, keeps missing GT in the comparison, and ignores phase.
- It rejects duplicate identical records, unsorted records and ploidy above two.
- It does not infer equivalence between differently normalized allele representations.

Success requires exit status 0, a completed JSON, and counts consistent with the independent `bcftools stats`.

Other checks:

- The Python comparison checks the header dictionary against the FASTA.
- `bcftools norm -c e` verifies every REF allele of the main VCF; its output is discarded.
- Both files are checked for their BGZF EOF block, exact size and SHA256.

Compile the helper against HTSlib:

```bash
cc -O3 -Wall -Wextra scripts/lib/compare_population_genotypes.c \
  $(pkg-config --cflags --libs htslib) -o compare_population_genotypes
```

## Independent bounded Python comparison

```bash
python scripts/30_audit_population_vcf.py --left joint_calls_merged.vcf.gz \
  --right Dtenuipes_population.vcf.gz --reference genome.fa --limit 100000 --out out/prefix_comparison
```

- `--limit` counts merged coordinate groups, not matching variants.
- QC statistics include one lookahead group, while the comparison counts include only processed groups.
- Such runs are marked `complete: false` and must not be extrapolated to the whole genome.
- `--limit 0` scans the whole files but is much slower than the C helper.

## Sample diagnostics

```bash
python scripts/31_diagnose_population_samples.py --vcf Dtenuipes_population.vcf.gz --fai genome.fa.fai --out out/diagnostics
```

Requires a matching VCF index.

**Window selection**: a fixed RNG seed (20260914) selects 1,000 positions uniformly across the concatenated reference, and 1-kb windows around them are merged to avoid double counting. Final coverage: 998 windows, 996,926 reference bases, 50,043 retained SNPs.

**Thresholds**:

- Diagnostic calls require a non-missing diploid GT, DP ≥5 and GQ ≥20.
- Allele-balance summaries additionally require heterozygosity, DP ≥10 and a two-element AD summing to ≥10.
- These are provisional diagnostic cutoffs, not final analysis filters.

**Limits**:

- Windows are not LD-pruned.
- Pairwise concordance is a screen, not a kinship estimate.
- Allele-balance summaries cannot diagnose contamination, hybrid ancestry or biological ploidy.

## Tests

```bash
DTE_AUDIT_BINARY=/absolute/path/compare_population_genotypes \
  python -m unittest discover -s tests -p 'test_population*.py'
```

The tests cover:

- allele and coordinate matching;
- missing-to-reference changes, including partial and haploid missing calls;
- phase-insensitive comparison and multi-record coordinates;
- rejection of duplicate records, bad order and ploidy above two;
- low-depth exclusion and allele-balance arithmetic.

Passing tests is an engineering check, not biological validation.
