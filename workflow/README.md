# Snakemake workflow

The `Snakefile` chains the numbered scripts in `scripts/` into one reproducible workflow. It does not modify the scripts. The default target covers the three population-genomics modules:

- filtering, then PCA, ADMIXTURE and FST (scripts 32–34);
- individual relationships (44–49);
- sample influence (53–55).

The historical assembly and annotation scripts (01–15) are available as an optional template module.

An end-to-end run reproduced the original genotype panels, PCA, ADMIXTURE and FST byte for byte. See [docs/WORKFLOW_VERIFICATION.md](../docs/WORKFLOW_VERIFICATION.md).

## Quick start

```bash
conda create -y -n snakemake -c conda-forge -c bioconda 'snakemake-minimal>=8'
conda activate snakemake
# Edit config/config.yaml: input paths, a new results_dir, and tool locations
snakemake -n                                     # dry run
snakemake --profile workflow/profiles/default    # run (16 cores, 26 GB memory cap)
```

To let Snakemake create the software environments instead:

```bash
snakemake --sdm conda --profile workflow/profiles/default \
  --configfile config/config.yaml config/config.conda.yaml
```

Targets:

- `all`: the modules listed in `config["targets"]`;
- `population_all`, `relationships_all`, `influence_all`;
- `assembly_all`: needs `--configfile config/config.yaml config/assembly.yaml`.

## What each rule runs

| Rule | Script | Main output (under `results_dir`) |
|---|---|---|
| `pop_filter_variants` | 32 | `population_analysis/{masked17,main14,strict14,complete14,supplement17}.bcf` |
| `pop_structure` | 33 | `population_analysis/structure/` (PCA, 36 ADMIXTURE runs, 6 FST) |
| `pop_report` | 34 | `population_analysis/report/` |
| `rel_prepare_panels` | 44 | `relationship_panels/` (19 panels) |
| `rel_relatedness` | 46 | `relationship_panels/relatedness/` (KING, PC-Relate) |
| `rel_stage` | (symlinks) | `relationships/data/` (22 panels) |
| `rel_primary17_flat_distance` | PLINK | `relationships/data/primary17_flat.mdist` |
| `rel_analyze` | 45 | `relationships/results/` (147 panels) |
| `rel_r_nj_network`, `rel_r_fixed_order` | 47, 47b | ape/phangorn cross-checks |
| `rel_network_cycle`, `rel_network_compare` | helper | `cycle.txt`, `network_cross_validation.json` |
| `rel_validate` | 49 | `relationships/VALIDATION.json` |
| `infl_pca_panels` | 53 | `influence/` (38 PCA panels) |
| `infl_fst` | 55 | `influence/fst/` (114 FST) |
| `infl_summarize` | 54 | `influence_summary/` |

Arguments, thread pools and PLINK memory settings come from the scripts; the rules only declare them, so the scheduler can respect the memory cap in the profile.

## Safety design

The numbered scripts refuse to write into an existing output directory. The workflow keeps that behaviour and adds three barriers:

1. **`frozen_dirs`** (`config/config.yaml`): the workflow stops at parse time if `results_dir` is one of these directories or lies inside one.
2. **`protected()` on every data output** of the computational rules: Snakemake makes them read-only and raises `ProtectedOutputException` *before deleting anything* if a job would rerun, including under `--forcerun` and `--forceall`.
3. **`fresh()`**: Snakemake pre-creates the parents of declared outputs, which would trip the scripts' `mkdir(exist_ok=False)`. Before each script, `fresh()` removes that tree only if it holds nothing but empty directories. Otherwise it stops with "Refusing to reuse non-empty …".

If a job fails midway, inspect the partial directory and move it aside before rerunning. Do not rely on `--rerun-incomplete`.

## Software environments

`workflow/envs/` pins the versions used in the original analysis:

- `popgen.yaml`: PLINK 1.90b6.21, ADMIXTURE 1.3.0, VCFtools 0.1.16. The full package list of the original environment is `popgen_production.explicit.txt`.
- `relationships_py.yaml`: numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.2, biopython 1.88, splitspy 0.0.10.
- `relationships_r.yaml`: R 4.5 with SNPRelate 1.44.0, GENESIS 2.40.0, GWASTools 1.56.0, ape 5.8-1, phangorn 2.12.1.

Known limitations:

- conda cannot install bcftools/htslib 1.13 alongside Python 3.12, so the conda environments resolve to a current bcftools. The original analysis used bcftools 1.13. To match it exactly, point `tools.extra_path` at a bcftools 1.13 build (for example the Ubuntu 22.04 package).
- splitspy 0.0.10 supplies only the Neighbor-Net circular ordering; weights are fitted with SciPy NNLS (see [docs/NETWORK_SOLVER_AUDIT.md](../docs/NETWORK_SOLVER_AUDIT.md)).

## Network cross-check helper

`workflow/scripts/network_crosscheck.py` writes the fixed circular ordering read by script 47b. It then compares the NNLS split weights with phangorn's same-order and automatic-order networks. The reported `maximum_weight_difference` is taken over the union of splits, with an absent split weighing 0.

## Assembly module (template)

`config/assembly.yaml` enables rules for scripts 01–15, from Guppy basecalling through BRAKER. The original runs used an unrecorded environment (`environment.lock.template.yml` is best-effort), so this module documents the chain rather than reproducing the original assembly. Replace the placeholder paths before use.

Not wrapped:

- organelle scripts 16–20 and 35–43,
- the historical variant-calling and audit scripts 21–31,
- the public-data probes 48–52.

## Lint

`snakemake --lint` reports findings that are accepted by design:

- output prefixes passed as directories (the scripts take directories and assume a shared filesystem),
- no conda env on the assembly rules,
- string concatenation in shell directives,
- no log on the small copy/link rules,
- a `run:` block in `rel_stage`.

`workflow/tests/dryrun.sh` fails only on shell commands that read global variables.

## CI

`.github/workflows/ci.yml` runs on each push:

- all unit tests in `tests/`, with bcftools/htslib and the compiled genotype-comparison helper;
- `workflow/tests/dryrun.sh`, which dry-runs both targets with placeholder inputs and checks the `frozen_dirs` guard.
