# Later-analysis workflow

The later scripts are grouped by their biological or technical purpose. The complete [script catalogue](../../docs/SCRIPT_CATALOG.md) describes every file; [ENTRYPOINTS.tsv](ENTRYPOINTS.tsv) lists command-line options; [SOURCE_MAP.csv](SOURCE_MAP.csv) records source stages and hashes.

## Scope

| Module | Purpose | Available execution entry |
|---|---|---|
| Nuclear markers | Build BUSCO amino-acid matrices, audit marker identity, infer locus trees, and count concordant/alternative splits | `build_busco_matrices.py`; `research_nuclear_all` for prepared loci → gene trees → concordance |
| Population sensitivity | Compare allele-balance filters, likelihood treatments, genomic regions and balanced sample sets | CLI scripts; `research_population_all` for prepared regional/balanced FST plans |
| Chloroplast comparisons | Prepare expanded panels, test source/topology alternatives and audit sequence recovery | CLI scripts; `research_chloroplast_trees_all` for prepared 95/94-taxon tree runs |
| Panthaica comparisons | Prepare matched markers and compare with/without-taxon and site/locus sensitivity arms | `matched27/prepare_matched27.py`, `run_matched27.py`, `sensitivity/run_sensitivity.py` |
| Nuclear coalescent models | Generate conditional gene histories and re-estimate their relationships | `simulator/simulate_panel.py`, `inference/reestimate_panel.py` |
| Chloroplast coalescent models | Simulate linked plastid genealogies and distance scenarios, then check them independently | `research_chloroplast_ils_all` |

These are distinct prepared-input workflows. The 173-locus nuclear panel, earlier marker sensitivity panels, CP95/94 and CP98/97 comparisons, and matched27 comparisons do not share one interchangeable input definition. Historical assertions about hashes, taxa, sites and seeds are retained. This release does not regenerate those input bundles or claim a complete raw-read-to-manuscript replay.

## Snakemake targets

Use the repository root as the working directory. Configure `config/research.yaml`, preferably through an ignored `config/research.local.yaml` overlay. The commands below only plan jobs:

```bash
snakemake research_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml

snakemake research_nuclear_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml

snakemake research_population_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml

snakemake research_chloroplast_ils_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml

# Optional bounded Linux driver; not included in research_all by default.
snakemake research_chloroplast_trees_all -n --cores 4 \
  --configfile config/config.yaml config/research.yaml
```

`research_all` uses `research.modules`, independently of the original `targets` list. Remove `-n` only when intentionally executing an analysis with the matching prepared inputs and tools. No scientific job was executed while preparing this public update.

## Input and output contracts

| Target | Prepared input contract | Main output under `research.outdir` |
|---|---|---|
| Nuclear | `metadata/LOCUS_MANIFEST.json`, the manifest's `loci/<gene>.faa` files, a hash-bound reference tree and IQ-TREE executable | Per-locus trees/receipts, branch concordance and per-locus classification tables |
| Population | Two JSON plans plus their PLINK BED/BIM/FAM files; original sample order, panel definitions and hashes retained | `spatial_fst/REGION_COMPONENTS.tsv.gz`, `balanced_fst/FST_ESTIMATORS.tsv`, receipts |
| Chloroplast ILS | `inputs/guide.json`, `observed.json`, `hime2_hachi_joint_acgt.fa`, `cp98_tree.nwk`, `cp97_tree.nwk` | `results_<phase>/` topology counts, distance scenarios and `VALIDATION.json` |
| Chloroplast trees | `PREPARED.json`, `CONFIG.json`, prepared panel FASTA/model files and the hash-bound IQ-TREE executable | Per-panel tree/support files and native run status |

Population plan paths are resolved from the invocation working directory, as in the original scripts. Use absolute paths or run from the repository root consistently. Source-map paths identify retained archive files; they are not files automatically downloaded by this workflow. Data bundles are supplied separately.

The staging helper copies only declared nuclear/conditional-plastid inputs into a separate output root. Source and output roots must not overlap. `frozen_dirs` applies to population, assembly and research output roots, including symlink aliases. Completed workflow receipts are protected; native scientific outputs and partial failure receipts are retained for inspection. Use a new output directory after a failed or completed attempt; do not force a rerun over retained evidence.

The chloroplast and other bounded historical drivers require Linux cgroup v2 memory/swap constraints and CPU affinity. Snakemake `mem_mb` is scheduling metadata, not an operating-system memory cap. These resource contracts are not removed to make a dry run appear executable.

## Standalone examples

```bash
python scripts/research/nuclear/build_busco_matrices.py --help
python scripts/research/population/prepare_ab_beds.py --help
python scripts/research/chloroplast/run_source_au.py --help
python scripts/research/panthaica/matched27/prepare_matched27.py --help
python scripts/research/nuclear_ils/simulator/simulate_panel.py --help
```

The 22 historical stages now support safe `--help`, configurable paths, read-only default preflight and explicit `--execute`. [Portable workflows](../portable/README.md) connect compatible stages through declared artifacts and preserve failures. Their fixed panel/evidence contracts remain active; they are separate optional workflows, not additions to the default research DAG.

## Dependencies and checks

[Software versions](../../software_versions.md) contains the stage-specific evidence. The research Python requirements describe the current engineering-check environment; external executables and Linux resource controls remain separate. Hash-preserving assertions are part of the input contract, not obstacles to bypass.

```bash
python -m unittest discover -s tests -v
bash workflow/tests/dryrun.sh
python workflow/tests/research_dryrun.py
```

The test inputs are tiny and artificial. Syntax, CLI, hash-preservation, dry-run planning and artificial-input checks do not establish scientific reproduction or biological validity. Conditional lineage-sorting simulations do not yield a calibrated probability of the proposed evolutionary mechanism.
