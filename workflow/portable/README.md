# Configured historical stages

The 22 historical entry points now accept a JSON configuration and expose safe
command-line help. Relative paths in a stage configuration resolve beside that
configuration, independent of the invocation directory. Default invocation
performs read-only preflight; `--execute` is required to enter the scientific
calculation. A preflight can read and hash existing inputs, but never launches
an analysis or creates its result directory.

```bash
python scripts/research/nuclear/historical/prepare_loci.py --help
python scripts/research/nuclear/historical/prepare_loci.py --config config/portable/nuclear_prepare_loci.json
```

The configuration examples contain generic locations. Supply matching input
bundles and tools, then keep local paths in ignored `*.local.json` files. Fixed
scientific counts, taxa, seeds, column geometry and evidence hashes still apply.
Changing a path does not make a different panel interchangeable with the one a
script was written to evaluate.

## Connected workflows

Use the repository root as the working directory. A workflow manifest declares
stage configurations, consumed files, produced files and upstream stage names.
Paths inside a manifest resolve beside that manifest; script paths resolve from
the repository root. The examples below are plans until `-n` is removed:

```bash
snakemake -s workflow/portable/Snakefile -n --cores 8 \
  --config manifest=config/portable/population_chain.json
```

Each upstream edge must name the producer of an input. External inputs must
exist before planning. A future output is represented by the producer's
completion receipt, so an initially empty output directory can still produce a
complete dependency graph. Completed input, code and output files are checked against recorded SHA-256
hashes before planning. Files that changed or disappeared cause a clear refusal.
This integrity check reads the files and may take time for large datasets; it
does not rerun calculations.

Only completion receipts are Snakemake outputs. Native output files remain in
place on failure, and a failed stage cannot create a completion receipt. A
successful exit must also produce every declared nonempty output file. Native
script checks remain responsible for scientific invariants and status checks.
Existing outputs and receipts are never silently overwritten; preserve a failed
run and configure a new output location for another attempt. Source files,
configurations, declared tools and all research helper scripts are tracked as
dependencies. `input_directories` recursively declares existing source bundles;
`generated_input_directories` consumes a directory owned by one upstream stage,
whose completion receipt supplies a stable dependency even before files exist.
`output_directories` seals all native files, including dynamically named loci.
Before execution, the script's read-only plan must agree with declared inputs,
tools, output locations and allocated CPU count. Input/code/configuration hashes
must remain unchanged across execution.

For explicit `bind_inputs` entries, a downstream configuration can bind a newly
produced receipt to its upstream completion hash. Only declared upstream outputs
are eligible. The workflow writes a separate execution configuration with
absolute paths and hashes; the supplied configuration remains unchanged. An
already specified, different expected hash is rejected. This binds newly
produced artifacts and does not replace frozen scientific identities.

`threads` is a scheduling request. The workflow refuses allocation below the
requested count. Linux cgroup, swap and affinity requirements in individual
drivers remain enforced by those drivers; Snakemake does not create those
operating-system limits. Some stages need a separately prepared scope receipt
before execution. Their preflight describes that requirement.

## What the checks establish

The release checks use syntax parsing, safe CLI help, deliberately tiny fixtures,
fake external tools and Snakemake dry runs. No original dataset, simulation or
phylogenetic reconstruction was rerun. These checks establish the configuration
and orchestration behavior exercised by the tests, not scientific reproduction.

The [stage catalogue](../research/PORTABLE_STAGES.json) lists the 22 entry points,
configuration examples and declared dependencies. The [script catalogue](../../docs/SCRIPT_CATALOG.md)
explains their scientific purposes. Recorded historical versions remain in the
[software table](../../software_versions.md); the current Python check environment
is listed in [research requirements](../envs/research-requirements.txt).

The [dependency table](DEPENDENCIES.md) covers external tools and platform guards.

## Available chain manifests

| Manifest in `config/portable/` | Connected purpose |
|---|---|
| `population_chain.json` | Full-panel fits → same-panel PLINK arm → descriptive covariance comparison. |
| `nuclear_preparation_chain.json` | Curated nuclear matrix/tree → exact fixed-panel locus preparation. |
| `chloroplast_edo_chain.json` | Published-plastome comparison → source/callable-site controls. |
| `nuclear_chloroplast_ils_chain.json` | Conditional nuclear history → linked-plastid guide/input preparation. |
| `panthaica_chain.json` | Separate QC, pilot and full RNA stages → marker recovery → matched alignment controls. |
| `all_stages.json` | Inventory DAG for all 22 scripts (24 phase configurations); distinct scientific panels remain separate branches. |

Each manifest needs its declared external bundles. The inventory is not a claim
that all analyses start from one common input. A configured null downstream hash
in a chained example is filled only by an explicit verified upstream binding;
standalone execution still requires a concrete reviewed hash.
