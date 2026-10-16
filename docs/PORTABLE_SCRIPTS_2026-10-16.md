# Historical-script portability update — 2026-10-16

This update converts all **22 historical reference scripts** from release
`58319f3` into configurable entry points. Every one now has safe `--help`, a JSON
configuration, a read-only default preflight and explicit `--execute`. All new
public documentation, configuration, help and workflow messages are in English.
No research analysis, read mapping, assembly, model fitting, tree inference or
simulation was rerun for this update.

## What changed

- Input/output/tool locations are explicit and resolve beside the configuration
  file. Paths inside immutable historical receipts can be relocated through
  explicit mappings without rewriting the evidence files.
- Directly required helpers are bundled, and the [dependency table](../workflow/portable/DEPENDENCIES.md)
  identifies external executables, Python packages, lineage data and resource
  controls. Historical versions remain evidence-based in [software_versions.md](../software_versions.md).
- [24 stage configurations](../config/portable/) cover the 22 scripts, including
  separate RNA QC, pilot and full-assembly phases. Five focused manifests connect
  population fit/comparison, nuclear matrix/locus preparation, published-plastome
  controls, nuclear-guide/chloroplast-input preparation, and panthaica stages.
- The [optional Snakefile](../workflow/portable/Snakefile) tracks explicit files,
  source bundles, generated directories and upstream receipts. It reconciles a
  script's preflight with the declared dependencies, outputs and CPU allocation.
  Only separate completion receipts belong to Snakemake; partial native outputs
  remain available after failure.
- The [script catalogue](SCRIPT_CATALOG.md), [source map](../workflow/research/SOURCE_MAP.csv)
  and [stage catalogue](../workflow/research/PORTABLE_STAGES.json) describe purposes,
  configurations, dependencies and fixed-panel limits.

These are compatible stage chains, not one interchangeable raw-data workflow.
CP95/94, CP98/97, the fixed nuclear panel and the matched panthaica controls retain
their distinct taxon, site, seed, sequence and evidence requirements.

## Problems found and corrected

The nuclear matrix builder writes a miniprot audit sidecar into its BUSCO inputs.
The configured driver now stages only consumed BUSCO files into its fresh output
root, verifies the copies and runs the unchanged builder there. Artificial tests
confirm that original source files stay unchanged.

Two previously sanitized metadata files no longer matched their original frozen
hashes. The original Git backup reproduced both original hashes. Comparison
identified only path/provenance strings, including an accidentally altered public
adapter-manual URL, which was restored. No scientific scalar, list structure,
marker identity, sequence hash or lineage hash changed. The [migration receipt](../config/portable/metadata/METADATA_MIGRATION.json)
binds original and portable file hashes. A configured replacement must explicitly
supply the original binding and a matching migration receipt; general hash checks
remain active.

Workflow review found and corrected four integrity problems: completion hashes
were not rechecked; configured dependencies could be omitted from the manifest;
input hashes were recorded only after execution; and scheduler CPU requests were
not reconciled with native requirements. Further checks cover generated-directory
ownership and undeclared sidecars. Inputs/code/configurations are now hashed before
and after execution, and only unchanged inputs can produce a completion receipt.
Completed file hashes are rechecked on subsequent plans. This costs file I/O for
large bundles but does not repeat calculations.

## Verification scope

Final local checks on 2026-10-16 (Asia/Tokyo):

| Check | Result |
|---|---|
| Combined unit/integration suite | **119 passed**, zero failures. |
| Portable entry-point checks | **22/22 passed**; help, explicit execution flags and JSON examples. |
| Portable workflow templates | **6/6 dry runs passed**, covering five focused chains and all 24 stage instances. |
| Existing workflow checks | Original population/assembly lint/dry run and later research DAG/guard checks passed. |
| Additional preserved panthaica checks | **39 passed** on artificial fixtures. |
| Independent review | 15 entry points and nested helpers reviewed; 19 workflow integrity tests passed with source hashes stable. |
| Public-content review | No unresolved privacy findings; no Chinese text in newly changed public files. Original data/manuscripts excluded. |


Tests use artificial sequences/tables, synthetic fixed-dimension matrices and fake
executables. A synthetic 173-locus fixture tests splitting and missing-row handling;
it does not reproduce the historical biological panel. Snakemake checks plan the
example graphs without starting any scientific stage. Fake-stage tests execute
only purpose-built fixture scripts, including failures and input mutations.

Original and distributed hashes remain in the source map. The [AST comparison](provenance/portable_code_audit.json)
identifies unchanged functions and those relocated or adapted. Structural
preservation is an engineering check, not scientific validation.

## Remaining limits

Real-data replay and execution under actual cgroup/systemd limits remain untested.
PCAngsd production fitting was not executed, and its four reviewed build-specific
module hashes are still required. Placeholder paths and all-zero example digests
must be replaced with matching local artifacts before execution; failures of
those checks are intentional. The publication does not supply raw reads, private
correspondence or the complete historical result archive.

A new derived receipt can be bound automatically only when an explicit manifest
binding names a verified completed upstream artifact. This writes a separate
execution configuration; it never changes the user's configuration or silently
replaces an already specified different hash. Fixed historical scientific
identities and admission checks are not relaxed by that mechanism.

## CI environment repair

The first GitHub run for `ac18c48` stopped before tests because Conda's msprime
1.4.4 required GSL 2.8 while available bcftools packages required GSL below 2.8.
The CI specification installs the same pinned msprime 1.4.4 and tskit 1.0.3 through
pip wheels, matching the locally tested Python environment and leaving bcftools
with its compatible Conda libraries. This changes environment installation only;
no scientific scripts or research results are changed by the repair.
