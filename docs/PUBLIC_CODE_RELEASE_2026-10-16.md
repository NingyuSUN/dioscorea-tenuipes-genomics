# Public code update — 2026-10-16

This release adds later-analysis source code, selected Snakemake wrappers, a complete script-purpose catalogue, recorded software versions, and three pipeline diagrams. No research dataset was reanalysed and no simulation or phylogenetic inference was rerun for the release.

## Content

- **70 later source files** across population sensitivity, nuclear markers, chloroplast comparisons, panthaica comparisons, nuclear coalescent models and conditional chloroplast models.
- Every numbered script (56 files, including 47b), later source file and shared helper has a purpose in [SCRIPT_CATALOG.md](SCRIPT_CATALOG.md).
- [SOURCE_MAP.csv](../workflow/research/SOURCE_MAP.csv) preserves original source identifiers/hashes, distributed-file hashes and descriptions of portability edits. Source paths are logical archive identifiers; the private evidence archive is not part of this public release.
- [software_versions.md](../software_versions.md) and its evidence table distinguish historical runs, environment snapshots, configuration constraints and current code checks. Unknown runtime versions remain unknown.
- Three vector figures in [figures/](figures/README.md) have a deterministic editable generator.

## Engineering verification

Checks were performed on 2026-10-16 (Asia/Tokyo), using the recorded [current environment](provenance/current_validation_environment.json).

| Check | Result and scope |
|---|---|
| Existing and new unit/integration tests | **66 passed, 0 failed**. Inputs were small artificial fixtures; five new integration cases used fake analysis drivers. |
| Existing population/assembly dry run | Passed, including the existing lint policy and frozen-directory guard. |
| Later-analysis dry runs | Passed for the default research modules and optional CP tree target, including paths containing spaces and frozen/symlink refusal. No scientific jobs executed. |
| Later CLI and source checks | 26 safe entry-point help checks passed; all 70 sources parsed and matched their distributed source-map hashes. |
| Original scientific helper checks | Nine concordance tests, six quartet fixtures, and artificial panthaica projection/site-mask/leave-one-locus preparation checks passed. |
| Source preservation | Original and distributed hashes were compared. The [preservation summary](provenance/code_preservation_summary.json) records unchanged bytes and the limited adapters/path normalizations. |
| Diagrams | All three SVGs were rendered and visually checked; regeneration was deterministic. |
| Public content | Checked for private identities/paths, operational metadata, credentials and broken local document links. Source correspondence, manuscript drafts and large data were excluded. |

The source comparison found 496 original functions: 493 retained identical parsed syntax trees. The remaining three were two command-line entry functions and a MAFFT executable-path substitution; after normalizing that tool path, the latter function also matched. The new panthaica adapter delegates projection/masking to the preserved scientific functions.

## Issues corrected

- Added directly imported helper modules that would otherwise have been missing from the later code release.
- Changed tree-workflow completion outputs to separate markers. A native driver that exits successfully but records a partial failure now fails the workflow while preserving its status and partial files.
- Declared individual prepared-input files, rather than depending only on directory timestamps. Changes inside nested nuclear/chloroplast input directories are now visible to the DAG.
- Checked the prepared chloroplast thread/affinity settings against allocated cores before launching the driver.
- Extended frozen-output checks to assembly and research roots, including symlink aliases and ancestor overlap.
- Removed old operational/temporary-directory labels from historical source paths and recorded those edits explicitly.

## Limits

The 22 `historical_reference` sources retain stage-specific data/layout contracts and are not advertised as automatic default workflow steps. Source restoration and dry-run checks do not establish end-to-end scientific reproduction. Native cgroup-restricted scientific execution was not tested in this update; those resource limits remain required. Completion markers record successful native driver status, not an independent revalidation of every sidecar output.

The earlier [population replay record](WORKFLOW_VERIFICATION.md) applies to its documented inputs and date. It is not evidence that all newly added research modules were replayed. Exact original environment versions that were not recovered are not replaced with guessed or current versions.
