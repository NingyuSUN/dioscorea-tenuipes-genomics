# Project index — Dioscorea tenuipes genomics

## Current state

- **Objective:** make the current assembly, population and marker QC easy to inspect in the public README.
- **Phase:** documentation update prepared from preserved results; no scientific rerun.
- **Last verified:** 2026-10-19T15:08:54+09:00 (Asia/Tokyo). Live public main equals baseline `2bc5ab2bf9ce3d7c0d673ff0b2320dd2b32daa4b`. The earlier CI environment repair succeeded in [run 38847206068](https://github.com/NingyuSUN/dioscorea-tenuipes-genomics/actions/runs/38847206068); it is no longer pending.
- **Entry points:** [README QC](README.md#quality-control), [aggregate QC snapshot](docs/provenance/qc_summary.json), [script purposes](docs/SCRIPT_CATALOG.md), [software versions](software_versions.md), [workflows](workflow/README.md).
- **Verified scope and limits:** existing h2 continuity, genotype/site thresholds, retained population-panel counts, PCA sensitivity and marker sequence correspondence. These do not establish genome-wide base accuracy, specimen identity or an evolutionary mechanism. Full raw inputs/result archives remain external. Historical engineering checks remain distinct from scientific validation.
- **Git / publication:** branch `public-update-20261016`, baseline above; README, QC snapshot and index are intended for one local commit followed by the requested public-main update. Remote state after publication is verified separately; query `git log -1 -- PROJECT_INDEX.md` for this record's commit.
- **Next:** use the published summary to review the existing QC and associated code; no analysis is scheduled by this documentation change.

## 2026-10-19 — Public QC summary

| Field | Record |
|---|---|
| Date | 2026-10-19T15:08:54+09:00, Asia/Tokyo. Source results retain their September–October 2026 dates. |
| Purpose | Add existing QC details to the public README in English, using aggregate evidence and explicit denominators. Documentation and source verification only. |
| Inputs and code | Preserved assembly/read audits, population result summaries, nuclear CDS-correspondence table and plastid sensitivity report. Source-record checksums are in the QC snapshot; public scripts 32–34 and the nuclear matrix/CDS-audit implementations were checked against the text. Known public baseline `2bc5ab2bf9ce3d7c0d673ff0b2320dd2b32daa4b`; no scientific code executed or large data rehashed. |
| Main results and limits | Current h2: 403,772,333 bp, 9,423 contigs, N50 103,611 bp. Main population panel: 14/16 study individuals, 5,564,311 SNPs, 37,566 PCA SNPs; strict/complete-call/one-per-contig comparisons preserve the recorded high distance correlations. Nuclear correspondence: 149/173 exact, all 24 other records retained. QC categories and limits are stated; no claim of blanket short-contig removal or whole-genome validation. |
| Deliverables and validation | README quality-control section and `docs/provenance/qc_summary.json`; numeric/source, relative-link, privacy and scoped-diff checks. Existing prior CI success verified live; this is not a new scientific validation. Original correspondence, private infrastructure details and biological mechanism claims are excluded. |
| Git / GitHub | https://github.com/NingyuSUN/dioscorea-tenuipes-genomics ; local branch `public-update-20261016`, publication target `main`. Baseline was live-verified equal to public main before editing; these three files will be committed and published. Publication readback is required before reporting synchronization. Raw data and full private evidence archives are outside this upload. |

## 2026-10-16 — CI environment compatibility repair

| Field | Record |
|---|---|
| Date | 2026-10-16, Asia/Tokyo; environment failure inspected after public push. |
| Purpose | Repair the GitHub Actions environment without changing scientific code or rerunning analysis. |
| Inputs and code | Published portability commit `ac18c48a19f8519fd4fb86f80a3a1b2d59c2ef3f`; Actions run `38846967458`. Its setup log reports Conda msprime 1.4.4 requiring GSL 2.8 while available bcftools requires GSL below 2.8. |
| Results and limits | Move the same pinned msprime/tskit versions to pip wheels, matching the locally checked Python environment. The 119 local tests remain passed; the first remote run never reached tests. Successful CI after this repair must be verified separately. No scientific calculation was executed. |
| Deliverables | `workflow/envs/ci.yaml`; updated portability verification note and this index. |
| Git / GitHub | Remote `main` was live-verified at `ac18c48a19f8519fd4fb86f80a3a1b2d59c2ef3f`. This environment-only repair is pending commit/push at index preparation. Code/configuration/docs are the synchronized scope; raw data remain external. |

## 2026-10-16 — Configurable historical stages

| Field | Record |
|---|---|
| Date | 2026-10-16T12:40:13+09:00; Asia/Tokyo. Engineering checks completed during this round. |
| Purpose | Port the 22 historical directory-bound scripts, complete explicit dependencies, and connect compatible stages. This was code/configuration maintenance and synthetic verification, with zero historical data analyses. |
| Inputs and code | Public baseline `58319f3a385ff764053160d7bc8d62b1ed02d876`; branch `public-update-20261016`. The 22 selected sources, original/distributed hashes and edits are in `workflow/research/SOURCE_MAP.csv`. Input scope: source files, reviewed metadata and small artificial fixtures; research sample/site counts are not applicable to this execution round. Reproducible checks: `python -m unittest discover -s tests -v`, `python workflow/tests/portable_checks.py`, `python workflow/tests/portable_dryrun.py`, and the existing workflow checks. |
| Results and limits | 22 configurable entry points; 24 phase configurations; five focused chains plus a combined inventory DAG. 119 combined tests, 22 CLI checks and six portable DAG plans passed. Source-preservation and privacy reviews completed. Failed native outputs are preserved; completion requires stable input hashes and required outputs. Real scientific execution, installed external build compatibility and live resource scopes were not tested. |
| Deliverables | `scripts/research/*/historical/`; `config/portable/`; `workflow/portable/`; `workflow/research/PORTABLE_STAGES.json`; `docs/PORTABLE_SCRIPTS_2026-10-16.md`; `docs/provenance/portable_code_audit.json`; updated source map, catalogue, dependencies and CI. |
| Git / GitHub | Public repository https://github.com/NingyuSUN/dioscorea-tenuipes-genomics ; target `main`. Known baseline `58319f3a385ff764053160d7bc8d62b1ed02d876`; this round is pending commit/push at index preparation. The containing commit is available with `git log -1 -- PROJECT_INDEX.md`; final synchronization is verified against the live remote after push. Coverage is public code/configuration/docs only; private evidence archives and raw data remain outside the upload. |

## 2026-10-16 — Later code and workflow release

| Field | Record |
|---|---|
| Date | 2026-10-16, Asia/Tokyo; release checks performed during this update. |
| Purpose | Restore later scientific code for method review, explain each script's purpose, document versions and visualize the pipeline. This was code preservation and engineering verification, not new analysis. |
| Inputs and code | Public baseline `774f4cf9f51a36dbe744e1243bff6d7db35e374d`; selected preserved research source snapshot `1ad0747e7b1c4631d05d2d7620863bb53407f6c7`; exact logical source paths and original/distributed SHA-256 values in `workflow/research/SOURCE_MAP.csv`. Working branch `public-update-20261016`. Research data sample/site range: not applicable to this release because no research payload was analysed. |
| Results and limits | 70 later source files, a per-script purpose catalogue, three diagrams, recorded software versions and selected research Snakemake targets. 66 artificial-input tests passed; existing and new DAG checks passed. Engineering checks do not establish new scientific or external validation. See release verification for coverage and unresolved limits. |
| Deliverables | `README.md`; `docs/SCRIPT_CATALOG.md`; `software_versions.md`; `docs/provenance/`; `docs/figures/`; `scripts/research/`; `config/research.yaml`; `workflow/rules/research.smk`; `workflow/research/`. |
| Git / GitHub | Public repository: https://github.com/NingyuSUN/dioscorea-tenuipes-genomics ; target branch `main`. At index preparation, remote main was live-verified at baseline `774f4cf9f51a36dbe744e1243bff6d7db35e374d`; this release was pending commit/push. Query `git log -1 -- PROJECT_INDEX.md` for the containing release commit. Final publication is verified separately against the remote branch, without repeatedly editing this index to name its own commit. Raw reads, private correspondence and manuscripts are not uploaded by this release. |

## 2026-09-24 — Original public code release

The original public tree contained 104 tracked files, including scripts 01–55 (with 47b), population workflows and an optional assembly template. Its latest baseline is `774f4cf9f51a36dbe744e1243bff6d7db35e374d`. The public Git history and [workflow verification](docs/WORKFLOW_VERIFICATION.md) provide the available record. Earlier analysis history is not reconstructed here; unrecorded dates and environments remain unverified.
