# Project index — Dioscorea tenuipes genomics

## Current state

- **Objective:** provide a readable public code repository covering the later analysis scripts, workflows, purposes and software versions.
- **Phase:** public code/documentation release; no new research analysis.
- **Last verified:** 2026-10-16T11:56:14+09:00 (Asia/Tokyo). Sources, engineering tests and workflow dry runs checked; details in [release verification](docs/PUBLIC_CODE_RELEASE_2026-10-16.md).
- **Entry points:** [README](README.md), [script catalogue](docs/SCRIPT_CATALOG.md), [software versions](software_versions.md), [workflow contracts](workflow/research/README.md).
- **Limits:** historical scripts retain stage-specific input contracts; full scientific replay and original unrecorded runtime versions remain unverified. Raw sequencing data and manuscript drafts are outside this public code release.
- **Next:** review the published code and purpose map; any scientific rerun requires a separate execution decision and matching inputs.

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
