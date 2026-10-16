# Portable stage configurations

The catalogue covers 22 scientific scripts. RNA QC, prefix pilot and full assembly use three configurations of one script, giving 24 stage instances. These templates require verified external inputs; they do not include research data or executable tools.

Paths resolve relative to each JSON file. Place immutable inputs under `data/portable`, install reviewed executables under `tools`, and write new runs under `results_portable`. Edit the concrete paths when using another layout. Keep raw reads outside the generated output directories.

Run a script with `--config CONFIG.json` for read-only preflight. Add `--execute` only to perform the declared stage. The workflow templates are consumed by `workflow/portable/Snakefile`; a dry run plans dependency edges and does not establish scientific validity.

## Required evidence

- Replace every all-zero SHA-256 value with a reviewed digest. Preserve pinned original hashes. A mismatched hash is an error, not permission to recompute the expected value.
- The full-population plan must bind all nine nested input files, the portable runner, the reference implementation and all four installed PCAngsd modules. Put the verified PCAngsd package below `tools/pcangsd` and configure the Python environment to import those exact files. The 37566-site, 14-sample panel and original fit settings remain fixed.
- The `nuclear_loci.json` null tree digest is filled only from the declared upstream tree completion receipt by the nuclear chain. For direct standalone use, supply the reviewed reference-tree digest explicitly.
- The nuclear taxa file must point to the complete original BUSCO run inputs. The example workflow tracks `data/portable/busco_runs`, the offline lineage and cache. Curate the full taxon panel before execution; candidate recovery does not automatically admit a taxon.
- `metadata/` contains byte-preserved public metadata templates and their reviewed migration proof. Its hashes remain bound in the RNA/marker configurations. Map typed evidence paths with `path_roots` or `path_map`; do not edit the evidence to make checks pass.
- CP97 keeps the archived `original_argv` as immutable provenance. The executable, alignment, model and output paths used for execution are configured separately. The expanded CP run requires a reviewed future deadline, program freeze, complete entry configuration and source receipts.

## Workflow boundaries

`population_chain.json`, `nuclear_preparation_chain.json`, `chloroplast_edo_chain.json`, `nuclear_chloroplast_ils_chain.json` and `panthaica_chain.json` describe the feasible dependencies. `all_stages.json` includes all 24 instances, including independent validation branches. Use fresh output and workflow-receipt directories for an independent run.

The workflow binds explicitly declared upstream artifacts for downstream checks and seals generated directories. Existing expected hashes remain authoritative. Frozen scientific dimensions, lineage restrictions, fixed seeds and execution resource guards still apply. Full scientific replay requires the matching evidence, software and resource allocation; configuration and synthetic tests alone do not prove that replay.

Nuclear simulation generation/inference and conditional chloroplast simulation validation remain separate stages with their existing interfaces. The historical nuclear audit recounts its completed fixed panel. The figure stage consumes validated summary tables and cannot substitute for simulation validation. RNA full-assembly eligibility does not prove marker orthology or independent specimen identity.
