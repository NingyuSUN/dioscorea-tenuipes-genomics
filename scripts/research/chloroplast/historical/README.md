# Portable historical chloroplast stages

These eight entry points expose the historical analyses through explicit JSON
configuration. They retain fixed panels, sequence geometry, receipt checks,
scientific assertions, and bounded execution policies. They are **not** a generic
pipeline for arbitrary new samples, and portability tests are not reproduction of
the historical results.

Every entry point supports:

```bash
python SCRIPT.py --help
python SCRIPT.py --config /path/to/stage.json
python SCRIPT.py --config /path/to/stage.json --execute
```

Without `--execute`, the command reads configuration and dependency metadata,
checks that input files and executable files exist, and prints its plan. It does
not launch tools, create result directories, or validate the biological results.
Paths, including executable paths, are relative to the configuration file's
parent directory, unless explicitly absolute. A bare tool name is a relative
filename, not a `PATH` search. Use a fresh `out` directory for every execution.
Existing results and failed partial outputs are preserved; there is no automatic
restart. Normal Python is required: `python -O` disables scientific assertions
and is rejected.

The Python API `plan(config_path)` returns `inputs`, `outputs`, `tools`,
`named_inputs`, `named_outputs`, `named_tools`, and `output_root`. Paths are
absolute. Inputs include code dependencies and files inside declared nested
sample directories. `plan()` is dependency discovery; the command-line default
also checks path availability. A successful path preflight does not imply that
hash, sequence, sample identity, memory, or cgroup checks will pass during execution.

## Configuration keys

All six simple configurations include `out` plus these keys:

| Entry point | Input paths | Tool paths | Fixed scope |
| --- | --- | --- | --- |
| `prepare_inputs.py` | `source_dir`, `p19_dir` | `iqtree` | Copy and bind 93/94/95-tip matrices; 124448 columns; historical SHA-256 values |
| `prepare_source_inputs.py` | `design_dir`, `source_dir`, `historical_dir` | `iqtree` | Paired 10-tip matrices; 122234 columns; design and source hashes |
| `audit_hachi_nuclear_coverage.py` | `binding`, `genome`, `cp`, `reads` | `minimap2`, `samtools` | Coverage of 143 bound coding loci; fixed first-million-pair source hashes |
| `compare_published_edo.py` | `genotypes`, `cp_snps` | None; NumPy is required on execution | Published eight-sample encoded calls; 1054866 nuclear columns and 829 CP SNP columns |
| `published_edo_controls.py` | `genotypes`, `comparison_receipt` | None; NumPy is required on execution | Same genotype digest as the preceding comparison; scalar control checks |
| `validate_denovo_cp_reads.py` | `contigs`, `assembly_receipt`, `screen`, `reads` | `minimap2`, `samtools`; pysam on execution | Three admitted CP-like contigs; original coordinates and masked unsupported bases |

`reads` is an ordered two-element list `[R1, R2]`. `design_dir` contains the
metadata, design, partitions, and constraint files. `historical_dir` contains
`candidates.nwk` and `CANDIDATE_MANIFEST.json`. For example:

```json
{
  "genotypes": "inputs/published_calls.gt.gz",
  "cp_snps": "inputs/published_cp_snps.fa",
  "out": "results/published_comparison"
}
```

The corresponding controls configuration uses the same `genotypes` and
`comparison_receipt: "results/published_comparison/RECEIPT.json"`.

## Frozen expanded-tree and CP97 execution

`run_expanded_trees.py` retains the full98/control95 design and the original
strict scheduling contract. It requires Linux cgroups, the exact 6 GiB/no-swap
scope, CPU IDs 4–7, the full startup memory gate, recovery receipts, a shared tree
lock, a configured future timezone-aware deadline, and a lingering systemd user.
The configured `systemd_user` defaults to the current account; no historical
account or hostname is assumed. `required_hostname` is optional when an operator
needs an additional execution-host guard. The fixed CPU IDs and memory gates are
retained scientific execution constraints, not silently relaxed for a new machine.

Start from a copy of the native bridge policy and resource budget. Set these
portable filesystem fields: `out`, `entry_config`, `sample_parent`,
`frozen_manifest`, `python`, `loginctl`, `tree_lock`, `required_receipts`,
`policy.iqtree`, and `fixed_model.fixed_model_file_remote`. `source_profile`
defaults to `local`. `entry_script` and `legacy_runner` default to the bundled
`expanded_entry.py` and `cp_tree_runner.py`. The explicit
`latest_new_stage_start` replaces the expired historical scheduling default.

The entry configuration keeps its original `default_sample_parent` as provenance
used by sample receipts. Set the wrapper's `sample_parent` to the actual relocated
sample directory. Paths in `entry_config.source_profiles` resolve relative to
that entry configuration. For external sources referenced by immutable ledger
receipts, `path_map` maps historical path prefixes to actual directories; mapping
values resolve relative to the wrapper configuration. Source receipt bytes and
all hash comparisons remain unchanged. A mapping cannot replace a mismatched file.

`run_cp97.py` requires `out`, `program_freeze`, `python`, `systemd_run`,
`systemctl`, `taskset`, `timeout`, `inputs`, `executable`, `guard`, and the original
`resources`, `original_argv`, `columns`, `expected_labels98`, `deleted_label`, and
`alphabet`. Each `inputs` record and the `executable`/`guard` records contain
`path` and `sha256`. A bundled `exec_verified_scope.py` is available for `guard`.
Preserve the original successful runtime receipt and `original_argv` verbatim;
only execution filesystem locators are relocated, and the copied matrix/model
hashes are checked. The 98-to-97 matrix must remain a byte-exact deletion, with
124448 columns. CPU IDs 0–3, 6 GiB/no swap, the startup checks, four-hour hard
limit, 1000 bootstrap trees, exact tip checks, and failure receipts are retained.
`--precheck` performs the additional read-only frozen-input and resource checks.
`--execute` waits for the final bounded supervisor exit, so its return code is
suitable for a workflow. Internal worker modes require a bound dispatch token.

Generate a new code/configuration freeze **outside** the fresh result directory:

```bash
python freeze_contract.py --stage run_cp97 --config cp97.json
python freeze_contract.py --stage run_cp97 --config cp97.json --execute
# Or use --stage run_expanded_trees for the expanded bridge.
```

The configuration names the manifest location but does not contain its hash.
The manifest binds the configuration, entry point, and required helper code;
there is no self-referential hash. Freeze generation refuses to overwrite an
existing manifest. It is a configuration-integrity step, not scientific validation.

## Supported connections and limits

- `compare_published_edo` → `published_edo_controls`: same genotype input and the
  comparison receipt bind the two stages.
- `prepare_inputs` → `../run_cp_trees.py`: the prepared matrices, `CONFIG.json`,
  and `PREPARED.json` are the runner's input bundle.
- `prepare_source_inputs` → `../run_source_au.py`: the paired matrices, constraints,
  candidate manifest, and prepared configuration form the source-sensitivity bundle.
- Expanded trees → CP97: use the main98 terminal runtime/tree, the prepared
  control97 matrix, entry receipt, and identical model as an explicitly frozen
  CP97 contract. This is a separate opt-in stage, not an automatic fallback.
- Read-back validation alone cannot feed the expanded entry: fixed-coordinate
  projection, sample binding, ledger, quality, and assembly receipts are also
  required. Coverage-only output does not create a nuclear genotype or tree tip.

The engineering tests use synthetic files and fake tools. They test read-only
plans, path relocation, preserved partial results, failed tools, frozen-code
mismatches, exact sequence deletion, and resource rejection. They do not run
real mapping, real IQ-TREE inference, systemd scopes, or the historical data sets.
