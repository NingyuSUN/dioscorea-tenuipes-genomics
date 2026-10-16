"""Resolve later-analysis DAGs and negative guards with empty temporary inputs.

Every Snakemake invocation uses --dry-run. No scientific job is launched.
"""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import yaml


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snakemake", default="snakemake")
    args = p.parse_args()
    repo = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="dt workflow ") as tmp:
        root = Path(tmp)
        source = root / "input"
        for name in ["nuclear/metadata/LOCUS_MANIFEST.json", "nuclear/loci/example.faa", "reference.treefile",
                     "population/panel.bed", "population/panel.bim", "population/panel.fam",
                     "chloroplast_ils/inputs/guide.json", "chloroplast_ils/inputs/observed.json",
                     "chloroplast_ils/inputs/hime2_hachi_joint_acgt.fa", "chloroplast_ils/inputs/cp98_tree.nwk",
                     "chloroplast_ils/inputs/cp97_tree.nwk", "chloroplast/PREPARED.json"]:
            path = source / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("{}\n")
        (source / "nuclear/metadata/LOCUS_MANIFEST.json").write_text(json.dumps({"loci": [{"gene": "example"}]}))
        cp = source / "chloroplast"
        cp_config = {"iqtree": str(cp / "fake-iqtree"), "threads": 4, "cpu_affinity_max_cpus": 4,
                     "model_file": "model.nex", "sequential_panels": ["main95"]}
        (cp / "CONFIG.json").write_text(json.dumps(cp_config))
        for name in ["fake-iqtree", "model.nex", "main95.fixed_columns.fa"]:
            (cp / name).write_text("synthetic dry-run placeholder\n")
        (cp / "PREPARED.json").write_text(json.dumps({"files": {name: {} for name in
            ["CONFIG.json", "model.nex", "main95.fixed_columns.fa"]}}))
        for filename, key in [("spatial_plan.json", "input_prefix"), ("balanced_plan.json", "source_prefix")]:
            (source / "population" / filename).write_text(json.dumps({key: str(source / "population/panel")}))
        override = {"research": {"outdir": str(root / "output"),
                    "nuclear_input_root": str(source / "nuclear"), "reference_tree": str(source / "reference.treefile"),
                    "spatial_plan": str(source / "population/spatial_plan.json"),
                    "balanced_plan": str(source / "population/balanced_plan.json"),
                    "chloroplast_ils_input_root": str(source / "chloroplast_ils"),
                    "cp_inputs": str(source / "chloroplast")}, "frozen_dirs": [str(root / "frozen")]}
        config = root / "override.yaml"

        def run(target, values=override, fail=False):
            config.write_text(yaml.safe_dump(values))
            cmd = [args.snakemake, target, "--dry-run", "--cores", "4", "--configfile",
                   "config/config.yaml", "config/research.yaml", str(config)]
            result = subprocess.run(cmd, cwd=repo, text=True, capture_output=True)
            text = result.stdout + result.stderr
            if fail:
                assert result.returncode != 0 and "completed evidence is read-only" in text, text
            else:
                assert result.returncode == 0, text
            return text

        plan = run("research_all")
        for rule in ["research_locus_trees", "research_gene_concordance", "research_spatial_fst",
                     "research_balanced_fst", "research_cp_ils_simulate", "research_cp_ils_validate"]:
            assert rule in plan, rule
        run("research_chloroplast_trees_all")
        # Independent output roots and symlink aliases cannot bypass freezing.
        (root / "frozen").mkdir()
        (root / "alias").symlink_to(root / "frozen", target_is_directory=True)
        for out in [root / "frozen/run", root / "alias/run"]:
            bad = {**override, "research": {**override["research"], "outdir": str(out)}}
            run("research_all", bad, fail=True)
        # --dry-run must not have created either research output directory.
        assert not (root / "output").exists()
        print("Research DAGs, space-containing paths and frozen-output guards: PASS; no analysis executed.")


if __name__ == "__main__":
    main()
