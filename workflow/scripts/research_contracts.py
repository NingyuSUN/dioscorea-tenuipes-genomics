"""Enumerate prepared research inputs without reading scientific payloads."""
import json
from pathlib import Path

CP_ILS_MEMBERS = ["inputs/guide.json", "inputs/observed.json", "inputs/hime2_hachi_joint_acgt.fa",
                  "inputs/cp98_tree.nwk", "inputs/cp97_tree.nwk"]


def read_json(path):
    with Path(path).open() as handle:
        return json.load(handle)


def member(root, name):
    root = Path(root).resolve()
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Prepared input member must be relative and contained: {name}")
    path = root / relative
    if not path.resolve().is_relative_to(root):
        raise ValueError(f"Prepared input member escapes its source root: {name}")
    if not path.is_file():
        raise ValueError(f"Required input is missing or not a file: {name}")
    return str(path)


def nuclear_inputs(root):
    manifest = member(root, "metadata/LOCUS_MANIFEST.json")
    rows = read_json(manifest)["loci"]
    if not isinstance(rows, list) or not rows:
        raise ValueError("Nuclear locus manifest must contain a nonempty loci list")
    names = [row["gene"] for row in rows]
    if any(not isinstance(gene, str) or not gene or Path(gene).name != gene or gene in (".", "..") for gene in names):
        raise ValueError("Each locus gene must be a plain nonempty filename stem")
    if len(set(names)) != len(names):
        raise ValueError("Duplicate locus names are not allowed")
    return [manifest] + [member(root, f"loci/{gene}.faa") for gene in names]


def cp_ils_inputs(root):
    return [member(root, name) for name in CP_ILS_MEMBERS]


def cp_tree_inputs(root):
    prepared = member(root, "PREPARED.json")
    configuration = member(root, "CONFIG.json")
    receipt, config = read_json(prepared), read_json(configuration)
    files = receipt["files"]
    if not isinstance(files, dict) or not files:
        raise ValueError("CP PREPARED.json must declare its input files")
    panels = config["sequential_panels"]
    if not isinstance(panels, list) or not panels or len(set(panels)) != len(panels):
        raise ValueError("CP sequential_panels must contain unique prepared panels")
    required = {"CONFIG.json", config["model_file"]} | {panel + ".fixed_columns.fa" for panel in panels}
    if not required.issubset(files):
        raise ValueError("CP PREPARED.json does not bind every configured scientific input")
    paths = [prepared, configuration] + [member(root, name) for name in files]
    # The native driver hashes this executable from the job's working directory.
    binary = Path(config["iqtree"])
    if not binary.is_file():
        raise ValueError("Prepared IQ-TREE executable is missing")
    return list(dict.fromkeys(paths + [str(binary)]))
