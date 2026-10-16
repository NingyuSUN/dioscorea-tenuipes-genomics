"""Copy declared research inputs to a fresh, separate run directory."""
import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from research_contracts import nuclear_inputs, cp_ils_inputs


def stage(source, destination, kind):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError("Source and output must be separate, non-nested directories")
    paths = {"nuclear": nuclear_inputs, "chloroplast_ils": cp_ils_inputs}[kind](source)
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("Output directory is not empty; inspect it before reuse")
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for path in paths:
        src = Path(path)
        name = src.relative_to(source)
        dst = destination / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        digest = hashlib.sha256(src.read_bytes()).hexdigest()
        if hashlib.sha256(dst.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Staged input copy differs from source: {name}")
        hashes[str(name)] = digest
    with (destination / "STAGED.json").open("x") as handle:
        json.dump({"kind": kind, "members": list(hashes), "sha256": hashes}, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--kind", choices=["nuclear", "chloroplast_ils"], required=True)
    args = parser.parse_args()
    stage(args.source, args.destination, args.kind)
