"""Validate native driver status and resource contracts without changing analyses.

Only a new workflow completion marker is declared as a Snakemake output. Native
receipts and partial scientific files stay in place when a driver fails.
"""
import argparse
import hashlib
import json
from pathlib import Path


def check_cp_resources(inputs, allocated_threads):
    config = json.loads((Path(inputs) / "CONFIG.json").read_text())
    if allocated_threads < 1:
        raise ValueError("Allocated threads must be positive")
    for key in ("threads", "cpu_affinity_max_cpus"):
        value = config[key]
        if type(value) is not int or not 1 <= value <= allocated_threads:
            raise ValueError(f"Prepared CP {key} exceeds the allocated threads or is invalid; allocate sufficient cores without changing scientific parameters")
    return {key: config[key] for key in ("threads", "cpu_affinity_max_cpus")}


def complete(native, completion, kind):
    native, completion = Path(native), Path(completion)
    data = native.read_bytes()
    receipt = json.loads(data)
    expected = {"loci": "COMPLETE_REVIEW_PENDING", "cp_trees": "FULL_COMPLETE_SCIENTIFIC_REVIEW_PENDING"}[kind]
    if receipt.get("status") != expected:
        raise ValueError(f"Native driver did not complete: {receipt.get('status')}")
    if kind == "loci":
        requested = receipt.get("requested_loci", 0)
        if requested < 1 or receipt.get("failed") or len(receipt.get("completed", [])) != requested:
            raise ValueError("Locus completion count or failure list disagrees with successful status")
    else:
        jobs = receipt.get("jobs", [])
        if not jobs or any(job.get("status") != "COMPLETE_REVIEW_PENDING" for job in jobs):
            raise ValueError("A CP tree job lacks a successful native completion status")
    result = {"status": "DRIVER_COMPLETE_SCIENTIFIC_REVIEW_PENDING", "kind": kind,
              "native_receipt": native.name, "native_sha256": hashlib.sha256(data).hexdigest(),
              "native_status": receipt["status"], "scientific_validation_claimed": False}
    completion.parent.mkdir(parents=True, exist_ok=True)
    with completion.open("x") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    resources = commands.add_parser("cp-resources")
    resources.add_argument("--inputs", required=True, type=Path)
    resources.add_argument("--allocated-threads", required=True, type=int)
    finish = commands.add_parser("complete")
    finish.add_argument("--native", required=True, type=Path)
    finish.add_argument("--completion", required=True, type=Path)
    finish.add_argument("--kind", required=True, choices=["loci", "cp_trees"])
    args = parser.parse_args()
    if args.command == "cp-resources":
        check_cp_resources(args.inputs, args.allocated_threads)
    else:
        complete(args.native, args.completion, args.kind)


if __name__ == "__main__":
    main()
