"""Small configuration and preflight helpers for the historical drivers."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def read_config(description, argv=None):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", type=Path, required=True,
                        help="JSON config; relative paths resolve beside this file")
    parser.add_argument("--execute", action="store_true",
                        help="Run the configured calculation; default is read-only preflight")
    args = parser.parse_args(argv)
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text())
    if not isinstance(config, dict):
        raise ValueError("Configuration must be a JSON object")
    return config, config_path.parent, args.execute


def path_value(config, key, base, default=None):
    value = config.get(key, default)
    if not isinstance(value, (str, Path)) or not str(value):
        raise ValueError("Missing path configuration: " + key)
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def require_file(path, executable=False):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing or empty input file: " + str(path))
    if executable and not os.access(path, os.X_OK):
        raise ValueError("File is not executable: " + str(path))
    return path


def require_directory(path):
    if not path.is_dir():
        raise ValueError("Missing input directory: " + str(path))
    return path


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def expected_hash(path, expected, label):
    if not isinstance(expected, str) or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
        raise ValueError("Expected lowercase SHA256 is required for " + label)
    actual = sha256(require_file(path))
    if actual != expected:
        raise ValueError("SHA256 mismatch: " + label)
    return actual


def fresh_output(path, inputs=()):
    if path.exists():
        raise ValueError("Output already exists; choose a fresh output directory: " + str(path))
    for source in inputs:
        source = Path(source).resolve()
        if path == source or path in source.parents or source in path.parents:
            raise ValueError("Input and output paths overlap: " + str(path))
    return path


def report_plan(script, inputs, outputs, dependencies=(), **details):
    tool_keys = {'plink', 'megahit', 'busco_python', 'busco_script', 'mafft', 'iqtree'}
    named_inputs = {key: str(value) for key, value in inputs.items()}
    named_outputs = {key: str(value) for key, value in outputs.items()}
    required_threads = details.pop("required_threads", details.get("threads", 1))
    if not isinstance(required_threads, int) or isinstance(required_threads, bool) or required_threads < 1:
        raise ValueError("required_threads must be a positive integer")
    return dict(status="PREFLIGHT_PASSED", script=script, required_threads=required_threads,
                inputs=list(named_inputs.values()), outputs=list(named_outputs.values()),
                tools=[value for key, value in named_inputs.items() if key in tool_keys],
                named_inputs=named_inputs, named_outputs=named_outputs,
                dependencies=list(dependencies), **details)
