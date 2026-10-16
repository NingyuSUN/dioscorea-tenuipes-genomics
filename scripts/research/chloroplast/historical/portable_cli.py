"""Read-only configuration and path contracts for historical chloroplast stages."""
import argparse
import json
import os
from pathlib import Path


def read_config(path):
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError('Configuration must be a JSON object')
    return path, data


def resolve(base, value):
    if not isinstance(value, str) or not value:
        raise ValueError('A nonempty path string is required')
    p = Path(value).expanduser()
    return (p if p.is_absolute() else Path(base).parent / p).resolve()


def configure(path, paths=(), lists=()):
    path, data = read_config(path)
    for key in paths:
        data[key] = resolve(path, data[key])
    for key in lists:
        if not isinstance(data[key], list) or not data[key]:
            raise ValueError(key + ' must be a nonempty list of paths')
        data[key] = [resolve(path, x) for x in data[key]]
    data['_config'] = path
    return data


def contract(config, inputs, outputs, tools=None, output_root=None, **details):
    tools = tools or {}
    named_inputs = {'configuration': str(config), 'portable_helper': str(Path(__file__).resolve()),
                    **{k: str(v) for k, v in inputs.items()}}
    named_outputs = {k: str(v) for k, v in outputs.items()}
    return {'status': 'PLAN_ONLY_NOT_EXECUTED',
            'inputs': list(dict.fromkeys(named_inputs.values())),
            'outputs': list(dict.fromkeys(named_outputs.values())),
            'tools': list(dict.fromkeys(str(v) for v in tools.values())),
            'named_inputs': named_inputs, 'named_outputs': named_outputs,
            'named_tools': {k: str(v) for k, v in tools.items()},
            'output_root': str(output_root) if output_root else None, 'required_threads': 1, **details}


def preflight(plan):
    """Check paths without starting tools, creating directories, or reading matrices."""
    if not __debug__:
        raise ValueError('Python optimization disables scientific assertions; use normal Python')
    for name in plan['inputs']:
        if not Path(name).is_file():
            raise ValueError('Missing input file: ' + name)
    for name in plan['tools']:
        if not Path(name).is_file() or not os.access(name, os.X_OK):
            raise ValueError('Missing or non-executable tool: ' + name)
    root = Path(plan['output_root']) if plan.get('output_root') else None
    if root:
        if root.exists():
            raise ValueError('Fresh output directory required; existing results are preserved: ' + str(root))
        for name in [*plan['inputs'], *plan['tools']]:
            if Path(name).is_relative_to(root):
                raise ValueError('Output directory would contain an input: ' + name)
    return {**plan, 'status': 'PREFLIGHT_ONLY_NOT_EXECUTED',
            'scientific_validation': 'Not performed by path preflight'}


def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--config', type=Path, required=True,
                   help='JSON configuration; paths are relative to this file, never the working directory')
    p.add_argument('--execute', action='store_true', help='Execute after preflight; default only checks the plan')
    return p


def emit(plan):
    print(json.dumps(plan, indent=2, allow_nan=False))


def existing_parent(path):
    """Find the existing filesystem ancestor for a read-only capacity check."""
    p=Path(path)
    while not p.exists():
        if p==p.parent:raise ValueError('No existing output filesystem ancestor')
        p=p.parent
    if not p.is_dir():raise ValueError('Output parent is not a directory: '+str(p))
    return p
