"""Dependency planning and receipt handling for configured historical stages.

Only explicit execution runs scientific scripts. Failed native outputs remain
in place; Snakemake owns only a separate completion receipt for each stage.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def resolve(base, value):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def overlaps(a, b):
    return a == b or a in b.parents or b in a.parents


def load_manifest(path):
    if not __debug__:
        raise ValueError('Python optimization disables scientific assertions; use normal Python')
    path = Path(path).resolve()
    data = json.loads(path.read_text())
    base = path.parent
    receipts = resolve(base, data['receipts'])
    frozen = [resolve(base, p) for p in data.get('frozen_dirs', [])]
    stages = data['stages']
    if not stages:
        raise ValueError('The workflow has no stages')
    producers = {}
    for name, stage in stages.items():
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', name):
            raise ValueError(f'Invalid stage name: {name}')
        stage['script'] = resolve(REPO, stage['script'])
        if not stage['script'].is_file() or not stage['script'].is_relative_to(REPO / 'scripts/research'):
            raise ValueError(f'Stage script must exist under scripts/research: {name}')
        stage['config'] = resolve(base, stage['config'])
        if not stage['config'].is_file():
            raise FileNotFoundError(f'Stage configuration is missing: {name}: {stage["config"]}')
        stage['inputs'] = [resolve(base, p) for p in stage['inputs']]
        stage['tools'] = [resolve(base, p) for p in stage.get('tools', [])]
        stage['inputs'].extend(stage['tools'])
        stage['generated_input_directories'] = [resolve(base, p) for p in stage.get('generated_input_directories', [])]
        stage['path_fields'] = stage.get('path_fields', [])
        stage['bind_fields'] = {key: resolve(base, value) for key, value in stage.get('bind_fields', {}).items()}
        stage['bind_inputs'] = {key: resolve(base, value) for key, value in stage.get('bind_inputs', {}).items()}
        stage['input_directories'] = [resolve(base, p) for p in stage.get('input_directories', [])]
        stage['output_directories'] = [resolve(base, p) for p in stage.get('output_directories', [])]
        for folder in stage['input_directories']:
            if not folder.is_dir():
                raise FileNotFoundError(f'External input directory is missing: {folder}')
            stage['inputs'].extend(p.resolve() for p in sorted(folder.rglob('*')) if p.is_file())
        stage['inputs'] = list(dict.fromkeys(stage['inputs']))
        for folder in stage['output_directories']:
            if overlaps(folder, receipts) or any(overlaps(folder, root) for root in frozen):
                raise ValueError(f'Output directory overlaps receipts or frozen evidence: {folder}')
            if any(overlaps(folder, p) for p in stage['inputs']):
                raise ValueError(f'Output directory overlaps source evidence: {folder}')
        stage['outputs'] = [resolve(base, p) for p in stage['outputs']]
        stage['after'] = stage.get('after', [])
        if not stage['outputs'] or len(stage['outputs']) != len(set(stage['outputs'])):
            raise ValueError(f'Stage needs unique explicit output files: {name}')
        stage['threads'] = int(stage.get('threads', 1))
        if stage['threads'] < 1:
            raise ValueError(f'Invalid thread count: {name}')
        for out in stage['outputs']:
            if out in producers:
                raise ValueError(f'Multiple producers for {out}')
            producers[out] = name
            if overlaps(out, receipts) or any(overlaps(out, root) for root in frozen):
                raise ValueError(f'Output overlaps receipts or frozen evidence: {out}')
    for root in frozen:
        if overlaps(receipts, root):
            raise ValueError('Receipt directory overlaps frozen evidence')
    for out in producers:
        for stage in stages.values():
            for src in [stage['config'], stage['script']]:
                if overlaps(out, src):
                    raise ValueError(f'Output overlaps code/configuration: {out}')
        for other in producers:
            if out != other and overlaps(out, other):
                raise ValueError(f'Nested output declarations: {out}, {other}')
    output_owners = [(name, folder) for name, stage in stages.items() for folder in stage['output_directories']]
    for owner, folder in output_owners:
        for other, other_folder in output_owners:
            if owner != other and overlaps(folder, other_folder):
                raise ValueError('Output directories of different stages must not overlap')
        for stage in stages.values():
            if any(overlaps(folder, source) for source in stage['input_directories']):
                raise ValueError('External input directory overlaps generated output; declare generated_input_directories and an upstream edge')
    def ancestors(name, active=()):
        if name in active:
            raise ValueError('Cyclic stage dependencies')
        found = set()
        for parent in stages[name]['after']:
            if parent not in stages:
                raise ValueError(f'Unknown upstream stage: {parent}')
            found.add(parent)
            found.update(ancestors(parent, (*active, name)))
        return found
    for name, stage in stages.items():
        upstream = ancestors(name)
        for folder in stage['generated_input_directories']:
            owners = [other for other in upstream for root in stages[other]['output_directories'] if folder == root or folder.is_relative_to(root)]
            if len(set(owners)) != 1:
                raise ValueError(f'Generated input directory needs one upstream output owner: {folder}')
        for role, src in {**stage['bind_inputs'], **stage['bind_fields']}.items():
            if producers.get(src) not in upstream or src not in stage['inputs']:
                raise ValueError(f'Bound input must be an explicitly consumed upstream artifact: {role}')
        for src in stage['inputs']:
            producer = producers.get(src)
            if producer and producer not in upstream:
                raise ValueError(f'{name} consumes {src} without an upstream dependency')
            if not producer:
                if any(src == folder or src.is_relative_to(folder) for owner, folder in output_owners):
                    raise ValueError('Input is an unlisted generated artifact; declare generated_input_directories or an explicit produced artifact')
                if any(overlaps(src, out) for out in producers) or overlaps(src, receipts):
                    raise ValueError(f'External input overlaps workflow outputs: {src}')
                if not src.is_file():
                    raise FileNotFoundError(f'External input is missing: {src}')
        stage['receipt'] = receipts / f'{name}.json'
    for stage in stages.values():
        if stage['receipt'].exists():
            record = json.loads(stage['receipt'].read_text())
            for source, expected in record.get('inputs', {}).items():
                if not Path(source).is_file() or sha256(Path(source)) != expected:
                    raise ValueError(f'Completed stage input/code/config changed; preserve evidence and use a fresh run: {source}')
            current_files = set(map(str, stage['outputs']))
            for folder in stage['output_directories']:
                current_files.update(str(p.resolve()) for p in folder.rglob('*') if p.is_file())
            if current_files != set(record.get('output_stats', {})):
                raise ValueError('Completed artifact set changed; preserve evidence and use a fresh run')
            for out in map(Path, current_files):
                saved = record.get('output_stats', {}).get(str(out))
                current = {'bytes': out.stat().st_size, 'mtime_ns': out.stat().st_mtime_ns} if out.is_file() else None
                if current != saved or sha256(out) != record.get('outputs', {}).get(str(out)):
                    raise ValueError(f'Completed artifact changed or disappeared; preserve evidence and use a fresh run: {out}')
    return {'path': path, 'receipts': receipts, 'stages': stages, 'producers': producers}


def tracked_inputs(workflow, name):
    stage = workflow['stages'][name]
    # Generated artifacts are bound by producer receipts; keep DAG inputs stable
    # before and after execution. External source changes remain direct inputs.
    paths = [workflow['path'], stage['config'], stage['script'], Path(__file__).resolve()]
    paths.extend(sorted((REPO / 'scripts/research').rglob('*.py')))
    paths.extend(src for src in stage['inputs'] if src not in workflow['producers'])
    paths.extend(workflow['stages'][parent]['receipt'] for parent in stage['after'])
    return list(dict.fromkeys(map(str, paths)))


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def plan_paths(values):
    if isinstance(values, dict):
        values = list(values.values())
    paths = []
    for value in values:
        if isinstance(value, dict):
            value = value['path']
        path = Path(value).resolve()
        if path.is_dir():
            paths.extend(p.resolve() for p in path.rglob('*') if p.is_file())
        else:
            paths.append(path)
    return set(paths)


def execution_config(workflow, stage):
    if not stage['bind_inputs'] and not stage['bind_fields']:
        return stage['config']
    base = stage['config'].parent
    data = json.loads(stage['config'].read_text())
    for field in ['paths', 'path_map', 'path_roots']:
        data[field] = {key: str(resolve(base, value)) for key, value in data.get(field, {}).items()}
    data['tools'] = {key: str(resolve(base, value)) if '/' in value or '\\' in value or value.startswith('.') else value for key, value in data.get('tools', {}).items()}
    for key in stage['path_fields']:
        data[key] = str(resolve(base, data[key]))
    for row in data.get('inputs', {}).values():
        row['path'] = str(resolve(base, row['path']))
        if row.get('metadata_migration'):
            row['metadata_migration']['path'] = str(resolve(base, row['metadata_migration']['path']))
    for role, source in stage['bind_inputs'].items():
        row = data['inputs'][role]
        if Path(row['path']) != source:
            raise ValueError(f'Bound artifact differs from the configured input: {role}')
        parent = workflow['stages'][workflow['producers'][source]]
        record = json.loads(parent['receipt'].read_text())
        actual = sha256(source)
        if record['outputs'].get(str(source)) != actual:
            raise ValueError(f'Upstream completion does not bind the input: {role}')
        if row.get('sha256') not in (None, actual):
            raise ValueError(f'Configured expected input SHA differs from upstream: {role}')
        row['sha256'] = actual
    for field, source in stage['bind_fields'].items():
        if not field.endswith('_sha256') or str(resolve(base, data[field[:-7]])) != str(source):
            raise ValueError(f'Hash field must have a matching configured input path: {field}')
        parent = workflow['stages'][workflow['producers'][source]]
        record = json.loads(parent['receipt'].read_text())
        actual = sha256(source)
        if record['outputs'].get(str(source)) != actual or data.get(field) not in (None, actual):
            raise ValueError(f'Hash field differs from verified upstream output: {field}')
        data[field] = actual
    target = stage['receipt'].with_suffix('.execution.json')
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x') as handle:
        json.dump(data, handle, indent=2)
        handle.write('\n')
    return target


def execute(workflow, name, python=sys.executable, allocated_threads=None):
    stage = workflow['stages'][name]
    receipt = stage['receipt']
    if receipt.exists():
        raise FileExistsError(f'Completed evidence is read-only: {receipt}')
    for src in stage['inputs']:
        if not src.is_file():
            raise FileNotFoundError(f'Required stage input is missing: {src}')
    # Refuse silent overwrite even if an individual historical script is permissive.
    for out in stage['outputs']:
        if out.exists():
            raise FileExistsError(f'Native output already exists; use a fresh run: {out}')
    generated_config = execution_config(workflow, stage)
    command = [python, str(stage['script']), '--config', str(generated_config)]
    watched = set(map(Path, tracked_inputs(workflow, name))) | set(stage['inputs'])
    watched.add(generated_config)
    for folder in stage['generated_input_directories']:
        if not folder.is_dir():
            raise FileNotFoundError(f'Generated input directory is missing: {folder}')
        watched.update(p.resolve() for p in folder.rglob('*') if p.is_file())
    before = {str(src): sha256(src) for src in watched}
    checked = subprocess.run(command, cwd=REPO, check=True, text=True, capture_output=True)
    contract = json.loads(checked.stdout)
    required = int(contract.get('required_threads', contract.get('threads', 1)))
    available = allocated_threads if allocated_threads is not None else stage['threads']
    if available < required:
        raise ValueError(f'Configured stage needs {required} threads; only {available} allocated')
    actual = plan_paths(contract.get('inputs', [])) | plan_paths(contract.get('tools', []))
    undeclared = actual - watched
    if undeclared:
        raise ValueError('Preflight found undeclared input/tool dependencies: ' + ', '.join(map(str, sorted(undeclared))))
    for out in contract.get('outputs', []):
        path = Path(out).resolve()
        if path not in stage['outputs'] and not any(path == folder or path.is_relative_to(folder) for folder in stage['output_directories']):
            raise ValueError(f'Preflight output is outside declared artifacts: {path}')
    subprocess.run([*command, '--execute'], cwd=REPO, check=True)
    after = {str(src): sha256(src) for src in watched}
    if after != before:
        raise RuntimeError('Input, code or configuration changed during execution; no completion receipt written')
    for out in stage['outputs']:
        if not out.is_file() or out.stat().st_size == 0:
            raise RuntimeError(f'Stage exited without a required nonempty output: {out}')
    artifacts = set(stage['outputs'])
    for folder in stage['output_directories']:
        if not folder.is_dir():
            raise RuntimeError(f'Required output directory is missing: {folder}')
        artifacts.update(p.resolve() for p in folder.rglob('*') if p.is_file())
    artifacts = sorted(artifacts)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    record = {'stage': name, 'status': 'complete', 'command': [*command, '--execute'],
              'inputs': before,
              'outputs': {str(out): sha256(out) for out in artifacts},
              'output_stats': {str(out): {'bytes': out.stat().st_size, 'mtime_ns': out.stat().st_mtime_ns} for out in artifacts},
              'config_sha256': sha256(stage['config']), 'execution_config_sha256': sha256(generated_config), 'script_sha256': sha256(stage['script'])}
    with receipt.open('x') as handle:
        json.dump(record, handle, indent=2)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--stage')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--allocated-threads', type=int)
    args = parser.parse_args()
    workflow = load_manifest(args.manifest)
    if args.execute:
        if args.stage not in workflow['stages']:
            parser.error('--execute requires a known --stage')
        execute(workflow, args.stage, args.python, args.allocated_threads)
    else:
        print(json.dumps({name: {key: list(map(str, value)) if isinstance(value, list) else str(value)
                                for key, value in stage.items()}
                          for name, stage in workflow['stages'].items()}, indent=2))


if __name__ == '__main__':
    main()
