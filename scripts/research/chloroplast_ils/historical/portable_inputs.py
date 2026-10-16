"""Explicit configuration and read-only preflight for archived scientific scripts."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil

class Config:
    def __init__(self, filename):
        if not __debug__:
            raise RuntimeError('Python optimization disables required scientific assertions; run without -O.')
        self.filename = Path(filename).resolve()
        self.base = self.filename.parent
        self.data = json.loads(self.filename.read_text())
        if not isinstance(self.data, dict):
            raise ValueError('Configuration must be a JSON object')
        for name in ('paths', 'inputs', 'tools', 'path_map', 'path_roots'):
            if not isinstance(self.data.get(name, {}), dict):
                raise ValueError(name + ' must be an object')
    def local(self, value):
        if not isinstance(value, (str, Path)) or not str(value):
            raise ValueError('A nonempty configured path is required')
        p = Path(value)
        return (self.base / p).resolve() if not p.is_absolute() else p.resolve()
    def resolve(self, value):
        """Relocate an evidence path only through explicit exact/prefix mappings."""
        value = str(value)
        mapping = self.data.get('path_map', {})
        if value in mapping:
            return self.local(mapping[value])
        candidates = []
        for old, new in self.data.get('path_roots', {}).items():
            old = old.rstrip('/')
            if value == old or value.startswith(old + '/'):
                candidates.append((len(old), self.local(new) / value[len(old):].lstrip('/')))
        if candidates:
            return max(candidates, key=lambda x:x[0])[1].resolve()
        return self.local(value)
    def path(self, key, default=None):
        value = self.data.get('paths', {}).get(key, default)
        if value is None:
            raise ValueError('Missing paths.' + key)
        return self.local(value)
    def input(self, key, expected=None, default=None):
        row = self.data.get('inputs', {}).get(key)
        if row is None:
            if default is None:raise ValueError('Missing inputs.' + key)
            row = {'path':str(default)}
        if not isinstance(row, dict) or 'path' not in row:
            raise ValueError('inputs.' + key + ' needs a path')
        h = row.get('sha256', expected)
        migration=None
        if expected is not None and h != expected:
            if row.get('original_sha256') != expected:
                raise ValueError('Changed frozen binding requires original_sha256 for inputs.' + key)
            migration=row.get('metadata_migration')
            if not isinstance(migration,dict) or not re.fullmatch('[0-9a-f]{64}',str(migration.get('sha256',''))):
                raise ValueError('Changed frozen binding requires a hash-bound metadata migration receipt')
            receipt_path=self.local(migration.get('path'))
            if digest(receipt_path)!=migration['sha256']:
                raise ValueError('Metadata migration receipt SHA-256 mismatch')
            receipt=json.loads(receipt_path.read_text())
            migration={'path':str(receipt_path),'sha256':migration['sha256']}
            matches=[x for x in receipt.get('bindings',[]) if x.get('original_sha256')==expected and x.get('current_sha256')==h]
            if receipt.get('schema')!='reviewed-metadata-migration.v1' or len(matches)!=1 or matches[0].get('metadata_only') is not True or matches[0].get('scientific_scalars_and_structures_unchanged') is not True:
                raise ValueError('No reviewed metadata-only binding matches this original/current SHA-256 pair')
        if h is not None and not re.fullmatch('[0-9a-f]{64}', h):
            raise ValueError('Invalid SHA-256 for inputs.' + key)
        return {'role':key, 'path':str(self.local(row['path'])), 'sha256':h, 'original_sha256':expected, 'metadata_migration':migration}
    def tool(self, key, default=None):
        value = self.data.get('tools', {}).get(key, default)
        if not isinstance(value,str) or not value:raise ValueError('Missing tools.' + key)
        if '/' in value or '\\' in value or value.startswith('.'):
            return str(self.local(value))
        return shutil.which(value) or value
    def relocated(self, value):
        """Resolve typed path fields without changing immutable evidence files."""
        if isinstance(value,list):return [self.relocated(x) for x in value]
        if isinstance(value,dict):
            return {k:(str(self.resolve(v)) if isinstance(v,str) and (k=='path' or k.endswith('_path') or k in ('baseline_root','resolved_path')) else self.relocated(v)) for k,v in value.items()}
        return value

def run_directory(config, root):
    if 'run_dir' not in config.data.get('paths',{}):return None
    result=config.path('run_dir')
    if result==root or not result.is_relative_to(root):
        raise ValueError('paths.run_dir must be a child of paths.root')
    return result

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def preflight(plan):
    """Check local files/hashes/dependencies; never execute tools or create files."""
    issues=[]
    for row in plan['inputs']:
        p=Path(row['path'])
        if not p.is_file():issues.append(row['role']+': missing file');continue
        if row.get('sha256') and digest(p)!=row['sha256']:issues.append(row['role']+': SHA-256 mismatch')
    for key,value in plan['tools'].items():
        p=Path(value)
        if not (p.is_file() and (key in plan.get('non_executable_tools',[]) or os.access(p,os.X_OK))):issues.append(key+': executable not found')
    for package in plan.get('python_packages',[]):
        if importlib.util.find_spec(package) is None:issues.append(package+': Python package not found')
    output_paths=[Path(p) for p in plan['outputs']]
    input_paths={Path(r['path']).resolve() for r in plan['inputs']}
    if any(p.resolve() in input_paths for p in output_paths):issues.append('Output overlaps an input file')
    result={**plan,'status':'PREFLIGHT_FAILED' if issues else ('PREFLIGHT_PASSED_DECLARED_INPUTS_ONLY' if plan.get('transitive_inputs') else 'PREFLIGHT_PASSED_NO_EXECUTION'),'issues':issues,'scientific_validation':'NOT_RUN'}
    return result

def report(plan):
    result=preflight(plan)
    print(json.dumps(result,indent=2))
    return 2 if result['issues'] else 0

def require_preflight(plan):
    result=preflight(plan)
    if result['issues']:raise ValueError('; '.join(result['issues']))
    return result
