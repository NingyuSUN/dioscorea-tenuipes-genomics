#!/usr/bin/env python3
"""Freeze configured code outside result directories; no analysis or dispatch.

The configuration names its manifest path but never contains that manifest's hash,
so the manifest can safely bind the configuration without a circular dependency.
"""
import argparse
import hashlib
import json
from pathlib import Path
import portable_cli as portable


def main():
    p=portable.parser(__doc__)
    p.add_argument('--stage', choices=['run_cp97','run_expanded_trees'], required=True)
    a=p.parse_args()
    if a.stage=='run_cp97':
        import run_cp97 as stage
        c=stage.configure(a.config);out=c['program_freeze']
        paths=[c['_config'],Path(stage.__file__).resolve(),Path(portable.__file__).resolve(),Path(c['guard']['path'])]
    else:
        import run_expanded_trees as stage
        c=stage.configure(a.config);out=c['frozen_manifest']
        paths=[c['_config'],Path(stage.__file__).resolve(),Path(portable.__file__).resolve(),
               c['entry_script'],c['entry_config'],c['legacy_runner']]
    if out.is_relative_to(c['out']):
        raise ValueError('The program freeze must live outside the fresh result directory')
    if out in paths or out.exists():
        raise ValueError('Require a fresh manifest outside the files being frozen')
    for path in paths:
        if not path.is_file():raise ValueError('Missing file to freeze: '+str(path))
    planned=portable.contract(c['_config'],{str(i):f for i,f in enumerate(paths)},
                              {'freeze':out},output_root=None)
    if not a.execute:portable.emit(planned);return
    if not __debug__:raise ValueError('Python optimization disables scientific assertions')
    records={str(path):{'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths}
    manifest={'status':'CODE_AND_CONFIGURATION_FROZEN_NOT_SCIENTIFIC_VALIDATION',
              'files':records,'sha256':{name:item['sha256'] for name,item in records.items()}}
    with out.open('x') as handle:json.dump(manifest,handle,indent=2);handle.write('\n')
    portable.emit({'status':'FROZEN_NO_ANALYSIS','manifest':str(out),'file_count':len(records)})


if __name__=='__main__':main()
