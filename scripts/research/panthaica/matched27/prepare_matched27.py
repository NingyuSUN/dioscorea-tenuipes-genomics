#!/usr/bin/env python3
"""Prepare matched panthaica controls using the preserved historical calculation functions."""
import argparse
import hashlib
import json
from pathlib import Path
from run_matched27 import fasta, project
from sensitivity.run_sensitivity import prepare_inputs


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_fasta(path, rows):
    with path.open('x') as handle:
        for label, sequence in rows.items():
            handle.write('>' + label + '\n' + sequence + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['project', 'sensitivity'], required=True)
    parser.add_argument('--with-alignment', type=Path, required=True)
    parser.add_argument('--without-alignment', type=Path, required=True)
    parser.add_argument('--partitions', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True,
                        help='JSON containing the expected taxon labels/order and panthaica label.')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    with_rows, without_rows = fasta(args.with_alignment), fasta(args.without_alignment)
    partitions = json.loads(args.partitions.read_text())
    if not isinstance(partitions, list):
        parser.error('--partitions must contain a JSON list of the 27 fixed gene intervals')
    config = json.loads(args.config.read_text())
    if args.mode == 'project':
        matrices, proof = project(with_rows, without_rows, config, partitions)
        extra = []
    else:
        matrices, extra, proof = prepare_inputs(with_rows, without_rows, partitions, config)
    # Validate the biological projections before creating the output directory.
    args.out.mkdir(parents=True, exist_ok=False)
    for name, rows in matrices.items():
        write_fasta(args.out / (name + '.faa'), rows)
    if extra:
        (args.out / 'loo').mkdir()
        for row in extra:
            write_fasta(args.out / 'loo' / (row['job_id'] + '.faa'), row['rows'])
        (args.out / 'LOO_JOBS.json').write_text(json.dumps(
            [{k: v for k, v in row.items() if k != 'rows'} for row in extra], indent=2) + '\n')
    partition_name = 'PARTITIONS.json' if args.mode == 'project' else 'SOURCE_PARTITIONS.json'
    (args.out / partition_name).write_text(json.dumps(partitions, indent=2) + '\n')
    (args.out / 'PROJECTION_PROOF.json').write_text(json.dumps(proof, indent=2) + '\n')
    receipt = {
        'status': 'PREPARED_NO_TREE_INFERENCE', 'mode': args.mode,
        'inputs_sha256': {name: sha(getattr(args, name)) for name in
                         ['with_alignment', 'without_alignment', 'partitions', 'config']},
        'scientific_source_sha256': {
            'run_matched27.py': sha(Path(__file__).with_name('run_matched27.py')),
            'sensitivity/run_sensitivity.py': sha(Path(__file__).parent / 'sensitivity/run_sensitivity.py')},
        'new_tree_inferences': 0,
        'outputs_sha256': {str(path.relative_to(args.out)): sha(path)
                           for path in sorted(args.out.rglob('*')) if path.is_file()}}
    (args.out / 'PREPARATION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
