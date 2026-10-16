"""Artificial orchestration checks. No scientific scripts or inputs are used."""
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('portable_stages', ROOT / 'workflow/portable/stages.py')
stages = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stages)


class PortableWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='portable workflow ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'scripts/research').mkdir(parents=True)
        (self.root / 'workflow/portable').mkdir(parents=True)
        for name in ['stages.py', 'Snakefile']:
            shutil.copyfile(ROOT / 'workflow/portable' / name, self.root / 'workflow/portable' / name)
        self.fake = self.root / 'scripts/research/fake.py'
        self.fake.write_text('''import argparse,json,pathlib,sys
p=argparse.ArgumentParser();p.add_argument('--config');p.add_argument('--execute',action='store_true');a=p.parse_args()
c=json.loads(pathlib.Path(a.config).read_text())
if a.execute:
 out=pathlib.Path(c['out']);out.parent.mkdir(parents=True,exist_ok=True);out.write_text('fixture')
 if c.get('mutate'): pathlib.Path(c['mutate']).write_text('mutated')
 if c.get('fail'): sys.exit(7)
else:
 print(json.dumps({'inputs':c.get('inputs',[]),'outputs':[c['out']],'required_threads':c.get('threads',1)}))
''')
        self.source = self.root / 'input.txt'
        self.source.write_text('source fixture')
        for name in ['a', 'b']:
            (self.root / f'{name}.json').write_text(json.dumps({'out': str(self.root / f'outputs/{name}.txt')}))
        self.data = {'receipts': 'receipts', 'frozen_dirs': ['frozen'], 'stages': {
            'first': {'script': 'scripts/research/fake.py', 'config': 'a.json', 'inputs': ['input.txt'],
                      'outputs': ['outputs/a.txt']},
            'second': {'script': 'scripts/research/fake.py', 'config': 'b.json', 'inputs': ['outputs/a.txt'],
                       'outputs': ['outputs/b.txt'], 'after': ['first']}}}
        self.path = self.root / 'workflow.json'
        self.override = patch.object(stages, 'REPO', self.root)
        self.override.start()
        self.addCleanup(self.override.stop)

    def load(self):
        self.path.write_text(json.dumps(self.data))
        return stages.load_manifest(self.path)

    def test_read_only_plan_with_future_inputs(self):
        workflow = self.load()
        inputs = stages.tracked_inputs(workflow, 'second')
        self.assertIn(str(self.root / 'receipts/first.json'), inputs)
        self.assertFalse((self.root / 'outputs').exists())
        self.assertFalse((self.root / 'receipts').exists())

    def test_missing_dependency_edge_is_rejected(self):
        self.data['stages']['second']['after'] = []
        with self.assertRaisesRegex(ValueError, 'without an upstream'):
            self.load()

    def test_cycles_are_rejected(self):
        self.data['stages']['first']['after'] = ['second']
        with self.assertRaisesRegex(ValueError, 'Cyclic'):
            self.load()

    def test_missing_external_input_is_rejected(self):
        self.source.unlink()
        with self.assertRaises(FileNotFoundError):
            self.load()

    def test_alias_to_frozen_output_is_rejected(self):
        (self.root / 'frozen').mkdir()
        (self.root / 'alias').symlink_to(self.root / 'frozen', target_is_directory=True)
        self.data['stages']['first']['outputs'] = ['alias/a.txt']
        with self.assertRaisesRegex(ValueError, 'frozen'):
            self.load()

    def test_failure_keeps_partial_artifact_without_completion(self):
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'outputs/a.txt'), 'fail': True}))
        workflow = self.load()
        with self.assertRaises(subprocess.CalledProcessError):
            stages.execute(workflow, 'first')
        self.assertTrue((self.root / 'outputs/a.txt').is_file())
        self.assertFalse((self.root / 'receipts/first.json').exists())
        with self.assertRaises(FileExistsError):
            stages.execute(workflow, 'first')

    def test_exit_zero_missing_output_is_not_completion(self):
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'wrong.txt')}))
        with self.assertRaisesRegex(ValueError, 'outside declared artifacts'):
            stages.execute(self.load(), 'first')
        self.assertFalse((self.root / 'receipts/first.json').exists())

    def test_connected_fake_stages_hash_receipts_and_refuse_overwrite(self):
        workflow = self.load()
        stages.execute(workflow, 'first')
        stages.execute(workflow, 'second')
        record = json.loads((self.root / 'receipts/second.json').read_text())
        self.assertEqual(record['status'], 'complete')
        self.assertIn(str(self.root / 'outputs/a.txt'), record['inputs'])
        with self.assertRaises(FileExistsError):
            stages.execute(workflow, 'second')

    def test_completed_artifact_change_is_rejected_on_next_plan(self):
        workflow = self.load()
        stages.execute(workflow, 'first')
        (self.root / 'outputs/a.txt').write_text('changed artificial artifact')
        with self.assertRaisesRegex(ValueError, 'changed or disappeared'):
            self.load()

    def test_same_size_mtime_byte_change_is_rejected(self):
        workflow = self.load()
        stages.execute(workflow, 'first')
        out = self.root / 'outputs/a.txt'
        stat = out.stat()
        out.write_text('changed')
        os.utime(out, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        with self.assertRaisesRegex(ValueError, 'changed or disappeared'):
            self.load()

    def test_undeclared_preflight_dependency_is_rejected(self):
        hidden = self.root / 'hidden.txt'
        hidden.write_text('hidden input')
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'outputs/a.txt'), 'inputs': [str(hidden)]}))
        with self.assertRaisesRegex(ValueError, 'undeclared input'):
            stages.execute(self.load(), 'first')
        self.assertFalse((self.root / 'outputs').exists())

    def test_insufficient_declared_cpus_are_rejected(self):
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'outputs/a.txt'), 'threads': 4}))
        with self.assertRaisesRegex(ValueError, 'needs 4 threads'):
            stages.execute(self.load(), 'first')
        self.assertFalse((self.root / 'outputs').exists())

    def test_input_mutation_is_not_sealed_as_success(self):
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'outputs/a.txt'), 'mutate': str(self.source)}))
        with self.assertRaisesRegex(RuntimeError, 'changed during execution'):
            stages.execute(self.load(), 'first')
        self.assertTrue((self.root / 'outputs/a.txt').exists())
        self.assertFalse((self.root / 'receipts/first.json').exists())

    def test_generated_directory_has_upstream_owner_and_is_hashed(self):
        self.data['stages']['first']['output_directories'] = ['outputs/first']
        self.data['stages']['first']['outputs'] = ['outputs/first/a.txt']
        self.data['stages']['second']['inputs'] = ['outputs/first/a.txt']
        self.data['stages']['second']['generated_input_directories'] = ['outputs/first']
        (self.root / 'a.json').write_text(json.dumps({'out': str(self.root / 'outputs/first/a.txt')}))
        (self.root / 'b.json').write_text(json.dumps({'out': str(self.root / 'outputs/b.txt'), 'inputs': [str(self.root / 'outputs/first')]}))
        workflow = self.load()
        stages.execute(workflow, 'first')
        stages.execute(workflow, 'second')
        saved = json.loads((self.root / 'receipts/second.json').read_text())
        self.assertIn(str(self.root / 'outputs/first/a.txt'), saved['inputs'])

    def test_generated_directory_without_owner_is_rejected(self):
        self.data['stages']['second']['generated_input_directories'] = ['unknown']
        with self.assertRaisesRegex(ValueError, 'upstream output owner'):
            self.load()

    def test_upstream_hash_binding_preserves_original_config(self):
        self.data['stages']['second']['bind_inputs'] = {'prior': 'outputs/a.txt'}
        original = {'inputs': {'prior': {'path': 'outputs/a.txt'}}, 'paths': {'root': 'outputs'}, 'tools': {}}
        (self.root / 'b.json').write_text(json.dumps(original))
        workflow = self.load()
        stages.execute(workflow, 'first')
        generated = stages.execution_config(workflow, workflow['stages']['second'])
        bound = json.loads(generated.read_text())
        self.assertEqual(bound['inputs']['prior']['sha256'], stages.sha256(self.root / 'outputs/a.txt'))
        self.assertEqual(bound['paths']['root'], str(self.root / 'outputs'))
        self.assertEqual(json.loads((self.root / 'b.json').read_text()), original)

    def test_existing_unlisted_sidecar_requires_explicit_owner(self):
        (self.root / 'outputs/first').mkdir(parents=True)
        (self.root / 'outputs/first/sidecar.txt').write_text('sidecar')
        self.data['stages']['first']['output_directories'] = ['outputs/first']
        self.data['stages']['second']['inputs'] = ['outputs/first/sidecar.txt']
        self.data['stages']['second']['after'] = []
        with self.assertRaisesRegex(ValueError, 'unlisted generated artifact'):
            self.load()

    def test_top_level_hash_field_binds_only_matching_upstream_path(self):
        self.data['stages']['second']['bind_fields'] = {'reference_tree_sha256': 'outputs/a.txt'}
        self.data['stages']['second']['path_fields'] = ['reference_tree', 'outdir']
        original = {'reference_tree': 'outputs/a.txt', 'reference_tree_sha256': None, 'outdir': 'outputs/loci'}
        (self.root / 'b.json').write_text(json.dumps(original))
        workflow = self.load()
        stages.execute(workflow, 'first')
        generated = stages.execution_config(workflow, workflow['stages']['second'])
        bound = json.loads(generated.read_text())
        self.assertEqual(bound['reference_tree_sha256'], stages.sha256(self.root / 'outputs/a.txt'))
        self.assertEqual(bound['reference_tree'], str(self.root / 'outputs/a.txt'))
        self.assertEqual(json.loads((self.root / 'b.json').read_text()), original)

    @unittest.skipUnless(shutil.which('snakemake'), 'Snakemake is not installed')
    def test_snakemake_plans_and_runs_only_fake_stages(self):
        self.load()
        command = ['snakemake', '-s', 'workflow/portable/Snakefile', '--cores', '1', '--config', f'manifest={self.path}']
        result = subprocess.run([*command, '-n'], cwd=self.root, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((self.root / 'outputs').exists())
        result = subprocess.run(command, cwd=self.root, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.root / 'receipts/second.json').is_file())
        result = subprocess.run([*command, '-n'], cwd=self.root, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('Nothing to be done', result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
