"""Synthetic portability contracts; no historical analysis or real inference runs."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/'scripts/research/chloroplast/historical'
NAMES=['prepare_inputs','prepare_source_inputs','audit_hachi_nuclear_coverage',
       'compare_published_edo','published_edo_controls','validate_denovo_cp_reads',
       'run_expanded_trees','run_cp97']
sys.path.insert(0,str(SCRIPTS))
import portable_cli
import run_cp97
import run_expanded_trees
import audit_hachi_nuclear_coverage
import cp_tree_runner


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def module(name):
    spec=importlib.util.spec_from_file_location('portable_cp_test_'+name,SCRIPTS/(name+'.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def cp97_resources():
    return {'cpus':[0,1,2,3],'memory_bytes':6442450944,'swap_bytes':0,
            'minimum_mem_available_bytes':17179869184,'soft_term_seconds':14385,'hard_seconds':14400,
            'startup_cpu_sample_seconds':2,'maximum_mean_cpu_busy_fraction':0.2}


def expanded_budget():
    return {'tree_cpu_ids':[4,5,6,7],'RNA_planned_cpu_ids':[0,1,2,3],
            'tree_memory_max_bytes':6*1024**3,'rna_coordinator_memory_max_bytes':2*1024**3,
            'rna_stage_max_bytes':16*1024**3,'combined_planned_caps_bytes':24*1024**3,
            'startup_headroom_bytes':2*1024**3,'anonymous_memory_credit_bytes':0,
            'minimum_MemTotal_bytes':30*1024**3,'minimum_MemAvailable_bytes':26*1024**3}


def fixture(directory,name):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    tool=directory/'fake tool';tool.write_text('#!/bin/sh\nexit 97\n');tool.chmod(0o755)
    c={'out':'results/'+name}
    if name=='prepare_inputs':c.update(source_dir='source',p19_dir='p19',iqtree='fake tool')
    elif name=='prepare_source_inputs':
        c.update(design_dir='design',source_dir='source',historical_dir='history',iqtree='fake tool')
        (directory/'design').mkdir();(directory/'design/INPUT_METADATA.json').write_text('{"source_file_hashes":{}}')
    elif name=='audit_hachi_nuclear_coverage':
        c.update(binding='binding.json',genome='genome.fa.gz',cp='cp.fa.gz',reads=['R1.fq.gz','R2.fq.gz'],minimap2='fake tool',samtools='fake tool')
    elif name=='compare_published_edo':c.update(genotypes='published.gt.gz',cp_snps='cp.fa')
    elif name=='published_edo_controls':c.update(genotypes='published.gt.gz',comparison_receipt='comparison/RECEIPT.json')
    elif name=='validate_denovo_cp_reads':
        c.update(contigs='contigs.fa',assembly_receipt='assembly.json',screen='screen.json',reads=['R1.fq.gz','R2.fq.gz'],minimap2='fake tool',samtools='fake tool')
        (directory/'screen.json').write_text('{"alignments":[{"contig":"ctg1","qualifying_cp_like":true}]}')
    elif name=='run_cp97':
        c.update(program_freeze='freeze.json',resources=cp97_resources(),python='fake tool',systemd_run='fake tool',systemctl='fake tool',taskset='fake tool',timeout='fake tool',
                 inputs={key:{'path':key+'.data','sha256':'0'*64} for key in ['main98','control97','models','old_runtime','old_tree','entry_receipt']},
                 executable={'path':'fake tool','sha256':digest(tool)},guard={'path':str(SCRIPTS/'exec_verified_scope.py'),'sha256':digest(SCRIPTS/'exec_verified_scope.py')},
                 original_argv=['old-iqtree','-s','old98.fa','-p','old.nex','-T','4','-o','a,b','--prefix','old/main'],columns=124448,expected_labels98=[],deleted_label='Hachi',alphabet='ACGTN-?')
    elif name=='run_expanded_trees':
        c.update(entry_config='entry.json',sample_parent='samples',frozen_manifest='freeze.json',python='fake tool',loginctl='fake tool',tree_lock='tree.lock',
                 policy={'iqtree':'fake tool','iqtree_sha256':digest(tool),'minimum_free_disk_bytes':2*1024**3},
                 required_receipts={'recovery03':'recovery.json'},fixed_model={'fixed_model_file_remote':'fixed_models.nex','fixed_model_file_sha256':'0'*64},
                 resource_budget=expanded_budget(),latest_new_stage_start='2099-01-01T00:00:00+00:00',path_map={})
        ec={'source_profiles':{'local':{'original93':'original.fa','kita':'kita.fa','hachi':'hachi.fa','kita_receipt':'kita.json','hachi_receipt':'hachi.json','parent_ledger':'ledger'}},
            'source_hashes':{x:'0'*64 for x in ['original93','kita','hachi','kita_receipt','hachi_receipt']},'geometry_hashes':{'regions.nex':'0'*64},
            'default_sample_parent':'historical-samples','samples':[{'run':'DRR351754','name':'Mie1'}]}
        (directory/'entry.json').write_text(json.dumps(ec))
    config=directory/'config.json';config.write_text(json.dumps(c))
    stage=module(name);plan=stage.plan(config)
    for path in plan['inputs']:
        path=Path(path)
        if not path.exists():path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{}\n')
    # A valid empty ledger source map is sufficient for read-only dependency discovery.
    if name=='run_expanded_trees':
        (directory/'samples/DRR351754/fixed_ledger_v2/CONFIG.json').write_text('{"sources":{}}')
    return config,stage


class PortableChloroplast(unittest.TestCase):
    def call(self,name,*args,cwd=None,optimize=False):
        cmd=[sys.executable]+(['-O'] if optimize else [])+[str(SCRIPTS/(name+'.py')),*map(str,args)]
        env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1'}
        return subprocess.run(cmd,cwd=cwd,text=True,capture_output=True,env=env,timeout=15)

    def test_all_help_without_optional_packages_or_writes(self):
        with tempfile.TemporaryDirectory() as td:
            for name in NAMES:
                r=subprocess.run([sys.executable,'-S',str(SCRIPTS/(name+'.py')),'--help'],cwd=td,text=True,capture_output=True,timeout=10)
                self.assertEqual(r.returncode,0,(name,r.stderr));self.assertIn('--execute',r.stdout)
            self.assertEqual(list(Path(td).iterdir()),[])

    def test_all_plans_are_config_relative_and_read_only(self):
        with tempfile.TemporaryDirectory() as td:
            for name in NAMES:
                config,stage=fixture(Path(td)/name,name)
                before={str(p):digest(p) for p in config.parent.rglob('*') if p.is_file()}
                result=self.call(name,'--config',config,cwd=td)
                self.assertEqual(result.returncode,0,(name,result.stderr))
                plan=json.loads(result.stdout);self.assertEqual(plan['status'],'PREFLIGHT_ONLY_NOT_EXECUTED')
                self.assertTrue(all(Path(x).is_absolute() for x in plan['inputs']+plan['outputs']+plan['tools']))
                self.assertFalse(Path(plan['output_root']).exists())
                self.assertEqual(before,{str(p):digest(p) for p in config.parent.rglob('*') if p.is_file()})
                self.assertEqual(stage.plan(config)['outputs'],plan['outputs'])

    def test_missing_input_and_existing_output_fail_without_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            config,stage=fixture(td,'compare_published_edo');plan=stage.plan(config)
            missing=Path(plan['named_inputs']['cp_snps']);missing.unlink()
            result=self.call('compare_published_edo','--config',config)
            self.assertNotEqual(result.returncode,0);self.assertIn('Missing input',result.stderr)
            missing.write_text('synthetic\n');out=Path(plan['output_root']);out.mkdir(parents=True);sentinel=out/'partial.txt';sentinel.write_text('keep')
            result=self.call('compare_published_edo','--config',config,'--execute')
            self.assertNotEqual(result.returncode,0);self.assertIn('Fresh output',result.stderr);self.assertEqual(sentinel.read_text(),'keep')

    def test_optimized_python_cannot_disable_scientific_assertions(self):
        with tempfile.TemporaryDirectory() as td:
            config,_=fixture(td,'prepare_inputs')
            result=self.call('prepare_inputs','--config',config,optimize=True)
            self.assertNotEqual(result.returncode,0);self.assertIn('optimization disables',result.stderr)

    def test_fixed_hash_mismatch_blocks_preparation_before_writing(self):
        with tempfile.TemporaryDirectory() as td:
            config,stage=fixture(td,'prepare_inputs')
            result=self.call('prepare_inputs','--config',config,'--execute')
            self.assertNotEqual(result.returncode,0);self.assertIn('AssertionError',result.stderr)
            self.assertFalse(Path(stage.plan(config)['output_root']).exists())

    def test_program_freeze_generated_once_and_detects_config_change(self):
        with tempfile.TemporaryDirectory() as td:
            config,stage=fixture(td,'run_cp97');freeze=config.parent/'freeze.json';freeze.unlink()
            planned=self.call('freeze_contract','--stage','run_cp97','--config',config)
            self.assertEqual(planned.returncode,0,planned.stderr);self.assertFalse(freeze.exists())
            result=self.call('freeze_contract','--stage','run_cp97','--config',config,'--execute')
            self.assertEqual(result.returncode,0,result.stderr)
            run_cp97.CONFIG_PATH=config;run_cp97.ROOT=config.parent/'results/run_cp97'
            self.assertEqual(run_cp97.verify_bundle()['resources'],cp97_resources())
            c=json.loads(config.read_text());c['relocation_note']='changed after freezing';config.write_text(json.dumps(c))
            with self.assertRaisesRegex(ValueError,'freeze mismatch'):run_cp97.verify_bundle()
            again=self.call('freeze_contract','--stage','run_cp97','--config',config,'--execute')
            self.assertNotEqual(again.returncode,0)

    def test_exact_deletion_rejects_changed_retained_sequence(self):
        raw=b'>a\nACG\n>drop\nNNN\n>b\nA-T\n';deleted=b'>a\nACG\n>b\nA-T\n'
        proof=run_cp97.validate_deletion(raw,deleted,3,['a','drop','b'],'drop','ACGTN-')
        self.assertTrue(proof['record_bytes_equal'])
        with self.assertRaisesRegex(ValueError,'byte-exact deletion'):
            run_cp97.validate_deletion(raw,deleted.replace(b'A-T',b'ACT'),3,['a','drop','b'],'drop','ACGTN-')

    def test_resource_guards_reject_weakened_contracts(self):
        with tempfile.TemporaryDirectory() as td:
            mem=Path(td)/'meminfo';mem.write_text('MemTotal: 33554432 kB\nMemAvailable: 1024 kB\n')
            c={'resource_budget':expanded_budget()}
            check=run_expanded_trees.resource_budget_snapshot(c,mem)
            self.assertFalse(check['passes'])
            c['resource_budget']['anonymous_memory_credit_bytes']=1
            with self.assertRaisesRegex(ValueError,'zero anonymous'):run_expanded_trees.resource_budget_snapshot(c,mem)
        with patch.object(run_cp97.os,'sched_getaffinity',return_value={8,9,10,11}):
            with self.assertRaisesRegex(ValueError,'not allowed'):run_cp97.resource_preflight({'resources':cp97_resources()})

    def test_nested_dependencies_include_new_and_relocated_files(self):
        with tempfile.TemporaryDirectory() as td:
            config,stage=fixture(td,'run_expanded_trees');base=config.parent
            nested=base/'samples/DRR351754/fixed_projection/new.evidence.tsv';nested.write_text('position\tbase\n1\tA\n')
            external=base/'relocated/reference.fa';external.parent.mkdir();external.write_text('>ref\nA\n')
            ledger=base/'samples/DRR351754/fixed_ledger_v2/CONFIG.json'
            ledger.write_text(json.dumps({'sources':{'external':{'path':'old-reference/reference.fa','sha256':digest(external)}}}))
            c=json.loads(config.read_text());c['path_map']={'old-reference':'relocated'};config.write_text(json.dumps(c))
            plan=stage.plan(config)
            self.assertIn(str(nested),plan['inputs']);self.assertIn(str(external),plan['inputs'])
            old=digest(nested);nested.write_text('position\tbase\n1\tC\n');self.assertNotEqual(old,digest(nested))
            self.assertIn(str(nested),stage.plan(config)['inputs'])

    def test_fake_iqtree_failure_preserves_partial_outputs_and_nonzero(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);tool=root/'fake_iqtree.py'
            tool.write_text('#!'+sys.executable+'\nimport pathlib,sys\np=pathlib.Path(sys.argv[sys.argv.index("--prefix")+1]);p.with_suffix(".checkpoint").write_text("partial")\nsys.exit(17)\n');tool.chmod(0o755)
            inputs=root/'inputs';inputs.mkdir();out=root/'trees';out.mkdir()
            c={'iqtree':str(tool),'model_file':'fixed.nex','threads':4,'seed':20260930,'outgroups':['a','b'],'full_bootstrap':1000,'full_alrt':1000}
            with self.assertRaisesRegex(RuntimeError,'exit 17'):
                cp_tree_runner.run_job('main98',inputs,out,'full',c,{'panels':{'main98':{'tips':['a','b']}}},None)
            self.assertEqual((out/'main98/main.checkpoint').read_text(),'partial')
            self.assertEqual(json.loads((out/'main98/RUNTIME.json').read_text())['status'],'FAILED_NO_AUTO_RESTART')
            with self.assertRaises(FileExistsError):cp_tree_runner.run_job('main98',inputs,out,'full',c,{},None)

    def test_template_depth_collapses_mates_and_excludes_disagreement(self):
        observations=[(('rg','one'),'A'),(('rg','one'),'A'),(('rg','two'),'A'),(('rg','two'),'C'),(('rg','three'),'N')]
        self.assertEqual(audit_hachi_nuclear_coverage.collapse(observations),1)


if __name__=='__main__':unittest.main()
