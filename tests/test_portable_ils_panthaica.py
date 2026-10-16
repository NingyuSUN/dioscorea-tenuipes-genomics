"""Synthetic CLI, relocation, and guard tests; no biological analysis or simulation."""
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=[
 'panthaica/historical/run_rna_pipeline_draft.py',
 'panthaica/historical/run_marker_recovery_draft.py',
 'panthaica/historical/prepare_alignment_controls.py',
 'nuclear_ils/historical/fit_history.py',
 'nuclear_ils/historical/audit_panel.py',
 'chloroplast_ils/historical/prepare_inputs.py',
 'chloroplast_ils/historical/make_report_figures.py']

def load(relative):
 p=ROOT/'scripts/research'/relative
 name='portable_test_'+p.stem
 spec=importlib.util.spec_from_file_location(name,p)
 module=importlib.util.module_from_spec(spec)
 sys.path.insert(0,str(p.parent))
 try:spec.loader.exec_module(module)
 finally:sys.path.pop(0)
 return module

class PortableTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.base=Path(self.tmp.name)
 def tearDown(self):self.tmp.cleanup()
 def config(self,data):
  p=self.base/'nested/config.json';p.parent.mkdir(exist_ok=True);p.write_text(json.dumps(data));return p
 def test_help_never_needs_data_or_writes(self):
  for name in SCRIPTS:
   with self.subTest(script=name):
    result=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/research'/name),'--help'],cwd=self.base,capture_output=True,text=True)
    self.assertEqual(result.returncode,0,result.stderr);self.assertIn('--config',result.stdout);self.assertIn('--execute',result.stdout)
    self.assertEqual(list(self.base.iterdir()),[])
 def test_import_has_no_file_writes_or_reads_of_data(self):
  for name in SCRIPTS:
   with self.subTest(script=name),patch.object(Path,'write_text',side_effect=AssertionError('import write')),patch.object(Path,'mkdir',side_effect=AssertionError('import mkdir')),patch.object(Path,'read_text',side_effect=AssertionError('import data read')):
    load(name)
 def test_optimized_python_rejects_configuration_but_help_stays_safe(self):
  cfg=self.config({});before=set(self.base.rglob('*'))
  for name in SCRIPTS:
   with self.subTest(script=name):
    command=[sys.executable,'-B','-O',str(ROOT/'scripts/research'/name)]
    help_result=subprocess.run(command+['--help'],cwd=self.base,capture_output=True,text=True)
    self.assertEqual(help_result.returncode,0,help_result.stderr)
    result=subprocess.run(command+['--config',str(cfg)],cwd=self.base,capture_output=True,text=True)
    self.assertNotEqual(result.returncode,0)
    self.assertIn('optimization disables required scientific assertions',result.stderr)
    self.assertEqual(set(self.base.rglob('*')),before)
 def test_explicit_relocation_uses_config_parent(self):
  m=load(SCRIPTS[0]);cfg=self.config({'path_map':{'old/exact':'../exact'},'path_roots':{'old/tree':'data'}});c=m.Config(cfg)
  self.assertEqual(c.resolve('old/exact'),self.base/'exact')
  self.assertEqual(c.resolve('old/tree/a.fa'),cfg.parent/'data/a.fa')
  self.assertEqual(c.resolve('unmapped.txt'),cfg.parent/'unmapped.txt')
  self.assertEqual(c.resolve('old/trees/a.fa'),cfg.parent/'old/trees/a.fa')
 def test_frozen_hash_cannot_be_overridden(self):
  m=load(SCRIPTS[0]);c=m.Config(self.config({'inputs':{'x':{'path':'x','sha256':'0'*64}}}))
  with self.assertRaisesRegex(ValueError,'frozen binding'):c.input('x','1'*64)
 def test_explicit_metadata_migration_binds_original_and_current_hashes(self):
  m=load(SCRIPTS[0]);payload=self.base/'metadata.json';payload.write_text('{"path":"portable/data"}')
  current=hashlib.sha256(payload.read_bytes()).hexdigest();old='1'*64
  receipt=self.base/'migration.json';receipt.write_text(json.dumps({'schema':'reviewed-metadata-migration.v1','bindings':[{'original_sha256':old,'current_sha256':current,'metadata_only':True,'scientific_scalars_and_structures_unchanged':True}]}))
  cfg=self.config({'inputs':{'metadata':{'path':'../metadata.json','sha256':current,'original_sha256':old,'metadata_migration':{'path':'../migration.json','sha256':hashlib.sha256(receipt.read_bytes()).hexdigest()}}}})
  c=m.Config(cfg);row=c.input('metadata',old);self.assertEqual(row['sha256'],current);self.assertEqual(row['original_sha256'],old)
  receipt.write_text('{}')
  with self.assertRaisesRegex(ValueError,'receipt SHA-256 mismatch'):c.input('metadata',old)
 def test_missing_inputs_fail_read_only(self):
  cfg=self.config({'paths':{'output_dir':'out'},'inputs':{k:{'path':k} for k in ['source_alignment','cp98_tree','cp97_tree','nuclear_S0']}})
  before=set(self.base.rglob('*'))
  r=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/research'/SCRIPTS[5]),'--config',str(cfg)],cwd=self.base,capture_output=True,text=True)
  self.assertEqual(r.returncode,2,r.stderr);self.assertEqual(json.loads(r.stdout)['status'],'PREFLIGHT_FAILED');self.assertEqual(set(self.base.rglob('*')),before)
 def test_figure_preflight_and_explicit_synthetic_render(self):
  cfg=self.config({'paths':{'output_dir':'figures'},'inputs':{'distance_grid':{'path':'grid.csv'},'topology_counts':{'path':'counts.csv'}}})
  (cfg.parent/'grid.csv').write_text('a,b,conditional_tail_probability_analytic\n0,0,0.5\n0,0.01,0.1\n0.01,0,0.01\n0.01,0.01,0.001\n')
  (cfg.parent/'counts.csv').write_text('multiplier,topology,frequency\n'+''.join(f'{m},{t},{v}\n' for m in [1,2,4] for t,v in [('HT',.8),('HG',.1),('HZ',.1)]))
  command=[sys.executable,'-B',str(ROOT/'scripts/research'/SCRIPTS[6]),'--config',str(cfg)]
  result=subprocess.run(command,cwd=self.base,capture_output=True,text=True)
  self.assertEqual(result.returncode,0,result.stderr);self.assertFalse((cfg.parent/'figures').exists())
  result=subprocess.run(command+['--execute'],cwd=self.base,capture_output=True,text=True)
  self.assertEqual(result.returncode,0,result.stderr)
  self.assertEqual({p.name for p in (cfg.parent/'figures').iterdir()},{'distance_compatibility.png','distance_compatibility.pdf','genealogy_frequencies.png','genealogy_frequencies.pdf'})
  again=subprocess.run(command+['--execute'],cwd=self.base,capture_output=True,text=True)
  self.assertNotEqual(again.returncode,0)
 def test_fit_guard_still_rejects_missing_resource_evidence(self):
  m=load(SCRIPTS[3])
  with self.assertRaisesRegex(ValueError,'scope-receipt'):m.guard(None)
  self.assertEqual(m.NUMERIC['CU_starts'],[.1,1.,3.]);self.assertEqual(m.PROFILE_POINTS,41)
 def test_panel_threshold_and_topology_unchanged(self):
  m=load(SCRIPTS[4]);self.assertEqual(m.threshold({'DT_TOK':950,'DT_GRA':25,'DT_ZIN':25,'unresolved':0}),'DT_TOK')
  self.assertIsNone(m.threshold({'DT_TOK':949,'DT_GRA':26,'DT_ZIN':25,'unresolved':0}))
  with self.assertRaisesRegex(ValueError,'denominator'):m.threshold({'DT_TOK':100})
 def _run_qc_fixture(self, external_raw=False):
  m=load(SCRIPTS[0]);root=self.base/'task';raw=self.base/'data/raw' if external_raw else root/'raw';raw.mkdir(parents=True);root.mkdir(exist_ok=True)
  files=[]
  for mate,name in [(1,'synthetic_R1.gz'),(2,'synthetic_R2.gz')]:
   payload=b''.join(b'@read'+str(i).encode()+b'/'+str(mate).encode()+b'\n'+b'A'*150+b'\n+\n'+b'I'*150+b'\n' for i in range(2))
   p=raw/name;p.write_bytes(gzip.compress(payload,mtime=0));files.append((name,p.stat().st_size,hashlib.md5(p.read_bytes()).hexdigest()))
  settings=self.base/'settings.json';settings.write_text(json.dumps({'diagnostic_prefix_pairs':1,'adapter_r1':'CGTACGTACGTACGTA','adapter_r2':'CGTACGTACGTACGTA'}))
  download=self.base/'download.json';download.write_text(json.dumps({'status':'COMPLETE','pairs_done':2,'paired_fastq_validation':{'pairs':2,'bases':600,'bases_per_mate_record':150,**{k:True for k in ['mate_ids_all_match','record_structure_all_pass','sequence_quality_lengths_all_pass','gzip_crc_and_complete_eof_both_pass']}},'files':[{'name':n,'bytes':z,'md5':h,'official_md5_match':True} for n,z,h in files]}))
  cfg=self.config({'phase':'qc','paths':{'root':'../task','run_dir':'../task/qc/configured',**({'raw_dir':'../data/raw'} if external_raw else {})},'inputs':{'workflow_settings':{'path':'../settings.json'},'download_receipt':{'path':'../download.json'}}})
  with patch.object(m,'FILES',files),patch.object(m,'RAW_PAIRS',2),patch.object(m,'CONFIG_EXPECTED_SHA',m.hashfile(settings)),patch.object(sys,'argv',[str(m.__file__),'--config',str(cfg),'--execute']),patch.object(m.signal,'signal'),contextlib.redirect_stdout(io.StringIO()):
   self.assertEqual(m.main(),0)
  receipts=list((root/'qc').glob('*/QC_RECEIPT.json'));self.assertEqual(len(receipts),1)
  result=json.loads(receipts[0].read_text());self.assertEqual(result['status'],'COMPLETE');self.assertEqual(result['source_raw_pairs'],2)
  self.assertEqual(result['download_receipt_sha256'],m.hashfile(download));self.assertEqual(result['portable_config_sha256'],m.hashfile(cfg));self.assertFalse(result['formal_marker_eligible'])
 def test_rna_qc_config_to_execution_with_synthetic_receipts(self):
  self._run_qc_fixture()
 def test_rna_qc_external_raw_directory_is_used_for_all_hash_checks(self):
  self._run_qc_fixture(external_raw=True)
  self.assertFalse((self.base/'task/raw').exists())
 def test_marker_config_accepts_external_sources_then_rejects_pilot(self):
  m=load(SCRIPTS[1]);root=self.base/'output';root.mkdir()
  evidence=self.base/'evidence.json';evidence.write_text(json.dumps({'baseline':{'genes':173,'markers':[{'marker':str(i)} for i in range(173)],'main_sha256':'982addd02b6a6d0a188aa0de9b8e688b95bec926703094ab7f9238a8a91cc01e'}}))
  transcript=self.base/'transcripts.fa';transcript.write_text('>synthetic\nACGT\n')
  assembly=self.base/'assembly.json';assembly.write_text(json.dumps({'status':'COMPLETE','input_scope':'NONRANDOM_PREFIX_1M','source_raw_pairs':28665951,'formal_marker_eligible':False}))
  download=self.base/'download.json';download.write_text('{}')
  cfg=self.config({'paths':{'root':'../output'},'inputs':{k:{'path':str(p),'sha256':m.sha(p)} for k,p in [('marker_evidence',evidence),('transcripts',transcript),('assembly_receipt',assembly),('download_receipt',download)]},'tools':{k:sys.executable for k in ['busco','busco_python','metaeuk','systemd_run','systemctl','time']}})
  with patch.object(m,'EVIDENCE_SHA',m.sha(evidence)),patch.object(sys,'argv',[str(m.__file__),'--config',str(cfg),'--execute']):
   with self.assertRaisesRegex(ValueError,'Pilot/prefix'):m.main()
  self.assertEqual(list(root.iterdir()),[])
 def test_fake_mafft_command_uses_configured_tool_and_preserves_rows(self):
  m=load(SCRIPTS[2]);tool=self.base/'fake_mafft'
  tool.write_text('#!'+sys.executable+'\nimport pathlib,sys\nq=pathlib.Path(sys.argv[sys.argv.index("--add")+1]);ref=pathlib.Path(sys.argv[-1])\nprint(ref.read_text()+q.read_text(),end="")\nseq="".join(q.read_text().splitlines()[1:]);q.with_name(q.name+".map").write_text(q.read_text().splitlines()[0]+"\\n"+"".join(f"{a},{i},{i}\\n" for i,a in enumerate(seq,1)))\n')
  tool.chmod(0o755);m.MAFFT=str(tool)
  added,mapping,receipt=m.mafft_add({'old':'ACDE'},'ACDE',self.base/'result')
  self.assertEqual(added['old'],'ACDE');self.assertEqual(len(mapping),4)
  self.assertEqual(receipt['command'][0],str(tool));self.assertIn('--keeplength',receipt['command']);self.assertIn('--mapout',receipt['command'])

if __name__=='__main__':unittest.main()
