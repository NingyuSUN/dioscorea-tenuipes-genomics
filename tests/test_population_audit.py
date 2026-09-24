import importlib.util,unittest,tempfile,pathlib,gzip,subprocess,json,os
SCRIPT_DIR=pathlib.Path(__file__).resolve().parent.parent/'scripts'
if not (SCRIPT_DIR/'30_audit_population_vcf.py').exists():SCRIPT_DIR=pathlib.Path(__file__).resolve().parent
AUDIT_BINARY=pathlib.Path(os.environ.get('DTE_AUDIT_BINARY',str(SCRIPT_DIR/'compare_population_genotypes'))).resolve()
spec=importlib.util.spec_from_file_location('audit',SCRIPT_DIR/'30_audit_population_vcf.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Tests(unittest.TestCase):
 def test_alias(self):
  self.assertEqual(m.alias('Dte_AIC-F20_S13_L001'),'AIC-F20-S13');self.assertEqual(m.alias('DRR351754'),'DRR351754')
 def test_missing_and_phase(self):
  self.assertTrue(m.missing((0,None)));self.assertTrue(m.missing((None,None)));self.assertFalse(m.missing((0,0)))
  self.assertEqual(m.canonical((1,0)),m.canonical((0,1)));self.assertNotEqual(m.canonical((0,None)),m.canonical((0,0)))
 def test_full_compare(self):
  with tempfile.TemporaryDirectory() as td:
   p=pathlib.Path(td)
   header='##fileformat=VCFv4.2\n##contig=<ID=c,length=8>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
   (p/'a.vcf').write_text(header+'Dte_AIC-F20_S13_L001\n'+'c\t1\t.\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t./.:0:.\n'+'c\t2\t.\tC\tT\t30\tPASS\t.\tGT:DP:GQ\t0/1:8:30\n'+'c\t4\t.\tT\tG\t30\tPASS\t.\tGT:DP:GQ\t0/.:3:0\n')
   (p/'b.vcf').write_text(header+'AIC-F20-S13\n'+'c\t1\t.\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t0/0:0:.\n'+'c\t3\t.\tG\tA\t30\tPASS\t.\tGT:DP:GQ\t1:5:30\n'+'c\t4\t.\tT\tG\t30\tPASS\t.\tGT:DP:GQ\t0/0:3:0\n')
   (p/'ref.fa').write_text('>c\nACGTACGT\n')
   r=m.run(str(p/'a.vcf'),str(p/'b.vcf'),str(p/'ref.fa'),str(p/'out'),0)
   self.assertEqual(r['comparison']['matched_records'],2);self.assertEqual(r['comparison']['missing_to_homref'],2)
   self.assertEqual(r['comparison']['left_only'],1);self.assertEqual(r['comparison']['right_only'],1)
   self.assertEqual(r['right']['samples'][0]['homref_dp_zero_or_missing'],1)
   self.assertEqual(r['right']['samples'][0]['non_diploid'],1)
   self.assertEqual(r['left']['ref_mismatch'],0);self.assertTrue(r['complete'])
   exe=AUDIT_BINARY
   env=dict(os.environ,LD_LIBRARY_PATH=str(exe.parent/'bcftools_runtime'))
   proc=subprocess.run([str(exe),str(p/'a.vcf'),str(p/'b.vcf'),str(p/'fast')],text=True,capture_output=True,env=env,check=True)
   fast=json.loads(proc.stdout)
   self.assertEqual(fast['matched_records'],2);self.assertEqual(fast['missing_to_homref'],2)
   self.assertEqual(fast['changed_gt'],2);self.assertEqual(fast['compared_genotypes'],2)
   self.assertEqual(fast['left_only'],1);self.assertEqual(fast['right_only'],1)
   # Duplicate exact variants must not be silently paired.
   with (p/'a.vcf').open('a') as f:f.write('c\t4\t.\tT\tG\t30\tPASS\t.\tGT:DP:GQ\t0/.:3:0\n')
   bad=subprocess.run([str(exe),str(p/'a.vcf'),str(p/'b.vcf'),str(p/'bad')],text=True,capture_output=True,env=env)
   self.assertNotEqual(bad.returncode,0);self.assertIn('Duplicate exact variant',bad.stderr)
 def test_phase_and_multirecord_coordinate(self):
  with tempfile.TemporaryDirectory() as td:
   p=pathlib.Path(td)
   h='##fileformat=VCFv4.2\n##contig=<ID=c,length=20>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
   a=['c\t1\t.\tA\tG\t30\tPASS\t.\tGT\t1|0','c\t1\t.\tA\tT\t30\tPASS\t.\tGT\t0/0','c\t2\t.\tC\tT\t30\tPASS\t.\tGT\t.']
   b=['c\t1\t.\tA\tT\t30\tPASS\t.\tGT\t0|0','c\t1\t.\tA\tG\t30\tPASS\t.\tGT\t0/1','c\t2\t.\tC\tT\t30\tPASS\t.\tGT\t0']
   (p/'a').write_text(h+'\n'.join(a)+'\n');(p/'b').write_text(h+'\n'.join(b)+'\n')
   exe=AUDIT_BINARY;env=dict(os.environ,LD_LIBRARY_PATH=str(exe.parent/'bcftools_runtime'))
   def call(tag):return subprocess.run([str(exe),str(p/'a'),str(p/'b'),str(p/tag)],text=True,capture_output=True,env=env)
   r=call('ok');self.assertEqual(r.returncode,0,r.stderr);d=json.loads(r.stdout)
   self.assertEqual(d['matched_records'],3);self.assertEqual(d['changed_gt'],1);self.assertEqual(d['missing_to_homref'],1)
   (p/'a').write_text(h+'\n'.join(reversed(a))+'\n');r=call('unsorted');self.assertNotEqual(r.returncode,0);self.assertIn('Unsorted',r.stderr)
   (p/'a').write_text(h+a[0].replace('1|0','0/0/1')+'\n');r=call('ploidy');self.assertNotEqual(r.returncode,0);self.assertIn('Ploidy',r.stderr)
if __name__=='__main__':unittest.main()
