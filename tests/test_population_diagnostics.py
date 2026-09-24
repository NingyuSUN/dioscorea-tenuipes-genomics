import unittest,tempfile,subprocess,json,sys
from pathlib import Path
import pysam
SCRIPT_DIR=Path(__file__).resolve().parent.parent/'scripts'
if not (SCRIPT_DIR/'31_diagnose_population_samples.py').exists():SCRIPT_DIR=Path(__file__).resolve().parent
class TestDiagnostic(unittest.TestCase):
 def test_low_depth_and_allele_balance(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'ref.fai').write_text('c\t10000\t0\t0\t0\n')
   h='##fileformat=VCFv4.2\n##contig=<ID=c,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n'
   (p/'x.vcf').write_text(h+'c\t5000\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:AD\t0/1:20:50:10,10\t0/0:0:0:0,0\n'+'c\t5001\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ:AD\t0/0:10:30:10,0\t0/1:20:50:5,15\n')
   pysam.tabix_compress(str(p/'x.vcf'),str(p/'x.vcf.gz'));pysam.tabix_index(str(p/'x.vcf.gz'),preset='vcf')
   subprocess.run([sys.executable,str(SCRIPT_DIR/'31_diagnose_population_samples.py'),'--vcf',str(p/'x.vcf.gz'),'--fai',str(p/'ref.fai'),'--out',str(p/'out')],check=True,capture_output=True)
   r=json.loads((p/'out/summary.json').read_text());a,b=r['samples'];self.assertEqual(r['records'],2)
   self.assertEqual(a['good_dp5_gq20'],2);self.assertEqual(b['good_dp5_gq20'],1)
   self.assertEqual(a['allele_balance_median'],.5);self.assertEqual(b['allele_balance_median'],.75)
   self.assertEqual(b['homref_dp0_or_missing'],1)
if __name__=='__main__':unittest.main()
