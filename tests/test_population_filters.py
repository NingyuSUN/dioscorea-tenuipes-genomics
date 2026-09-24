import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
if not (HERE / "32_filter_population_variants.py").exists():
    HERE = HERE.parent / "scripts"


class FilterTests(unittest.TestCase):
    def test_mask_is_per_sample_and_never_fills_reference(self):
        spec = importlib.util.spec_from_file_location('pipeline', HERE / '32_filter_population_variants.py')
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            header = ('##fileformat=VCFv4.2\n##contig=<ID=contig_1,length=100>\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\tE\tF\n')
            (p/'in.vcf').write_text(header + 'contig_1\t1\t.\tA\tG\t100\tPASS\t.\tGT:DP:GQ\t0/1:5:20\t1/1:4:50\t0/0:20:.\t0/1:.:40\t./.:30:50\t0/1:101:60\n' +
                'contig_1\t2\t.\tA\tG\t100\tPASS\t.\tGT:DP:GQ\t0/1:100:20\t0/0:20:50\t0/0:20:50\t0/0:20:40\t0/0:30:50\t0/0:20:60\n')
            m.run_pipe(m.mask_commands(p/'in.vcf', 5, 20, 100), p/'mask.bcf', p/'mask.log')
            q = subprocess.check_output(['bcftools','query','-f','[%GT\t]\n',str(p/'mask.bcf')],text=True)
            self.assertEqual(q.splitlines()[0].split(), ['0/1','./.','./.','./.','./.','./.'])
            self.assertEqual(q.splitlines()[1].split()[0], '0/1')
            (p/'keep.txt').write_text('B\nC\nD\nE\nF\n')
            m.run_pipe(m.panel_commands(p/'mask.bcf',p/'keep.txt',.1),p/'panel.bcf',p/'panel.log')
            q = subprocess.check_output(['bcftools','query','-f','%POS\n',str(p/'panel.bcf')],text=True)
            self.assertEqual(q,'')  # ALT disappears after subset; no stale AC survives.

    def test_sample_contract(self):
        spec = importlib.util.spec_from_file_location('pipeline', HERE / '32_filter_population_variants.py')
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        self.assertEqual(m.EXCLUDED, {'AIC-M6-S11','AIC-F20-S13'})


if __name__ == '__main__':
    unittest.main()
