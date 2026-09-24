import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
HERE=Path(__file__).resolve().parent
if not (HERE / "32_filter_population_variants.py").exists():
    HERE = HERE.parent / "scripts"
s=importlib.util.spec_from_file_location('summary',HERE/'34_summarize_population_results.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class SummaryTests(unittest.TestCase):
 def test_bed_decoding_and_padding(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'fixture'
   p.with_suffix('.fam').write_text('a a 0 0 0 -9\nb b 0 0 0 -9\nc c 0 0 0 -9\n')
   p.with_suffix('.bim').write_text('ctg v1 0 1 A G\nctg v2 0 2 A G\n')
   # sample codes [00,01,10], then [11,10,00], with ignored fourth-sample bits.
   p.with_suffix('.bed').write_bytes(bytes([108,27,1,36,11]))
   np.testing.assert_equal(m.read_bed(p),[[2,np.nan,1],[0,1,2]])
 def test_component_alignment(self):
  q=np.array([[.1,.9],[.8,.2]])
  np.testing.assert_equal(m.align_q(q[:,::-1],q),q)
 def test_constant_qc_is_undefined(self):
  self.assertIsNone(m.correlation([0,0,0],[1,2,3]))
if __name__=='__main__':unittest.main()
