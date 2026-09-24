import sys,io,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
class ProbeTests(unittest.TestCase):
 def test_valid_fastq_and_mate_names(self):
  from wgs_probe_core import read_record,read_id
  x=read_record(io.BytesIO(b'@SRR33602354.1 1 length=4\nACGN\n+\nI?5!\n'))
  self.assertEqual(x[1],b'ACGN');self.assertEqual(read_id(x[0]),b'SRR33602354.1')
  self.assertEqual(read_id(b'@sample/1'),read_id(b'@sample/2'))
 def test_truncated_record_rejected(self):
  from wgs_probe_core import read_record
  with self.assertRaises(ValueError):read_record(io.BytesIO(b'@x\nAC\n+\n'))
 def test_mismatched_quality_rejected(self):
  from wgs_probe_core import read_record
  with self.assertRaises(ValueError):read_record(io.BytesIO(b'@x\nAC\n+\nI\n'))
 def test_prefix_byte_budget_enforced(self):
  from wgs_probe_core import LimitedReader
  r=LimitedReader(io.BytesIO(b'abcdef'),4);self.assertEqual(r.read(3),b'abc');self.assertEqual(r.read(3),b'd')
  with self.assertRaises(RuntimeError):r.read(1)
if __name__=='__main__':unittest.main()
