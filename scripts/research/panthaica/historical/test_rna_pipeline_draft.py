"""Tiny synthetic validation of RNA QC/assembly preparation; no public data or jobs."""
import gzip,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import run_rna_pipeline_draft as d

def rec(name,mate,seq,quality=None):
    return b'@'+name+b'/'+str(mate).encode()+b'\n'+seq+b'\n+\n'+(quality or b'I'*len(seq))+b'\n'
class RNATests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name);self.cfg={'adapter_r1':'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA','adapter_r2':'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'}
    def tearDown(self):self.tmp.cleanup()
    def files(self,a,b):
        p=[self.p/'1.gz',self.p/'2.gz']
        for x,y in zip(p,[a,b]):x.write_bytes(gzip.compress(y,mtime=0))
        return p
    def test_diagnostic_known_motif_and_native_tail(self):
        a=b'ACGT'*20+self.cfg['adapter_r1'].encode();b=b'ACGT'*20+b'A'*30
        r=d.diagnostic(self.files(rec(b'x',1,a),rec(b'x',2,b)),1,self.cfg)
        self.assertEqual(r['mates'][0]['motif_read_counts']['NEB_R1_prefix16'],1)
        self.assertEqual(r['mates'][1]['terminal_homopolymer_ge10']['A'],1)
        self.assertEqual(sum(r['mates'][1]['motif_read_counts'].values()),0)
        self.assertFalse(r['random_sample'])
    def test_quality_numerators(self):
        r=d.diagnostic(self.files(rec(b'x',1,b'AC',b'I+'),rec(b'x',2,b'AC',b'I+')),1,self.cfg)
        self.assertEqual(r['mates'][0]['q20_bases'],1);self.assertEqual(r['mates'][0]['bases'],2)
        self.assertEqual(r['mates'][0]['mean_phred33'],25)
    def test_pair_id_mismatch(self):
        with self.assertRaisesRegex(ValueError,'Pair ID'):list(d.pairs(self.files(rec(b'x',1,b'AC'),rec(b'y',2,b'AC')),1))
    def test_prefix_does_not_include_later_records(self):
        p=self.files(rec(b'x',1,b'AC')+rec(b'z',1,b'GT'),rec(b'x',2,b'AC')+rec(b'z',2,b'GT'))
        self.assertEqual(len(list(d.pairs(p,1))),1)
    def test_incomplete_prefix_rejected(self):
        with self.assertRaisesRegex(ValueError,'Insufficient'):list(d.pairs(self.files(rec(b'x',1,b'AC'),rec(b'x',2,b'AC')),2))
    def test_trimmed_variable_lengths_and_exact_count(self):
        p=self.files(rec(b'x',1,b'A'*50)+rec(b'y',1,b'C'*70),rec(b'x',2,b'A'*60)+rec(b'y',2,b'G'*80))
        r=d.validate_trimmed(p,2,lambda **x:None);self.assertEqual(r['bases_per_mate'],[120,140])
        with self.assertRaisesRegex(ValueError,'count mismatch'):d.validate_trimmed(p,3,lambda **x:None)
    def test_corrupt_gzip_rejected(self):
        p=self.files(rec(b'x',1,b'A'*50),rec(b'x',2,b'A'*50));b=bytearray(p[1].read_bytes());b[-8]^=0xff;p[1].write_bytes(b)
        with self.assertRaises((gzip.BadGzipFile,EOFError)):d.validate_trimmed(p,1,lambda **x:None)
    def test_transcript_summary(self):
        p=self.p/'tx.fasta';p.write_text('>a\n'+'A'*100+'\n>b\n'+'C'*50+'\n')
        r=d.summarize_fasta(p);self.assertEqual((r['records'],r['bases'],r['N50']),(2,150,100))
    def test_invalid_transcript_fasta_rejected(self):
        p=self.p/'tx.fasta';p.write_text('>a\nAZ\n')
        with self.assertRaisesRegex(ValueError,'Malformed'):d.summarize_fasta(p)
    def test_failed_download_receipt_cannot_run(self):
        (self.p/'DOWNLOAD_RECEIPT.json').write_text(json.dumps({'status':'FAILED','pairs_done':d.RAW_PAIRS}))
        with patch.object(d,'ROOT',self.p),self.assertRaisesRegex(ValueError,'COMPLETE'):d.check_download()
    def test_missing_complete_pair_validation_cannot_run(self):
        (self.p/'DOWNLOAD_RECEIPT.json').write_text(json.dumps({'status':'COMPLETE','pairs_done':d.RAW_PAIRS}))
        with patch.object(d,'ROOT',self.p),self.assertRaisesRegex(ValueError,'paired FASTQ'):d.check_download()
    def test_config_hash_remains_declared(self):self.assertEqual(d.CONFIG_EXPECTED_SHA,'ed19abb14b85a9d10a626f2b43146a1a55bf565da038b4b07b0f0c4996a1611c')
if __name__=='__main__':unittest.main(verbosity=2)
