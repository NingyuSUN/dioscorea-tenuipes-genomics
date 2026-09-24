"""Hand-calculated fixtures catch exon skipping, wrong strand and false support."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import subprocess
import json
import sys

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/29_prepare_rna_evidence.py'

class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.is_file(), 'RNA evidence builder is not implemented')
        spec=importlib.util.spec_from_file_location('rna',SCRIPT)
        self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)

    def test_middle_exon_is_never_swallowed_and_models_are_deduplicated(self):
        rows=[('c','+',1,30,'t1','transcript'),('c','+',1,5,'t1','exon'),
              ('c','+',11,15,'t1','exon'),('c','+',26,30,'t1','exon'),
              ('c','+',1,30,'t2','transcript'),('c','+',1,5,'t2','exon'),
              ('c','+',11,15,'t2','exon'),('c','+',26,30,'t2','exon')]
        valid,issues=self.parse(rows,{'c':'A'*30})
        introns=self.m.derive_introns(valid)
        self.assertEqual(set(introns),{('c',6,10,'+'),('c',16,25,'+')})
        self.assertEqual(set(introns[('c',6,10,'+')]),{'t1','t2'})
        self.assertEqual(issues,{})

    def parse(self,rows,genome):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.gtf'
            p.write_text(''.join(f'{c}\tStringTie\t{kind}\t{s}\t{e}\t1000\t{st}\t.\ttranscript_id "{t}"; gene_id "g";\n' for c,st,s,e,t,kind in rows))
            return self.m.parse_gtf(p,genome)

    def test_bad_exon_quarantines_whole_transcript_without_reconnecting_neighbors(self):
        rows=[('c','+',1,30,'t1','transcript'),('c','+',1,5,'t1','exon'),
              ('c','+',8,40,'t1','exon'),('c','+',26,30,'t1','exon')]
        valid,issues=self.parse(rows,{'c':'A'*30})
        self.assertEqual(valid,{})
        self.assertIn('coordinate_out_of_bounds',issues['t1'])

    def test_overlap_or_missing_transcript_is_quarantined(self):
        rows=[('c','+',1,30,'over','transcript'),('c','+',1,15,'over','exon'),
              ('c','+',12,30,'over','exon'),('c','+',1,5,'orphan','exon')]
        valid,issues=self.parse(rows,{'c':'A'*30})
        self.assertEqual(valid,{})
        self.assertIn('overlapping_or_adjacent_exons',issues['over'])
        self.assertIn('missing_or_multiple_transcript_record',issues['orphan'])

    def test_reverse_lift_uses_component_offset_and_flips_strand(self):
        placement=('chr1',101,130,11,40,'-')
        self.assertEqual(self.m.lift(('c',15,20,'+'),{'c':placement}),('chr1',121,126,'-'))
        self.assertEqual(self.m.lift(('c',15,20,'-'),{'c':placement}),('chr1',121,126,'+'))
        plus=('chr1',101,130,11,40,'+')
        self.assertEqual(self.m.lift(('c',15,20,'+'),{'c':plus}),('chr1',105,110,'+'))

    def test_strand_aware_motif_and_ambiguous_sequence(self):
        genome={'p':'AAAAAGTCCCCAGAAAA','m':'AAAAACTCCCCACAAAA','n':'AAAAAGTNNNNAGAAAA'}
        self.assertEqual(self.m.classify(('p',6,13,'+'),genome,4),('GT-AG','canonical'))
        self.assertEqual(self.m.classify(('m',6,13,'-'),genome,4),('GT-AG','canonical'))
        self.assertEqual(self.m.classify(('n',6,13,'+'),genome,4)[1],'ambiguous_sequence')
        self.assertEqual(self.m.classify(('p',6,13,'+'),genome,20)[1],'short_intron')

    def test_agp_sequence_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.agp';p.write_text('scaf\t1\t4\t1\tW\tc\t1\t4\t-\n')
            with self.assertRaisesRegex(ValueError,'sequence mismatch'):
                self.m.load_agp(p,{'c':'AAAA'},{'scaf':'AAAA'})
            self.assertEqual(self.m.load_agp(p,{'c':'AAAA'},{'scaf':'TTTT'})['c'],('scaf',1,4,1,4,'-'))

    def test_gff_does_not_invent_read_counts(self):
        row=self.m.gff_row('rna0001',('c',6,30,'+'),'GT-AG',2)
        self.assertNotIn('mult=',row)
        self.assertIn('transcript_models=2',row)
        self.assertEqual(row.split('\t')[5],'.')

    def test_unknown_strand_remains_unknown_on_reverse_component(self):
        self.assertEqual(self.m.lift(('c',15,20,'.'),{'c':('chr1',101,130,11,40,'-')}),('chr1',121,126,'.'))

    def test_cli_removes_false_bridge_and_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'genome.fa').write_text('>c\n'+'AAAAA'+'GT'+'C'*16+'AG'+'AAAAA'+'GT'+'C'*16+'AG'+'AAAAA\n')
            (p/'scaffold.fa').write_text('>scaf\n'+'TTTTT'+'CT'+'G'*16+'AC'+'TTTTT'+'CT'+'G'*16+'AC'+'TTTTT\n')
            (p/'x.agp').write_text('scaf\t1\t55\t1\tW\tc\t1\t55\t-\n')
            (p/'x.gtf').write_text(''.join(f'c\tStringTie\t{kind}\t{s}\t{e}\t1000\t+\t.\ttranscript_id "t";\n' for kind,s,e in [('transcript',1,55),('exon',1,5),('exon',26,30),('exon',51,55)]))
            (p/'old.gff').write_text('c\tStringTie\tintron\t6\t50\t.\t+\t.\tmult=100;src=E\n')
            cmd=[sys.executable,str(SCRIPT),'--gtf',str(p/'x.gtf'),'--genome',str(p/'genome.fa'),'--agp',str(p/'x.agp'),'--scaffold',str(p/'scaffold.fa'),'--old-introns',str(p/'old.gff'),'--out',str(p/'out')]
            run=subprocess.run(cmd,text=True,capture_output=True)
            self.assertEqual(run.returncode,0,run.stderr)
            summary=json.loads((p/'out/summary.json').read_text())
            self.assertEqual(summary['unique_introns'],2)
            self.assertEqual(summary['old_intron_comparison'],{'absent_from_valid_full_transcripts':1})
            gff=p/'out/introns.canonical.Kita1.gff3';before=gff.read_bytes()
            coords=[(x.split('\t')[3],x.split('\t')[4],x.split('\t')[6]) for x in gff.read_text().splitlines() if not x.startswith('#')]
            self.assertEqual(coords,[('31','50','-'),('6','25','-')])
            repeat=subprocess.run(cmd,text=True,capture_output=True)
            self.assertNotEqual(repeat.returncode,0)
            self.assertEqual(gff.read_bytes(),before)

if __name__=='__main__':unittest.main()
