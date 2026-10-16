"""Synthetic profile/map/control tests. Does not invoke MAFFT or read server data."""
import tempfile,unittest
from pathlib import Path
import prepare_alignment_controls as m

class ProfileControlTests(unittest.TestCase):
    def setUp(self):
        self.old={'old1':'ACD-EFGH','old2':'A-D-E-GH'}
        self.added={**self.old,m.NEW:'ACD-EFG-'}
        self.query='ACDPQEFG'
        self.mapping='>'+m.NEW+'\n# letter, original, reference\nA,1,1\nC,2,2\nD,3,3\nP,4,-\nQ,5,-\nE,6,5\nF,7,6\nG,8,7\n'
    def test_reference_rows_frozen(self):m.check_old_rows(self.old,self.added)
    def test_changed_old_residue_rejected(self):
        bad={**self.added,'old1':'ACD-EFGA'}
        with self.assertRaisesRegex(ValueError,'changed old row'):m.check_old_rows(self.old,bad)
    def test_deleted_gap_column_rejected(self):
        bad={k:v.replace('-','') for k,v in self.added.items()}
        with self.assertRaisesRegex(ValueError,'changed old row'):m.check_old_rows(self.old,bad)
    def test_missing_old_taxon_rejected(self):
        bad={k:v for k,v in self.added.items() if k!='old2'}
        with self.assertRaisesRegex(ValueError,'taxon set'):m.check_old_rows(self.old,bad)
    def test_insertion_and_deletion_map(self):
        rows=m.validate_map(self.mapping,self.query,self.added[m.NEW])
        self.assertEqual(sum(r[2] is None for r in rows),2)
        self.assertEqual([r[1] for r in rows if r[2] is None],[4,5])
    def test_map_residue_mutation_rejected(self):
        with self.assertRaisesRegex(ValueError,'residue differs'):m.validate_map(self.mapping.replace('P,4,-','A,4,-'),self.query,self.added[m.NEW])
    def test_map_missing_source_position_rejected(self):
        with self.assertRaisesRegex(ValueError,'exact candidate positions'):m.validate_map(self.mapping.replace('P,4,-\n',''),self.query,self.added[m.NEW])
    def test_map_duplicate_target_rejected(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):m.validate_map(self.mapping.replace('F,7,6','F,7,5'),self.query,self.added[m.NEW])
    def test_map_output_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError,'reconstruct'):m.validate_map(self.mapping,self.query,'ACD-EFGA')
    def test_project_mask_missing_taxon_and_frozen_columns(self):
        p=m.project(self.added,['old1','old2','missing_old',m.NEW],[0,2,4,5,7])
        self.assertEqual(p,{'old1':'ADEFH','old2':'ADE-H','missing_old':'-----',m.NEW:'ADEF-'})
    def test_unknown_aa_masked_consistently(self):
        p=m.project({'old':'AXD-'},['old'],[0,1,2,3]);self.assertEqual(p['old'],'A-D-')
    def test_nonmonotonic_mask_rejected(self):
        with self.assertRaisesRegex(ValueError,'non-monotonic'):m.project(self.old,['old1'],[0,2,1])
    def test_same_locus_old_rows_and_byte_prefix(self):
        loci={'g1':{'a':'AC-','b':'A--',m.NEW:'ACD'},'g2':{'a':'EF','b':'-F',m.NEW:'--'},'g3':{'a':'GH','b':'-H',m.NEW:'GH'}}
        full=m.concatenate(loci,['a','b',m.NEW],['g1','g2','g3'])
        matched=m.concatenate(loci,['a','b',m.NEW],['g1','g3']);control={t:matched[t] for t in ['a','b']}
        positions=[0,1,2,5,6]
        for t in control:self.assertEqual(control[t],''.join(full[t][i] for i in positions))
        self.assertTrue(m.fasta_bytes(matched).startswith(m.fasta_bytes(control)))
    def test_fasta_duplicate_id_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.faa';p.write_text('>x\nAC\n>x\nGT\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):m.read_fasta(p)
if __name__=='__main__':unittest.main(verbosity=2)
