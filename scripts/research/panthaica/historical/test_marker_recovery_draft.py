"""Preparation-only tests; synthetic strings and temporary files, no biological run."""
import tempfile,unittest
from pathlib import Path
import run_marker_recovery_draft as m
class PreparedGateTests(unittest.TestCase):
    def test_complete_normal_length(self):self.assertIsNone(m.quality_reason('A'*100,100))
    def test_internal_stop(self):self.assertEqual(m.quality_reason('A'*50+'*'+'A'*49,100),'internal_stop')
    def test_single_terminal_stop(self):self.assertIsNone(m.quality_reason('A'*100+'*',100))
    def test_multiple_terminal_stops_rejected(self):self.assertEqual(m.quality_reason('A'*100+'**',100),'internal_stop')
    def test_short_sequence(self):self.assertEqual(m.quality_reason('A'*49,60),'shorter_than_50_aa')
    def test_excess_noncanonical(self):self.assertEqual(m.quality_reason('A'*94+'X'*6,100),'canonical_aa_fraction_below_0.95')
    def test_length_lower_boundary(self):
        self.assertIsNone(m.quality_reason('A'*70,100));self.assertEqual(m.quality_reason('A'*69,100),'length_ratio_outside_0.70_to_1.50')
    def test_length_upper_boundary(self):
        self.assertIsNone(m.quality_reason('A'*150,100));self.assertEqual(m.quality_reason('A'*151,100),'length_ratio_outside_0.70_to_1.50')
    def test_duplicate_fasta_id_fails(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.faa';p.write_text('>a\nAA\n>a\nAA\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):m.fasta(p)
    def test_preheader_sequence_fails(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.faa';p.write_text('AA\n>a\nAA\n')
            with self.assertRaisesRegex(ValueError,'before'):m.fasta(p)
    def test_pilot_cannot_enter_marker_recovery(self):
        with self.assertRaisesRegex(ValueError,'Pilot/prefix'):
            m.validate_assembly_scope({'input_scope':'NONRANDOM_PREFIX_1M','source_raw_pairs':28665951,'formal_marker_eligible':False})
    def test_full_input_scope_can_proceed_to_other_checks(self):
        m.validate_assembly_scope({'input_scope':'FULL_CRR515167','source_raw_pairs':28665951,'formal_marker_eligible':True})
    def test_fixed_evidence_has_173(self):
        e=m.evidence();self.assertEqual(len(e['baseline']['markers']),173)
        self.assertNotIn('2661at3193',{x['marker'] for x in e['baseline']['markers']})
if __name__=='__main__':unittest.main()
