"""Regression checks for coordinate/strand/allele logic, not biological proof."""
import importlib.util
import unittest
from pathlib import Path


def module(name, file):
    path = Path(__file__).parent/file
    if not path.exists():
        path = Path(__file__).parent.parent/'scripts'/file
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


v = module('indels', '41_validate_plastome_indels.py')


class IndelTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(v.rc('ACGTA'), 'TACGT')

    def test_circular_origin(self):
        self.assertEqual(v.circular('ACGT', -2, 3), 'GTACG')

    def test_homopolymer_shift(self):
        self.assertEqual(v.shift_range('CGAAAATC',3,4,''),(-1,2))
        self.assertEqual(v.shift_range('CGAAATC',3,3,'A'),(-1,2))

    def test_tandem_repeat_shift(self):
        lo,hi=v.shift_range('GCCATATATGC',3,5,'')
        self.assertEqual((lo,hi),(0,4))

    def test_no_marker_is_unresolved(self):
        e={label+'_diagnostic_kmers':0 for label in ['own_left','own_right','public_left','public_right']}
        e['own_marker_min_genome_copies']=0
        from collections import defaultdict
        self.assertEqual(v.classify(e,defaultdict(set)), 'unresolved_no_discriminating_junction')
        self.assertEqual(e['repeat_copy_resolution'],'no_marker')

    def test_strong_contrary_evidence_not_passed(self):
        e={label+'_diagnostic_kmers':1 for label in ['own_left','own_right','public_left','public_right']}
        e['own_marker_min_genome_copies']=1
        from collections import defaultdict
        s=defaultdict(set,own_left={'a'},own_right={'a'},public_left={'b','c','d'},public_right={'b','c','d'})
        self.assertEqual(v.classify(e,s),'conflicting_or_public_allele_evidence')


class AnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=module('annotation','40_annotate_plastome.py')

    def test_gap_projection(self):
        m,ev=self.m.alignment_map('AC-GTAA','ACTG-AA')
        self.assertEqual(m,{0:0,1:1,2:3,3:None,4:4,5:5})
        self.assertEqual([e['type'] for e in ev],['insertion','deletion'])

    def test_inserted_bases_inside_exon(self):
        from Bio.SeqFeature import FeatureLocation
        m,_=self.m.alignment_map('AC-GT','ACTGT')
        loc,_=self.m.lift_location(FeatureLocation(0,4,strand=-1),m)
        from Bio.Seq import Seq
        self.assertEqual(str(loc.extract(Seq('ACTGT'))),'ACAGT')

    def test_trans_spliced_part_order(self):
        from Bio.SeqFeature import FeatureLocation,CompoundLocation
        loc=CompoundLocation([FeatureLocation(6,9,strand=-1),FeatureLocation(0,3,strand=-1)])
        new,_=self.m.lift_location(loc,{i:i+2 for i in range(12)})
        self.assertEqual([(int(p.start),int(p.end),p.strand) for p in new.parts],[(8,11,-1),(2,5,-1)])

    def test_deleted_boundary_flag(self):
        from Bio.SeqFeature import FeatureLocation
        new,issues=self.m.lift_location(FeatureLocation(0,3),{0:None,1:0,2:1})
        self.assertIn('deleted_feature_boundary',issues)
        self.assertEqual((int(new.start),int(new.end)),(0,2))


if __name__=='__main__':
    unittest.main()
