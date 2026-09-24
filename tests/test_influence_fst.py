import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from influence_fst_core import population_counts, wc84_pair, WeightedFstAccumulator, panel_site_mask, validate_genotypes

class FstTests(unittest.TestCase):
    def test_separation_monomorphism_and_missing(self):
        g=np.array([[0,0,0,0,2,2,2,2], [0]*8, [-1]*4+[0,1,1,2], [1]*8],dtype=np.int8)
        c=wc84_pair(population_counts(g,range(4)),population_counts(g,range(4,8)))
        self.assertAlmostEqual(c.fst[0],1.0)
        self.assertTrue(np.isnan(c.fst[1]))
        self.assertTrue(np.isnan(c.fst[2]))
        self.assertAlmostEqual(c.fst[3],0.0)
        a=WeightedFstAccumulator();a.update(c)
        self.assertEqual(a.selected_sites,4); self.assertEqual(a.defined_sites,2)
        self.assertAlmostEqual(a.weighted_fst,float(c.numerator[[0,3]].sum()/c.denominator[[0,3]].sum()))

    def test_mask_subset_refilter_and_complete(self):
        # A rare variant carried only by removed sample is excluded after refilter.
        g=np.zeros((4,14),dtype=np.int8);g[:,0]=1;g[1,1]=-1;g[2,1:3]=-1;g[3,1]=1
        np.testing.assert_array_equal(panel_site_mask(g,range(14)),[True,True,False,True])
        np.testing.assert_array_equal(panel_site_mask(g,range(1,14)),[False,False,False,True])
        np.testing.assert_array_equal(panel_site_mask(g,range(14),max_missing=0),[True,False,False,True])

    def test_negative_retained_and_chunk_invariance(self):
        # Same intermediate frequency with homozygotes in both groups gives negative WC84 estimate.
        g=np.array([[0,0,2,2,0,0,2,2],[0,0,0,0,2,2,2,2]],dtype=np.int8)
        c=wc84_pair(population_counts(g,range(4)),population_counts(g,range(4,8)))
        self.assertLess(c.fst[0],0)
        a=WeightedFstAccumulator(); a.update(c)
        b=WeightedFstAccumulator()
        for i in range(2):b.update(wc84_pair(population_counts(g[i:i+1],range(4)),population_counts(g[i:i+1],range(4,8))))
        self.assertAlmostEqual(a.weighted_fst,b.weighted_fst)
        self.assertEqual(a.negative_sites,1)

    def test_invalid_encoding_rejected(self):
        for bad in [np.array([[0,3]],dtype=np.int8),np.array([[0.,1.5]])]:
            with self.assertRaises(ValueError):validate_genotypes(bad)

if __name__=='__main__':unittest.main()
