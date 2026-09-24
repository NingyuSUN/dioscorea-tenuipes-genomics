import sys, unittest, tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
    from relationship_core import read_bed, ibs_distance, nj_tree, tree_splits, canonical_split
except ImportError:
    read_bed=ibs_distance=nj_tree=tree_splits=canonical_split=None

class RelationshipsTest(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(ibs_distance,'relationship implementation not present yet')
    def test_bed_little_endian_missing_and_padding(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'x.bed';p.write_bytes(bytes([0x6c,0x1b,1,56,33]))
            np.testing.assert_array_equal(read_bed(p,3,2),[[2,1,0],[-1,2,1]])
    def test_missing_is_pairwise_excluded_not_reference(self):
        g=np.array([[2,1,0],[-1,2,1]])
        d,n=ibs_distance(g)
        np.testing.assert_allclose(d,[[0,.5,1],[.5,0,.5],[1,.5,0]])
        self.assertEqual(n[0,1],1);self.assertEqual(n[1,2],2)
    def test_heterozygote_identity_and_opposite_homozygotes(self):
        d,n=ibs_distance(np.array([[1,1,0,2],[1,1,0,2]]))
        self.assertEqual(d[0,1],0);self.assertEqual(d[2,3],1);self.assertEqual(d[0,2],.5)
    def test_allele_reversal_leaves_ibs_unchanged(self):
        g=np.array([[0,1,2],[2,0,1],[-1,2,0]])
        rev=np.where(g<0,-1,2-g)
        np.testing.assert_allclose(ibs_distance(g)[0],ibs_distance(rev)[0])
    def test_pair_without_comparable_sites_rejected(self):
        with self.assertRaises(ValueError):ibs_distance(np.array([[0,-1],[-1,2]]))
    def test_additive_quartet_topology(self):
        names=['A','B','C','D'];d=np.array([[0,2,6,6],[2,0,6,6],[6,6,0,2],[6,6,2,0.]])
        tree=nj_tree(names,d)
        self.assertEqual(tree_splits(tree),{('A','B')})
        for i,a in enumerate(names):
            for j,b in enumerate(names):self.assertAlmostEqual(tree.distance(a,b),d[i,j])
    def test_split_canonical_complement(self):
        self.assertEqual(canonical_split({'C','D'},set('ABCD')),('A','B'))
        self.assertIsNone(canonical_split({'A'},set('ABCD')))


class NetworkChecks(unittest.TestCase):
    def test_network_recovers_additive_quartet_distances(self):
        from splitspy.nnet.nnet_algo import neighbor_net
        from splitspy.splits.basic_split import split_dist
        d=np.array([[0,2,6,6],[2,0,6,6],[6,6,0,2],[6,6,2,0.]])
        cycle,splits=neighbor_net(list('ABCD'),d.tolist(),cutoff=1e-10)
        np.testing.assert_allclose(np.array(split_dist(4,splits))[1:,1:],d,atol=1e-8)
    def test_residual_fit_detects_nonzero_residual(self):
        from relationship_core import residual_fit
        self.assertAlmostEqual(residual_fit([[0,2],[2,0]],[[0,1],[1,0]]),75)
    def test_invalid_bed_size_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'bad.bed';p.write_bytes(b'\x6c\x1b\x01')
            with self.assertRaises(ValueError):read_bed(p,3,2)

class CheckedNetworkTest(unittest.TestCase):
    def test_repeatable_nonnegative_weights_and_optimality(self):
        from relationship_core import checked_neighbor_net
        d=np.array([[0,.2,.31,.28],[.2,0,.29,.32],[.31,.29,0,.18],[.28,.32,.18,0]])
        a=checked_neighbor_net(list('ABCD'),d);b=checked_neighbor_net(list('ABCD'),d)
        self.assertEqual(a[0],b[0]);np.testing.assert_array_equal([s.weight for s in a[1]],[s.weight for s in b[1]])
        self.assertLess(a[2]['kkt_violation'],1e-9)

if __name__=='__main__':unittest.main()
