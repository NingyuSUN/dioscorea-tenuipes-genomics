import unittest
from io import StringIO
from Bio import Phylo
from gene_concordance import canonical,gene_splits,reference_edges,classify
T=lambda s:Phylo.read(StringIO(s),'newick')
class ConcordanceTests(unittest.TestCase):
 def setUp(self):self.U,self.b=reference_edges(T('((A,B),(C,D));'));self.groups=next(iter(self.b.values()))
 def test_root_suppression(self):self.assertEqual(len(self.b),1)
 def test_concordant_high(self):
  u,e=gene_splits(T('((A,B)90/99,(C,D)90/99);'));self.assertEqual(classify(self.groups,u,e),('concordant',True))
 def test_concordant_low(self):
  u,e=gene_splits(T('((A,B)79/99,(C,D)90/99);'));self.assertEqual(classify(self.groups,u,e),('concordant',False))
 def test_two_alternatives(self):
  outcomes=[]
  for s in ['((A,C)90/99,(B,D)90/99);','((A,D)90/99,(B,C)90/99);']:
   u,e=gene_splits(T(s));outcomes.append(classify(self.groups,u,e)[0])
  self.assertEqual(set(outcomes),{'alternative_1','alternative_2'})
 def test_missing_quadrant(self):
  u,e=gene_splits(T('(A,B,C);'));self.assertEqual(classify(self.groups,u,e),('not_decisive',False))
 def test_polytomy(self):
  u,e=gene_splits(T('(A,B,C,D);'));self.assertEqual(classify(self.groups,u,e),('other_or_unresolved',False))
 def test_many_tips_other(self):
  U,b=reference_edges(T('(((A,E),B),(C,D));'));groups=b[canonical({'A','E','B'},U)];u,e=gene_splits(T('((A,C),(B,E),D);'));self.assertEqual(classify(groups,u,e)[0],'other_or_unresolved')
 def test_root_invariance(self):
  u,e=gene_splits(T('(A,(B,(C,D)90/99)90/99);'));self.assertEqual(classify(self.groups,u,e),('concordant',True))
 def test_invalid_support(self):
  with self.assertRaises(ValueError):gene_splits(T('((A,B)90/101,C,D);'))
if __name__=='__main__':unittest.main()
