"""Independent split-based gene-concordance counts, including explicit low-support resolution."""
import argparse,csv,json,hashlib,math
from collections import Counter
from pathlib import Path
from Bio import Phylo

def canonical(side,universe):
 a,b=frozenset(side),frozenset(universe)-frozenset(side)
 return min((a,b),key=lambda s:(len(s),tuple(sorted(s))))
def gene_splits(tree):
 names=[t.name for t in tree.get_terminals()]
 if len(names)!=len(set(names)):raise ValueError('duplicate tips')
 U=frozenset(names);edges={}
 for c in tree.find_clades():
  if c is tree.root:continue
  side=frozenset(t.name for t in c.get_terminals())
  if min(len(side),len(U-side))<2:continue
  key=canonical(side,U);raw=c.name if c.name is not None else c.confidence;high=False
  if isinstance(raw,str) and '/' in raw:
   values=raw.split('/')
   if len(values)!=2:raise ValueError('Unexpected gene support label')
   sh,uf=map(float,values)
   if not all(math.isfinite(x) and 0<=x<=100 for x in [sh,uf]):raise ValueError('Invalid support')
   high=sh>=80 and uf>=95
  edges[key]=edges.get(key,True) and high
 return U,edges

def reference_edges(tree):
 nodes=list(tree.find_clades());adj={n:set() for n in nodes};labels={n:n.name for n in tree.get_terminals()};U=frozenset(labels.values())
 if len(U)!=len(labels):raise ValueError('duplicate reference tips')
 for n in nodes:
  for ch in n.clades:adj[n].add(ch);adj[ch].add(n)
 for n in list(nodes):
  if n not in labels and len(adj[n])==2:
   a,b=adj[n];adj[a].remove(n);adj[b].remove(n);adj[a].add(b);adj[b].add(a);del adj[n]
 if any(len(adj[n])!=3 for n in adj if n not in labels):raise ValueError('Reference must be unrooted binary after root suppression')
 def leaves(n,parent):
  if n in labels:return {labels[n]}
  return set().union(*(leaves(ch,n) for ch in adj[n] if ch is not parent))
 out={}
 for a in adj:
  if a in labels:continue
  for b in adj[a]:
   if b in labels:continue
   side=leaves(a,b);key=canonical(side,U)
   if key in out:continue
   groups=[frozenset(leaves(ch,a)) for ch in adj[a] if ch is not b]+[frozenset(leaves(ch,b)) for ch in adj[b] if ch is not a]
   # Stable ordering preserves the reference pairing but makes reporting deterministic.
   first=sorted(groups[:2],key=lambda x:tuple(sorted(x)));second=sorted(groups[2:],key=lambda x:tuple(sorted(x)))
   if canonical(first[0]|first[1],U)!=first[0]|first[1]:first,second=second,first
   out[key]=first+second
 return U,out

def classify(groups,universe,edges):
 g=[s&universe for s in groups]
 if any(not s for s in g):return 'not_decisive',False
 pairings=[g[0]|g[1],g[0]|g[2],g[0]|g[3]]
 matches=[(i,canonical(s,universe)) for i,s in enumerate(pairings) if canonical(s,universe) in edges]
 if len(matches)>1:raise ValueError('Incompatible splits coexist in gene tree')
 if not matches:return 'other_or_unresolved',False
 i,key=matches[0];return ['concordant','alternative_1','alternative_2'][i],edges[key]

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--reference-tree',type=Path,required=True);parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
 R=args.root.resolve();manifest=json.loads((R/'metadata/LOCUS_MANIFEST.json').read_text());ref=args.reference_tree;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();assert sha(ref)==manifest['source_tree_sha256'];U,branches=reference_edges(Phylo.read(ref,'newick'));counts={s:Counter() for s in branches};details=[];trees=[]
 if args.pilot:manifest['loci']=[r for r in manifest['loci'] if r['gene'] in manifest['pilot_genes']]
 for row in manifest['loci']:
  g=row['gene'];p=R/'gene_trees'/g;receipt=json.loads((p/'RECEIPT.json').read_text());treepath=p/'gene.treefile';assert receipt['status']=='COMPLETE' and sha(treepath)==receipt['tree_sha256'] and receipt['input_sha256']==row['input_sha256'];t=Phylo.read(treepath,'newick');names,edges=gene_splits(t);assert names==set(row['retained_taxa']);trees.append(treepath.read_text().strip())
  for branch,groups in branches.items():
   result,high=classify(groups,names,edges);counts[branch][result]+=1
   if high:counts[branch]['high_'+result]+=1
   details.append({'gene':g,'reference_split':json.dumps(sorted(branch)),'classification':result,'resolving_edge_dual_high_support':high})
 out=R/('pilot_concordance_independent' if args.pilot else 'concordance_independent');out.mkdir(exist_ok=False);rows=[]
 for branch,c in sorted(counts.items(),key=lambda x:tuple(sorted(x[0]))):
  den=len(manifest['loci'])-c['not_decisive'];rows.append({'reference_split':json.dumps(sorted(branch)),'total_loci':len(manifest['loci']),'decisive_loci':den,'not_decisive':c['not_decisive'],'concordant':c['concordant'],'alternative_1':c['alternative_1'],'alternative_2':c['alternative_2'],'other_or_unresolved':c['other_or_unresolved'],'gCF_percent':100*c['concordant']/den if den else None,'high_concordant':c['high_concordant'],'high_alternative_1':c['high_alternative_1'],'high_alternative_2':c['high_alternative_2'],'decisive_without_dual_high_pairing':den-sum(c['high_'+s] for s in ['concordant','alternative_1','alternative_2'])})
  assert den==sum(c[s] for s in ['concordant','alternative_1','alternative_2','other_or_unresolved'])
 for name,rs in [('BRANCH_CONCORDANCE.tsv',rows),('LOCUS_BRANCH_CLASSIFICATIONS.tsv',details)]:
  with (out/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
 (out/'all_gene_trees.nwk').write_text('\n'.join(trees)+'\n');(out/'SUMMARY.json').write_text(json.dumps({'status':'INDEPENDENT_COUNTS_COMPLETE_IQTREE_CROSSCHECK_PENDING','loci':len(trees),'reference_internal_edges':len(branches),'reference_tree_sha256':sha(ref),'branches':rows,'limits':['Raw gCF includes estimated low-support gene resolutions','Dual-high resolved counts are descriptive and retain all decisive loci in denominator','Gene-tree errors, linkage and missing data prevent interpreting this as hybridization probability']},indent=2)+'\n');print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
