"""Compare integer per-branch counts with IQ-TREE, allowing the two NNI alternatives to swap labels."""
import argparse,csv,json
from pathlib import Path
from Bio import Phylo
from gene_concordance import canonical
p=argparse.ArgumentParser();p.add_argument('--independent',type=Path,required=True);p.add_argument('--iqtree-prefix',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();ind=json.loads((a.independent/'SUMMARY.json').read_text());by={frozenset(json.loads(r['reference_split'])):r for r in ind['branches']};tree=Phylo.read(str(a.iqtree_prefix)+'.cf.branch','newick');U=frozenset(t.name for t in tree.get_terminals());ids={}
for c in tree.get_nonterminals():
 if c is tree.root:continue
 raw=c.name if c.name is not None else c.confidence;side=frozenset(t.name for t in c.get_terminals())
 if min(len(side),len(U-side))<2:continue
 ids[str(int(float(raw)))]=canonical(side,U)
lines=[s for s in Path(str(a.iqtree_prefix)+'.cf.stat').read_text().splitlines() if not s.startswith('#')];rows=list(csv.DictReader(lines,delimiter='\t'));assert set(ids)=={r['ID'] for r in rows};assert set(ids.values())==set(by);checks=[]
for r in rows:
 ours=by[ids[r['ID']]];n=ours['decisive_loci'];fields=['gCF_N','gDF1_N','gDF2_N','gDFP_N','gN']
 if n==0:assert all(r[k]=='NA' for k in fields),r;checks.append({'branch':r['ID'],'decisive_loci':0,'status':'PASS_NA_DENOMINATOR'});continue
 assert int(r['gN'])==n
 assert int(r['gCF_N'])==ours['concordant']
 assert int(r['gDFP_N'])==ours['other_or_unresolved']
 assert sorted([int(r['gDF1_N']),int(r['gDF2_N'])])==sorted([ours['alternative_1'],ours['alternative_2']])
 assert abs(float(r['gCF'])-ours['gCF_percent'])<=.00501
 checks.append({'branch':r['ID'],'reference_split':json.loads(ours['reference_split']),'decisive_loci':n,'status':'PASS_INTEGER_COUNTS_AND_ROUNDED_PERCENT'})
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps({'status':'PASS_INDEPENDENT_IQTREE_CROSSCHECK','loci':ind['loci'],'branches_checked':len(checks),'checks':checks,'limits':'Two implementations agree on estimated gene trees; biological and gene-tree-error interpretations remain separate.'},indent=2)+'\n');print(a.out.read_text())
