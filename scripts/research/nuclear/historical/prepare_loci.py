"""Preserve all retained nuclear loci; omit only completely missing taxon rows per locus."""
import hashlib,json,subprocess
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
from Bio import SeqIO
R=Path(__file__).resolve().parents[1];SRC=Path('workspace/wgs/dt_research_p18_20260928');D=SRC/'matrices/genome_main_no2661';AA=set('ACDEFGHIKLMNPQRSTVWY');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();S=json.loads((D/'summary.json').read_text());assert S['genes']==173 and S['columns']==68796 and S['main_sha256']==sha(D/'main.faa');assert '2661at3193' not in {p['gene'] for p in S['partitions']};assert sha(SRC/'results/genome_main_no2661/main.treefile')=='276c1f42bfabdfe3965eb87c131a066a90e19a3141c3288428aa6ce77d3d98f0'
assert not (R/'metadata/LOCUS_MANIFEST.json').exists();(R/'loci').mkdir(exist_ok=True);whole={s.id:str(s.seq) for s in SeqIO.parse(D/'main.faa','fasta')};rows=[]
for i,p in enumerate(S['partitions']):
 g=p['gene'];f=D/'trimmed'/(g+'.faa');assert sha(f)==p['sha256'];recs=list(SeqIO.parse(f,'fasta'));assert len(recs)==len({r.id for r in recs});seq={r.id:str(r.seq) for r in recs};assert set(seq)==set(whole);assert all(seq[t]==whole[t][p['start_1based']-1:p['end_1based']] for t in whole)
 keep={t:s for t,s in seq.items() if any(c in AA for c in s)};removed=sorted(set(seq)-set(keep));assert len(keep)>=4;assert {'Dioscorea_tenuipes_historical','Dioscorea_gracillima'}<=set(keep);assert all(set(seq[t])<=set('-?X') for t in removed)
 dest=R/'loci'/(g+'.faa');dest.write_text(''.join('>'+t+'\n'+s+'\n' for t,s in keep.items()));info=0;variable=0
 for chars in zip(*keep.values()):
  cnt=Counter(c for c in chars if c in AA);variable+=len(cnt)>1;info+=sum(n>=2 for n in cnt.values())>=2
 rows.append({'gene':g,'length':p['length'],'source_sha256':sha(f),'input_sha256':sha(dest),'retained_taxa':sorted(keep),'all_missing_taxa_omitted':removed,'variable_columns':variable,'parsimony_informative_columns':info,'seed':20260928+i})
ordered=sorted(rows,key=lambda p:(p['length'],p['gene']));pilot=[ordered[0]['gene'],ordered[len(ordered)//2]['gene'],ordered[-1]['gene']];out={'utc':datetime.now(timezone.utc).isoformat(),'status':'FROZEN_BEFORE_LOCUS_INFERENCE','loci':rows,'locus_count':173,'columns':68796,'source_matrix_sha256':sha(D/'main.faa'),'source_tree_sha256':sha(SRC/'results/genome_main_no2661/main.treefile'),'excluded_orthogroups':['2661at3193'],'pilot_genes':pilot,'pilot_selection':'minimum, median and maximum aligned length; no topology-based selection','command_policy':{'model_selection':'MFP among JTT,WAG,LG with F and G4','bootstrap':1000,'bnni':True,'SH_aLRT':1000,'threads_per_locus':1,'parallel_loci':3},'matrix_slice_checks':173*len(whole),'limits':['Same loci and sequence evidence as concatenated tree; not independent biological replication','Missing rows omitted, residues and columns unchanged','Gene-tree error and short loci can lower concordance']};(R/'metadata/LOCUS_MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='loci'},indent=2))
