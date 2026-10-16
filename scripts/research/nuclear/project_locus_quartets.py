"""Fixed quartet projection of existing gene bootstrap trees; no new tree inference."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import argparse,hashlib,json,re,math,time,os,traceback
R=Path(__file__).resolve().parents[1]
SRC=Path('workspace/wgs/dt_research_p24_20260928')
TAXA=['Dioscorea_tenuipes_historical','Dioscorea_tokoro','Dioscorea_gracillima','Dioscorea_zingiberensis']
LABEL={3:'DT_TOK',5:'DT_GRA',6:'DT_ZIN'}
KEYS=list(LABEL.values())+['unresolved']
AA=set('ACDEFGHIKLMNPQRSTVWY')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def project(s,expected):
    tok=re.findall(r'[(),:;]|[^(),:;\s]+',s);i=0;leaves=[];splits=set();edge_lengths=[]
    def node():
        nonlocal i
        if tok[i]=='(':
            i+=1;bits=node()
            while tok[i]==',':i+=1;bits|=node()
            assert tok[i]==')';i+=1
            if tok[i] not in [',',')',':',';']:i+=1
        else:
            name=tok[i].strip("'");i+=1;leaves.append(name)
            bits=(1<<TAXA.index(name)) if name in TAXA else 0
        length=None
        if tok[i]==':':
            i+=1;length=float(tok[i]);i+=1
            assert math.isfinite(length) and length>=0
        if bits.bit_count()==2:
            splits.add(LABEL[min(bits,15^bits)]);edge_lengths.append(length)
        return bits
    assert node()==15
    assert tok[i]==';' and i==len(tok)-1
    assert len(leaves)==len(set(leaves)) and set(leaves)==set(expected)
    assert len(splits)<=1,splits
    return next(iter(splits)) if splits else 'unresolved',edge_lengths

def signal(p):
    seq={};name=None
    for line in p.read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];seq[name]=''
        elif line.strip():seq[name]+=line.strip()
    assert set(TAXA)<=set(seq)
    lengths={len(seq[t]) for t in TAXA};assert len(lengths)==1
    counts=Counter({k:0 for k in LABEL.values()});common=variable=0
    for column in zip(*(seq[t] for t in TAXA)):
        if not all(c in AA for c in column):continue
        common+=1;c=Counter(column)
        if len(c)>1:variable+=1
        if sorted(c.values())==[2,2]:
            mask=sum(1<<j for j,x in enumerate(column) if x==column[0]);counts[LABEL[min(mask,15^mask)]]+=1
    return {'alignment_columns':next(iter(lengths)),'common_canonical_columns':common,'variable_columns':variable,'informative_2_2':dict(counts)}

def fixtures():
    d,t,g,z=TAXA
    cases=[(f'(({d},{t}),({g},{z}));','DT_TOK'),(f'({d},({g},({t},{z})));','DT_GRA'),(f'(({d},{z}),({t},{g}));','DT_ZIN'),(f'({d},{t},{g},{z});','unresolved'),(f'((({d}:0,{t}:0)99/100:0,X:1):1,({g}:1,{z}:1):1);','DT_TOK'),(f'({d}:1,(({g}:1,{t}:1):1,{z}:1):1);','DT_ZIN')]
    for tree,target in cases:
        expected=TAXA+(['X'] if 'X:' in tree else [])
        actual,_=project(tree,expected);assert actual==target,(tree,target,actual)
    return len(cases)

def main():
    start=time.time();out=R/'results';val=R/'validation';out.mkdir(exist_ok=True);val.mkdir(exist_ok=True)
    dest=out/'NUCLEAR_QUARTETS.json';assert not dest.exists(),'Existing output; do not overwrite or rerun blindly'
    plan=json.loads((R/'PLAN.json').read_text());assert plan['fixed_taxa']==TAXA
    for name,h in plan['input_hashes'].items():assert sha(R/'inputs'/name)==h,name
    manifest=json.loads((R/'inputs/LOCUS_MANIFEST.json').read_text());aud=json.loads((R/'inputs/ALL_LOCUS_AUDIT.json').read_text())
    checks=aud['records']
    original={x['gene']:x for x in checks};assert len(original)==173
    warnings=json.loads((R/'inputs/WARNING_SENSITIVITY.json').read_text())['excluded_loci']
    refs=(R/'inputs/all_gene_trees.nwk').read_text().splitlines();assert len(refs)==173
    state={'status':'RUNNING','pid':os.getpid(),'started_utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':sha(R/'PLAN.json'),'completed_loci':0,'new_tree_inferences':False}
    status=R/'RUN_STATUS.json'
    status.write_text(json.dumps(state,indent=2)+'\n')
    rows=[]
    for idx,m in enumerate(manifest['loci']):
        gene=m['gene'];assert gene!='2661at3193';p=SRC/'gene_trees'/gene
        tree=p/'gene.treefile';boots=p/'gene.ufboot';fasta=SRC/'loci'/(gene+'.faa')
        assert sha(tree)==original[gene]['tree_sha256'] and sha(boots)==original[gene]['bootstrap_sha256'],gene
        assert sha(fasta)==m['input_sha256'],gene
        assert refs[idx].strip()==tree.read_text().strip(),gene
        absent=sorted(set(TAXA)-set(m['retained_taxa']));eligible=not absent
        row={'gene':gene,'eligible':eligible,'expected_taxa':m['retained_taxa'],'missing_taxa':absent,'point_split':'missing_taxa','bootstrap_counts':None,'log_warning':gene in warnings,'source_sha256':{'tree':sha(tree),'bootstrap':sha(boots),'alignment':sha(fasta)}}
        if eligible:
            row['point_split'],edges=project(tree.read_text(),m['retained_taxa']);row['projecting_point_edge_lengths']=edges;counts=Counter({k:0 for k in KEYS});n=0
            for line in boots.read_text().splitlines():
                if not line.strip():continue
                split,_=project(line,m['retained_taxa']);counts[split]+=1;n+=1
            assert n==1000
            row['bootstrap_counts']=dict(counts);row['signal']=signal(fasta)
        rows.append(row);state['completed_loci']=idx+1;state['elapsed_seconds']=time.time()-start;status.write_text(json.dumps(state,indent=2)+'\n')
    summary={'status':'PROJECTED_AUDIT_PENDING','fixed_taxa':TAXA,'topologies':plan['outcomes'],'loci':rows,'fixtures_passed':fixtures(),'elapsed_seconds':time.time()-start,'locus_count':len(rows),'eligible_loci':sum(x['eligible'] for x in rows),'limits':plan['limits']}
    dest.write_text(json.dumps(summary,indent=2)+'\n')
    state.update(status='COMPLETE_INDEPENDENT_AUDIT_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),result_sha256=sha(dest));status.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='loci'},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--source-root',type=Path,required=True)
    args=parser.parse_args()
    R=args.root.resolve();SRC=args.source_root.resolve()
    fixtures()
    try:main()
    except Exception:
        (R/'FAILED.txt').write_text(traceback.format_exc());raise
