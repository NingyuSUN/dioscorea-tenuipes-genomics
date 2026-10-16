#!/usr/bin/env python3
"""Two estimators on frozen original14 SNP sets and all fifteen IZU4 subsets."""
import argparse,csv,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from influence_fst_core import population_counts,mask_from_counts,wc84_pair

PAIRS=[('AIC','HCO'),('AIC','IZU'),('HCO','IZU')]
CORE_SHA='58d033294bd33d0cfb4616eba7481e8c6dfbc685c6df75326248a91dd00661d8'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def decode(data,n):
    a=np.frombuffer(data,dtype=np.uint8).reshape(-1,(n+3)//4)
    return np.array([0,-1,1,2],dtype=np.int8)[(a[:,np.arange(n)//4]>>(2*(np.arange(n)%4)))&3]

def hudson_pair(x,y):
    """n counts called chromosome copies, not diploid people; no ratio averaging."""
    n1=2*x.called.astype(np.float64);n2=2*y.called.astype(np.float64)
    with np.errstate(divide='ignore',invalid='ignore'):
        p1=x.allele_count/n1;p2=y.allele_count/n2
        numerator=(p1-p2)**2-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1)
        denominator=p1*(1-p2)+p2*(1-p1)
    valid=(n1>1)&(n2>1)&np.isfinite(numerator)&np.isfinite(denominator)&(denominator>0)
    return numerator,denominator,valid

def components(g,plan):
    """Yield config,pair,two estimators and keep per-site definition explicit."""
    samples=plan['sample_order'];indices={s:i for i,s in enumerate(samples)}
    base={group:population_counts(g,[indices[s] for s in plan['populations'][group]]) for group in ['AIC','HCO','IZU']}
    cache={}
    for cfg in plan['configurations']:
        counts=dict(base)
        if cfg['name']!='original_4_4_6':counts['IZU']=population_counts(g,[indices[s] for s in cfg['IZU']])
        for pair in PAIRS:
            if pair==('AIC','HCO') and pair in cache:
                yield cfg,pair,cache[pair];continue
            x,y=[counts[p] for p in pair];wc=wc84_pair(x,y)
            if np.isinf(wc.fst).any():raise ValueError('Infinite WC site FST')
            wvalid=np.isfinite(wc.fst)&np.isfinite(wc.numerator)&(wc.denominator>0)
            ha,hb,hvalid=hudson_pair(x,y)
            if not np.array_equal(wvalid,hvalid):raise ValueError('Estimator-defined masks differ; review before interpretation')
            values={'WC84':(wc.numerator,wc.denominator,wvalid),'Hudson':(ha,hb,hvalid)}
            if pair==('AIC','HCO'):cache[pair]=values
            yield cfg,pair,values

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--chunk-sites',type=int,default=100000);args=ap.parse_args();t0=time.monotonic()
    p=json.loads(args.plan.read_text());files={e:Path(p['source_prefix']+'.'+e) for e in ['bed','bim','fam']}
    if sha(Path(__file__).with_name('influence_fst_core.py'))!=CORE_SHA:raise ValueError('Frozen WC core mismatch')
    before={e:sha(f) for e,f in files.items()}
    if before!=p['source_sha256']:raise ValueError('Source changed')
    samples=[line.split()[1] for line in files['fam'].read_text().splitlines()]
    if samples!=p['sample_order'] or len(set(samples))!=len(samples):raise ValueError('Samples/order mismatch')
    n=len(samples);width=(n+3)//4
    if files['bed'].stat().st_size!=3+width*p['source_sites']:raise ValueError('Source BED dimensions')
    if not 1000<=args.chunk_sites<=500000:raise ValueError('Chunk out of bounds')
    args.out.mkdir(parents=True,exist_ok=False)
    # entries: A, B, defined, negative numerator, selected positions
    sums={(cfg['name'],panel,'_'.join(pair),method):np.zeros(5) for cfg in p['configurations'] for panel in p['panels'] for pair in PAIRS for method in ['WC84','Hudson']}
    selected={panel:0 for panel in p['panels']};processed=0
    with files['bed'].open('rb') as f:
        if f.read(3)!=b'\x6c\x1b\x01':raise ValueError('Not SNP-major BED')
        while data:=f.read(width*args.chunk_sites):
            g0=decode(data,n);processed+=len(g0)
            allcounts=population_counts(g0,range(n))
            mainmask=mask_from_counts(allcounts.called,allcounts.allele_count,n,p['panels']['main']['max_missing'],1)
            g=g0[mainmask];called=allcounts.called[mainmask];allele=allcounts.allele_count[mainmask]
            masks={panel:mask_from_counts(called,allele,n,v['max_missing'],v['min_mac']) for panel,v in p['panels'].items()}
            for panel,m in masks.items():selected[panel]+=int(m.sum())
            for cfg,pair,methods in components(g,p):
                for method,(num,den,valid) in methods.items():
                    for panel,mask in masks.items():
                        use=mask&valid;a=sums[(cfg['name'],panel,'_'.join(pair),method)]
                        a+=np.array([num[use].sum(dtype=np.float64),den[use].sum(dtype=np.float64),int(use.sum()),int((num[use]<0).sum()),int(mask.sum())])
            if processed%(args.chunk_sites*10)==0:print(json.dumps({'processed':processed,'total':p['source_sites'],'elapsed_s':round(time.monotonic()-t0,2)}),flush=True)
    assert processed==p['source_sites']
    for panel,v in p['panels'].items():assert selected[panel]==v['expected_sites']
    rows=[]
    for (cfg,panel,pair,method),a in sums.items():
        if not np.isfinite(a[:2]).all() or a[1]<=0:raise ValueError('Global not estimable')
        rows.append({'configuration':cfg,'panel':panel,'pair':pair,'estimator':method,'fst':float(a[0]/a[1]),'numerator':float(a[0]),'denominator':float(a[1]),'defined_sites':int(a[2]),'negative_numerator_sites':int(a[3]),'selected_sites':int(a[4]),'undefined_sites':int(a[4]-a[2])})
    expected={'main':[.3312913257919642,.332169705556553,.195661227449393],'complete':[.33216814878392154,.3335648649456691,.19871837422413705]}
    defined={'main':[4538382,4776008,4300810],'complete':[2756842,2897559,2609873]}
    for panel in p['panels']:
        for i,pair in enumerate(PAIRS):
            a=sums[('original_4_4_6',panel,'_'.join(pair),'WC84')]
            assert abs(a[0]/a[1]-expected[panel][i])<1e-11
            assert int(a[2])==defined[panel][i]
        for cfg in p['configurations']:
            for method in ['WC84','Hudson']:
                assert np.array_equal(sums[(cfg['name'],panel,'AIC_HCO',method)],sums[('original_4_4_6',panel,'AIC_HCO',method)])
    with (args.out/'FST_ESTIMATORS.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    after={e:sha(f) for e,f in files.items()}
    if after!=before:raise ValueError('Source mutated')
    receipt={'status':'COMPLETE','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-t0,
       'source_sha256':before,'source_unchanged':True,'sample_order':samples,'source_sites':processed,'selected_sites':selected,
       'configurations':p['configurations'],'rows':len(rows),'defined_masks_identical':True,'all14_WC_reproduced':True,'AIC_HCO_constant_across_subsets':True,
       'script_sha256':sha(__file__),'core_sha256':CORE_SHA,'plan_sha256':sha(args.plan),'results_sha256':sha(args.out/'FST_ESTIMATORS.tsv'),
       'python':sys.version,'numpy':np.__version__,'limitations':p['limitations']}
    (args.out/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'COMPLETE','elapsed_seconds':receipt['elapsed_seconds'],'rows':len(rows),'selected_sites':selected}),flush=True)

if __name__=='__main__':main()
