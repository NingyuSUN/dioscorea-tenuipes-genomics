#!/usr/bin/env python3
"""Stream frozen BED/BIM into WC84 components by physical region.

This computes descriptive sensitivity, never independent-block uncertainty.
Output excludes individual genotypes. Original files remain read-only.
"""
import argparse, csv, datetime, gzip, hashlib, itertools, json, sys, time
from pathlib import Path
import numpy as np
from influence_fst_core import population_counts, mask_from_counts, wc84_pair

PAIRS = ['AIC_HCO', 'AIC_IZU', 'HCO_IZU']
CORE_SHA = '58d033294bd33d0cfb4616eba7481e8c6dfbc685c6df75326248a91dd00661d8'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()

def decode(data, n):
    width = (n+3)//4
    packed = np.frombuffer(data, dtype=np.uint8).reshape(-1, width)
    return np.array([0,-1,1,2], dtype=np.int8)[(packed[:,np.arange(n)//4] >> (2*(np.arange(n)%4))) & 3]

def region_keys(contig_codes, positions, scale):
    if np.any(positions < 1): raise ValueError('Positions must be 1-based positive integers')
    bins = np.zeros(len(positions), dtype=np.int64) if scale == 0 else (positions-1)//scale
    if np.any(bins >= 2**32): raise ValueError('Coordinate exceeds key capacity')
    return (contig_codes.astype(np.int64) << 32) + bins

def component_matrix(genotypes, plan):
    """Columns: selected, then numerator/denominator/defined for each pair."""
    columns = {g:[plan['sample_order'].index(s) for s in ids] for g,ids in plan['populations'].items()}
    counts = {g:population_counts(genotypes, c) for g,c in columns.items()}
    called = sum(c.called for c in counts.values())
    alleles = sum(c.allele_count for c in counts.values())
    arrays = {}
    base = np.zeros((len(genotypes), 10), dtype=np.float64)
    for i, pair in enumerate(PAIRS):
        a,b = pair.split('_')
        comp = wc84_pair(counts[a], counts[b])
        if np.isinf(comp.fst).any(): raise ValueError('Infinite site FST')
        keep = ~np.isnan(comp.fst)
        if not np.isfinite(comp.numerator[keep]).all() or not np.isfinite(comp.denominator[keep]).all():
            raise ValueError('Nonfinite defined components')
        base[keep,1+3*i] = comp.numerator[keep]
        base[keep,2+3*i] = comp.denominator[keep]
        base[keep,3+3*i] = 1
    for name,p in plan['panels'].items():
        mask = mask_from_counts(called, alleles, len(plan['sample_order']), p['max_missing'], p['min_mac'])
        a = base.copy(); a[~mask] = 0; a[:,0] = mask
        arrays[name] = a
    return arrays

def group_components(keys, values):
    uniq, inverse = np.unique(keys, return_inverse=True)
    sums = np.column_stack([np.bincount(inverse, weights=values[:,j], minlength=len(uniq)) for j in range(values.shape[1])])
    return uniq, sums

def ratio_and_contrasts(values):
    values = np.asarray(values)
    den = values[..., [2,5,8]]
    num = values[..., [1,4,7]]
    fst = np.full(den.shape, np.nan)
    np.divide(num, den, out=fst, where=den>0)
    return fst, fst[..., :2]-fst[..., 2,None]

def deletion_stats(total, regions):
    return ratio_and_contrasts(np.asarray(total)-np.asarray(regions))

def additive_q(total, regions):
    d = np.asarray(total)[[2,5,8]]
    if np.any(d<=0): raise ValueError('Global denominator must be positive')
    return regions[:,[1,4]]/d[:2] - regions[:,8-1,None]/d[2]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--plan',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True); ap.add_argument('--chunk-sites',type=int,default=200000)
    args=ap.parse_args(); plan=json.loads(args.plan.read_text()); t0=time.monotonic()
    if not 1000<=args.chunk_sites<=1000000: raise ValueError('Invalid chunk size')
    if sha(Path(__file__).with_name('influence_fst_core.py')) != CORE_SHA: raise ValueError('FST core changed')
    files={e:Path(plan['input_prefix']+'.'+e) for e in ['bed','bim','fam']}
    hashes={e:sha(p) for e,p in files.items()}
    if hashes != plan['input_sha256']: raise ValueError('Input hash mismatch')
    samples=[line.split()[1] for line in files['fam'].read_text().splitlines()]
    if samples != plan['sample_order'] or len(set(samples))!=len(samples): raise ValueError('FAM order mismatch')
    flat=[s for ids in plan['populations'].values() for s in ids]
    if sorted(flat)!=sorted(samples): raise ValueError('Populations must partition samples')
    n=len(samples); width=(n+3)//4
    if files['bed'].stat().st_size!=3+width*plan['source_sites']: raise ValueError('BED size mismatch')
    args.out.mkdir(parents=True,exist_ok=False)
    scales=[0]+plan['scales_bp']; acc={(p,s):{} for p in plan['panels'] for s in scales}
    totals={p:np.zeros(10) for p in plan['panels']}
    contigs=[]; lookup={}; previous_contig=None; previous_pos=0; processed=0
    with files['bed'].open('rb') as bed,files['bim'].open() as bim:
        if bed.read(3)!=b'\x6c\x1b\x01': raise ValueError('Not SNP-major BED')
        while data:=bed.read(width*args.chunk_sites):
            g=decode(data,n); cid=np.empty(len(g),dtype=np.int64); pos=np.empty(len(g),dtype=np.int64)
            for k in range(len(g)):
                f=bim.readline().split()
                if len(f)!=6: raise ValueError('Truncated/malformed BIM')
                c=f[0]; v=int(f[3])
                if not f[1].startswith(c+':'+str(v)+':'): raise ValueError('Marker ID versus coordinate mismatch')
                if c!=previous_contig:
                    if c in lookup: raise ValueError('Disjoint contig segments')
                    lookup[c]=len(contigs);contigs.append(c);previous_contig=c;previous_pos=0
                if v<=previous_pos: raise ValueError('Non-increasing or duplicated coordinate')
                previous_pos=v;cid[k]=lookup[c];pos[k]=v
            # All retained panels are subsets of main's calls>=13/MAC>=1.
            called=(g>=0).sum(1); allele=np.maximum(g,0).sum(1)
            possible=mask_from_counts(called,allele,n,plan['panels']['main']['max_missing'],1)
            g=g[possible];cid=cid[possible];pos=pos[possible]
            arrays=component_matrix(g,plan)
            for p,values in arrays.items():
                totals[p]+=values.sum(0)
                for scale in scales:
                    keys,sums=group_components(region_keys(cid,pos,scale),values)
                    target=acc[(p,scale)]
                    for key,val in zip(keys,sums):
                        if val[0]==0:continue
                        key=int(key)
                        if key in target:target[key]+=val
                        else:target[key]=val.copy()
            processed+=len(data)//width
            if processed%(args.chunk_sites*5)==0:
                print(json.dumps({'sites_processed':processed,'total_sites':plan['source_sites'],'elapsed_s':round(time.monotonic()-t0,2)}),flush=True)
        if bim.readline(): raise ValueError('BIM longer than BED')
    if processed!=plan['source_sites']:raise ValueError('Source count mismatch')
    expected={'main':[.3312913257919642,.332169705556553,.195661227449393],
              'complete':[.33216814878392154,.3335648649456691,.19871837422413705]}
    expected_defined={'main':[4538382,4776008,4300810],'complete':[2756842,2897559,2609873]}
    rows=[];checks=[]
    columns=['panel','scale_bp','contig','window_index','nominal_start_1based','nominal_end_1based','selected_sites']
    for pair in PAIRS:columns.extend([pair+'_numerator',pair+'_denominator',pair+'_defined_sites'])
    for p in plan['panels']:
        total=totals[p];fst,contrast=ratio_and_contrasts(total)
        assert int(total[0])==plan['panels'][p]['expected_sites']
        assert np.array_equal(total[[3,6,9]],expected_defined[p])
        assert np.allclose(fst,expected[p],atol=1e-11,rtol=0)
        rows.append({'panel':p,'selected_sites':int(total[0]),'components':total.tolist(),'fst':fst.tolist(),'contrasts':contrast.tolist()})
        for scale in scales:
            values=np.array(list(acc[(p,scale)].values()))
            error=np.abs(values.sum(0)-total)
            assert np.allclose(values.sum(0),total,rtol=1e-11,atol=1e-7)
            checks.append({'panel':p,'scale_bp':scale,'occupied_regions':len(values),'maximum_sum_error':float(error.max())})
    with gzip.open(args.out/'REGION_COMPONENTS.tsv.gz','wt') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(columns)
        for (p,scale),regions in acc.items():
            for key,val in sorted(regions.items()):
                c=contigs[key>>32];idx=key&((1<<32)-1)
                bounds=['',''] if scale==0 else [idx*scale+1,(idx+1)*scale]
                w.writerow([p,scale,c,idx,*bounds,*val.tolist()])
    after={e:sha(p) for e,p in files.items()}
    if after!=hashes:raise ValueError('Input changed during run')
    receipt={'status':'COMPLETE','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-t0,'input_sha256':hashes,'input_unchanged':True,'plan_sha256':sha(args.plan),
        'script_sha256':sha(__file__),'core_sha256':CORE_SHA,'python':sys.version,'numpy':np.__version__,
        'source_sites':processed,'source_contigs':len(contigs),'sample_order':samples,'global_results':rows,
        'sum_checks':checks,'region_components_sha256':sha(args.out/'REGION_COMPONENTS.tsv.gz'),
        'coordinate_contract':'All BIM marker IDs agree with CHROM/POS; strictly increasing positions in contiguous original-contig segments. Nominal window ends may exceed unknown terminal contig lengths.'}
    (args.out/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'COMPLETE','elapsed_seconds':receipt['elapsed_seconds'],'global_results':rows}),flush=True)

if __name__=='__main__':main()
