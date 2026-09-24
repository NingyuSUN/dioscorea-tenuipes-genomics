#!/usr/bin/env python3
"""Bounded chunked WC84 sensitivity analysis from an unfiltered masked BED.

PCA and FST use different site filters. This program uses MAC>=1 and never
applies spacing or LD pruning. FIXED retains the baseline all-14 MAC1/missing
site set after removing individuals. REFILTER recomputes MAC and missingness
using retained individuals. COMPLETE requires no missing genotypes.

Input JSON: sample_order, populations (group -> IDs), panels (objects with
name, mode fixed/refilter/complete, samples list); optional expected_baseline_sites
and expected_complete_sites. No genotype masking is performed here: input must
be derived from the already DP/GQ-masked biallelic diploid source, without site
MAC/missing filtering before BED import. BED dosage counts allele A2, which
is FST-equivalent to consistently reversing the allele at any whole site.
"""
import argparse,csv,hashlib,itertools,json,sys,time
from pathlib import Path
import numpy as np
from influence_fst_core import (population_counts,subtract_counts,mask_from_counts,
                                wc84_pair,WeightedFstAccumulator,validate_genotypes)


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def bed_chunks(prefix,n_samples,n_sites,chunk_size):
    width=(n_samples+3)//4
    path=Path(str(prefix)+'.bed')
    if path.stat().st_size!=3+width*n_sites:
        raise ValueError('BED size disagrees with FAM sample count and BIM site count')
    lookup=np.array([0,-1,1,2],dtype=np.int8)
    cols=np.arange(n_samples)//4; shifts=2*(np.arange(n_samples)%4)
    with path.open('rb') as f:
        if f.read(3)!=b'\x6c\x1b\x01':raise ValueError('Require SNP-major PLINK BED')
        while data:=f.read(width*chunk_size):
            packed=np.frombuffer(data,dtype=np.uint8).reshape(-1,width)
            dosage=lookup[(packed[:,cols]>>shifts)&3]
            yield dosage


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bed-prefix',type=Path,required=True)
    ap.add_argument('--panels',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--chunk-sites',type=int,default=200000)
    a=ap.parse_args(); t0=time.monotonic()
    if not 1000<=a.chunk_sites<=1000000:raise ValueError('Chunk must be1000..1000000sites')
    a.out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(a.panels.read_text())
    files={ext:Path(str(a.bed_prefix)+'.'+ext) for ext in ['bed','bim','fam']}
    hashes={ext:sha256(path) for ext,path in files.items()}
    samples=[line.split()[1] for line in files['fam'].read_text().splitlines() if line.strip()]
    if samples!=cfg['sample_order'] or len(samples)!=len(set(samples)):raise ValueError('FAM sample order mismatch')
    sample_idx={s:i for i,s in enumerate(samples)}
    n_sites=sum(1 for line in files['bim'].open() if line.strip())
    groups=list(cfg['populations']); pairs=list(itertools.combinations(groups,2))
    popcols={g:[sample_idx[s] for s in ids] for g,ids in cfg['populations'].items()}
    flat=[i for cols in popcols.values() for i in cols]
    if sorted(flat)!=list(range(len(samples))):raise ValueError('Populations must partition all input samples')
    panels=cfg['panels']; names=[p['name'] for p in panels]
    if len(names)!=len(set(names)):raise ValueError('Duplicate panel name')
    retain_groups={}
    for p in panels:
        if p['mode'] not in ('fixed','refilter','complete'):raise ValueError('Unsupported site mode')
        if len(p['samples'])!=len(set(p['samples'])):raise ValueError('Duplicate retained sample')
        keep=tuple(sorted(sample_idx[s] for s in p['samples']))
        for g in groups:
            if len(set(keep)&set(popcols[g]))<2:raise ValueError('Each population requires at least two retained individuals')
        retain_groups.setdefault(keep,[]).append(p)
    min_n=min(len(keep) for keep in retain_groups)
    # Any admissible site must be polymorphic overall and have at least this
    # many called samples. This only discards sites impossible in every panel.
    min_calls=int(np.ceil(.9*min_n-1e-12))
    acc={(p['name'],pair):WeightedFstAccumulator() for p in panels for pair in pairs}
    processed=potential_sites=baseline_sites=complete_sites=0
    for chunk_i,g0 in enumerate(bed_chunks(a.bed_prefix,len(samples),n_sites,a.chunk_sites),1):
        validate_genotypes(g0);processed+=len(g0)
        all0=population_counts(g0,range(len(samples)))
        possible=(all0.called>=min_calls)&(np.minimum(all0.allele_count,2*all0.called-all0.allele_count)>=1)
        g=g0[possible];potential_sites+=len(g)
        if len(g)==0:continue
        base_counts={pop:population_counts(g,cols) for pop,cols in popcols.items()}
        called=sum(c.called for c in base_counts.values());allele_count=sum(c.allele_count for c in base_counts.values())
        base_mask=mask_from_counts(called,allele_count,len(samples),.1,1)
        comp_mask=mask_from_counts(called,allele_count,len(samples),0,1)
        baseline_sites+=int(base_mask.sum());complete_sites+=int(comp_mask.sum())
        for keep,subpanels in retain_groups.items():
            keep_set=set(keep)
            counts={pop:subtract_counts(base_counts[pop],g,[i for i in cols if i not in keep_set]) for pop,cols in popcols.items()}
            subcalled=sum(c.called for c in counts.values());suballele=sum(c.allele_count for c in counts.values())
            selections={'fixed':base_mask,
                        'refilter':mask_from_counts(subcalled,suballele,len(keep),.1,1),
                        'complete':mask_from_counts(subcalled,suballele,len(keep),0,1)}
            for pair in pairs:
                components=wc84_pair(counts[pair[0]],counts[pair[1]])
                for p in subpanels:acc[(p['name'],pair)].update(components,selections[p['mode']])
        if chunk_i%5==0: print(f'processed={processed}/{n_sites}; potential={potential_sites}; elapsed_s={time.monotonic()-t0:.1f}',flush=True)
    if processed!=n_sites:raise AssertionError('Source site count mismatch')
    for field,value in [('expected_baseline_sites',baseline_sites),('expected_complete_sites',complete_sites)]:
        if field in cfg and cfg[field]!=value:raise AssertionError(f'{field}: expected {cfg[field]}, observed {value}')
    if any(x.infinite_sites for x in acc.values()):raise AssertionError('Infinite site FST present; inspect before interpreting')
    rows=[]
    for p in panels:
        for pair in pairs:
            rows.append(dict(panel=p['name'],pair='_'.join(pair),mode=p['mode'],**acc[(p['name'],pair)].as_dict()))
    fields=['panel','pair','weighted_fst','numerator','denominator','selected_sites','defined_sites','mode','undefined_sites','mean_fst','negative_sites','infinite_sites']
    with (a.out/'fst_summary.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    after={ext:sha256(path) for ext,path in files.items()}
    if after!=hashes:raise AssertionError('Input file changed during run')
    audit=dict(status='COMPLETE',source_prefix=str(a.bed_prefix),source_sha256=hashes,
        input_files_unchanged=True,panels_sha256=sha256(a.panels),sample_order=samples,
        raw_sites=n_sites,potential_sites=potential_sites,baseline_sites=baseline_sites,
        complete_sites=complete_sites,panels=len(panels),rows=len(rows),chunk_sites=a.chunk_sites,
        elapsed_seconds=time.monotonic()-t0,python=sys.version,numpy=np.__version__,
        script_sha256=sha256(__file__),core_sha256=sha256(Path(__file__).with_name('influence_fst_core.py')),
        estimator='VCFtools v0.1.16 WC84 global ratio of summed both-allele variance components; no clipping',
        source_contract='upstream DP/GQ masked; diploid biallelic; missing -1; dosage A2; source unfiltered for MAC/missing',
        fixed_site_policy='all input individuals: MAC>=1, missing<=0.1; retain positions after removal',
        refilter_policy='retained individuals: MAC>=1, missing<=0.1',
        complete_policy='retained individuals: MAC>=1, zero missing')
    (a.out/'fst_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    (a.out/'FST_COMPLETE').write_text('All requested FST estimates completed; expected counts and input hashes verified.\n')
    print(json.dumps(audit),flush=True)

if __name__=='__main__':main()
