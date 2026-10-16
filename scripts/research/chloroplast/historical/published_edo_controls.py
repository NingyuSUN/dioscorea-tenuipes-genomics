#!/usr/bin/env python3
"""Controls for published encoded calls; no genotype re-calling or ancestry inference."""
from pathlib import Path
import collections, datetime, gzip, hashlib, json, resource, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'inputs/published/Fig_4_nighbor-net__vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.gt.gz'
OUT=ROOT/'published_controls'
N=1054866
NAMES=['dt_kita1','dt_waka1','dt_mzawa','ed_hachi1','tn_yuga','tn_utsu','tn_mie1','dq_mie']

def run():
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    started=time.monotonic()
    previous=json.loads((ROOT/'published_comparison/RECEIPT.json').read_text())
    digest=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert digest==previous['GT_file_sha256']
    OUT.mkdir(exist_ok=False)
    A=np.empty((N,8),dtype=np.int8)
    with gzip.open(SOURCE,'rt') as f:
        for i,line in enumerate(f):
            values=line.split(); assert len(values)==10 and i<N
            vals=list(map(int,values[2:])); assert set(vals)<=set([0,1,2])
            A[i]=vals
    assert i+1==N
    all_het=np.all(A==1,axis=1); all_alt=np.all(A==2,axis=1)
    assert int(all_het.sum())==32 and int(all_alt.sum())==566
    keep=~(all_het|all_alt)
    per_sample=[]
    for j,name in enumerate(NAMES):
        for scope,mask in [('all_published_sites',np.ones(N,dtype=bool)),('exclude_all_het_and_all_alt',keep)]:
            c=np.bincount(A[mask,j],minlength=3)
            per_sample.append(dict(sample=name,scope=scope,denominator=int(mask.sum()),homozygous_reference_code=int(c[0]),heterozygous_code=int(c[1]),homozygous_alternate_code=int(c[2])))
    controls=[]
    independent=[]
    for omitted in [None,0,1,2,4,5,6]:
        T=[x for x in [0,1,2] if x!=omitted]; U=[x for x in [4,5,6] if x!=omitted]
        mask=np.all(A[:,T]==A[:,[T[0]]],axis=1)&np.isin(A[:,T[0]],[0,2])&np.all(A[:,U]==A[:,[U[0]]],axis=1)&np.isin(A[:,U[0]],[0,2])&(A[:,T[0]]!=A[:,U[0]])
        assert not np.any(mask&~keep)
        row={'omitted_comparator':None if omitted is None else NAMES[omitted],'defining_tokoro':[NAMES[x] for x in T],'defining_tenuipes':[NAMES[x] for x in U],'diagnostic_sites':int(mask.sum()),'selected_without_Hachi':True,'exclude_all_het_and_all_alt_changes_sites':False,'calls':[]}
        for j in ([3] if omitted is None else [omitted,3]):
            counts={'same_tokoro_homozygous_code':int((A[mask,j]==A[mask,T[0]]).sum()),'heterozygous_code':int((A[mask,j]==1).sum()),'same_tenuipes_homozygous_code':int((A[mask,j]==A[mask,U[0]]).sum())}
            assert sum(counts.values())==row['diagnostic_sites']
            row['calls'].append({'sample':NAMES[j],**counts})
        controls.append(row)
        independent.append(collections.Counter())
    # Independently rescan original rows with scalar logic, not the array masks.
    with gzip.open(SOURCE,'rt') as f:
        for line in f:
            g=list(map(int,line.split()[2:]))
            for k,omitted in enumerate([None,0,1,2,4,5,6]):
                T=[x for x in [0,1,2] if x!=omitted]; U=[x for x in [4,5,6] if x!=omitted]
                ts={g[x] for x in T}; us={g[x] for x in U}
                if len(ts)==len(us)==1 and ts<=set([0,2]) and us<=set([0,2]) and ts!=us:
                    independent[k]['n']+=1
                    for j in ([3] if omitted is None else [omitted,3]):
                        label='heterozygous_code' if g[j]==1 else 'same_tokoro_homozygous_code' if g[j] in ts else 'same_tenuipes_homozygous_code'
                        independent[k][(j,label)]+=1
    for row,ctr in zip(controls,independent):
        assert ctr['n']==row['diagnostic_sites']
        for calls in row['calls']:
            j=NAMES.index(calls['sample'])
            for key,value in calls.items():
                if key!='sample': assert ctr[(j,key)]==value
    def save(name,value): (OUT/name).write_text(json.dumps(value,indent=2)+'\n')
    save('ALL_SAMPLE_STATES.json',per_sample); save('LEAVE_ONE_COMPARATOR_OUT.json',controls)
    save('RECEIPT.json',{'status':'COMPLETE_DESCRIPTIVE_CONTROLS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_sha256':digest,'input_sites':N,'all_heterozygous_sites':32,'all_alternate_homozygous_sites':566,'independent_scalar_rescan_all_controls_passed':True,'elapsed_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limits':['Repository species labels and genotype encoding retained; source calls not revalidated.','Controls characterize the archived selected-SNP panel, not external validation.','Leaving one comparator out changes ascertainment; its comparator and Hachi use exactly the same reduced panel.','No ancestry proportions, molecular ploidy validation, hybrid generation or mechanism inferred.','Per-site observations are linked and not independent biological replicates.']})
    print(json.dumps({'status':'complete','controls':controls,'elapsed_seconds':time.monotonic()-started},indent=2))
if __name__=='__main__':run()
