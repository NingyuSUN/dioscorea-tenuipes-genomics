#!/usr/bin/env python3
"""Fixed full-panel paired comparison. Requires a root-reviewed, hash-bound plan.

No remote dispatch, no site replacement/filter, no actual PL in official kernel.
All result arrays and source calls are private and output-directory gitignored.
"""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
import contextlib,csv,gzip,hashlib,importlib.util,io,json,re,sys,time,traceback
from pathlib import Path
import sys
if not __debug__:
    raise RuntimeError('Optimized Python disables required scientific assertions; run without -O.')
_HELPER_DIR = str(Path(__file__).resolve().parent)
if _HELPER_DIR not in sys.path:
    sys.path.insert(0, _HELPER_DIR)
import os
# Bound native numerical libraries before importing them in the CLI process.
for _thread_env in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_thread_env] = '1'
import numpy as np

HERE=Path(__file__).resolve().parent
REFERENCE=HERE.parent/'explicit3.py'
REFERENCE_SHA='e01d701ae77893533eb73025bd19834755d200a684b5961c4ece3bf6dcd4db55'
COMPILED={'shared_cy':'b6839ac55c5d2c3813377daea7bd8ea939fe20be8589ad02dd54c5a6094a9c55','covariance_cy':'a0664e6bca8ed9174a841e9c96bebaf7f7e177f11b2ef5cdd0205b693fc5a41f'}
SETTINGS={'selected_sites':37566,'samples':14,'panel_sites':37566,'ranks':[2,1,3],'primary_rank':2,'pca_max_iter':1000,'pca_tolerance':1e-5,'reference_em_max_iter':10000,'reference_em_tolerance':1e-10,'official_em_max_iter':1000,'official_em_tolerance':1e-5,'epsilon':1e-6,'epsilon_sensitivity':1e-4,'seed':20260924}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj):
    with Path(p).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def bind_actual_alleles(source_marker):
    m=dict(source_marker)
    if not m.get('source_REF') or not m.get('source_ALT'):raise ValueError('Missing authoritative source REF/ALT')
    for key in ['REF','ALT']:
        if key in m and m[key]!=m['source_'+key]:raise ValueError('Conflicting source/canonical allele field')
        m[key]=m['source_'+key]
    return m

def canonical_key(m):
    return json.dumps([m['CHROM'],int(m['POS']),m['REF'],m['ALT']],ensure_ascii=False,separators=(',',':'))
def full_panel_indices(markers):
    """Every source row, in original order. No hashing/subsampling/filtering."""
    if [int(m['index']) for m in markers]!=list(range(len(markers))):
        raise ValueError('Full panel indices must be contiguous original order')
    return list(range(len(markers)))
def parse_pl(raw):
    if raw in ('','.'):return None,'whole_missing'
    fields=raw.split(',')
    if len(fields)!=3:raise ValueError('Illegal PL cardinality')
    if any(x in ('','.') for x in fields):return None,'partial_missing'
    if any(not re.fullmatch(r'-?\d+',x) for x in fields):raise ValueError('Illegal PL integer')
    values=tuple(map(int,fields))
    if min(values)<0 or max(values)>np.iinfo(np.int64).max:raise ValueError('PL outside nonnegative int64')
    if min(values)!=0:raise ValueError('Nonzero-minimum PL requires root review')
    return values,'valid'
def source_gt(raw):
    fields=raw.split(':')
    if len(fields)!=5:raise ValueError('Expected GT:AD:DP:GQ:PL')
    gt,_,dp,gq,_=fields;a=re.split(r'[/|]',gt)
    if any(x in ('','.') for x in a):return None,None
    if len(a)!=2 or any(x not in ('0','1') for x in a):raise ValueError('Unsupported complete GT')
    dosage=sum(map(int,a))
    if dp=='.':return dosage,None
    if not re.fullmatch(r'\d+',dp):raise ValueError('Illegal DP')
    if not 5<=int(dp)<=100:return dosage,None
    if gq=='.':return dosage,None
    if not re.fullmatch(r'\d+',gq):raise ValueError('Illegal GQ')
    return dosage,dosage if int(gq)>=20 else None

def verify_cell(raw,stored_gt,baseline,common):
    _,masked=source_gt(raw);values,status=parse_pl(raw.rsplit(':',1)[1])
    if masked!=(int(stored_gt) if baseline else None):raise ValueError('Source quality-masked GT differs from BED')
    if bool(common)!=bool(baseline and values is not None):raise ValueError('Common mask differs from source PL eligibility')
    return values,status

def decode_bed(data,markers,n):
    width=(n+3)//4
    if data[:3]!=b'\x6c\x1b\x01' or len(data)!=3+len(markers)*width:raise ValueError('BED header/dimensions')
    gt=np.full((len(markers),n),-1,dtype=np.int8)
    for i,m in enumerate(markers):
        block=int.from_bytes(data[3+i*width:3+(i+1)*width],'little')
        if block>>(2*n):raise ValueError('BED nonzero padding')
        for j in range(n):
            code=(block>>(2*j))&3
            if code==1:continue
            if code==2:gt[i,j]=1
            else:
                allele=m['A1'] if code==0 else m['A2']
                if allele not in (m['REF'],m['ALT']):raise ValueError('BIM/VCF allele mismatch')
                gt[i,j]=0 if allele==m['REF'] else 2
    return gt

def encode_bed(gt,markers):
    """Pack a fresh common-mask BED; independent decode uses whole-row integers."""
    n=gt.shape[1];width=(n+3)//4;out=bytearray(b'\x6c\x1b\x01')
    for calls,m in zip(gt,markers):
        packed=bytearray(width)
        for j,g in enumerate(calls):
            if g==-1:code=1
            elif g==1:code=2
            elif g in (0,2):
                allele=m['REF'] if g==0 else m['ALT'];code=0 if allele==m['A1'] else 3
            else:raise ValueError('Invalid ALT genotype dosage')
            packed[j//4]|=code<<(2*(j%4))
        out.extend(packed)
    return bytes(out)

def hc_likelihood(gt,common,epsilon):
    if not 0<epsilon<1/3:raise ValueError('epsilon outside (0,1/3)')
    if np.shape(common)!=np.shape(gt) or np.any(common & ~np.isin(gt,[0,1,2])):raise ValueError('HC genotype/common mask mismatch')
    x=np.full((*gt.shape,3),epsilon,dtype=np.float64)
    for g in (0,1,2):x[...,g][common & (gt==g)]=1-2*epsilon
    x[~common,:]=1/3
    return x

def pl_likelihood(pl,common):
    if pl.shape!=(*common.shape,3) or pl.dtype!=np.int64:raise ValueError('PL shape/type')
    active=pl[common]
    if (active<0).any() or (active.min(axis=1)!=0).any():raise ValueError('Invalid active PL')
    delta=active-active.min(axis=1,keepdims=True)
    with np.errstate(under='ignore'):weights=np.power(10.,-.1*delta.astype(np.float64))
    under=(weights==0)&(delta>0)
    weights/=weights.sum(axis=1,keepdims=True)
    x=np.full((*common.shape,3),1/3,dtype=np.float64);x[common]=weights
    return x,{'underflow_nonminimum_entries':int(under.sum()),'underflow_cells':int(under.any(axis=1).sum()),'called_cells':int(common.sum()),'minimum_subtraction_performed':True,'epsilon_or_PL_clipping_applied':False}

def load_panel(plan):
    files={k:Path(v['path']) for k,v in plan['inputs'].items()}
    for k,p in files.items():
        if sha(p)!=plan['inputs'][k]['sha256']:raise ValueError('Input hash mismatch: '+k)
    receipt=json.loads(files['extraction_receipt'].read_bytes())
    if receipt['status']!='EXTRACTED_PANEL_IDENTITY_VERIFIED' or receipt.get('script_exit_code')!=0:raise ValueError('Extraction not successful')
    binding=plan['source_binding']
    if sha(files['source_header'])!=receipt['source_header_sha256'] or receipt['source_header_sha256']!=binding['header_sha256']:raise ValueError('Source header identity')
    if receipt['VCF_path']!=binding['VCF_path'] or receipt['historical_source_VCF_sha256']!=binding['historical_body_sha256']:raise ValueError('Source VCF provenance mismatch')
    if not receipt['original_inputs_unchanged'] or receipt['VCF_stat_before']!=receipt['VCF_stat_after'] or receipt['indexes_before']!=receipt['indexes_after']:raise ValueError('Extraction source changed')
    if receipt['FAI']['sha256']!=binding['FAI_sha256']:raise ValueError('Reference dictionary provenance')
    header=files['source_header'].read_text();source_samples=next(x for x in header.splitlines() if x.startswith('#CHROM\t')).split('\t')[9:]
    if source_samples!=binding['source_samples'] or [x for x in source_samples if x in plan['samples']]!=plan['samples']:raise ValueError('Header sample identities/order')
    if not re.search(r'^##FORMAT=<ID=PL,Number=G,Type=Integer(?:,|>)',header,re.M):raise ValueError('Header PL definition')
    for label in ['markers','calls','baseline_mask','common_mask','source_header']:
        if sha(files[label])!=receipt['output_files'][files[label].name]['sha256']:raise ValueError('Extraction output binding differs: '+label)
    for label in ['bed','bim','fam']:
        if sha(files[label])!=receipt['panel_inputs']['original'][label]['sha256']:raise ValueError('Original input binding differs: '+label)
    with files['markers'].open() as f:markers=[bind_actual_alleles(m) for m in csv.DictReader(f,delimiter='\t')]
    for i,m in enumerate(markers):
        if int(m['index'])!=i:raise ValueError('Marker index/order')
        m['index']=i;m['POS']=int(m['POS'])
        if m['REF']==m['ALT'] or any(x not in 'ACGT' or len(x)!=1 for x in (m['REF'],m['ALT'])):raise ValueError('Not distinct biallelic SNP')
        if {m['A1'],m['A2']}!={m['REF'],m['ALT']}:raise ValueError('BIM/VCF allele set mismatch')
    if len(markers)!=SETTINGS['panel_sites'] or len({(m['CHROM'],m['POS']) for m in markers})!=len(markers):raise ValueError('Panel dimension/duplicate coordinate')
    fam=files['fam'].read_text().splitlines();bim=files['bim'].read_text().splitlines()
    if any(len(x.split())!=6 for x in fam+bim):raise ValueError('Malformed BIM/FAM')
    samples=[x.split()[1] for x in fam]
    if samples!=plan['samples'] or len(samples)!=14 or len(set(samples))!=14 or len(bim)!=len(markers):raise ValueError('Samples/BIM dimensions')
    for m,line in zip(markers,bim):
        f=line.split()
        if (f[0],int(f[3]),f[4],f[5],f[1])!=(m['CHROM'],m['POS'],m['A1'],m['A2'],m['original_ID']):raise ValueError('BIM row identity')
    baseline=np.frombuffer(files['baseline_mask'].read_bytes(),dtype=np.uint8)
    common=np.frombuffer(files['common_mask'].read_bytes(),dtype=np.uint8)
    if baseline.size!=len(markers)*14 or common.size!=baseline.size or (baseline>1).any() or (common>1).any():raise ValueError('Mask dimensions/values')
    baseline=baseline.reshape(-1,14).astype(bool);common=common.reshape(-1,14).astype(bool)
    if np.any(common & ~baseline):raise ValueError('Mask restored missing baseline GT')
    gt=decode_bed(files['bed'].read_bytes(),markers,14)
    if not np.array_equal(gt>=0,baseline):raise ValueError('BED/mask missingness mismatch')
    selected=full_panel_indices(markers);where={idx:k for k,idx in enumerate(selected)}
    pl=np.full((len(markers),14,3),-1,dtype=np.int64);raw_selected=[];counts={'whole_missing':0,'partial_missing':0,'valid':0}
    with gzip.open(files['calls'],'rt') as f:
        header=f.readline()
        expected=['#CHROM','POS','REF','ALT','FILTER']+[s+'[GT:AD:DP:GQ:PL]' for s in samples]
        if header.rstrip('\n').split('\t')!=expected:raise ValueError('Raw header/sample order')
        i=-1
        for i,line in enumerate(f):
            if i>=len(markers):raise ValueError('Extra raw row')
            fields=line.rstrip('\n').split('\t');m=markers[i]
            if len(fields)!=19 or fields[:5]!=[m['CHROM'],str(m['POS']),m['REF'],m['ALT'],'PASS']:raise ValueError('Raw VCF authoritative key differs from bound marker table')
            for j,cell in enumerate(fields[5:]):
                values,status=verify_cell(cell,gt[i,j],baseline[i,j],common[i,j]);counts[status]+=1
                if i in where and values is not None:pl[where[i],j]=values
            if i in where:raw_selected.append(line)
        if i+1!=len(markers):raise ValueError('Missing raw rows')
    picked=[markers[i] for i in selected];called=common[selected];pgt=gt[selected].copy();pgt[~called]=-1
    if np.any(called.sum(axis=1)==0):raise ValueError('Full panel all-missing site requires review')
    ac=np.maximum(pgt,0).sum(axis=1);an=2*called.sum(axis=1)
    if np.any((ac==0)|(ac==an)):raise ValueError('Full panel HC monomorphic site requires review; no replacement')
    gl,tail=pl_likelihood(pl,called)
    if np.any(np.all(gl==gl[...,:1],axis=(1,2))):raise ValueError('Full panel all-flat PL site requires review')
    return dict(markers=picked,selected=selected,samples=samples,bim_lines=[bim[i] for i in selected],fam_bytes=files['fam'].read_bytes(),header=header,raw_rows=raw_selected,gt=pgt,common=called,pl=pl,gl=gl,tail=tail,pl_counts_entire_panel=counts,input_files=files)

def describe_cov(c):
    c=np.asarray(c,dtype=np.float64)
    if c.shape!=(14,14) or not np.isfinite(c).all():raise ValueError('Covariance dimensions/nonfinite')
    symmetry=float(np.max(np.abs(c-c.T)))
    if symmetry>1e-10:raise ValueError('Covariance asymmetric')
    values,vectors=np.linalg.eigh(c);values=values[::-1];vectors=vectors[:,::-1]
    if values[-1]<-1e-8:raise ValueError('Covariance negative eigenvalue requires review')
    return {'all_14_eigenvalues':values.tolist(),'diagonal':np.diag(c).tolist(),'all_finite':True,'max_asymmetry':symmetry,'minimum_eigenvalue':float(values[-1])},vectors

PYTHON_SOURCES={'shared':'ea89ab7a0924d4611d2ecb93091361e1d770b256ce055701e45605252aa23b51','covariance':'abaf0c0e9538c73230299efece68faef128465815c5cd07d20f3f671d20dbf3e'}
DEFAULT_MODULE_SHA256={**COMPILED,**PYTHON_SOURCES}
REPRESENTATIONS=['B_HC_e1e6','C_PL','B_HC_e1e4','B_official_HC_e1e6']

def array_sha(x):return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
def check_inputs(plan):
    for label,item in plan['inputs'].items():
        if sha(item['path'])!=item['sha256']:raise ValueError('Input identity changed: '+label)

def arm_schedule():
    return [(name,rank) for name in REPRESENTATIONS for rank in SETTINGS['ranks']]

def progress(out,state,event):
    """Immutable event receipts plus replaceable current state, only in this fresh output."""
    state['progress_event']=event
    state['elapsed_seconds']=time.monotonic()-state['_monotonic_start']
    state['utc']=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()
    public={k:v for k,v in state.items() if not k.startswith('_')}
    dump(out/'progress'/('%03d.json'%state['_event_index']),public)
    temporary=out/'LIVE_STATE.json.tmp'
    with temporary.open('x') as f:json.dump(public,f,indent=2,allow_nan=False);f.write('\n')
    temporary.replace(out/'LIVE_STATE.json')
    state['_event_index']+=1
    print(json.dumps({'event':event,'utc':state['utc'],'status':state['status'],'elapsed_seconds':state['elapsed_seconds']},ensure_ascii=False),flush=True)

def save_arrays(path,result):
    np.savez_compressed(path,**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})

def covariance_metadata(c):
    description,vectors=describe_cov(c)
    return description,vectors

def run_reference_ranks(ref,gl,frequency,name,out,state,plan):
    f=frequency['f'];fhash=array_sha(f)
    for rank in SETTINGS['ranks']:
        arm=f'{name}_r{rank}';check_inputs(plan)
        state['arms'][arm]={'status':'RUNNING','representation':name,'rank':rank,'frequency_array_sha256':fhash}
        progress(out,state,arm+'_START');t=time.monotonic()
        try:
            result=ref.fit(gl,rank=rank,f=f,pca_max_iter=SETTINGS['pca_max_iter'],pca_tolerance=SETTINGS['pca_tolerance'])
            save_arrays(out/(arm+'.npz'),result)
            if array_sha(f)!=fhash or array_sha(result['f'])!=fhash:raise ValueError('Fixed population frequency altered')
            description,vectors=covariance_metadata(result['covariance']);np.save(out/(arm+'.eigenvectors.npy'),vectors)
            state['arms'][arm].update(description)
            state['arms'][arm].update(status='CONVERGED' if result['pca_converged'] else 'NONCONVERGED_RETAINED',pca_converged=bool(result['pca_converged']),pca_iterations=int(result['pca_iterations']),pca_v_rmse_trace=result['v_rmse_trace'],f_min=float(f.min()),f_max=float(f.max()),all_singular_values=result['all_singular_values'].tolist(),rank_cut_gap=result['rank_cut_singular_value_gap'],rank_cut_relative_gap=result['rank_cut_gap_relative_to_last_retained'],all_sites_flat=result['all_sites_flat'],all_flat_site_count=int(result['all_cells_flat_site_mask'].sum()),singular_scope='Full singular spectrum of final iteration centered dosage, returned by frozen reference; pre-final full dosage matrix not separately saved')
        except Exception as error:
            state['arms'][arm].update(status='FAILED_RETAINED',error=type(error).__name__+': '+str(error))
            (out/(arm+'.failure.txt')).write_text(traceback.format_exc())
        finally:
            state['arms'][arm]['seconds']=time.monotonic()-t
            progress(out,state,arm+'_END')
        check_inputs(plan)
        if array_sha(f)!=fhash:raise ValueError('Population frequency mutated; subsequent ranks unsafe')

def run_official_ranks(ref,covariance,legacy,represented_hc,f,out,state,plan):
    name='B_official_HC_e1e6';fhash=array_sha(f)
    for rank in SETTINGS['ranks']:
        arm=f'{name}_r{rank}';check_inputs(plan)
        state['arms'][arm]={'status':'RUNNING','representation':name,'rank':rank,'frequency_array_sha256':fhash,'seed':SETTINGS['seed']}
        progress(out,state,arm+'_START');t=time.monotonic()
        try:
            np.random.seed(SETTINGS['seed'])
            with (out/(arm+'.log')).open('x') as log,contextlib.redirect_stdout(log):
                c,pi,k,it,converged=covariance.emPCA(legacy,f,rank,SETTINGS['pca_max_iter'],SETTINGS['pca_tolerance'])
            np.savez_compressed(out/(arm+'.npz'),covariance=c,pi=pi,f=f)
            if int(k)!=rank or array_sha(f)!=fhash:raise ValueError('Official rank/frequency changed')
            description,vectors=covariance_metadata(c);np.save(out/(arm+'.eigenvectors.npy'),vectors)
            e=ref.dosage(ref.posterior(represented_hc,pi))-2*f[:,None]
            s=np.linalg.svd(e,compute_uv=False)
            state['arms'][arm].update(description)
            state['arms'][arm].update(status='CONVERGED' if converged else 'NONCONVERGED_RETAINED',pca_converged=bool(converged),pca_updates=int(it)+1,f_min=float(f.min()),f_max=float(f.max()),all_singular_values_diagnostic=s.tolist(),rank_cut_gap_diagnostic=float(s[rank-1]-s[rank]),rank_cut_relative_gap_diagnostic=float((s[rank-1]-s[rank])/s[rank-1]) if s[rank-1]>0 else None,singular_scope='Recomputed from saved float32 final pi; not intercepted final iterative SVD or pre-cast pi',input='Only HC epsilon1e-6 in two-float32/double residual representation; actual PL excluded')
        except Exception as error:
            state['arms'][arm].update(status='FAILED_RETAINED',error=type(error).__name__+': '+str(error))
            (out/(arm+'.failure.txt')).write_text(traceback.format_exc())
        finally:
            state['arms'][arm]['seconds']=time.monotonic()-t
            progress(out,state,arm+'_END')
        check_inputs(plan)
        if array_sha(f)!=fhash:raise ValueError('Official population frequency mutated; subsequent ranks unsafe')

def final_status(arms):
    complete=len(arms)==len(arm_schedule()) and all(x.get('status')=='CONVERGED' for x in arms.values())
    return 'COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT' if complete else 'PARTIAL_PENDING_ROOT_SCIENTIFIC_AUDIT'

def run(plan_path,out):
    plan_path=Path(plan_path);out=Path(out);plan=load_execution_plan(plan_path)
    if plan['status']!='ROOT_APPROVED_FIXED_FULL_PANEL' or plan['settings']!=SETTINGS:raise ValueError('Plan not root-reviewed/fixed settings differ')
    if sha(__file__)!=plan['expected_runner_sha256']:raise ValueError('Reviewed runner changed')
    if sha(REFERENCE)!=REFERENCE_SHA or plan['expected_reference_sha256']!=REFERENCE_SHA:raise ValueError('Reference kernel changed')
    check_inputs(plan)
    from pcangsd import shared,covariance,shared_cy,covariance_cy
    modules=dict(shared=shared,covariance=covariance,shared_cy=shared_cy,covariance_cy=covariance_cy)
    for name,expected in {**COMPILED,**PYTHON_SOURCES}.items():
        if sha(modules[name].__file__)!=expected:raise ValueError('Official module changed: '+name)
    out.mkdir(parents=True,exist_ok=False);(out/'.gitignore').write_text('*\n!.gitignore\n');(out/'progress').mkdir()
    state={'status':'RUNNING','plan_sha256':sha(plan_path),'reference_sha256':REFERENCE_SHA,'runner_sha256':sha(__file__),'settings':SETTINGS,'frequencies':{},'arms':{},'PLINK_not_run':True,'private_data_exclude_Git':True,'python':sys.version,'python_executable':sys.executable,'numpy_version':np.__version__,'official_modules':{name:{'path':m.__file__,'sha256':sha(m.__file__)} for name,m in modules.items()},'source_binding':plan['source_binding'],'historical_VCF_body_hash_not_fresh_rehash':True,'_monotonic_start':time.monotonic(),'_event_index':0}
    progress(out,state,'INPUT_CHECK_START')
    try:
        panel=load_panel(plan)
        if panel['selected']!=list(range(SETTINGS['panel_sites'])):raise ValueError('Full-panel order/coverage mismatch')
        bed=encode_bed(panel['gt'],panel['markers'])
        if not np.array_equal(decode_bed(bed,panel['markers'],14),panel['gt']):raise ValueError('Independent full BED decode mismatch')
        (out/'paired.bed').write_bytes(bed);(out/'paired.bim').write_text('\n'.join(panel['bim_lines'])+'\n');(out/'paired.fam').write_bytes(panel['fam_bytes'])
        with (out/'panel_markers_with_source_fields.tsv').open('x') as f:
            w=csv.DictWriter(f,fieldnames=list(panel['markers'][0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(panel['markers'])
        np.savez_compressed(out/'private_inputs.npz',integer_pl=panel['pl'],common_mask=panel['common'],reported_GT_ALT_dosage=panel['gt'],PL_explicit3=panel['gl'])
        state.update(samples=panel['samples'],selected_sites=SETTINGS['panel_sites'],called_cells=int(panel['common'].sum()),missing_cells=int((~panel['common']).sum()),mask_sha256=array_sha(panel['common'].astype(np.uint8)),PL_tails=panel['tail'],PL_counts_entire_panel=panel['pl_counts_entire_panel'],BED_all_525924_cells_independently_decoded=True,selection='All original panel rows in original BIM order; no sampling/filtering/replacement',original_input_hashes={k:sha(p) for k,p in panel['input_files'].items()})
        progress(out,state,'INPUT_CHECK_COMPLETE')
        spec=importlib.util.spec_from_file_location('explicit3_fullpanel_reference',REFERENCE);ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
        for name in REPRESENTATIONS[:3]:
            check_inputs(plan);t=time.monotonic();progress(out,state,name+'_FREQUENCY_START')
            gl=panel['gl'] if name=='C_PL' else hc_likelihood(panel['gt'],panel['common'],SETTINGS['epsilon'] if name=='B_HC_e1e6' else SETTINGS['epsilon_sensitivity'])
            frequency=ref.em_maf(gl,max_iter=SETTINGS['reference_em_max_iter'],tolerance=SETTINGS['reference_em_tolerance'])
            np.savez_compressed(out/(name+'.frequency.npz'),f=frequency['f'],max_abs_update_trace=np.array(frequency['max_abs_update_trace']))
            f=frequency['f'];valid=bool(np.isfinite(f).all() and np.all((f>0)&(f<1)))
            state['frequencies'][name]={'converged':bool(frequency['converged']),'iterations':int(frequency['iterations']),'max_abs_update_trace':frequency['max_abs_update_trace'],'f_min':float(f.min()),'f_max':float(f.max()),'strict_interior':valid,'frequency_array_sha256':array_sha(f),'seconds':time.monotonic()-t,'reused_without_reestimation_at_ranks':SETTINGS['ranks']}
            progress(out,state,name+'_FREQUENCY_END')
            if not frequency['converged'] or not valid:
                for rank in SETTINGS['ranks']:state['arms'][f'{name}_r{rank}']={'status':'SKIPPED_FREQUENCY_INVALID','rank':rank,'representation':name}
                progress(out,state,name+'_RANKS_SKIPPED');continue
            run_reference_ranks(ref,gl,frequency,name,out,state,plan)
        check_inputs(plan)
        hc=hc_likelihood(panel['gt'],panel['common'],SETTINGS['epsilon'])
        legacy=np.ascontiguousarray(hc[...,:2].reshape(SETTINGS['panel_sites'],28),dtype=np.float32)
        residual=1-legacy[:,0::2].astype(np.float64)-legacy[:,1::2].astype(np.float64)
        if (residual<0).any():raise ValueError('Official HC residual negative')
        represented_hc=np.stack((legacy[:,0::2].astype(np.float64),legacy[:,1::2].astype(np.float64),residual),axis=2)
        name='B_official_HC_e1e6';progress(out,state,name+'_FREQUENCY_START');t=time.monotonic()
        with (out/(name+'.frequency.log')).open('x') as log,contextlib.redirect_stdout(log):
            f=shared.emMAF(legacy,SETTINGS['official_em_max_iter'],SETTINGS['official_em_tolerance'])
        np.savez_compressed(out/(name+'.frequency.npz'),f=f)
        text=(out/(name+'.frequency.log')).read_text();maf_ok='EM (MAF) converged' in text and 'EM (MAF) failed to converge' not in text
        valid=bool(np.isfinite(f).all() and np.all((f>0)&(f<1)))
        state['frequencies'][name]={'converged_log':maf_ok,'strict_interior':valid,'f_min':float(f.min()),'f_max':float(f.max()),'frequency_array_sha256':array_sha(f),'seconds':time.monotonic()-t,'residual_GL2_min':float(residual.min()),'reused_without_reestimation_at_ranks':SETTINGS['ranks'],'frequency_log':name+'.frequency.log'}
        progress(out,state,name+'_FREQUENCY_END')
        if maf_ok and valid:run_official_ranks(ref,covariance,legacy,represented_hc,f,out,state,plan)
        else:
            for rank in SETTINGS['ranks']:state['arms'][f'{name}_r{rank}']={'status':'SKIPPED_FREQUENCY_INVALID','rank':rank,'representation':name}
            progress(out,state,name+'_RANKS_SKIPPED')
        check_inputs(plan);state['status']=final_status(state['arms'])
    except Exception as error:
        state['status']='FAILED_PENDING_ROOT_REVIEW';state['error']=type(error).__name__+': '+str(error)
        (out/'failure.txt').write_text(traceback.format_exc())
        raise
    finally:
        progress(out,state,'FINAL')
        state['output_files']=[{'path':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(out.rglob('*')) if p.is_file()]
        dump(out/'RUN_RECEIPT.json',{k:v for k,v in state.items() if not k.startswith('_')})
    return 0 if state['status']=='COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT' else 2


def load_execution_plan(plan_path):
    """Resolve declared artifact paths without rewriting their provenance fields."""
    plan_path = Path(plan_path).resolve()
    plan = json.loads(plan_path.read_text())
    for item in plan['inputs'].values():
        value = Path(item['path']).expanduser()
        item['path'] = str((value if value.is_absolute() else plan_path.parent / value).resolve())
    return plan


def configure(config, base):
    from _population_cli import path_value, require_file, fresh_output
    global REFERENCE, COMPILED, PYTHON_SOURCES
    REFERENCE = require_file(path_value(config, 'reference', base, HERE.parent / 'explicit3.py'))
    declared = config.get('official_module_sha256', DEFAULT_MODULE_SHA256)
    if set(declared) != {'shared', 'covariance', 'shared_cy', 'covariance_cy'}:
        raise ValueError('All four PCAngsd module SHA256 values must be supplied together')
    COMPILED = {k: declared[k] for k in ('shared_cy', 'covariance_cy')}
    PYTHON_SOURCES = {k: declared[k] for k in ('shared', 'covariance')}
    plan_path = require_file(path_value(config, 'plan', base))
    plan = load_execution_plan(plan_path)
    out = fresh_output(path_value(config, 'outdir', base), [plan_path, REFERENCE, *[Path(x['path']) for x in plan['inputs'].values()]])
    return plan_path, plan, out


def preflight(plan_path, plan, out):
    from _population_cli import expected_hash, report_plan
    if plan['status'] != 'ROOT_APPROVED_FIXED_FULL_PANEL' or plan['settings'] != SETTINGS:
        raise ValueError('Plan is not reviewed or fixed scientific settings differ')
    expected_hash(Path(__file__), plan['expected_runner_sha256'], 'reviewed runner')
    expected_hash(REFERENCE, REFERENCE_SHA, 'frozen explicit-three-genotype reference')
    if plan['expected_reference_sha256'] != REFERENCE_SHA:
        raise ValueError('Plan reference identity differs')
    required = {'markers', 'calls', 'baseline_mask', 'common_mask', 'source_header', 'extraction_receipt', 'bed', 'bim', 'fam'}
    if not required <= set(plan['inputs']):
        raise ValueError('Plan lacks required input bindings: ' + ', '.join(sorted(required - set(plan['inputs']))))
    if len(plan.get('samples', [])) != 14 or len(set(plan['samples'])) != 14:
        raise ValueError('Plan requires 14 unique samples in their frozen order')
    binding_keys = {'header_sha256', 'VCF_path', 'historical_body_sha256', 'FAI_sha256', 'source_samples'}
    if not binding_keys <= set(plan.get('source_binding', {})):
        raise ValueError('Plan lacks source-identity fields')
    check_inputs(plan)
    from pcangsd import shared, covariance, shared_cy, covariance_cy
    modules = dict(shared=shared, covariance=covariance, shared_cy=shared_cy, covariance_cy=covariance_cy)
    for name, expected in {**COMPILED, **PYTHON_SOURCES}.items():
        expected_hash(Path(modules[name].__file__), expected, 'PCAngsd ' + name)
    declared_inputs = {'plan': plan_path, 'reference': REFERENCE, **{key: Path(item['path']) for key, item in plan['inputs'].items()}, **{'pcangsd_' + name: Path(module.__file__) for name, module in modules.items()}}
    declared_outputs = {'receipt': out / 'RUN_RECEIPT.json', **{name: out / name for name in ['paired.bed', 'paired.bim', 'paired.fam']}, **{name + '_r' + str(rank): out / (name + '_r' + str(rank) + '.npz') for name, rank in arm_schedule()}}
    return report_plan('run_full_panel.py', declared_inputs,
                       declared_outputs,
                       ['numpy', 'pcangsd'], fixed_settings=SETTINGS,
                       note='Preflight verifies bindings only; full sample/site/likelihood checks run before fitting.')


def plan(config_path):
    """Return a read-only, JSON-compatible artifact and dependency plan."""
    config_path=Path(config_path).resolve();config=json.loads(config_path.read_text());p,data,out=configure(config,config_path.parent);return preflight(p,data,out)


def main(argv=None):
    from _population_cli import read_config
    config, base, execute = read_config(__doc__, argv)
    plan_path, plan, out = configure(config, base)
    print(json.dumps(preflight(plan_path, plan, out), indent=2))
    return run(plan_path, out) if execute else 0


if __name__ == '__main__':
    raise SystemExit(main())
