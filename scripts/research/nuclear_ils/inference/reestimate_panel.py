#!/usr/bin/env python3
"""Original gene inference on one simulated panel; no hypothesis-test P value.

Default is a local plan validation. Actual execution requires a newly named
Linux output, verified frozen program/input hashes, live resources and bounded
systemd scopes. Failed genes are never silently retried or settings shortened.
"""
from pathlib import Path
from collections import Counter
import argparse,concurrent.futures,datetime,hashlib,importlib.util,json,math,os,re,signal,subprocess,sys,threading,time,traceback

P=Path(__file__).resolve().parent
FIRST=('209041at3193','140075at3193','49273at3193')
IQ_SHA='a5a11554ca60b5682fd50af81301ea247c887d93a7d7014fcd61c0fced033832'
AA=set('ARNDCQEGHILKMFPSTWYV')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,indent=2)+'\n');q.replace(p)
def require(ok,message):
    if not ok:raise ValueError(message)
def fasta(p):
    out={};n=None
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):
            n=line[1:].split()[0];require(n not in out,'Duplicate FASTA label');out[n]=''
        elif line.strip():require(n is not None,'FASTA before header');out[n]+=line.strip()
    return out
def load_vendor(name):
    spec=importlib.util.spec_from_file_location(name,P/'vendor'/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def prepare(jobs_path,out):
    data=json.loads(jobs_path.read_text());jobs=data['jobs'];order=data['gene_order']
    require(len(jobs)==len(order)==len(set(order))==173,'Require one original 173-locus panel')
    require(all(re.fullmatch(r'[0-9]+at3193',g) for g in order),'Unexpected locus identifier/path component')
    require([j['gene'] for j in jobs]==order and '2661at3193' not in order,'Gene order/exclusion mismatch')
    require(set(FIRST)<=set(order),'Missing fixed engineering connection loci')
    require(len(data['taxa12'])==len(set(data['taxa12']))==12,'Full species sample changed')
    sources={str(jobs_path):sha(jobs_path)};planned=[]
    for j in jobs:
        path=Path(j['input_fasta']);require(sha(path)==j['input_sha256'],'Masked source changed')
        rows=fasta(path);expected=j['expected_inference_taxa']
        require(list(rows)==expected and set(expected)<=set(data['taxa12']),'Retained leaf order/set changed')
        require({len(v) for v in rows.values()}=={j['original_full_length']} and all(set(v)&AA for v in rows.values()),'Alignment length/all-missing row policy changed')
        cmd=list(j['argv']);prefix=out/'genes'/j['gene']/'inferred'
        require(cmd.count('-s')==cmd.count('--prefix')==cmd.count('-seed')==1,'Ambiguous original command binding')
        require(Path(cmd[cmd.index('-s')+1])==path,'Input command binding mismatch')
        require(int(cmd[cmd.index('-seed')+1])==j['seed'],'Seed binding mismatch')
        original=list(j['original_argv']);expected_cmd=list(original);expected_cmd[0]=cmd[0]
        for key in ['-s','-seed','--prefix']:
            require(original.count(key)==1,'Original command binding is ambiguous')
            expected_cmd[original.index(key)+1]=cmd[cmd.index(key)+1]
        require(expected_cmd==cmd,'Command differs beyond declared binary/input/seed/prefix replacements')
        settings={'-m':'MFP','-mset':'JTT,WAG,LG','-mfreq':'F','-mrate':'G4','-B':'1000','-alrt':'1000','-T':'1'}
        require(all(cmd.count(k)==1 and cmd[cmd.index(k)+1]==v for k,v in settings.items()),'Original inference settings changed')
        require(all(k in cmd for k in ['-safe','-bnni','--boot-trees']),'Original support settings missing')
        require(not any(k in cmd for k in ['-redo','-n','-nstop','--fast','-fast']), 'Shortened/restart search not admitted')
        cmd[cmd.index('--prefix')+1]=str(prefix)
        true_tree=Path(j['true_tree_path']);require(sha(true_tree)==j['true_tree_sha256'],'True genealogy changed')
        bindings=j['source_bindings']
        if isinstance(bindings,dict):bindings=list(bindings.values())
        require(bool(bindings),'No simulation provenance binding')
        for binding in bindings:
            require(sha(Path(binding['path']))==binding['sha256'],'Generation source seal changed');sources[binding['path']]=binding['sha256']
        q=dict(j);q.update(argv=cmd,output_prefix=str(prefix),full12_taxa=data['taxa12']);planned.append(q);sources[str(path)]=sha(path);sources[str(true_tree)]=sha(true_tree)
    return data,planned,sources

def summarize_gene(job,projection,tree_parser):
    prefix=Path(job['output_prefix']);tree=prefix.with_suffix('.treefile');boot=prefix.with_suffix('.ufboot')
    tips,splits=tree_parser.full_splits(tree.read_text());require(tips==set(job['expected_inference_taxa']),'Estimated tips changed')
    lines=[s for s in boot.read_text().splitlines() if s.strip()];require(len(lines)==1000,'Expected 1000 saved bootstrap trees')
    true_tips,true_splits=tree_parser.full_splits(Path(job['true_tree_path']).read_text())
    require(true_tips==set(job['full12_taxa']) and tips<=true_tips,'True genealogy does not include original full12 sample')
    induced=set()
    for split in true_splits:
        a=set(split)&tips;b=tips-a
        if min(len(a),len(b))>=2:induced.add(min(tuple(sorted(a)),tuple(sorted(b))))
    true_point,_=projection.project(Path(job['true_tree_path']).read_text(),true_tips)
    eligible=set(projection.TAXA)<=tips;counts=Counter({k:0 for k in projection.KEYS});point=None
    if eligible:
        point,_=projection.project(tree.read_text(),tips)
        for line in lines:counts[projection.project(line,tips)[0]]+=1
    else:
        # Validate every bootstrap's full retained label set even without focal quartet.
        for line in lines:require(tree_parser.full_splits(line)[0]==tips,'Bootstrap tips changed')
    supported=next((k for k in projection.LABEL.values() if counts[k]>=950),None) if eligible else None
    result=dict(gene=job['gene'],taxa=len(tips),AA_columns=job['original_full_length'],focal_eligible=eligible,true_focal_split=true_point,point_split=point,full_unrooted_RF_on_retained_taxa=len(induced^splits),bootstrap_counts=dict(counts) if eligible else None,bootstrap_denominator=1000 if eligible else None,saved_bootstrap_trees=1000,bootstrap_supported_split_at_950=supported,true_tree_sha256=sha(Path(job['true_tree_path'])),output_tree_sha256=sha(tree),output_bootstrap_sha256=sha(boot),scope='One conditional simulated locus, not independent evidence for ILS')
    save(prefix.parent/'GENE_RESULT.json',result);return result

def resource_preflight(planned):
    binaries={j['argv'][0] for j in planned};require(len(binaries)==1,'Multiple IQ-TREE binaries')
    binary=Path(next(iter(binaries)));require(sha(binary)==IQ_SHA,'Pinned IQ-TREE binary differs')
    mem={k:int(v.split()[0])*1024 for k,v in (s.split(':',1) for s in Path('/proc/meminfo').read_text().splitlines())}
    require(mem['MemAvailable']>=18*1024**3,'Need >=18GiB currently available for eight 2GiB caps plus margin')
    require(set(range(8))<=os.sched_getaffinity(0),'CPUs0-7 unavailable')
    return dict(iqtree=str(binary),iqtree_sha256=sha(binary),MemAvailable_bytes=mem['MemAvailable'],CPU_ids=list(range(8)),workers=8,memory_each_bytes=2*1024**3,swap_each_bytes=0)

def execute(jobs_path,out,freeze_path,deadline_utc):
    # This is a technical guard; it never substitutes for tool/user authorization.
    require(str(out).startswith('workspace/wgs/dt_research_p57_20260930/'),'Must use declared Linux project output')
    require(not out.exists(),'Output already exists; inspect state, do not restart blindly')
    frozen=json.loads(freeze_path.read_text())
    require({'reestimate_panel.py','vendor/original_quartet_projector.py','vendor/engineering_parser.py','vendor/exec_verified_scope.py'}<=set(frozen['files']),'Incomplete inference program freeze')
    for rel,h in frozen['files'].items():require(sha(P/rel)==h,'Frozen inference implementation changed')
    data,planned,sources=prepare(jobs_path,out);resources=resource_preflight(planned)
    dt=datetime.datetime.fromisoformat(deadline_utc);require(dt.tzinfo is not None,'Deadline must include timezone offset');deadline=dt.timestamp()
    require(deadline>time.time(),'Explicit run deadline has passed')
    projection=load_vendor('original_quartet_projector');tree_parser=load_vendor('engineering_parser')
    require(projection.fixtures()==6,'Original quartet fixture failure')
    out.mkdir(parents=True);save(out/'INPUT_BINDING.json',dict(source_hashes=sources,jobs=planned,original_jobs_sha256=sha(jobs_path),driver_freeze_sha256=sha(freeze_path)))
    state=dict(status='RUNNING_CONDITIONAL_PANEL_INFERENCE',started_epoch=time.time(),pid=os.getpid(),completed_genes=[],jobs=173,resources=resources,deadline_utc=deadline_utc,not_calibrated_ILS=True)
    save(out/'STATE.json',state);stop=threading.Event();lock=threading.Lock();results=[];active_units=set()
    def interrupt(sig,frame):
        stop.set()
        with lock:units=list(active_units)
        for unit in units:subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=TERM',unit],capture_output=True)
        raise InterruptedError('Interrupted; results preserved and no next gene started')
    signal.signal(signal.SIGTERM,interrupt);signal.signal(signal.SIGINT,interrupt)
    def run(job,cpu):
        require(not stop.is_set(),'Sibling failed; no further launch')
        remaining=deadline-time.time();require(remaining>0,'Deadline reached; unstarted genes remain pending')
        prefix=Path(job['output_prefix']);prefix.parent.mkdir(parents=True,exist_ok=False)
        unit='dt-msc03-'+str(os.getpid())+'-'+job['gene'];guard=prefix.parent/'SCOPE_RECEIPT.json'
        argv=['systemd-run','--user','--scope','--quiet','--unit',unit,'-p','MemoryMax=2147483648','-p','MemorySwapMax=0','taskset','-c',str(cpu),sys.executable,str(P/'vendor/exec_verified_scope.py'),'--memory-bytes','2147483648','--cpus','1','--receipt',str(guard),'--','/usr/bin/time','-v','-o',str(prefix.parent/'TIME.txt'),*job['argv']]
        receipt=dict(status='RUNNING',argv=argv,gene=job['gene'],started_epoch=time.time(),timeout_seconds=min(480,remaining),unit=unit+'.scope')
        save(prefix.parent/'RUN_RECEIPT.json',receipt)
        try:
            with (prefix.parent/'stdout.txt').open('x') as stdout,(prefix.parent/'stderr.txt').open('x') as stderr:
                env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
                with lock:active_units.add(unit+'.scope')
                child=subprocess.Popen(argv,stdout=stdout,stderr=stderr,env=env,start_new_session=True);receipt['pid']=child.pid
                save(prefix.parent/'RUN_RECEIPT.json',receipt)
                try:code=child.wait(timeout=receipt['timeout_seconds'])
                except subprocess.TimeoutExpired:
                    subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=TERM',unit+'.scope'],capture_output=True)
                    try:child.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=KILL',unit+'.scope'],capture_output=True);child.wait(timeout=10)
                    raise TimeoutError('Inference incomplete; retained without settings reduction')
            require(code==0,'IQ-TREE failed '+str(code));g=json.loads(guard.read_text())
            require(g['memory_max_bytes']==2147483648 and g['swap_max_bytes']==0 and g['logical_cpu_ids']==[cpu],'Actual resource scope mismatch')
            result=summarize_gene(job,projection,tree_parser);receipt.update(status='COMPLETE',result_sha256=sha(prefix.parent/'GENE_RESULT.json'))
            with lock:
                results.append(result);state['completed_genes'].append(job['gene']);save(out/'STATE.json',state)
            return result
        except BaseException:
            stop.set();receipt.update(status='FAILED_OR_INCOMPLETE_NO_AUTORETRY',error=traceback.format_exc());raise
        finally:
            with lock:active_units.discard(unit+'.scope')
            receipt.update(finished_epoch=time.time(),elapsed_seconds=time.time()-receipt['started_epoch']);save(prefix.parent/'RUN_RECEIPT.json',receipt)
    by_gene={j['gene']:j for j in planned}
    try:
        for gene in FIRST:run(by_gene[gene],0)
        remaining=[j for j in planned if j['gene'] not in FIRST]
        def lane(cpu):
            for job in remaining[cpu::8]:run(job,cpu)
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            futures=[pool.submit(lane,i) for i in range(8)]
            for future in concurrent.futures.as_completed(futures):future.result()
        require(len(results)==173,'Incomplete panel')
        ordered=sorted(results,key=lambda v:data['gene_order'].index(v['gene']))
        eligible=[r for r in ordered if r['focal_eligible']];require(len(eligible)==139,'Original focal occupancy changed')
        cells=Counter((r['true_focal_split'],r['point_split']) for r in eligible)
        summary=dict(status='ONE_CONDITIONAL_PANEL_INFERRED_SCIENTIFIC_AUDIT_PENDING',results=ordered,completed_genes=173,focal_eligible_genes=len(eligible),true_counts_eligible=dict(Counter(r['true_focal_split'] for r in eligible)),point_counts=dict(Counter(r['point_split'] for r in eligible)),truth_estimate_cross_table=[dict(true=k[0],estimated=k[1],genes=v) for k,v in sorted(cells.items())],bootstrap_950_counts=dict(Counter(r['bootstrap_supported_split_at_950'] or 'none' for r in eligible)),calibrated_P_value=None,limitations=['One source-fitted conditional panel, not a calibrated test','Full-tree RF uses only each original retained sample; focal cross-table uses139 original eligible genes','No empirical independence, physical Ne, or gene-flow inference follows'])
        save(out/'PANEL_RESULTS.json',summary);state['status']=summary['status']
    except BaseException:
        stop.set();state.update(status='FAILED_OR_PARTIAL_NO_AUTORETRY',error=traceback.format_exc());raise
    finally:
        state.update(finished_epoch=time.time(),completed_count=len(results));save(out/'STATE.json',state)
        save(out/'OUTPUT_MANIFEST.json',{'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'}})
    return state

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--jobs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--execute',action='store_true');p.add_argument('--freeze',type=Path);p.add_argument('--deadline-utc');a=p.parse_args()
    if a.execute:
        require(a.freeze is not None and a.deadline_utc is not None,'Execution needs frozen program and explicit UTC deadline')
        print(json.dumps(execute(a.jobs.resolve(),a.out.resolve(),a.freeze.resolve(),a.deadline_utc),indent=2))
    else:
        data,jobs,sources=prepare(a.jobs.resolve(),a.out.resolve());print(json.dumps({'status':'PLAN_ONLY_NO_INFERENCE','jobs':len(jobs),'first':FIRST,'remaining_order':[g for g in data['gene_order'] if g not in FIRST],'workers_after_first_three':8,'GiB_per_job':2,'maximum_seconds_per_job':480,'scope':'One conditional panel; not a calibrated ILS test'},indent=2))
