#!/usr/bin/env python3
"""Frozen one-shot CP97 execution. Default prints plan; no implicit dispatch."""
import argparse,collections,datetime,hashlib,json,math,os,re,signal,subprocess,sys,time,traceback,uuid
from pathlib import Path
ROOT=Path('workspace/wgs/dt_research_p57_20260930/cp_expansion_entry/tree_runs_recovery03/without_hachi97_01')
CPUS=[0,1,2,3]
def require(ok,msg):
    if not ok: raise ValueError(msg)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()
def save(path,obj):
    with Path(path).open('x') as f: json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def parse_fasta(raw,columns,alphabet):
    require(raw.startswith(b'>'),'FASTA must begin with header')
    out=[];seen=set()
    for part in raw.split(b'>')[1:]:
        lines=part.splitlines();label=lines[0].decode('ascii')
        require(bool(re.fullmatch(r'[A-Za-z0-9_.-]+',label)) and label not in seen,'Invalid or duplicate exact FASTA label')
        seq=b''.join(lines[1:]);require(len(seq)==columns,'FASTA column count changed')
        require(set(seq)<=set(alphabet.encode('ascii')),'FASTA contains undeclared character')
        seen.add(label);out.append((label,seq,b'>'+part))
    return out
def validate_deletion(raw98,raw97,columns,labels,drop,alphabet):
    a=parse_fasta(raw98,columns,alphabet);b=parse_fasta(raw97,columns,alphabet)
    require([x[0] for x in a]==labels,'Original exact labels/order changed')
    require(labels.count(drop)==1,'Deletion target must occur exactly once')
    expect=[x for x in a if x[0]!=drop]
    require([x[0] for x in b]==[x[0] for x in expect],'97 labels/order differ from exact deletion')
    require(raw97==b''.join(x[2] for x in expect),'97 matrix is not byte-exact deletion of the Hachi record')
    require(all(x[1]==y[1] for x,y in zip(expect,b)),'Retained AA/NT row changed')
    return {'rows_before':len(a),'rows_after':len(b),'columns':columns,'deleted_label':drop,'retained_order_equal':True,'retained_sequences_equal':True,'record_bytes_equal':True,'alphabet_counts':dict(collections.Counter(b''.join(x[1] for x in b).decode('ascii'))),'retained_row_sha256':{x[0]:hashlib.sha256(x[1]).hexdigest() for x in b}}
def verify_bundle():
    require(Path(__file__).resolve()==ROOT/'run_cp97.py','Run only from the frozen remote root')
    freeze=json.loads((ROOT/'PROGRAM_FREEZE.json').read_text())
    require(set(freeze['files'])=={'run_cp97.py','CONFIG.json'},'Unexpected program freeze')
    for n,v in freeze['files'].items():
        p=ROOT/n;require(p.stat().st_size==v['bytes'] and sha(p)==v['sha256'],'Program freeze mismatch: '+n)
    c=json.loads((ROOT/'CONFIG.json').read_text());require(c['remote_root']==str(ROOT),'Output root changed')
    require(c['resources']=={'cpus':CPUS,'memory_bytes':6442450944,'swap_bytes':0,'minimum_mem_available_bytes':17179869184,'soft_term_seconds':14385,'hard_seconds':14400,'startup_cpu_sample_seconds':2,'maximum_mean_cpu_busy_fraction':0.2},'Resource contract changed')
    return c
def cpu_counters():
    result={}
    for line in Path('/proc/stat').read_text().splitlines():
        parts=line.split()
        if parts and parts[0] in ['cpu'+str(i) for i in CPUS]:
            v=list(map(int,parts[1:9]));result[int(parts[0][3:])]=(sum(v),v[3]+v[4])
    require(set(result)==set(CPUS),'CPU 0-3 not present');return result
def resource_preflight(c):
    require(set(CPUS)<=os.sched_getaffinity(0),'CPU 0-3 not allowed by current affinity')
    mem={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith(('MemAvailable:','MemTotal:'))}
    require(mem['MemAvailable']>=c['resources']['minimum_mem_available_bytes'],'MemAvailable below 16 GiB startup gate')
    a=cpu_counters();time.sleep(c['resources']['startup_cpu_sample_seconds']);b=cpu_counters();busy={}
    for cpu in CPUS:
        total=b[cpu][0]-a[cpu][0];idle=b[cpu][1]-a[cpu][1];require(total>0,'CPU online/clock sample failed');busy[str(cpu)]=1-idle/total
    require(sum(busy.values())/4<=c['resources']['maximum_mean_cpu_busy_fraction'],'CPU0-3 mean busy exceeds declared 20% startup gate')
    return {'mem_bytes':mem,'cpu_busy_fractions':busy,'sample_seconds':2,'gate_is_not_a_resource_reservation':True}
def frozen_inputs(c):
    hashes={}
    for role,v in c['inputs'].items():
        p=Path(v['path']);require(p.is_file() and sha(p)==v['sha256'],'Frozen input SHA mismatch: '+role);hashes[role]=v['sha256']
    exe=c['executable'];require(os.access(exe['path'],os.X_OK) and sha(exe['path'])==exe['sha256'],'IQ-TREE executable mismatch')
    guard=c['guard'];require(sha(guard['path'])==guard['sha256'],'Resource guard mismatch')
    runtime=json.loads(Path(c['inputs']['old_runtime']['path']).read_text())
    require(runtime['exit_code']==0 and runtime['status']=='COMPLETE_REVIEW_PENDING' and runtime['tree_tip_count']==98 and runtime['consensus_tip_count']==98,'Original main98 not successful terminal')
    require(runtime['argv']==c['original_argv'],'Original successful argv changed')
    require(runtime['files']['main.treefile']['sha256']==c['inputs']['old_tree']['sha256'],'Original tree binding mismatch')
    oldpid=Path('/proc')/str(runtime['pid'])/'cmdline';old_running=False
    try: oldcmd=oldpid.read_bytes();old_running=c['original_argv'][c['original_argv'].index('--prefix')+1].encode() in oldcmd and Path(exe['path']).name.encode() in oldcmd
    except FileNotFoundError: pass
    require(not old_running,'Original main98 process is still running')
    a=Path(c['inputs']['main98']['path']).read_bytes();b=Path(c['inputs']['control97']['path']).read_bytes()
    matrix=validate_deletion(a,b,c['columns'],c['expected_labels98'],c['deleted_label'],c['alphabet'])
    require(matrix['rows_before']==98 and matrix['rows_after']==97,'Require exactly 98 -> 97 rows')
    entry=json.loads(Path(c['inputs']['entry_receipt']['path']).read_text());require(entry['matrix_hashes']['control97_without_Hachi']==c['inputs']['control97']['sha256'],'Prebuilt control97 receipt mismatch')
    argv=list(c['original_argv']);argv[argv.index('-s')+1]=c['inputs']['control97']['path'];argv[argv.index('--prefix')+1]=str(ROOT/'tree97/main')
    require(argv[argv.index('-p')+1]==c['inputs']['models']['path'],'Model path mismatch')
    require(all(x in matrix['retained_row_sha256'] for x in argv[argv.index('-o')+1].split(',')),'Original outgroups missing')
    return {'utc':utc(),'source_sha256':hashes,'iqtree_sha256':exe['sha256'],'guard_sha256':guard['sha256'],'matrix':matrix,'old_main98_exit_code':0,'old_main98_running':False,'argv':argv,'argv_changes_only':['-s','--prefix'],'new_tree_executed':False}
def process_identity(pid):
    t=(Path('/proc')/str(pid)/'stat').read_text().split(') ',1)[1].split()
    return {'pid':pid,'start_ticks':int(t[19]),'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip()}
def authorized_internal(token):
    d=json.loads((ROOT/'DISPATCH.json').read_text());require(token==d['token'],'Internal invocation not bound to dispatch')
def tree_tip_labels(s): return re.findall(r'(?<=[(,])([^():;,\s]+)\s*:',s)
def complete_outputs(labels):
    d=ROOT/'tree97';required=['main.treefile','main.contree','main.iqtree','main.log','main.ufboot','main.best_model.nex']
    require(all((d/n).is_file() and (d/n).stat().st_size for n in required),'Required IQ-TREE output missing or empty')
    for n in ['main.treefile','main.contree']:
        tips=tree_tip_labels((d/n).read_text());require(len(tips)==97 and set(tips)==set(labels),'Output exact tree tips mismatch: '+n)
    trees=[x for x in (d/'main.ufboot').read_text().splitlines() if x.strip()]
    require(len(trees)==1000,'Require 1000 bootstrap trees')
    require(all(len(tree_tip_labels(x))==97 and set(tree_tip_labels(x))==set(labels) for x in trees),'Bootstrap tree labels changed')
    require('Alignment has 97 sequences with 124448 columns' in (d/'main.log').read_text(),'IQ-TREE alignment dimensions unconfirmed')
    return {'tree_tips':97,'bootstrap_trees':1000,'scientific_comparison_audit':'NOT_PERFORMED'}
def worker(token):
    c=verify_bundle();authorized_internal(token);control=ROOT/'control';scope=json.loads((control/'SCOPE.json').read_text())
    require(scope['memory_max_bytes']==6442450944 and scope['swap_max_bytes']==0 and scope['logical_cpu_ids']==CPUS and sorted(os.sched_getaffinity(0))==CPUS,'Effective resource scope mismatch')
    actual='/sys/fs/cgroup/'+next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::')).lstrip('/')
    require(scope['cgroup']==actual,'Worker not inside verified cgroup')
    child=None;stop=[];start=time.monotonic();rc=None;error=None;validated=None;before=None
    def term(signum,frame):
        stop.append(signum)
        if child is not None and child.poll() is None: child.terminate()
    signal.signal(signal.SIGTERM,term);signal.signal(signal.SIGINT,term)
    try:
        before=frozen_inputs(c);before['resource_preflight']=resource_preflight(c);save(control/'PREFLIGHT.json',before)
        require(not stop,'Termination before tree launch');(ROOT/'tree97').mkdir(exist_ok=False)
        save(control/'TREE_COMMAND.json',{'argv':before['argv'],'utc':utc(),'program_freeze_sha256':sha(ROOT/'PROGRAM_FREEZE.json'),'one_tree_only':True})
        env=dict(os.environ);env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
        with (ROOT/'tree97/stdout.log').open('x') as out,(ROOT/'tree97/stderr.log').open('x') as err:
            child=subprocess.Popen(before['argv'],stdout=out,stderr=err,env=env)
            save(control/'IQTREE_PROCESS.json',dict(process_identity(child.pid),utc=utc()))
            rc=child.wait()
        if rc==0 and not stop: validated=complete_outputs([x for x in c['expected_labels98'] if x!=c['deleted_label']])
        require(rc==0 and not stop,'IQ-TREE nonzero/terminated; retain partial files without retry')
        after=frozen_inputs(c);require(after['source_sha256']==before['source_sha256'],'Input changed during run')
    except BaseException: error=traceback.format_exc()
    finally:
        if child is not None and child.poll() is None:
            child.terminate()
            try: child.wait(timeout=5)
            except subprocess.TimeoutExpired: child.kill();child.wait()
            rc=child.returncode
        save(control/'TREE_EXIT.json',{'status':'COMPLETE_TECHNICAL_ONLY' if error is None and rc==0 else 'FAILED_OR_TERMINATED_PARTIAL_PRESERVED','utc':utc(),'iqtree_exit_code':rc,'signals':stop,'error':error,'elapsed_seconds':time.monotonic()-start,'output_validation':validated,'scientific_pass':None,'automatic_retry':False})
    return 0 if error is None and rc==0 else 1

def supervise(token):
    c=verify_bundle();authorized_internal(token);control=ROOT/'control';control.mkdir(exist_ok=False)
    proc=None;stops=[];error=None;rc=1;start=time.monotonic();unit='dt-cp97-without-hachi-01.scope'
    def stop_owned():
        subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=TERM',unit],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=5)
    def term(signum,frame):
        stops.append(signum)
        if proc is not None and proc.poll() is None:
            try: stop_owned()
            except Exception: pass
    signal.signal(signal.SIGTERM,term);signal.signal(signal.SIGINT,term)
    try:
        cmd=['systemd-run','--user','--scope','--quiet','--unit='+unit,'-p','MemoryMax=6442450944','-p','MemorySwapMax=0','taskset','-c','0-3','/usr/bin/timeout','--signal=TERM','--kill-after=15','14385',c['python'],c['guard']['path'],'--memory-bytes','6442450944','--cpus','4','--receipt',str(control/'SCOPE.json'),'--',c['python'],str(ROOT/'run_cp97.py'),'--worker',token]
        save(control/'GUARD_COMMAND.json',{'argv':cmd,'utc':utc(),'hard_seconds_including_kill_grace':14400})
        require(not stops,'Termination before guarded launch')
        with (control/'stdout.log').open('x') as out,(control/'stderr.log').open('x') as err:
            proc=subprocess.Popen(cmd,stdout=out,stderr=err,start_new_session=True)
            save(control/'GUARD_PROCESS.json',dict(process_identity(proc.pid),utc=utc()))
            try: rc=proc.wait(timeout=14430)
            except subprocess.TimeoutExpired:
                stop_owned()
                try: rc=proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=KILL',unit],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=5)
                    os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
                raise TimeoutError('Outer wrapper timed out; no retry')
        require(rc==0 and not stops,'Guarded command did not complete successfully')
        tr=json.loads((control/'TREE_EXIT.json').read_text());require(tr['status']=='COMPLETE_TECHNICAL_ONLY','Tree result not complete')
    except BaseException: error=traceback.format_exc()
    finally:
        # If the guard failed, terminate only this exact owned scope before sealing partial files.
        if rc!=0 or error or stops:
            try: subprocess.run(['systemctl','--user','stop',unit],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=20)
            except Exception: pass
        active=[];scopepath=control/'SCOPE.json'
        if scopepath.exists():
            cg=Path(json.loads(scopepath.read_text())['cgroup'])/'cgroup.procs'
            if cg.exists(): active=cg.read_text().split()
        if active: error=(error or '')+'\nOwned cgroup still active; output manifest intentionally not sealed.'
        save(control/'EXIT_RECEIPT.json',{'status':'COMPLETE_TECHNICAL_ONLY' if rc==0 and error is None and not active else 'FAILED_OR_TERMINATED_NO_RETRY','utc':utc(),'guard_exit_code':rc,'error':error,'signals':stops,'elapsed_seconds':time.monotonic()-start,'owned_cgroup_remaining_pids':active,'automatic_retry':False,'scientific_pass':None})
        if not active:
            files=[p for folder in [control,ROOT/'tree97'] if folder.exists() for p in folder.iterdir() if p.is_file()]
            save(ROOT/'OUTPUT_MANIFEST.json',{'status':'OUTPUT_BYTES_SEALED_NOT_SCIENTIFIC_ACCEPTANCE','files':{str(p.relative_to(ROOT)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(files)}})
    return 0 if rc==0 and error is None and not stops else 1

def main():
    a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group();g.add_argument('--start',action='store_true');g.add_argument('--precheck',action='store_true');g.add_argument('--supervise');g.add_argument('--worker');x=a.parse_args()
    if not any([x.start,x.precheck,x.supervise,x.worker]):
        print(json.dumps({'status':'PLAN_ONLY_NOT_EXECUTED','remote_root':str(ROOT),'one_tree_only':'98 minus Hachi, fixed 124448 columns','cpus':CPUS,'memory_GiB':6,'swap':0,'hard_seconds':14400,'expected_hours':[2,3],'start':['workspace/miniconda3/bin/python',str(ROOT/'run_cp97.py'),'--start']}));return
    if x.worker:sys.exit(worker(x.worker))
    if x.supervise:sys.exit(supervise(x.supervise))
    c=verify_bundle();require(not (ROOT/'DISPATCH.json').exists(),'Dispatch already attempted: no restart/retry')
    info=frozen_inputs(c);info['resources']=resource_preflight(c)
    if x.precheck:print(json.dumps(dict(info,status='PRECHECK_ONLY_NO_DISPATCH'),indent=2));return
    require(not (ROOT/'tree97').exists() and not (ROOT/'control').exists(),'Output/control directory already exists; refuse checkpoint overwrite')
    token=uuid.uuid4().hex;save(ROOT/'DISPATCH.json',{'utc':utc(),'token':token,'program_freeze_sha256':sha(ROOT/'PROGRAM_FREEZE.json'),'precheck':info,'automatic_retry':False})
    try:
        with (ROOT/'launcher.stdout.log').open('x') as out,(ROOT/'launcher.stderr.log').open('x') as err:
            p=subprocess.Popen([c['python'],str(ROOT/'run_cp97.py'),'--supervise',token],stdout=out,stderr=err,start_new_session=True,close_fds=True)
        identity=dict(process_identity(p.pid),utc=utc(),status='DETACHED_SUPERVISOR_STARTED_NOT_COMPLETED')
        save(ROOT/'DISPATCH_PROCESS.json',identity);print(json.dumps(identity))
    except BaseException:
        save(ROOT/'DISPATCH_ERROR.json',{'error':traceback.format_exc(),'automatic_retry':False});raise
if __name__=='__main__':main()
