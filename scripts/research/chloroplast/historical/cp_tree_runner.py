#!/usr/bin/env python3
"""Fresh sequential CP 95/94 runs; dry-run default; mandatory bounded cgroup."""
from pathlib import Path
import argparse, collections, hashlib, io, json, os, re, shutil, signal, subprocess, time, traceback

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def verify_inputs(p):
    m=json.loads((p/'PREPARED.json').read_text());c=json.loads((p/'CONFIG.json').read_text())
    for n,r in m['files'].items():assert sha(p/n)==r['sha256'],n
    assert sha(c['iqtree'])==c['iqtree_sha256'],'IQ-TREE binary changed'
    assert m['same_columns']==124448 and m['original93_unchanged_in_both_panels'] and m['all94_rows_identical_in_95']
    return m,c
def argv(panel,inputs,out,mode,c):
    cmd=[c['iqtree'],'-s',str(inputs/(panel+'.fixed_columns.fa')),'-p',str(inputs/c['model_file']),'-T',str(c['threads']),'-safe','-keep-ident','-seed',str(c['seed']),'-o',','.join(c['outgroups']),'--prefix',str(out/panel/'main')]
    if mode=='pilot':cmd+=['-n',str(c['pilot_iterations'])]
    else:cmd+=['-B',str(c['full_bootstrap']),'-alrt',str(c['full_alrt']),'-bnni','--boot-trees']
    assert '-redo' not in cmd and '-mem' not in cmd
    return cmd
def cgroup_check(c):
    lines=Path('/proc/self/cgroup').read_text().splitlines();rel=next(l[3:] for l in lines if l.startswith('0::'))
    root=Path('/sys/fs/cgroup')/rel.lstrip('/')
    assert all((root/n).is_file() for n in ['memory.max','memory.swap.max','memory.events']),'Required cgroup memory controls missing; use the bounded launcher'
    mem=(root/'memory.max').read_text().strip();swap=(root/'memory.swap.max').read_text().strip()
    assert mem!='max' and 0<int(mem)<=c['cgroup_memory_max_bytes'],'Must run inside <=6GiB memory cgroup'
    assert swap=='0','Must disable swap in run cgroup'
    return {'path':str(root),'memory_max_bytes':int(mem),'memory_swap_max_bytes':int(swap),'cpu_max':(root/'cpu.max').read_text().strip() if (root/'cpu.max').exists() else None,'memory_events_before':(root/'memory.events').read_text()}
def limit_cpu_affinity(c):
    available=sorted(os.sched_getaffinity(0));chosen=available[:c['cpu_affinity_max_cpus']]
    assert chosen;os.sched_setaffinity(0,chosen)
    observed=sorted(os.sched_getaffinity(0));assert observed==chosen and len(observed)<=4
    return {'method':'OS sched_setaffinity inherited by all child processes/threads','logical_cpu_ids':observed,'iqtree_threads':c['threads']}
def terminate(child):
    if child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
        try:child.wait(timeout=10)
        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
def check_tips(path,expected):
    from Bio import Phylo
    tree=Phylo.read(str(path),'newick');observed=[t.name for t in tree.get_terminals()]
    assert len(observed)==len(expected) and set(observed)==set(expected),(path,'tip loss/duplication')
    return len(observed)
def run_job(panel,inputs,out,mode,c,m,timeout):
    p=out/panel;p.mkdir(exist_ok=False);cmd=argv(panel,inputs,out,mode,c)
    env=os.environ.copy();env.update(OMP_NUM_THREADS=str(c['threads']),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    r={'status':'RUNNING','panel':panel,'argv':cmd,'started_epoch':time.time(),'timeout_seconds':timeout};save(p/'RUNTIME.json',r)
    child=None
    try:
        with (p/'stdout.log').open('x') as stdout,(p/'stderr.log').open('x') as stderr:
            child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            r['pid']=child.pid;save(p/'RUNTIME.json',r)
            try:code=child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                terminate(child);code=child.returncode;r['status']='PILOT_TIME_LIMIT_NO_INFERENCE_CLAIM'
        r.update(exit_code=code,elapsed_seconds=time.time()-r['started_epoch'])
        if r['status']=='RUNNING':
            if code!=0:raise RuntimeError(f'{panel} IQ-TREE exit {code}')
            r['tree_tip_count']=check_tips(p/'main.treefile',m['panels'][panel]['tips'])
            if mode=='full':r['consensus_tip_count']=check_tips(p/'main.contree',m['panels'][panel]['tips'])
            text=(p/'main.iqtree').read_text();match=re.search(r'Log-likelihood of the tree:\s*([-+0-9.eE]+)',text)
            assert match,'Missing final likelihood';r['logL']=float(match[1])
            r['status']='PILOT_COMPLETE_TIMING_ONLY' if mode=='pilot' else 'COMPLETE_REVIEW_PENDING'
    except BaseException:
        if child is not None:terminate(child)
        r.update(status='FAILED_NO_AUTO_RESTART',error=traceback.format_exc());raise
    finally:
        r['finished_epoch']=time.time();r['files']={x.name:{'sha256':sha(x),'bytes':x.stat().st_size} for x in p.iterdir() if x.is_file() and x.name!='RUNTIME.json'}
        save(p/'RUNTIME.json',r)
    return r
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--mode',choices=['pilot','full'],default='pilot');ap.add_argument('--execute',action='store_true');a=ap.parse_args()
    a.inputs=a.inputs.resolve();a.out=a.out.resolve();m,c=verify_inputs(a.inputs)
    plan={'status':'DRY_RUN_NO_INFERENCE','mode':a.mode,'fresh_output':str(a.out),'sequential':True,'cgroup_required':{'MemoryMax':'6G','MemorySwapMax':'0'},'cpu_bound':'4 IQ-TREE threads; driver and children inherit affinity to 4 available logical CPUs','argvs':[argv(p,a.inputs,a.out,a.mode,c) for p in c['sequential_panels']]}
    if not a.execute:print(json.dumps(plan,indent=2));return
    def interrupted(signum,frame):raise InterruptedError(f'Received signal {signum}; stop child and preserve receipts')
    signal.signal(signal.SIGTERM,interrupted)
    cg=cgroup_check(c)
    cg['cpu_affinity']=limit_cpu_affinity(c)
    assert not a.out.exists(),'Existing run directory; inspect receipts, never auto-resume or overwrite'
    assert shutil.disk_usage(a.inputs).free>=c['minimum_free_disk_bytes'],'Insufficient disk'
    a.out.mkdir(parents=True,exist_ok=False)
    state={'status':'RUNNING','mode':a.mode,'started_epoch':time.time(),'input_receipt_sha256':sha(a.inputs/'PREPARED.json'),'driver_sha256':sha(__file__),'cgroup':cg,'jobs':[]};save(a.out/'RUN_STATUS.json',state)
    try:
        budget=c['pilot_total_budget_seconds']/len(c['sequential_panels']) if a.mode=='pilot' else None
        for panel in c['sequential_panels']:
            r=run_job(panel,a.inputs,a.out,a.mode,c,m,budget);state['jobs'].append(r);save(a.out/'RUN_STATUS.json',state)
        state['status']='PILOT_DONE_TIMING_ONLY_NO_AUTO_FULL' if a.mode=='pilot' else 'FULL_COMPLETE_SCIENTIFIC_REVIEW_PENDING'
    except BaseException:state.update(status='FAILED_NO_AUTO_RESTART',error=traceback.format_exc());raise
    finally:
        state['finished_epoch']=time.time();state['elapsed_seconds']=state['finished_epoch']-state['started_epoch'];cgroot=Path(cg['path'])
        state['cgroup']['memory_events_after']=(cgroot/'memory.events').read_text()
        if (cgroot/'memory.peak').exists():state['cgroup']['memory_peak_bytes']=(cgroot/'memory.peak').read_text().strip()
        save(a.out/'RUN_STATUS.json',state)
    print(json.dumps(state,indent=2))
if __name__=='__main__':main()
