#!/usr/bin/env python3
"""Bridge the sealed sample-admission entry to the existing full IQ-TREE job function.
Dry plan by default; --execute must run on Linux in the declared 6GiB scope.
"""
from pathlib import Path
import argparse, datetime, fcntl, hashlib, importlib.util, json, os, re, shutil, signal, socket, subprocess, sys, tarfile, time, traceback

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(p,obj):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def require(v,message):
 if not v:raise ValueError(message)
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def selected(subset):
 require(subset is None,'This recovery03 scheduling entry admits only full98/control95; partial requires a separate explicit design')
 return ['Mie1','Yuga','Utsu']
def before_deadline(deadline):
 if deadline:require(time.time()<datetime.datetime.fromisoformat(deadline).timestamp(),'Latest new-stage start reached; running stages were allowed to finish')
def prerequisite_check(c,subset):
 selected(subset);p=Path(c['required_receipts']['recovery03']);r=read(p)
 require(r.get('status')=='COMPLETE_REVIEW_PENDING','Recovery03 must complete successfully; failed/running never admitted')
 expected=['DRR351754','DRR351755','DRR351756']
 require(r.get('samples')==expected and type(r.get('expected_samples')) is int and r['expected_samples']==3,'Recovery03 sample identity/count mismatch')
 require(r.get('safe_to_start_rna') is True and r.get('heavy_lock_released') is True,'Recovery03 heavy process/lock release not confirmed')
 require(r.get('newly_executed_stage_count')==9 and r.get('reused_stage_count')==6 and r.get('total_sample_stage_count')==15,'Recovery03 completed/reused counts mismatch')
 new=[('DRR351756',x) for x in ['assembly','screen','readback','projection']]+[('DRR351755',x) for x in ['quality','assembly','screen','readback','projection']]
 reused=[('DRR351754',x) for x in ['quality','assembly','screen','readback','projection']]+[('DRR351756','quality')]
 require([(x.get('run'),x.get('phase')) for x in r.get('completed',[])]==new and all(x.get('exit_code')==0 for x in r['completed']),'Recovery03 nine actual successful stage rows required')
 require([(x.get('run'),x.get('phase')) for x in r.get('reused_stages',[])]==reused,'Recovery03 six reused stage identities mismatch')
 return {'recovery03':{'path':str(p),'sha256':sha(p),'samples':expected,'new_stage_count':9,'reused_stage_count':6,'heavy_lock_released':True}}

def resource_budget_snapshot(c,meminfo_path=Path('/proc/meminfo')):
 budget=c['resource_budget'];require(budget['tree_cpu_ids']==[4,5,6,7] and budget['RNA_planned_cpu_ids']==[0,1,2,3],'Frozen CPU allocation changed')
 require(budget['tree_memory_max_bytes']==6*1024**3 and budget['rna_coordinator_memory_max_bytes']==2*1024**3 and budget['rna_stage_max_bytes']==16*1024**3,'Declared scheduling caps changed')
 require(sum(budget[x] for x in ['tree_memory_max_bytes','rna_coordinator_memory_max_bytes','rna_stage_max_bytes'])==budget['combined_planned_caps_bytes']==24*1024**3,'Planned concurrent cap sum changed')
 require(budget['startup_headroom_bytes']==2*1024**3 and budget['anonymous_memory_credit_bytes']==0,'Require2GiB startup headroom and zero anonymous-memory credit')
 require(budget['minimum_MemTotal_bytes']==30*1024**3,'Declared total-memory gate changed')
 require(budget['minimum_MemAvailable_bytes']==budget['combined_planned_caps_bytes']+budget['startup_headroom_bytes']-budget['anonymous_memory_credit_bytes']==26*1024**3,'Require exact26GiB future-cap-plus-headroom availability gate; no lowered declaration')
 actual={}
 for line in Path(meminfo_path).read_text().splitlines():
  fields=line.split()
  if fields[0] in ['MemTotal:','MemAvailable:']:
   require(len(fields)==3 and fields[2]=='kB','Unexpected meminfo unit');actual[fields[0][:-1]+'_bytes']=int(fields[1])*1024
 require(len(actual)==2,'MemTotal and MemAvailable both required')
 actual['logical_cpus_reported']=os.cpu_count();actual['driver_cpu_ids_at_gate']=sorted(os.sched_getaffinity(0))
 passes=actual['MemTotal_bytes']>=budget['minimum_MemTotal_bytes'] and actual['MemAvailable_bytes']>=budget['minimum_MemAvailable_bytes']
 return {'status':'PASS_STARTUP_BUDGET' if passes else 'BLOCKED_STARTUP_BUDGET_NO_ENTRY_OR_TREES','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'actual':actual,'planned_budget':budget,'passes':passes,'limits':'Startup check only, not an actual reservation: tree6GiB + RNA future total18GiB + headroom2GiB =26GiB; anonymous_memory_credit_bytes=0, no current RNA memory deduction. May conservatively reject feasible overlap when RNA already occupies memory. External future load is not guaranteed; no automatic cap reduction.'}

def pin_tree_cpus(c,legacy):
 expected=c['resource_budget']['tree_cpu_ids'];available=sorted(os.sched_getaffinity(0));require(set(expected)<=set(available),'Required logical CPU4..7 not available in inherited affinity')
 os.sched_setaffinity(0,expected);require(sorted(os.sched_getaffinity(0))==expected,'Could not apply exact CPU4..7')
 proof=legacy.limit_cpu_affinity(c['policy']);require(proof['logical_cpu_ids']==expected,'Legacy affinity did not preserve CPU4..7');proof['inherited_available_cpu_ids']=available
 return proof

def verify_entry(entry,config,out,subset):
 r=read(out/'ENTRY_RECEIPT.json');partial=bool(subset);names=selected(subset)
 require(r['status']==('PARTIAL_MATRICES_PREPARED_REVIEW_PENDING' if partial else 'MATRICES98_95_97_PREPARED_REVIEW_PENDING'),'Entry did not finish complete output set')
 require(r['config_sha256']==sha(config),'Entry configuration changed')
 require(r['selected_samples']==names and r['omitted_samples']==(['Yuga'] if partial else []) and r['explicit_partial_mode']==partial,'Partial mode or selected samples changed')
 require(r['old93_bytes_unchanged'] is True and r['columns']==124448 and r['matrices_written'] is True and r['tree_run'] is False,'Entry matrix invariant failed')
 c=read(config);panels=[r['panels']['main'],r['panels']['without_new']]
 require(panels==(['main97_partial_Mie1_Utsu','control95_without_new_tenuipes'] if partial else ['main98','control95_without_new_tenuipes']),'Unexpected main/control naming')
 parsed={};proof={}
 for panel in r['matrix_hashes']:
  f=out/(panel+'.fixed_columns.fa');require(f.name in r['matrix_files_completed'],'Matrix not marked fully written');require(sha(f)==r['matrix_hashes'][panel],'Matrix hash changed: '+panel)
 for panel,n in zip(panels,[95+len(names),95]):
  order,seqs,blocks=entry.fasta_bytes((out/(panel+'.fixed_columns.fa')).read_bytes(),124448)
  require(len(order)==n,'Wrong taxon count: '+panel);parsed[panel]=(order,seqs,blocks);proof[panel]={'taxa':n,'columns':124448,'tips':order,'sha256':r['matrix_hashes'][panel]}
 main,control=parsed[panels[0]],parsed[panels[1]]
 require(main[0][:95]==control[0] and all(main[2][x]==control[2][x] for x in control[0]),'Shared95 rows/bytes differ')
 require(r['control95_sha256']==c['frozen95_sha256']==proof[panels[1]]['sha256'],'Frozen95 control changed')
 require(sha(out/'fixed_regions.nex')==c['geometry_hashes']['regions.nex'],'Fixed column partition changed')
 return r,panels,proof

def collect(out):
 files={str(f.relative_to(out)):{'bytes':f.stat().st_size,'sha256':sha(f)} for f in sorted(out.rglob('*')) if f.is_file() and f.name not in ['FULL_RESULTS.tar.gz','COLLECTION_RECEIPT.json','ALL_OUTPUT_FILES.json']}
 save(out/'ALL_OUTPUT_FILES.json',{'status':'OUTPUTS_HASHED','files':files,'tree_and_bootstrap_outputs_included':True,'scientific_review_passed':False})
 archive=out/'FULL_RESULTS.tar.gz'
 with tarfile.open(archive,'x:gz') as tar:
  for name in [*files,'ALL_OUTPUT_FILES.json']:tar.add(out/name,arcname=name,recursive=False)
 save(out/'COLLECTION_RECEIPT.json',{'status':'COMPLETE_COLLECTION','archive':str(archive),'bytes':archive.stat().st_size,'sha256':sha(archive),'file_count':len(files)+1,'archive_includes_final_run_status':True})

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=Path(__file__).with_name('BRIDGE_CONFIG.json'));p.add_argument('--run-name',required=True);p.add_argument('--subset',choices=['Mie1,Utsu']);p.add_argument('--latest-start',default='2026-10-01T00:00:00+00:00',help='Timezone-aware deadline for new stages only; default09JST');p.add_argument('--execute',action='store_true');a=p.parse_args()
 require(re.fullmatch(r'[a-z0-9][a-z0-9_-]{0,40}',a.run_name),'Invalid run name');selected(a.subset)
 if a.latest_start:require(datetime.datetime.fromisoformat(a.latest_start).tzinfo is not None,'latest-start needs timezone')
 c=read(a.config);out=Path(c['remote_root'])/'cp_expansion_entry/tree_runs_recovery03'/a.run_name;entrydir=Path(c['entry_dir']);entrycmd=[c['python'],str(entrydir/'prepare_expanded_cp.py'),'--config',str(entrydir/'ENTRY_CONFIG.json'),'--source-profile','remote','--execute','--out',str(out/'matrices')]
 if a.subset:entrycmd+=['--subset',a.subset]
 plan={'status':'PREPARED_NO_DISPATCH','output':str(out),'selected_samples':selected(a.subset),'entry_command_once':entrycmd,'tree_panels':['main97_partial_Mie1_Utsu' if a.subset else 'main98','control95_without_new_tenuipes'],'no_Hachi_removed_tree':True,'resources':{'memory_bytes':6442450944,'swap_bytes':0,'cpus':4},'tree_lock':c['tree_lock'],'scheduling_revision':c['scheduling_revision'],'planned_resource_budget':c['resource_budget'],'latest_new_stage_start':a.latest_start,'formal_tree_hours':'Two-tree planning range2-8hours, uncalibrated; no guarantee before09JST; previous150s pilots do not calibrate bootstrap/BNNI.'}
 if not a.execute:print(json.dumps(plan,indent=2));return
 require(socket.gethostname()=='compute-host-1ccbf7da4d','Linux compute-host-1ccbf7da4d execution required')
 for path,h in read(Path(__file__).with_name('FROZEN_MANIFEST.json'))['sha256'].items():require(sha(path)==h,'Frozen new bridge changed: '+path)
 require(a.latest_start==c['latest_new_stage_start'],'Approved09JST new-stage deadline changed')
 require(not out.exists(),'Fresh output required; no resume/overwrite')
 for f,key in [(entrydir/'prepare_expanded_cp.py','entry_code_sha256'),(entrydir/'ENTRY_CONFIG.json','entry_config_sha256'),(Path(c['legacy_runner']),'legacy_runner_sha256')]:require(sha(f)==c[key],'Sealed code/config changed: '+str(f))
 require(sha(c['policy']['iqtree'])==c['policy']['iqtree_sha256'],'IQ-TREE binary changed')
 legacy=load('frozen_cp_tree_runner',c['legacy_runner']);entry=load('sealed_cp_entry',entrydir/'prepare_expanded_cp.py');resources=legacy.cgroup_check(c['policy']);require(resources['memory_max_bytes']==6442450944,'Exactly6GiB required');resources['affinity']=pin_tree_cpus(c,legacy)
 require(shutil.disk_usage(Path(c['remote_root'])).free>=c['policy']['minimum_free_disk_bytes'],'Insufficient disk')
 out.mkdir(parents=True);state={'status':'WAITING_TREE_LOCK','started_epoch':time.time(),'pid':os.getpid(),'start_ticks':Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],'plan':plan,'resources':resources,'jobs':[],'bridge_code_sha256':sha(__file__),'bridge_config_sha256':sha(a.config)};save(out/'RUN_STATUS.json',state)
 def interrupted(signum,frame):raise InterruptedError('Signal '+str(signum)+'; preserve outputs')
 signal.signal(signal.SIGTERM,interrupted)
 try:
  with Path(c['tree_lock']).open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX);before_deadline(a.latest_start);state['prerequisites']=prerequisite_check(c,a.subset)
   state['resource_budget']=resource_budget_snapshot(c);save(out/'RUN_STATUS.json',state);save(out/'RESOURCE_BUDGET.json',state['resource_budget']);require(state['resource_budget']['passes'],'Startup memory budget insufficient; preserve plan, no entry/tree and no cap change')
   linger=subprocess.run(['loginctl','show-user','private-organization','-p','Linger','--value'],capture_output=True,text=True,check=True,timeout=20).stdout.strip();require(linger=='yes','Linger must be enabled before detached run');state['linger']=linger;state['status']='RUNNING_ENTRY';save(out/'RUN_STATUS.json',state)
   with (out/'entry.stdout').open('x') as stdout,(out/'entry.stderr').open('x') as stderr:subprocess.run(entrycmd,stdout=stdout,stderr=stderr,check=True)
   receipt,panels,proof=verify_entry(entry,entrydir/'ENTRY_CONFIG.json',out/'matrices',a.subset);state['entry_receipt_sha256']=sha(out/'matrices/ENTRY_RECEIPT.json');state['matrix_proof']=proof
   model=Path(c['fixed_model']['fixed_model_file_remote']);require(sha(model)==c['fixed_model']['fixed_model_file_sha256'],'Frozen model file changed');shutil.copyfile(model,out/'matrices/fixed_models.nex');require(sha(out/'matrices/fixed_models.nex')==c['fixed_model']['fixed_model_file_sha256'],'Model copy mismatch')
   provenance=out/'code_provenance';provenance.mkdir()
   for f in [Path(__file__),a.config,entrydir/'prepare_expanded_cp.py',entrydir/'ENTRY_CONFIG.json',Path(c['legacy_runner'])]:shutil.copyfile(f,provenance/f.name)
   (out/'trees').mkdir();state['commands']=[legacy.argv(n,out/'matrices',out/'trees','full',c['policy']) for n in panels];state['status']='RUNNING_TREES';save(out/'RUN_STATUS.json',state)
   for panel in panels:
    before_deadline(a.latest_start);require(sha(out/'matrices'/(panel+'.fixed_columns.fa'))==proof[panel]['sha256'],'Frozen matrix changed before tree');require(sha(out/'matrices/fixed_models.nex')==c['fixed_model']['fixed_model_file_sha256'],'Model changed before tree')
    job=legacy.run_job(panel,out/'matrices',out/'trees','full',c['policy'],{'panels':proof},None)
    boot=out/'trees'/panel/'main.ufboot';require(boot.is_file() and boot.stat().st_size>0,'Bootstrap trees missing');job['bootstrap_tree_count']=boot.read_text().count(';');require(job['bootstrap_tree_count']==c['policy']['full_bootstrap'],'Bootstrap tree count mismatch');state['jobs'].append(job);save(out/'RUN_STATUS.json',state)
   require(len(state['jobs'])==2 and all(x['status']=='COMPLETE_REVIEW_PENDING' for x in state['jobs']),'Not both full trees completed');state['status']='TWO_TREES_COMPLETE_SCIENTIFIC_REVIEW_PENDING'
 except BaseException:state.update(status='FAILED_PRESERVED_NO_AUTO_RETRY',error=traceback.format_exc());raise
 finally:
  state['finished_epoch']=time.time();state['elapsed_seconds']=state['finished_epoch']-state['started_epoch'];cg=Path(resources['path']);state['resources']['memory_events_after']=(cg/'memory.events').read_text();state['resources']['memory_peak_bytes']=(cg/'memory.peak').read_text().strip() if (cg/'memory.peak').exists() else None;save(out/'RUN_STATUS.json',state)
  try:collect(out)
  except BaseException as error:save(out/'COLLECTION_FAILED.json',{'status':'COLLECTION_FAILED_OUTPUTS_PRESERVED','error':repr(error)})
 print(json.dumps(state,indent=2))
if __name__=='__main__':main()
