#!/usr/bin/env python3
"""Bridge the sealed sample-admission entry to the existing full IQ-TREE job function.
Dry plan by default; --execute must run on Linux in the declared 6GiB scope.
"""
from pathlib import Path
import getpass
import portable_cli as portable
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

def configure(config_path):
 c=portable.configure(config_path, ('out','entry_config','sample_parent','frozen_manifest','python','loginctl','tree_lock'))
 for key,default in [('entry_script','expanded_entry.py'),('legacy_runner','cp_tree_runner.py')]:
  c[key]=portable.resolve(c['_config'],c[key]) if key in c else Path(__file__).with_name(default).resolve()
 c['policy']['iqtree']=str(portable.resolve(c['_config'],c['policy']['iqtree']))
 c['fixed_model']['fixed_model_file_remote']=str(portable.resolve(c['_config'],c['fixed_model']['fixed_model_file_remote']))
 c['required_receipts']={k:str(portable.resolve(c['_config'],v)) for k,v in c['required_receipts'].items()}
 c['path_map']={k:str(portable.resolve(c['_config'],v)) for k,v in c.get('path_map',{}).items()}
 return c

def entry_config(c):
 ec=read(c['entry_config'])
 ec['source_profiles']={profile:{k:str(portable.resolve(c['entry_config'],v)) for k,v in paths.items()}
                        for profile,paths in ec['source_profiles'].items()}
 ec['portable_path_map']=c['path_map']
 return ec

def entry_dependencies(c):
 ec=entry_config(c);profile=c.get('source_profile','local');paths=ec['source_profiles'][profile]
 inputs={k:Path(paths[k]) for k in ec['source_hashes']}
 ledger=Path(paths['parent_ledger']);inputs['parent_ledger_receipt']=ledger/'LEDGER_RECEIPT.json'
 inputs.update({'geometry/'+n:ledger/n for n in ec['geometry_hashes']})
 for spec in ec['samples']:
  root=c['sample_parent']/spec['run'];orig=Path(ec['default_sample_parent'])/spec['run']
  mandatory=['fixed_projection/SAMPLE_BINDING.json','fixed_projection/PROJECTION_RECEIPT.json',
             'fixed_ledger_v2/LEDGER_RECEIPT.json','fixed_ledger_v2/CONFIG.json',
             'denovo_readback/STATUS.json','pilot_1m_qualityonly/STATUS.json','denovo_1m/STATUS.json']
  for name in mandatory:inputs[spec['run']+'/'+name]=root/name
  # Enumerate actual nested artifacts, not a directory mtime, for workflow invalidation.
  if root.is_dir():
   for f in sorted(root.rglob('*')):
    if f.is_file():inputs[spec['run']+'/'+str(f.relative_to(root))]=f
  ledger_config=root/'fixed_ledger_v2/CONFIG.json'
  if ledger_config.is_file():
   for key,record in read(ledger_config)['sources'].items():
    q=Path(record['path'])
    if q.is_relative_to(orig):q=root/q.relative_to(orig)
    else:
     for old,new in sorted(c['path_map'].items(),key=lambda kv:len(kv[0]),reverse=True):
      if q.is_relative_to(Path(old)):q=Path(new)/q.relative_to(Path(old));break
    inputs[spec['run']+'/external/'+key]=q
 return inputs

def plan(config_path):
 c=configure(config_path)
 inputs={'entry_script':c['entry_script'],'entry_config':c['entry_config'],'legacy_runner':c['legacy_runner'],
         'fixed_model':c['fixed_model']['fixed_model_file_remote'],'frozen_manifest':c['frozen_manifest'],
         'script':Path(__file__),'helper':Path(portable.__file__),**c['required_receipts']}
 if c['entry_config'].is_file():inputs.update(entry_dependencies(c))
 names=['RUN_STATUS.json','RESOURCE_BUDGET.json','FULL_RESULTS.tar.gz','COLLECTION_RECEIPT.json',
        'ALL_OUTPUT_FILES.json','matrices/ENTRY_RECEIPT.json','matrices/main98.fixed_columns.fa',
        'matrices/control95_without_new_tenuipes.fixed_columns.fa','matrices/control97_without_Hachi.fixed_columns.fa',
        'matrices/fixed_regions.nex','matrices/fixed_models.nex',
        'trees/main98/main.treefile','trees/control95_without_new_tenuipes/main.treefile']
 return portable.contract(c['_config'],inputs,{n:c['out']/n for n in names},
                          {'python':c['python'],'loginctl':c['loginctl'],'iqtree':c['policy']['iqtree']},c['out'],
                          resources=c['resource_budget'],latest_new_stage_start=c['latest_new_stage_start'],required_threads=4)

def prepare_entry(entry,c,out):
 # Invoke the packaged entry's explicit algorithm API, keeping provenance receipts immutable.
 ec=entry_config(c);profile=c.get('source_profile','local');samples,partial=entry.choose_samples(ec,None)
 paths,base,base95=entry.validate_known(ec,profile);extras=[];provenance=[]
 for spec in samples:
  raw,proof=entry.validate_sample(ec,paths,c['sample_parent'],spec,base);extras.append(raw);provenance.append(proof)
 matrices=entry.build_matrices(base95,extras,ec['hachi_tip'],ec['columns']);names=entry.matrix_names(samples,partial)
 require(all(data.startswith(base) for data in matrices.values()),'Original93 bytes changed')
 result={'status':'ALL_SELECTED_INPUTS_VALIDATED_MATRICES_NOT_WRITTEN','tree_run':False,'scientific_review_passed':False,
         'config_sha256':sha(c['entry_config']),'selected_samples':[x['name'] for x in samples],
         'omitted_samples':[],'explicit_partial_mode':False,'panels':names,'columns':ec['columns'],
         'unknown_N_is_never_reference_filled':True,'all3_required_by_default':True,'source_profile':profile,
         'sample_parent':str(c['sample_parent']),'old93_bytes_unchanged':True,'control95_sha256':entry.digest(base95),
         'sample_provenance':provenance,'matrix_hashes':{names[k]:entry.digest(v) for k,v in matrices.items()}}
 focal=[{'tip':'DT_local_path2','name':'hime2','source':'frozen original93','row_from':paths['original93']},
        {'tip':ec['hachi_tip'],'name':'Hachi-1','source':'DRR351751 / SAMD00444457','row_from':paths['hachi']},
        {'tip':ec['kita_tip'],'name':'Kita1','source':'Natsume2022 author-declared public CP','row_from':paths['kita']},
        {'tip':'OQ525999.1','name':'tokoro OQ525999.1','source':'Hu2023 public CP','row_from':paths['original93']}]
 focal += [{'tip':x['tip'],'name':x['name'],'source':x['run']+' / '+x['sample_accession'],'included':True,
            'row_from':next(p['projection_path'] for p in provenance if p['name']==x['name']),'author_alias':x['author_alias']} for x in samples]
 result['focal_tip_provenance']=focal
 entry.write_outputs(out,matrices,names,(Path(paths['parent_ledger'])/'regions.nex').read_bytes(),focal,result,partial)

def main():
 p=portable.parser(__doc__);p.add_argument('--subset',choices=['Mie1,Utsu']);p.add_argument('--latest-start',help='Override must exactly equal the configured timezone-aware deadline');a=p.parse_args()
 selected(a.subset);c=configure(a.config);a.latest_start=a.latest_start or c['latest_new_stage_start']
 require(datetime.datetime.fromisoformat(a.latest_start).tzinfo is not None,'latest-start needs timezone')
 out=c['out'];checked=portable.preflight(plan(a.config))
 require(c['tree_lock'].parent.is_dir(),'Tree lock parent directory must exist')
 require(not c['tree_lock'].is_relative_to(out) and str(c['tree_lock']) not in checked['inputs']+checked['tools'],'Tree lock must be separate from inputs and results')
 if not a.execute:portable.emit(checked);return
 require(sys.platform.startswith('linux'),'Linux cgroup execution required')
 if c.get('required_hostname'):require(socket.gethostname()==c['required_hostname'],'Configured execution hostname differs')
 freeze=read(c['frozen_manifest']);bound=set()
 for path,h in freeze['sha256'].items():
  f=portable.resolve(c['frozen_manifest'],path);bound.add(f);require(sha(f)==h,'Frozen bridge changed: '+str(f))
 require({Path(__file__).resolve(),Path(portable.__file__).resolve(),c['_config'],c['entry_script'],c['entry_config'],c['legacy_runner']}<=bound,'Freeze must bind bridge, helper, configuration, entry and legacy implementations')
 require(a.latest_start==c['latest_new_stage_start'],'Configured new-stage deadline changed')
 require(not out.exists(),'Fresh output required; no resume/overwrite')
 # The explicit freeze above binds all code and entry configuration before imports.
 require(sha(c['policy']['iqtree'])==c['policy']['iqtree_sha256'],'IQ-TREE binary changed')
 legacy=load('frozen_cp_tree_runner',c['legacy_runner']);entry=load('sealed_cp_entry',c['entry_script']);resources=legacy.cgroup_check(c['policy']);require(resources['memory_max_bytes']==6442450944,'Exactly6GiB required');resources['affinity']=pin_tree_cpus(c,legacy)
 require(shutil.disk_usage(portable.existing_parent(out.parent)).free>=c['policy']['minimum_free_disk_bytes'],'Insufficient disk')
 out.mkdir(parents=True,exist_ok=False);state={'status':'WAITING_TREE_LOCK','started_epoch':time.time(),'pid':os.getpid(),'start_ticks':Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],'plan':checked,'resources':resources,'jobs':[],'bridge_code_sha256':sha(__file__),'bridge_config_sha256':sha(a.config)};save(out/'RUN_STATUS.json',state)
 def interrupted(signum,frame):raise InterruptedError('Signal '+str(signum)+'; preserve outputs')
 signal.signal(signal.SIGTERM,interrupted)
 try:
  with Path(c['tree_lock']).open('a') as lock:
   fcntl.flock(lock,fcntl.LOCK_EX);before_deadline(a.latest_start);state['prerequisites']=prerequisite_check(c,a.subset)
   state['resource_budget']=resource_budget_snapshot(c);save(out/'RUN_STATUS.json',state);save(out/'RESOURCE_BUDGET.json',state['resource_budget']);require(state['resource_budget']['passes'],'Startup memory budget insufficient; preserve plan, no entry/tree and no cap change')
   linger=subprocess.run([str(c['loginctl']),'show-user',c.get('systemd_user',getpass.getuser()),'-p','Linger','--value'],capture_output=True,text=True,check=True,timeout=20).stdout.strip();require(linger=='yes','Linger must be enabled before detached run');state['linger']=linger;state['status']='RUNNING_ENTRY';save(out/'RUN_STATUS.json',state)
   prepare_entry(entry,c,out/'matrices')
   receipt,panels,proof=verify_entry(entry,c['entry_config'],out/'matrices',a.subset);state['entry_receipt_sha256']=sha(out/'matrices/ENTRY_RECEIPT.json');state['matrix_proof']=proof
   model=Path(c['fixed_model']['fixed_model_file_remote']);require(sha(model)==c['fixed_model']['fixed_model_file_sha256'],'Frozen model file changed');shutil.copyfile(model,out/'matrices/fixed_models.nex');require(sha(out/'matrices/fixed_models.nex')==c['fixed_model']['fixed_model_file_sha256'],'Model copy mismatch')
   provenance=out/'code_provenance';provenance.mkdir()
   for i,f in enumerate([Path(__file__),a.config,c['entry_script'],c['entry_config'],Path(c['legacy_runner']),Path(portable.__file__)]):shutil.copyfile(f,provenance/(str(i)+'_'+f.name))
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
