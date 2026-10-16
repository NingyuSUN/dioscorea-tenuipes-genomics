#!/usr/bin/env python3
"""Independent matched27 supplement; default only prints a plan, never executes."""
from pathlib import Path
import argparse, concurrent.futures, datetime, hashlib, json, os, re, signal, subprocess, sys, threading, time, traceback
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def need(ok,why):
 if not ok:raise ValueError(why)
def fresh(p,obj):
 with Path(p).open('x') as out:json.dump(obj,out,indent=2,allow_nan=False);out.write('\n')
def update(p,obj):
 p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');q.replace(p)
def fasta(path):
 rows={};name=None
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):
   name=line[1:].strip();need(bool(name) and not any(c.isspace() for c in name) and name not in rows,'Duplicate/unsafe FASTA label');rows[name]=''
  elif line.strip():need(name is not None,'Sequence before label');rows[name]+=line.strip()
 need(rows and len({len(s) for s in rows.values()})==1,'Empty/ragged alignment');return rows

def project(with_rows,without_rows,cfg,partitions):
 need(list(without_rows)==cfg['original_without_taxa'] and list(with_rows)==cfg['original_with_taxa'],'Source row labels/order differ')
 need(len(with_rows)==13 and len(without_rows)==12,'Original13/12 rows required')
 need({len(s) for s in with_rows.values()}=={12000} and {len(s) for s in without_rows.values()}=={12000},'Fixed12000 columns required')
 need(all(with_rows[t]==seq for t,seq in without_rows.items()),'Shared original12 rows differ')
 need(set(with_rows)-set(without_rows)=={cfg['panthaica_taxon']},'Added taxon differs')
 need(len(partitions)==27 and len({g['gene'] for g in partitions})==27,'Exactly original27 genes required');next_col=1
 for g in partitions:
  need(g['start_1based']==next_col and g['end_1based']-g['start_1based']+1==g['length'],'Partition gaps/overlaps/length mismatch');next_col=g['end_1based']+1
 need(next_col==12001,'Partition does not exhaust12000 columns')
 drop=set(cfg['drop_taxa']);need(drop=={'Arabidopsis_thaliana','Oryza_sativa'} and drop<=set(without_rows),'Only fixed two-row deletion permitted')
 projected={'with':{t:s for t,s in with_rows.items() if t not in drop},'without':{t:s for t,s in without_rows.items() if t not in drop}}
 need(len(projected['with'])==11 and len(projected['without'])==10,'Projected11/10 rows required')
 need(all(projected['with'][t]==seq for t,seq in projected['without'].items()),'Old10 rows differ')
 proof={'status':'EXACT_ROWS_ONLY_PROJECTION','original_rows':{'with':13,'without':12},'projected_rows':{'with':11,'without':10},'columns':12000,'genes':27,'dropped_taxa':cfg['drop_taxa'],'column_reselection':False,'old10_rows_identical':True,'partitions':partitions,'row_sha256':{kind:{t:hashlib.sha256(s.encode()).hexdigest() for t,s in rows.items()} for kind,rows in projected.items()},'scope':'Independent matched27 supplement; original full173 admission threshold and nuclear reference unchanged.'}
 return projected,proof

def iq_argv(cfg,kind,phase,input_path,prefix):
 cmd=[cfg['iqtree'],'-s',str(input_path),'-m','JTT+F+G4','-T','4','-seed',str(cfg['seeds'][phase]),'-safe','-keep-ident','--prefix',str(prefix)]
 if phase=='main':cmd+=['-B','1000','-bnni','--alrt','1000','-wbt']
 return cmd

def canonical_splits(tree,taxa):
 need(len(tree.get_terminals())==len(taxa) and {t.name for t in tree.get_terminals()}==set(taxa),'Tree leaf set differs')
 all_taxa=set(taxa);splits=set()
 for node in tree.find_clades():
  side={t.name for t in node.get_terminals()};other=all_taxa-side
  if min(len(side),len(other))>=2:splits.add(min((tuple(sorted(side)),tuple(sorted(other))),key=lambda x:(len(x),x)))
 return splits

def parse_result(prefix,expected_taxa,main):
 from Bio import Phylo
 tree_path=Path(str(prefix)+'.treefile');report_path=Path(str(prefix)+'.iqtree')
 tree=Phylo.read(str(tree_path),'newick');splits=canonical_splits(tree,expected_taxa)
 hits=re.findall(r'Log-likelihood of the tree:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)',report_path.read_text())
 need(len(hits)==1,'Unique reported tree log likelihood required')
 result={'reported_logL':float(hits[0]),'reported_logL_text':hits[0],'tree_path':str(tree_path),'tree_sha256':sha(tree_path),'report_sha256':sha(report_path),'unrooted_splits':[list(x) for x in sorted(splits)],'rooting':'Arbitrary Newick orientation; no designated outgroup'}
 if main:
  boot=Path(str(prefix)+'.ufboot');count=0
  for bt in Phylo.parse(str(boot),'newick'):canonical_splits(bt,expected_taxa);count+=1
  need(count==1000,'Require exactly1000 full-tip bootstrap trees');result.update(bootstrap_trees=count,bootstrap_sha256=sha(boot),support_attached_only_to_this_main_tree=True)
 return result

def execute(cfg,config_path):
 import Bio
 from Bio import Phylo
 root=Path(cfg['remote_root']);need(sys.platform=='linux' and Path(__file__).resolve().parent==root,'Fixed Linux deployment root required')
 need(sha(cfg['iqtree'])==cfg['iqtree_sha256'] and sha(root/'exec_verified_scope.py')==cfg['guard_sha256'],'Binary/guard changed')
 for name,item in read(root/'PROGRAM_FREEZE.json')['files'].items():need(sha(root/name)==item['sha256'],'Program/configuration changed: '+name)
 mem=next(int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))
 need(mem>=10*1024**3,'Require at least10GiB currently available, not a reservation')
 need(set(range(8))<=set(os.sched_getaffinity(0)),'CPU0-7 unavailable')
 source={}
 for k,item in cfg['sources'].items():need(sha(item['path'])==item['sha256'],'Source SHA changed: '+k);source[k]=Path(item['path'])
 projected,proof=project(fasta(source['with']),fasta(source['without']),cfg,read(source['partitions']))
 out=root/'run_01';need(not out.exists(),'Existing attempt: inspect, never retry or overwrite');out.mkdir()
 start=time.monotonic();deadline=start+cfg['total_timeout_seconds'];stop=threading.Event();lock=threading.Lock();units=set();state={'status':'RUNNING_MATCHED27_SUPPLEMENT','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'coordinator_pid':os.getpid(),'Biopython_version':Bio.__version__,'completed_jobs':[],'available_memory_bytes_before_launch':mem,'not_full173_admission':True,'no_AU_test':True,'no_auto_retry':True};update(out/'STATE.json',state)
 inputs={}
 for kind,rows in projected.items():
  path=out/(kind+'.faa')
  with path.open('x') as handle:
   for t,seq in rows.items():handle.write('>'+t+'\n'+seq+'\n')
  inputs[kind]=path
 proof['source_bindings']=cfg['sources'];proof['projected_matrix_sha256']={k:sha(v) for k,v in inputs.items()};fresh(out/'PROJECTION_PROOF.json',proof)
 fresh(out/'INPUT_AND_COMMAND_FREEZE.json',{'config_sha256':sha(config_path),'input_proof_sha256':sha(out/'PROJECTION_PROOF.json'),'jobs':[{'matrix':kind,'phase':phase,'argv':iq_argv(cfg,kind,phase,inputs[kind],out/kind/phase/'tree'),'cpus':cfg['cpus'][kind],'max_seconds':cfg['phase_timeout_seconds'][phase]} for kind in ['with','without'] for phase in ['main','ml_check']]})
 def kill_units():
  with lock:active=sorted(units)
  if active:subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=KILL',*active],capture_output=True,timeout=10)
 def interrupted(signum,frame):stop.set();kill_units();raise RuntimeError('Interrupted; preserve partial outputs; no retry')
 signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
 def job(kind,phase):
  need(not stop.is_set(),'Another job failed');remaining=deadline-time.monotonic()-15;need(remaining>0,'Global60min limit reached')
  dest=out/kind/phase;dest.mkdir(parents=True,exist_ok=False);prefix=dest/'tree';argv=iq_argv(cfg,kind,phase,inputs[kind],prefix);cpus=cfg['cpus'][kind];unit='dt-panth27-'+str(os.getpid())+'-'+kind+'-'+phase.replace('_','-');guard=dest/'SCOPE_RECEIPT.json'
  wrapped=['systemd-run','--user','--scope','--quiet','--unit='+unit,'-p','MemoryMax=4294967296','-p','MemorySwapMax=0','taskset','-c',','.join(map(str,cpus)),'/usr/bin/time','-v','-o',str(dest/'TIME.txt'),cfg['python'],str(root/'exec_verified_scope.py'),'--memory-bytes','4294967296','--cpus','4','--receipt',str(guard),'--',*argv]
  record={'status':'FAILED_OR_PARTIAL_NO_RETRY','argv':argv,'scoped_argv':wrapped,'unit':unit+'.scope','matrix':kind,'phase':phase,'timeout_seconds':min(cfg['phase_timeout_seconds'][phase],remaining),'started_epoch':time.time()};child=None
  try:
   with lock:need(not stop.is_set(),'Sibling failed before launch');units.add(unit+'.scope')
   with (dest/'stdout.txt').open('x') as stdout,(dest/'stderr.txt').open('x') as stderr:
    child=subprocess.Popen(wrapped,stdin=subprocess.DEVNULL,stdout=stdout,stderr=stderr,start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'});record['pid']=child.pid;update(dest/'RUN_RECEIPT.json',record)
    code=child.wait(timeout=record['timeout_seconds']);record['exit_code']=code;need(code==0,'IQ-TREE nonzero exit')
   g=read(guard);need(g['status']=='VERIFIED_BEFORE_EXEC' and g['logical_cpu_ids']==cpus and g['memory_max_bytes']==4294967296 and g['swap_max_bytes']==0 and g['command']==argv,'Actual scope/argv mismatch')
   result=parse_result(prefix,list(projected[kind]),phase=='main');fresh(dest/'RESULT.json',result);record.update(status='COMPLETE',result_sha256=sha(dest/'RESULT.json'))
   with lock:state['completed_jobs'].append(kind+'/'+phase);update(out/'STATE.json',state)
   return result
  except BaseException:
   stop.set();record['error']=traceback.format_exc();kill_units()
   if child is not None and child.poll() is None:
    try:child.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait(timeout=5)
   raise
  finally:
   record['elapsed_seconds']=time.time()-record['started_epoch'];update(dest/'RUN_RECEIPT.json',record)
   with lock:units.discard(unit+'.scope')
 def lane(kind):
  main=job(kind,'main');check=job(kind,'ml_check');delta=check['reported_logL']-main['reported_logL'];return {'matrix':kind,'main':main,'ml_check':check,'within_matrix_reported_delta_logL_check_minus_main':delta,'check_has_higher_reported_logL':delta>0,'main_and_check_unrooted_splits_equal':main['unrooted_splits']==check['unrooted_splits'],'support_policy':'Main tree and its support preserved. ML-only tree has no transferred support, even if its reported likelihood is higher.'}
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   futures={pool.submit(lane,k):k for k in ['with','without']};results={}
   for future in concurrent.futures.as_completed(futures):results[futures[future]]=future.result()
  from Bio import Phylo
  a=Phylo.read(results['with']['main']['tree_path'],'newick');a.prune(cfg['panthaica_taxon']);taxa=list(projected['without']);sa=canonical_splits(a,taxa);b=Phylo.read(results['without']['main']['tree_path'],'newick');sb=canonical_splits(b,taxa)
  summary={'status':'COMPLETE_MATCHED27_SUPPLEMENT_SCIENTIFIC_REVIEW_PENDING','matrices':results,'common10_main_unrooted_split_comparison':{'with_tree_panthaica_pruned_splits':len(sa),'without_tree_splits':len(sb),'shared_splits':len(sa&sb),'with_only_splits':[list(x) for x in sorted(sa-sb)],'without_only_splits':[list(x) for x in sorted(sb-sa)],'RF_symmetric_difference':len(sa^sb),'pruned_support_transferred':False},'cross_matrix_logL_comparison':False,'AU_performed':False,'old_full173_coverage_gate_changed':False,'limitations':['27 recovered transcript markers only; not representative full173 coverage','No new reference or original main-tree replacement','No designated outgroup or root inference','Different-tip likelihoods are not compared','ML-only higher-likelihood result does not inherit primary support']}
  fresh(out/'RESULTS.json',summary);state['status']=summary['status']
 except BaseException:
  stop.set();kill_units();state.update(status='FAILED_OR_PARTIAL_NO_AUTORETRY',error=traceback.format_exc());raise
 finally:
  state['elapsed_seconds']=time.monotonic()-start;update(out/'STATE.json',state)
  fresh(out/'OUTPUT_MANIFEST.json',{'status':state['status'],'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'}})

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=Path(__file__).with_name('CONFIG.json'));p.add_argument('--execute',action='store_true');a=p.parse_args();cfg=read(a.config)
 if a.execute:execute(cfg,a.config)
 else:print(json.dumps({'status':'PLAN_ONLY_NO_EXECUTION','source_matrices':cfg['sources'],'drop_only':cfg['drop_taxa'],'projected_tips':[11,10],'columns':12000,'genes':27,'model':'JTT+F+G4','cpus':cfg['cpus'],'memory_GiB_per_job':4,'phase_caps_seconds':cfg['phase_timeout_seconds'],'total_timeout_seconds':cfg['total_timeout_seconds'],'main_seed':cfg['seeds']['main'],'ml_check_seed':cfg['seeds']['ml_check'],'no_new_outgroup':True,'no_AU_across_matrices':True,'no_retry':True},indent=2))
if __name__=='__main__':main()
