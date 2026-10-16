#!/usr/bin/env python3
"""Prespecified Pan-callable paired sites and original-column27 LOO; plan by default."""
from pathlib import Path
import argparse, concurrent.futures, datetime, hashlib, json, os, re, signal, subprocess, sys, threading, time, traceback
AA=set('ARNDCQEGHILKMFPSTWYV')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def need(ok,why):
 if not ok:raise ValueError(why)
def fresh(p,obj):
 with Path(p).open('x') as out:json.dump(obj,out,indent=2,allow_nan=False);out.write('\n')
def update(p,obj):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def fasta(path):
 rows={};name=None
 for line in Path(path).read_text().splitlines():
  if line.startswith('>'):
   name=line[1:].strip();need(name and not any(c.isspace() for c in name) and name not in rows,'Duplicate/unsafe FASTA label');rows[name]=''
  elif line.strip():need(name is not None,'Sequence before label');rows[name]+=line.strip()
 need(rows and len({len(s) for s in rows.values()})==1,'Empty/ragged FASTA');return rows

def prepare_inputs(with_rows,without_rows,parts,cfg):
 need(list(with_rows)==cfg['with_taxa'] and list(without_rows)==cfg['without_taxa'],'Wrong source taxon labels/order')
 need(len(with_rows)==11 and len(without_rows)==10 and {len(x) for x in with_rows.values()}=={12000} and {len(x) for x in without_rows.values()}=={12000},'Original projected11/10 rows x12000 required')
 need(set(with_rows)-set(without_rows)=={cfg['panthaica']} and all(with_rows[t]==s for t,s in without_rows.items()),'Shared10 source rows differ')
 need(len(parts)==27 and len({g['gene'] for g in parts})==27,'Exactly27 unique gene intervals required');next_col=1
 for g in parts:
  need(re.fullmatch(r'[A-Za-z0-9_.-]+',g['gene']) is not None,'Unsafe gene name')
  need(g['start_1based']==next_col and g['end_1based']-g['start_1based']+1==g['length'],'Noncontiguous/incorrect gene partition');next_col=g['end_1based']+1
 need(next_col==12001,'Partition does not exhaust12000 columns')
 keep=[i for i,a in enumerate(with_rows[cfg['panthaica']]) if a in AA];remove=[i for i,a in enumerate(with_rows[cfg['panthaica']]) if a not in AA]
 need(len(remove)==233 and len(keep)==11767,'Expected exactly233 Pan noncanonical columns')
 masked={k:{t:''.join(s[i] for i in keep) for t,s in rows.items()} for k,rows in [('with',with_rows),('without',without_rows)]}
 need(all(masked['with'][t]==s for t,s in masked['without'].items()),'Paired masked shared10 differs')
 loo=[]
 for g in parts:
  a=g['start_1based']-1;b=g['end_1based'];rows={t:s[:a]+s[b:] for t,s in with_rows.items()}
  need({len(s) for s in rows.values()}=={12000-g['length']},'LOO length error')
  loo.append({'job_id':'loo_'+g['gene'],'omitted_gene':g['gene'],'removed_start_1based':a+1,'removed_end_1based':b,'removed_length':g['length'],'rows':rows})
 proof={'status':'EXACT_PRESPECIFIED_SENSITIVITY_INPUTS','canonical_AA_order':'ARNDCQEGHILKMFPSTWYV','source_columns':12000,'masked_removed_columns_1based':[i+1 for i in remove],'masked_retained_columns_1based':[i+1 for i in keep],'masked_columns':11767,'masked_old10_rows_identical':True,'original_partitions':parts,'masked_partition_counts':[{'gene':g['gene'],'original_length':g['length'],'kept_columns':sum(g['start_1based']<=i+1<=g['end_1based'] for i in keep)} for g in parts],'LOO_basis':'Original WITH Pan12000 columns, independent of the separate233-site masking sensitivity','LOO_gene_order':[g['gene'] for g in parts],'LOO_runs':27,'not_independent_gene_votes':True}
 return masked,loo,proof

def splits(tree,taxa):
 all_taxa=set(taxa);leaves=tree.get_terminals();need(len(leaves)==len(taxa) and {x.name for x in leaves}==all_taxa,'Tree leaves differ');result=set()
 for node in tree.find_clades():
  a={x.name for x in node.get_terminals()};b=all_taxa-a
  if min(len(a),len(b))>=2:result.add(min((tuple(sorted(a)),tuple(sorted(b))),key=lambda x:(len(x),x)))
 return result

def iq_argv(cfg,input_path,prefix,paired):
 cmd=[cfg['iqtree'],'-s',str(input_path),'-m','JTT+F+G4','-T','4' if paired else '1','-seed',str(20261003 if paired else 20261004),'-safe','-keep-ident','--prefix',str(prefix)]
 if paired:cmd+=['-B','1000','-bnni','--alrt','1000','-wbt']
 return cmd

def parse_result(prefix,taxa,cfg,paired,baseline):
 from Bio import Phylo
 tree_path=Path(str(prefix)+'.treefile');report=Path(str(prefix)+'.iqtree');tree=Phylo.read(str(tree_path),'newick');sp=splits(tree,taxa)
 ll=re.findall(r'Log-likelihood of the tree:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)',report.read_text());need(len(ll)==1,'Unique tree log likelihood required')
 pairs={name:(tuple(sorted(labels)) in sp if set(labels)<=set(taxa) else None) for name,labels in cfg['focal_pairs'].items()}
 result={'reported_logL':float(ll[0]),'reported_logL_text':ll[0],'tree_sha256':sha(tree_path),'tree_path':str(tree_path),'unrooted_splits':[list(s) for s in sorted(sp)],'focal_full_tree_two_leaf_splits':pairs,'pair_definition':'Exact two-leaf side of full unrooted bipartition; not quartet projection or a root claim','RF_to_primary_same_taxa':len(sp^baseline),'splits_lost_from_primary':[list(s) for s in sorted(baseline-sp)],'splits_gained_vs_primary':[list(s) for s in sorted(sp-baseline)],'likelihood_comparison_across_changed_data':False,'support_computed':paired}
 if paired:
  count=0;hits={k:0 if set(v)<=set(taxa) else None for k,v in cfg['focal_pairs'].items()}
  for bt in Phylo.parse(str(prefix)+'.ufboot','newick'):
   bsp=splits(bt,taxa);count+=1
   for key,labels in cfg['focal_pairs'].items():
    if hits[key] is not None:hits[key]+=int(tuple(sorted(labels)) in bsp)
  need(count==1000,'Exactly1000 saved bootstrap trees required');result.update(bootstrap_trees=1000,focal_pair_bootstrap_counts=hits,bootstrap_sha256=sha(str(prefix)+'.ufboot'),support_transferred_from_primary=False)
 return result

def execute(cfg,config_path):
 import Bio
 from Bio import Phylo
 root=Path(cfg['remote_root']);need(sys.platform=='linux' and Path(__file__).resolve().parent==root,'Fixed Linux root required')
 for n,item in read(root/'PROGRAM_FREEZE.json')['files'].items():need(sha(root/n)==item['sha256'],'Program/config changed: '+n)
 need(sha(cfg['iqtree'])==cfg['iqtree_sha256'] and sha(root/'exec_verified_scope.py')==cfg['guard_sha256'],'IQTREE/guard changed')
 need(set(range(8))<=set(os.sched_getaffinity(0)),'CPU0-7 unavailable')
 for k,item in cfg['sources'].items():need(sha(item['path'])==item['sha256'],'Source hash differs: '+k)
 rows={k:fasta(cfg['sources'][k]['path']) for k in ['with','without']};parts=read(cfg['sources']['partitions']['path']);masked,loo,proof=prepare_inputs(rows['with'],rows['without'],parts,cfg)
 base={k:splits(Phylo.read(cfg['sources']['primary_'+k+'_tree']['path'],'newick'),list(rows[k])) for k in ['with','without']}
 out=root/'run_01';need(not out.exists(),'Refuse existing output; inspect instead of retry');out.mkdir();start=time.monotonic();deadline=start+5400;lock=threading.Lock();stop=threading.Event();active=set();jobs=[]
 state={'status':'PREPARING_SENSITIVITY','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'coordinator_pid':os.getpid(),'Biopython_version':Bio.__version__,'completed_jobs':[],'no_auto_retry':True,'formal_P_value':None};update(out/'STATE.json',state)
 def input_job(name,seqs,paired,matrix,extra=None):
  dest=out/'inputs';dest.mkdir(exist_ok=True);path=dest/(name+'.faa')
  with path.open('x') as f:
   for t,s in seqs.items():f.write('>'+t+'\n'+s+'\n')
  job={'id':name,'input':str(path),'input_sha256':sha(path),'taxa':list(seqs),'columns':len(next(iter(seqs.values()))),'paired':paired,'matrix':matrix,'omitted_gene':None}
  if extra:job.update(extra)
  jobs.append(job)
 for k in ['with','without']:input_job('masked_'+k,masked[k],True,k)
 for item in loo:input_job(item['job_id'],item['rows'],False,'with',{k:v for k,v in item.items() if k not in ['rows','job_id']})
 proof['source_bindings']=cfg['sources'];proof['input_jobs']=jobs;fresh(out/'INPUT_PROOF.json',proof)
 fresh(out/'COMMAND_FREEZE.json',{'config_sha256':sha(config_path),'input_proof_sha256':sha(out/'INPUT_PROOF.json'),'jobs':[{'id':j['id'],'argv':iq_argv(cfg,j['input'],out/'jobs'/j['id']/'tree',j['paired'])} for j in jobs]})
 def available_gib():return next(int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:'))/1024**3
 def kill_owned():
  with lock:units=sorted(active)
  if units:subprocess.run(['systemctl','--user','kill','--kill-whom=all','--signal=KILL',*units],capture_output=True,timeout=10)
 def interrupt(sig,frame):stop.set();kill_owned();raise RuntimeError('Interrupted; no retry')
 signal.signal(signal.SIGTERM,interrupt);signal.signal(signal.SIGINT,interrupt)
 def run_job(j,cpus):
  need(not stop.is_set(),'Sibling failed');remaining=deadline-time.monotonic()-20;need(remaining>0,'Whole90min deadline reached')
  need(sha(j['input'])==j['input_sha256'],'Prepared input changed');dest=out/'jobs'/j['id'];dest.mkdir(parents=True,exist_ok=False);prefix=dest/'tree';argv=iq_argv(cfg,j['input'],prefix,j['paired']);mem=4294967296 if j['paired'] else 2147483648;unit='dt-pan27sens-'+str(os.getpid())+'-'+j['id'].replace('_','-');guard=dest/'SCOPE_RECEIPT.json'
  cmd=['systemd-run','--user','--scope','--quiet','--unit='+unit,'-p','MemoryMax='+str(mem),'-p','MemorySwapMax=0','taskset','-c',','.join(map(str,cpus)),'/usr/bin/time','-v','-o',str(dest/'TIME.txt'),cfg['python'],str(root/'exec_verified_scope.py'),'--memory-bytes',str(mem),'--cpus',str(len(cpus)),'--receipt',str(guard),'--',*argv]
  rec={'status':'FAILED_OR_PARTIAL_NO_RETRY','job':j['id'],'argv':argv,'scope_argv':cmd,'unit':unit+'.scope','started_epoch':time.time(),'max_seconds':min(1800 if j['paired'] else 900,remaining),'cpus':cpus,'memory_max_bytes':mem};child=None
  try:
   with lock:need(not stop.is_set(),'Sibling failed before launch');active.add(unit+'.scope')
   with (dest/'stdout.txt').open('x') as o,(dest/'stderr.txt').open('x') as e:
    child=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=o,stderr=e,start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':str(len(cpus)),'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'});rec['pid']=child.pid;update(dest/'RUN_RECEIPT.json',rec);code=child.wait(timeout=rec['max_seconds']);rec['exit_code']=code;need(code==0,'IQTREE nonzero')
   g=read(guard);need(g['status']=='VERIFIED_BEFORE_EXEC' and g['logical_cpu_ids']==cpus and g['memory_max_bytes']==mem and g['swap_max_bytes']==0 and g['command']==argv,'Actual scope/argv differs')
   result=parse_result(prefix,j['taxa'],cfg,j['paired'],base[j['matrix']]);result.update(job=j['id'],columns=j['columns'],omitted_gene=j['omitted_gene']);fresh(dest/'RESULT.json',result);rec.update(status='COMPLETE',result_sha256=sha(dest/'RESULT.json'))
   with lock:state['completed_jobs'].append(j['id']);update(out/'STATE.json',state)
   return result
  except BaseException:
   stop.set();rec['error']=traceback.format_exc();kill_owned()
   if child is not None and child.poll() is None:
    try:child.wait(timeout=5)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait(timeout=5)
   raise
  finally:
   rec['elapsed_seconds']=time.time()-rec['started_epoch'];update(dest/'RUN_RECEIPT.json',rec)
   with lock:active.discard(unit+'.scope')
 results={}
 try:
  need(available_gib()>=10,'Paired phase requires current available memory>=10GiB');state['status']='RUNNING_MASKED_PAIR';update(out/'STATE.json',state)
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   futures={pool.submit(run_job,jobs[i],list(range(i*4,i*4+4))):jobs[i]['id'] for i in range(2)}
   for f in concurrent.futures.as_completed(futures):results[futures[f]]=f.result()
  need(available_gib()>=18,'LOO phase requires current available memory>=18GiB');state['status']='RUNNING_ORIGINAL_COLUMN_LOO';update(out/'STATE.json',state)
  def lane(cpu):return [run_job(j,[cpu]) for j in jobs[2:][cpu::8]]
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
   futures=[pool.submit(lane,cpu) for cpu in range(8)]
   for f in concurrent.futures.as_completed(futures):
    for result in f.result():results[result['job']]=result
  ordered=[results[j['id']] for j in jobs];loo_results=ordered[2:];need(len(loo_results)==27,'Complete27 LOO required')
  summary={'status':'COMPLETE_SENSITIVITY_SCIENTIFIC_REVIEW_PENDING','paired_results':ordered[:2],'LOO_results':loo_results,'LOO_pair_retention_counts_descriptive':{name:{'retained':sum(x['focal_full_tree_two_leaf_splits'][name] is True for x in loo_results),'denominator':27} for name in cfg['focal_pairs']},'LOO_distinct_unrooted_topologies':len({json.dumps(x['unrooted_splits'],separators=(',',':')) for x in loo_results}),'formal_P_value':None,'independent_gene_trials_assumed':False,'support_in_LOO':False,'cross_matrix_or_length_logL_comparison':False,'old_full173_gate_changed':False,'limitations':['Pan-callable-site conditioning is a sensitivity, not an unbiased locus sample','LOO uses original12000-column WITH-Pan input; not combined with masking','Pair presence denotes full unrooted cherry, not taxonomy or quartet-only affinity','No independence, calibration, ILS/gene-flow mechanism or paralog resolution follows']}
  fresh(out/'RESULTS.json',summary);state['status']=summary['status']
 except BaseException:
  stop.set();kill_owned();state.update(status='FAILED_OR_PARTIAL_NO_AUTORETRY',error=traceback.format_exc());raise
 finally:
  state['elapsed_seconds']=time.monotonic()-start;update(out/'STATE.json',state);fresh(out/'OUTPUT_MANIFEST.json',{'status':state['status'],'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'}})

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=Path(__file__).with_name('CONFIG.json'));p.add_argument('--execute',action='store_true');a=p.parse_args();cfg=read(a.config)
 if a.execute:execute(cfg,a.config)
 else:print(json.dumps({'status':'PLAN_ONLY_NO_EXECUTION','paired_masked_columns':11767,'paired_deleted_columns':233,'paired_seed':20261003,'LOO_genes':27,'LOO_basis':'Original WITH-Pan12000 columns','LOO_seed':20261004,'pair_cpus':[[0,1,2,3],[4,5,6,7]],'pair_memory_GiB_each':4,'LOO_lanes':8,'LOO_CPU_each':1,'LOO_memory_GiB_each':2,'max_simultaneous_CPU':8,'minimum_available_memory_GiB_by_phase':[10,18],'seconds_caps':{'pair_each':1800,'LOO_each':900,'whole':5400},'estimated_minutes':[10,30],'estimate_uncalibrated':True,'no_retry':True},indent=2))
if __name__=='__main__':main()
