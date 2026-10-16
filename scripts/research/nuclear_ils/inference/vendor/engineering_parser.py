#!/usr/bin/env python3
"""Exactly three frozen AA replay jobs, serial; no MSC history, retries or settings changes."""
from pathlib import Path
from collections import Counter
import datetime,hashlib,importlib.util,json,math,os,re,signal,subprocess,sys,time,traceback
P=Path(__file__).resolve().parent;AA=set('ARNDCQEGHILKMFPSTWYV');sys.dont_write_bytecode=True
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,o):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(o,indent=2)+'\n');q.replace(p)
def require(ok,msg):
 if not ok:raise ValueError(msg)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def fasta(p):
 d={};name=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];require(name not in d,'Duplicate FASTA label');d[name]=''
  elif line.strip():require(name is not None,'Sequence before header');d[name]+=line.strip()
 return d
def writefa(p,d):p.write_text(''.join('>'+t+'\n'+s+'\n' for t,s in d.items()))
def mask_digest(d):return hashlib.sha256(json.dumps([(t,''.join('1' if c in AA else c for c in s)) for t,s in d.items()],separators=(',',':')).encode()).hexdigest()
def full_splits(tree):
 tokens=re.findall(r'[(),:;]|[^(),:;\s]+',tree);i=0;clades=[];leaves=[]
 def node():
  nonlocal i
  if tokens[i]=='(':
   i+=1;s=set(node())
   while tokens[i]==',':i+=1;s.update(node())
   require(tokens[i]==')','Newick closing bracket');i+=1
   if tokens[i] not in [',',')',':',';']:i+=1
  else:s={tokens[i].strip("'")};leaves.extend(s);i+=1
  if tokens[i]==':':
   i+=1;b=float(tokens[i]);require(math.isfinite(b) and b>=0,'Invalid AA branch length');i+=1
  clades.append(frozenset(s));return s
 tips=frozenset(node());require(tokens[i]==';' and i==len(tokens)-1 and len(leaves)==len(tips),'Tree labels/syntax')
 result=set()
 for c in clades:
  other=tips-c
  if min(len(c),len(other))>=2:result.add(min(tuple(sorted(c)),tuple(sorted(other))))
 return tips,result

def cgroup_metrics(cg):
 d={}
 for name in ['memory.max','memory.swap.max','memory.current','memory.peak','memory.events','cpu.stat']:
  p=cg/name
  if p.exists():d[name]=p.read_text().strip()
 d['affinity']=sorted(os.sched_getaffinity(0));return d

def run_command(argv,folder,tag,timeout,state,deadline,cg):
 remaining=deadline-time.monotonic();require(remaining>1,'Total pilot deadline reached');limit=min(timeout,remaining)
 command=['/usr/bin/time','-v','-o',str(folder/(tag+'.time.txt'))]+argv
 record={'argv':argv,'wrapped_command':command,'timeout_seconds':limit,'started_utc':now()};started=time.monotonic()
 with (folder/(tag+'.stdout')).open('x') as out,(folder/(tag+'.stderr')).open('x') as err:
  child=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
  record['pid']=child.pid;state.update(stage=tag,active_gene=folder.name,active_pid=child.pid,active_command=argv);save(P/'STATE.json',state)
  try:code=child.wait(timeout=limit);record['exit_code']=code;require(code==0,tag+' command exit '+str(code))
  except BaseException as e:
   if child.poll() is None:
    os.killpg(child.pid,signal.SIGTERM)
    try:child.wait(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
   record.update(error=type(e).__name__+': '+str(e),exit_code=child.returncode);raise
  finally:
   record.update(finished_utc=now(),elapsed_seconds=time.monotonic()-started,cgroup=cgroup_metrics(cg));save(folder/(tag+'.RECEIPT.json'),record)
 return record

def main():
 started=time.monotonic();config=json.loads((P/'CONFIG.json').read_text());require(str(P)==config['remote_root'],'Must run in frozen Linux destination')
 frozen=json.loads((P/'EXECUTION_FREEZE.json').read_text())
 for rel,v in frozen['files'].items():require(sha(P/rel)==v['sha256'],'Frozen package SHA mismatch '+rel)
 require(sha(config['IQTREE'])==config['IQTREE_sha256'],'IQTREE binary changed')
 guard=json.loads((P/'SCOPE_RECEIPT.json').read_text());cg=Path(guard['cgroup']);metrics=cgroup_metrics(cg)
 require(metrics['memory.max']=='4294967296' and metrics['memory.swap.max']=='0' and len(metrics['affinity'])==1,'Effective resource guard mismatch')
 for item in config['original_remote_sources']:require(sha(item['path'])==item['sha256'],'Original server source changed '+item['path'])
 spec=importlib.util.spec_from_file_location('original_quartet_projector',P/'original_quartet_projector.py');proj=importlib.util.module_from_spec(spec);spec.loader.exec_module(proj);require(proj.fixtures()==6,'Original quartet fixtures')
 def interrupted(sig,frame):raise TimeoutError('Pilot interrupted by signal '+str(sig))
 signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
 output=P/'results';require(not output.exists(),'Refuse existing result directory');output.mkdir()
 state={'status':'RUNNING','pilot_kind':'AA_REPLAY_ENGINEERING','started_utc':now(),'driver_pid':os.getpid(),'scope':guard,'config_sha256':sha(P/'CONFIG.json'),'freeze_sha256':sha(P/'EXECUTION_FREEZE.json'),'completed_genes':0,'not_MSC':True,'no_auto_retry':True};save(P/'STATE.json',state);save(P/'START_RECEIPT.json',state)
 results=[];deadline=started+config['total_timeout_seconds']
 try:
  for gene in config['genes']:
   g=gene['gene'];src=P/'inputs'/g;dest=output/g;dest.mkdir();n=gene['length'];tree=src/'generation_tree.nwk';original=fasta(src/(g+'.faa'));full=fasta(src/'original_full12_mask_source.faa');tips,splits=full_splits(tree.read_text());require(tips==set(original)==set(gene['original_taxa']),'Original source tree/AA tips differ')
   require(sha(tree)==gene['tree_sha256'] and sha(src/'original_full12_mask_source.faa')==gene['full12_mask_source_sha256'],'Bound source inputs changed')
   cmd=[config['IQTREE'],'--alisim',str(dest/'generated'),'-t',str(tree),'--length',str(n),'--seqtype','AA','-m',gene['alisim_model'],'--seed',str(gene['simulation_seed']),'-T','1','-af','fasta']
   generated_receipt=run_command(cmd,dest,'generation',config['generation_timeout_seconds'],state,deadline,cg)
   generated=fasta(dest/'generated.fa');require(set(generated)==tips and {len(s) for s in generated.values()}=={n},'Generated dimensions or labels');require(all(set(s)<=AA for s in generated.values()),'Generated unexpected AA or gaps')
   masked={t:''.join(generated[t][i] if c in AA else c for i,c in enumerate(s)) if t in generated else s for t,s in full.items()}
   require(mask_digest(masked)==mask_digest(full),'Full12 missing-cell mask changed');kept={t:s for t,s in masked.items() if set(s)&AA};require(set(kept)==tips and list(kept)==list(original),'Original omitted-row policy changed')
   require(sum(sum(c not in AA for c in s) for s in masked.values())==gene['full12_missing_cells'],'Missing-cell denominator changed');writefa(dest/'masked.full12.faa',masked);writefa(dest/'masked.inference.faa',kept)
   mr={'full12_taxa':len(masked),'inference_taxa':len(kept),'AA_columns':n,'full12_total_cells':len(masked)*n,'full12_missing_cells':gene['full12_missing_cells'],'original_and_restored_mask_sha256':mask_digest(full),'generation_sha256':sha(dest/'generated.fa'),'masked_full12_sha256':sha(dest/'masked.full12.faa'),'inference_AA_sha256':sha(dest/'masked.inference.faa'),'original_tree_sha256':sha(tree),'no_branches_added_for_omitted_taxa':True};save(dest/'MASK_RECEIPT.json',mr)
   infer=list(gene['original_locus_command']);infer[0]=config['IQTREE'];infer[infer.index('-s')+1]=str(dest/'masked.inference.faa');infer[infer.index('-seed')+1]=str(gene['inference_seed']);infer[infer.index('--prefix')+1]=str(dest/'inferred')
   require(infer[infer.index('-B')+1]=='1000' and infer[infer.index('-alrt')+1]=='1000' and '-bnni' in infer and '--boot-trees' in infer and infer[infer.index('-T')+1]=='1' and infer[infer.index('-m')+1]=='MFP','Original full settings changed')
   inferred_receipt=run_command(infer,dest,'inference',config['inference_timeout_seconds'],state,deadline,cg)
   inferred_tree=dest/'inferred.treefile';inferred_tips,inferred_splits=full_splits(inferred_tree.read_text());require(inferred_tips==tips,'Inference taxon set changed')
   boot=dest/'inferred.ufboot';lines=[s for s in boot.read_text().splitlines() if s.strip()];require(len(lines)==1000,'Expected exactly1000 saved bootstrap trees')
   missing=sorted(set(proj.TAXA)-tips);point='missing_taxa';generation_point='missing_taxa';counts=None;edges=None
   if not missing:
    generation_point,_=proj.project(tree.read_text(),tips);point,edges=proj.project(inferred_tree.read_text(),tips);counts=Counter({k:0 for k in proj.KEYS})
    for line in lines:category,_=proj.project(line,tips);counts[category]+=1
   text=(dest/'inferred.iqtree').read_text();selected=re.search(r'Best-fit model according to BIC: (\S+)',text).group(1)
   result={'gene':g,'generation_model':gene['alisim_model'],'simulation_seed':gene['simulation_seed'],'inference_seed':gene['inference_seed'],'inference_selected_model':selected,'taxa':len(tips),'AA_length':n,'mask':mr,'generation_splits':[list(x) for x in sorted(splits)],'inferred_splits':[list(x) for x in sorted(inferred_splits)],'full_unrooted_RF_symmetric_difference':len(splits^inferred_splits),'generating_tree_exact_unrooted_topology_recovered':splits==inferred_splits,'focal_missing':missing,'generating_focal_projection':generation_point,'inferred_focal_projection':point,'inferred_projecting_edge_lengths':edges,'bootstrap_counts':counts,'bootstrap_denominator':1000 if not missing else None,'saved_bootstrap_trees':1000,'inference_tree_sha256':sha(inferred_tree),'bootstrap_sha256':sha(boot),'generation_elapsed_seconds':generated_receipt['elapsed_seconds'],'inference_elapsed_seconds':inferred_receipt['elapsed_seconds'],'interpretation':'Synthetic generation tree is an engineering input, not a known biological truth. Topological disagreement is recorded, not a run failure.'}
   save(dest/'RESULT.json',result);results.append(result);state.update(completed_genes=len(results),stage='GENE_COMPLETE',active_pid=None);save(P/'STATE.json',state)
  save(P/'RESULTS.json',{'pilot_kind':'AA_REPLAY_ENGINEERING','completed':len(results),'planned':3,'results':results,'not_MSC_or_ILS_test':True,'limitations':['One replicate per locus; no stable error-rate or tail-probability calibration.','No assembly/indel/alignment/trim-selection simulation; conditional on original masks and reported parameter precision.','WAG generation and frequency-special-case loci not covered.','No expanded panthaica nuclear tip; no fitting to observed116/15/8.']})
  state.update(status='COMPLETE_ENGINEERING_ONLY',result_sha256=sha(P/'RESULTS.json'))
 except BaseException as e:
  state.update(status='STOPPED_NO_RETRY',error=type(e).__name__+': '+str(e));(P/'FAILED.txt').write_text(traceback.format_exc());raise
 finally:
  state.update(finished_utc=now(),elapsed_seconds=time.monotonic()-started,final_cgroup_metrics=cgroup_metrics(cg));save(P/'STATE.json',state);save(P/'FINAL_RECEIPT.json',state)
  save(P/'RESULT_MANIFEST.json',{'files':{str(p.relative_to(P)):{'sha256':sha(p),'bytes':p.stat().st_size} for p in P.rglob('*') if p.is_file() and (output in p.parents or p.name in ['RESULTS.json','STATE.json','START_RECEIPT.json','FINAL_RECEIPT.json','SCOPE_RECEIPT.json','FAILED.txt'])}})
if __name__=='__main__':main()
