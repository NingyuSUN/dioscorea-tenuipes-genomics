#!/usr/bin/env python3
"""Independent completed173-panel recount. Default is plan only; never infer/simulate."""
import sys
from pathlib import Path
from collections import Counter
from decimal import Decimal
import argparse,csv,datetime,hashlib,json,math,os,re,shlex,time,traceback
try:
    from .portable_inputs import Config, report, require_preflight
except ImportError:
    from portable_inputs import Config, report, require_preflight
PORTABLE=None

ROOT=None
AUDIT_CPU=0
IQ_BINARY=None
AA=set('ARNDCQEGHILKMFPSTWYV')
FOCAL=('Dioscorea_tenuipes_historical','Dioscorea_tokoro','Dioscorea_gracillima','Dioscorea_zingiberensis')
KEYS=('DT_TOK','DT_GRA','DT_ZIN','unresolved')
FIRST=('209041at3193','140075at3193','49273at3193')
MODELS={'JTT+F+G4','LG+F+G4','WAG+F+G4'}
FIT_SHA='0d3ca650ae998a87e036e21c9a9131d38a9a638bf23cd3d9fc67b5a25d59f716'
GEN_SEAL_SHA='1b61b0a00399650b103ea2764868a0c012c1a42a7da96cb427492f051e84b326'
JOBS_SHA='066f48ad3beea88cce8bb87027ad3dd82081585cd9f3de5ca39415bfec0decff'
INPUT_FREEZE_SHA='862a1bf89e68993d5d47216cac93912fc6f395fc9864230db651094d87165ace'
PROGRAM_FREEZE_SHA='738caf940d8574714762c93147deeca17207c8143cb2eacb75834db6a44e7c26'
IQ_SHA='a5a11554ca60b5682fd50af81301ea247c887d93a7d7014fcd61c0fced033832'
SUPPORT_SHA='b564ff3be9e7be9fbc260188182c8b95ab5e03030726f0fa5479952cd23a1f18'
def need(ok,msg):
 if not ok:raise ValueError(msg)
def sha(p):
 h=hashlib.sha256()
 with (PORTABLE.resolve(p) if PORTABLE else Path(p)).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads((PORTABLE.resolve(p) if PORTABLE else Path(p)).read_text())
def save(p,o):
 with Path(p).open('x') as f:json.dump(o,f,indent=2,allow_nan=False);f.write('\n')
def fasta(p):
 d={};k=None
 for s in (PORTABLE.resolve(p) if PORTABLE else Path(p)).read_text().splitlines():
  if s.startswith('>'):k=s[1:].split()[0];need(k not in d,'Duplicate FASTA label');d[k]=''
  elif s.strip():need(k is not None,'Sequence without header');d[k]+=s.strip()
 return d

def parse_tree(s):
 """Independent descendant-set parser; zero/nearzero edges remain structural."""
 tok=re.findall(r"'(?:[^']|'')*'|[(),:;]|[^\s(),:;]+",s);i=0;nodes=[];leaves=[]
 def peek():return tok[i] if i<len(tok) else None
 def take():
  nonlocal i
  need(i<len(tok),'Truncated Newick');x=tok[i];i+=1;return x
 def node():
  kids=[];label=None
  if peek()=='(':
   take();kids.append(node())
   while peek()==',':take();kids.append(node())
   need(take()==')' and len(kids)>=2,'Invalid internal node')
   if peek() not in (None,',',')',':',';'):label=take()
   tips=set()
   for c in kids:need(not tips&nodes[c]['tips'],'Repeated tip');tips|=nodes[c]['tips']
  else:
   label=take();need(label not in ['(',')',',',':',';'],'Invalid tip');label=label[1:-1].replace("''", "'") if label.startswith("'") else label;need(label not in leaves,'Repeated tip');leaves.append(label);tips={label}
  length=None;bound=0.
  if peek()==':':
   take();token=take();length=float(token);need(math.isfinite(length) and length>=0,'Invalid length');bound=float(Decimal(5).scaleb(Decimal(token).as_tuple().exponent-1))
  n=len(nodes);nodes.append(dict(tips=tips,children=kids,label=label,length=length,bound=bound,parent=None))
  for c in kids:nodes[c]['parent']=n
  return n
 root=node();need(take()==';' and i==len(tok),'Trailing Newick');return dict(tips=set(leaves),nodes=nodes,root=root)
def splits(t,keep=None):
 keep=t['tips'] if keep is None else set(keep);need(keep<=t['tips'],'Absent projection taxon');out=set()
 for n in t['nodes']:
  a=n['tips']&keep;b=keep-a
  if min(len(a),len(b))>=2:out.add(min(tuple(sorted(a)),tuple(sorted(b))))
 return out
def project(t):
 need(set(FOCAL)<=t['tips'],'Missing focal tips');found=set();edges=[];H,T,G,Z=FOCAL
 for n in t['nodes']:
  side=n['tips']&set(FOCAL)
  if len(side)==2:
   if H not in side:side=set(FOCAL)-side
   found.add({T:'DT_TOK',G:'DT_GRA',Z:'DT_ZIN'}[next(iter(side-{H}))]);edges.append(dict(side=sorted(n['tips']),support_label=n['label'],branch_length=n['length']))
 need(len(found)<=1,'Conflicting quartet');return (next(iter(found)) if found else 'unresolved'),edges
def threshold(counts):
 need(sum(counts.values())==1000 and all(isinstance(x,int) and x>=0 for x in counts.values()),'Invalid bootstrap denominator')
 winners=[k for k in KEYS[:3] if counts.get(k,0)>=950];need(len(winners)<=1,'Multiple950 winners');return winners[0] if winners else None

def distance(t,a,b):
 ids={next(iter(n['tips'])):i for i,n in enumerate(t['nodes']) if not n['children']}
 def up(i):
  d=e=0.;out={i:(d,e)}
  while t['nodes'][i]['parent'] is not None:
   n=t['nodes'][i];need(n['length'] is not None,'Missing distance length');d+=n['length'];e+=n['bound'];i=n['parent'];out[i]=(d,e)
  return out
 x,y=up(ids[a]),up(ids[b]);c=next(i for i in x if i in y);return x[c][0]+y[c][0],x[c][1]+y[c][1]
def compare_report_tree(r,t,aln):
 need(r['tips']<=t['tips'],'Extra report tip');reps={}
 for m in sorted(t['tips']-r['tips']):
  same=sorted(x for x in r['tips'] if aln[x]==aln[m]);need(same,'Nonidentical omitted report sequence');reps[m]=same
 need(splits(r)==splits(t,r['tips']),'Report common-tip topology differs');tips=sorted(r['tips']);mx=0.;n=0
 for i,a in enumerate(tips):
  for b in tips[i+1:]:
   x,e=distance(r,a,b);y,f=distance(t,a,b);need(abs(x-y)<=e+f+1e-12,'Report/treefile distances differ');mx=max(mx,abs(x-y));n+=1
 return dict(missing_to_identical_candidates=reps,pairs_checked=n,max_distance_difference=mx)
def parse_report(text,aln,input_path,seed,tree):
 need(text.startswith('IQ-TREE 2.4.0'),'Report version')
 need(re.search(r'^Input file name:\s*(.+)$',text,re.M)[1]==str(input_path),'Report input')
 need(int(re.search(r'Random seed number:\s*(\d+)',text)[1])==seed,'Report seed')
 need('ultrafast bootstrap (1000 replicates)' in text,'Report bootstrap setting')
 d=re.search(r'Input data:\s*(\d+) sequences with (\d+) amino-acid sites',text);need(d and (int(d[1]),int(d[2]))==(len(aln),len(next(iter(aln.values())))),'Report dimensions')
 m=re.search(r'Best-fit model according to BIC:\s*(\S+)',text)[1];need(m in MODELS,'Selected model outside original set')
 table=text.split('List of models sorted by BIC scores:',1)[1].split('AIC, w-AIC',1)[0];bic={}
 for line in table.splitlines():
  z=line.split()
  if z and z[0] in MODELS:need(len(z)==11,'Model table format');bic[z[0]]=float(z[-3])
 need(set(bic)==MODELS and all(math.isfinite(x) for x in bic.values()) and bic[m]<=min(bic.values())+.0011,'Reported BIC selection')
 need(re.search(r'Model of substitution:\s*(\S+)',text)[1]==m,'Model report disagreement')
 ll=float(re.search(r'Log-likelihood of the tree:\s*([\d.eE+-]+)',text)[1]);alpha=float(re.search(r'Gamma shape alpha:\s*([\d.eE+-]+)',text)[1]);need(math.isfinite(ll) and math.isfinite(alpha) and alpha>0,'Nonfinite likelihood/alpha')
 region=text.split('Tree in newick format:',1)[1].split('CONSENSUS TREE',1)[0];rt=parse_tree(region[:region.index(';')+1].strip());binding=compare_report_tree(rt,tree,aln)
 return dict(selected_model=m,reported_BIC=bic,logL=ll,alpha=alpha,report_tree_binding=binding)
def scope_check(p):
 g=read(p);need(sys.platform=='linux' and g['status']=='VERIFIED_BEFORE_EXEC' and g['memory_max_bytes']==2147483648 and g['swap_max_bytes']==0 and g['logical_cpu_ids']==[AUDIT_CPU],'Audit receipt limits')
 rel=next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'));cg=Path('/sys/fs/cgroup')/rel.lstrip('/')
 need(cg.resolve()==Path(g['cgroup']).resolve() and sorted(os.sched_getaffinity(0))==[AUDIT_CPU],'Audit actual scope')
 need((cg/'memory.max').read_text().strip()=='2147483648' and (cg/'memory.swap.max').read_text().strip()=='0','Audit effective memory')
 a=list(g['command']);b=[sys.executable,sys.argv[0],*sys.argv[1:]];a[:2]=[str(Path(x).resolve()) for x in a[:2]];b[:2]=[str(Path(x).resolve()) for x in b[:2]];need(a==b,'Audit guard command');return g

def execute(args):
 start=time.monotonic();end=start+600;out=args.out.resolve();need(not out.exists() and out!=ROOT and not ROOT.is_relative_to(out),'Fresh audit output outside source ancestry required');guard=scope_check(args.scope_receipt)
 contract_path=Path(PORTABLE.input('audit_contract')['path']);contract=read(contract_path);need(sha(PORTABLE.input('original_audit_script')['path'])==contract['audit_script_sha256'],'Original audit source does not match historical contract')
 support=Path(PORTABLE.input('support_definition',SUPPORT_SHA)['path']);need(sha(support)==SUPPORT_SHA,'Support-definition source changed')
 G=ROOT/'generation';I=ROOT/'inference_results';hashes={}
 def tick():need(time.monotonic()<end,'Audit deadline; incomplete is not PASS')
 def check(p,h):
  p=PORTABLE.resolve(p);k=str(p)
  if k not in hashes:hashes[k]=sha(p)
  need(hashes[k]==h,'SHA mismatch '+k)
 chain=read(ROOT/'control/CHAIN_EXIT.json');outer=read(ROOT/'inference/EXIT_RECEIPT.json');state=read(I/'STATE.json')
 need(chain['status']=='ONE_CONDITIONAL_PANEL_INFERRED_SCIENTIFIC_AUDIT_PENDING' and chain['completed_genes']==173,'Chain incomplete')
 need(outer['status']=='INFERENCE_WRAPPER_EXITED' and outer['exit_code']==0,'Inference wrapper incomplete')
 need(state['status']=='ONE_CONDITIONAL_PANEL_INFERRED_SCIENTIFIC_AUDIT_PENDING' and state['completed_count']==173 and len(set(state['completed_genes']))==173,'Incomplete173 state')
 check(I/'PANEL_RESULTS.json',chain['panel_results_sha256']);producer=read(I/'PANEL_RESULTS.json')
 for p,h in [(G/'INFERENCE_JOBS.json',JOBS_SHA),(G/'GENERATION_OUTPUT_SEAL.json',GEN_SEAL_SHA),(G/'PANEL_INPUT_FREEZE.json',INPUT_FREEZE_SHA),(ROOT/'inference/PROGRAM_FREEZE.json',PROGRAM_FREEZE_SHA)]:check(p,h)
 jobs=read(G/'INFERENCE_JOBS.json');gf=read(G/'PANEL_INPUT_FREEZE.json');gs=read(G/'GENERATION_OUTPUT_SEAL.json');check(gf['fitted_path'],FIT_SHA);fit=read(gf['fitted_path']);order=[g['gene'] for g in fit['gene_entries']];taxa=fit['taxa'];need(len(order)==173 and len(set(order))==173 and len(taxa)==12 and order==jobs['gene_order']==[j['gene'] for j in jobs['jobs']],'Panel/order')
 for n,h in read(ROOT/'inference/PROGRAM_FREEZE.json')['files'].items():check(ROOT/'inference'/n,h)
 mf=read(I/'OUTPUT_MANIFEST.json');sealed={};need(len(mf['files'])<10000,'Manifest size')
 for n,h in mf['files'].items():
  p=(I/n).resolve();need(not Path(n).is_absolute() and p.is_relative_to(I) and p!=I/'OUTPUT_MANIFEST.json','Manifest path');check(p,h);sealed[str(p)]=h;tick()
 binding=read(I/'INPUT_BINDING.json');need(binding['original_jobs_sha256']==JOBS_SHA and binding['driver_freeze_sha256']==PROGRAM_FREEZE_SHA and [j['gene'] for j in binding['jobs']]==order,'Executed input binding');executed={j['gene']:j for j in binding['jobs']}
 for p,h in binding['source_hashes'].items():check(p,h)
 src=fit['gene_entries'][0]['source_mask'];check(src['resolved_path'],src['sha256']);original=fasta(src['resolved_path']);need(list(original)==taxa and {len(s) for s in original.values()}=={68796},'Original full12 mask')
 need(sum(g['length'] for g in fit['gene_entries'])==68796,'Columns')
 import tskit
 remaining=[g for g in order if g not in FIRST];rows=[];details=[];warnings=[];boots_total=0;missing_total=0;reported={r['gene']:r for r in producer['results']};need(len(reported)==173,'Producer duplicate/missing genes')
 out.mkdir()
 try:
  for source,j in zip(fit['gene_entries'],jobs['jobs']):
   tick();gene=j['gene'];gd=G/'genes'/gene;D=I/'genes'/gene;prefix=D/'inferred';length=j['original_full_length'];retained=j['expected_inference_taxa'];cpu=0 if gene in FIRST else remaining.index(gene)%8
   need(all(str(D/n) in sealed for n in ['RUN_RECEIPT.json','SCOPE_RECEIPT.json','TIME.txt','inferred.log','inferred.iqtree','inferred.treefile','inferred.ufboot','GENE_RESULT.json','stdout.txt','stderr.txt']),'Required output absent from seal '+gene)
   for x in j['source_bindings']:check(x['path'],x['sha256'])
   for n in ['masked.full12.faa','masked.inference.faa','true_genealogy.trees','GENEALOGY_RECEIPT.json','generation.AA.nwk']:check(gd/n,gs['files']['genes/'+gene+'/'+n])
   check(j['input_fasta'],j['input_sha256']);aln=fasta(j['input_fasta']);masked=fasta(gd/'masked.full12.faa');a=source['source_mask']['matrix_start_1based'];b=source['source_mask']['matrix_end_1based'];old={t:s[a-1:b] for t,s in original.items()}
   need(list(masked)==taxa and list(aln)==retained==source['retained_taxa'] and {len(s) for s in masked.values()}=={length},'Taxa/length '+gene)
   need(all((x==y if x not in AA else y in AA) for t in taxa for x,y in zip(old[t],masked[t])),'Cell mask changed '+gene);need(aln=={t:s for t,s in masked.items() if set(s)&AA},'Retained-row policy');missing=sum(c not in AA for s in old.values() for c in s);missing_total+=missing
   seed=j['seed'];argv=list(source['original_inference_argv']);argv[0]=j['argv'][0]
   for flag,value in [('-s',j['input_fasta']),('-seed',str(seed)),('--prefix',executed[gene]['argv'][executed[gene]['argv'].index('--prefix')+1])]:need(argv.count(flag)==1,'Ambiguous command');argv[argv.index(flag)+1]=value
   need(PORTABLE.resolve(argv[argv.index('--prefix')+1])==prefix,'Relocated inference prefix mismatch');check(IQ_BINARY,IQ_SHA);need(executed[gene]['argv']==argv,'Original inference settings altered')
   run=read(D/'RUN_RECEIPT.json');sc=read(D/'SCOPE_RECEIPT.json');unit='dt-msc03-'+str(state['pid'])+'-'+gene+'.scope';timecmd=list(sc['command']);need(len(timecmd)==4+len(argv) and timecmd[1:3]==['-v','-o'] and PORTABLE.resolve(timecmd[3])==D/'TIME.txt' and timecmd[4:]==argv,'Recorded timing command differs')
   need(run['status']=='COMPLETE' and run['unit']==unit and run['gene']==gene and 0<run['timeout_seconds']<=480 and run['elapsed_seconds']<=run['timeout_seconds']+20,'Run incomplete/bound')
   need(sc['status']=='VERIFIED_BEFORE_EXEC' and sc['command']==timecmd and sc['memory_max_bytes']==2147483648 and sc['swap_max_bytes']==0 and sc['logical_cpu_ids']==[cpu] and Path(sc['cgroup']).name==unit,'Real command/scope differs')
   need(run['argv'][run['argv'].index('--')+1:]==timecmd,'Run wrapper differs');tm=(D/'TIME.txt').read_text();need(re.search(r'^\s*Exit status:\s*0\s*$',tm,re.M),'Real exit nonzero/missing');peak=int(re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',tm)[1]);need(peak<=2097152,'PeakRSS')
   log=(D/'inferred.log').read_text();stdout=(D/'stdout.txt').read_text();stderr=(D/'stderr.txt').read_text()
   need(log.startswith('IQ-TREE multicore version 2.4.0') and 'Analysis results written to:' in log and 'UFBoot trees:' in log and 'Date and Time:' in log and '1 threads (1 CPU cores detected)' in log,'Incomplete actual log')
   for txt in [log,stdout]:
    m=re.search(r'^Command:\s*(.*)$',txt,re.M);need(m and shlex.split(m[1])==argv,'Real log command');need(int(re.search(r'^Seed:\s*(\d+)',txt,re.M)[1])==seed,'Real log seed')
   tree=parse_tree((D/'inferred.treefile').read_text());need(tree['tips']==set(retained),'Estimated taxa');text=(D/'inferred.iqtree').read_text();mi=parse_report(text,aln,j['input_fasta'],seed,tree)
   for stream,txt in [('log',log),('report',text),('stderr',stderr)]:
    lines=[s.strip() for s in txt.splitlines() if re.match(r'^\s*(WARNING|ERROR|FATAL):',s)];need(not any(re.match(r'^(ERROR|FATAL):',s) for s in lines),'Fatal log finding')
    if lines or (stream=='stderr' and txt.strip()):warnings.append(dict(gene=gene,stream=stream,lines=lines,stderr=txt if stream=='stderr' else None))
   ts=tskit.load(gd/'true_genealogy.trees');need(ts.num_samples==12 and ts.num_trees==1 and ts.sequence_length==1,'Truth file dimensions');tt=ts.first();need(tt.num_roots==1,'Truth uncoalesced');gr=read(gd/'GENEALOGY_RECEIPT.json');labels={int(k):v for k,v in gr['sample_node_labels'].items()};need(set(labels)==set(ts.samples()) and set(labels.values())==set(taxa),'Truth labels')
   truth=dict(tips=set(taxa),nodes=[dict(tips={labels[v] for v in tt.samples(u)},label=None,length=None) for u in tt.nodes()]);need(splits(truth)==splits(parse_tree((gd/'generation.AA.nwk').read_text())),'Actual genealogy/Newick topology differs')
   true_point,_=project(truth);rf=len(splits(truth,retained)^splits(tree));eligible=set(FOCAL)<=set(retained);point,edges=project(tree) if eligible else (None,[]);counts=Counter({k:0 for k in KEYS});n=0
   with (D/'inferred.ufboot').open() as f:
    for line in f:
     if not line.strip():continue
     bt=parse_tree(line);need(bt['tips']==set(retained),'Bootstrap taxa differ '+gene)
     if eligible:counts[project(bt)[0]]+=1
     n+=1
     if n%100==0:tick()
   need(n==1000,'Bootstrap count '+gene);boots_total+=n;cd=dict(counts) if eligible else None;high=threshold(counts) if eligible else None
   actual=dict(gene=gene,taxa=len(retained),AA_columns=length,focal_eligible=eligible,true_focal_split=true_point,point_split=point,full_unrooted_RF_on_retained_taxa=rf,bootstrap_counts=cd,bootstrap_denominator=1000 if eligible else None,saved_bootstrap_trees=n,bootstrap_supported_split_at_950=high,true_tree_sha256=sha(gd/'generation.AA.nwk'),output_tree_sha256=sha(D/'inferred.treefile'),output_bootstrap_sha256=sha(D/'inferred.ufboot'))
   check(D/'GENE_RESULT.json',run['result_sha256']);pr=read(D/'GENE_RESULT.json')
   for k,v in actual.items():need(pr.get(k)==v and reported[gene].get(k)==v,'Producer disagreement '+gene+'/'+k)
   row=dict(actual);row.update(selected_model=mi['selected_model'],selected_logL=mi['logL'],selected_alpha=mi['alpha'],missing_cells=missing,CPU=cpu,peak_RSS_KiB=peak,elapsed_seconds=run['elapsed_seconds']);rows.append(row);details.append(dict(gene=gene,model_report=mi,point_projecting_edges=edges))
  need(len(rows)==173 and boots_total==173000 and missing_total==98546,'Panel totals');eligible=[r for r in rows if r['focal_eligible']];need(len(eligible)==139,'Focal denominator');cross=Counter((r['true_focal_split'],r['point_split']) for r in eligible)
  summary=dict(true_counts_eligible=dict(Counter(r['true_focal_split'] for r in eligible)),point_counts=dict(Counter(r['point_split'] for r in eligible)),bootstrap_950_counts=dict(Counter(r['bootstrap_supported_split_at_950'] or 'none' for r in eligible)),truth_estimate_cross_table=[dict(true=a,estimated=b,genes=n) for (a,b),n in sorted(cross.items())])
  for k,v in summary.items():need(producer[k]==v,'Panel producer summary mismatch '+k)
  need(producer['calibrated_P_value'] is None,'Unexpected formalP')
  binding_paths=[I/'OUTPUT_MANIFEST.json',I/'PANEL_RESULTS.json',I/'INPUT_BINDING.json',I/'STATE.json',ROOT/'inference/EXIT_RECEIPT.json',ROOT/'control/CHAIN_EXIT.json',G/'INFERENCE_JOBS.json',G/'GENERATION_OUTPUT_SEAL.json',G/'PANEL_INPUT_FREEZE.json',Path(gf['fitted_path']),contract_path,Path(__file__),support]
  result=dict(status='PASS_ONE_CONDITIONAL_PANEL_INDEPENDENT_RECOUNT',role='Design-participant result consistency; not new independent scientific triad',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,audit_guard=guard,genes=173,AA_columns=68796,original_missing_cells=missing_total,bootstrap_trees_recounted=boots_total,focal_eligible_genes=139,summary=summary,selected_models=dict(Counter(r['selected_model'] for r in rows)),RF_histogram=dict(Counter(r['full_unrooted_RF_on_retained_taxa'] for r in rows)),warnings_recorded=len(warnings),manifest_files_verified=len(sealed),sources_hashed=len(hashes),formal_P_value=None,research_jobs_started=0,portable_implementation_sha256=sha(__file__),portable_config_sha256=sha(PORTABLE.filename),source_bindings={str(p):dict(sha256=sha(p),bytes=PORTABLE.resolve(p).stat().st_size) for p in binding_paths},support_definition='Any resolved projected split with >=950 of all1000 bootstrap trees; no SH gate or point-match gate. Zero/nearzero branches retained structurally.',limitations=['One conditional S0 realization, not a calibrated general ILS test.','True genealogy is simulation input, not biological truth.','Parameter/root/clock/rate/linkage/refitting calibration remains outstanding.','Warnings retained, no outcome-based exclusion.'])
  save(out/'RESULT.json',result);save(out/'MODEL_AND_EDGE_DETAILS.json',details);save(out/'WARNINGS.json',warnings)
  with (out/'PER_GENE.tsv').open('x',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader()
   for r in rows:r=dict(r);r['bootstrap_counts']=json.dumps(r['bootstrap_counts'],sort_keys=True,separators=(',',':'));w.writerow(r)
  save(out/'OUTPUT_SEAL.json',dict(status=result['status'],files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}));print(json.dumps({k:result[k] for k in ['status','genes','bootstrap_trees_recounted','focal_eligible_genes','summary','warnings_recorded','elapsed_seconds']},indent=2));return result
 except BaseException:
  save(out/'FAILED_OR_INCOMPLETE.json',dict(status='FAILED_OR_INCOMPLETE_NOT_ACCEPTED',elapsed_seconds=time.monotonic()-start,completed_gene_checks=len(rows),bootstrap_trees_checked=boots_total,error=traceback.format_exc(),research_outputs_unchanged=True));raise


def plan(config_path):
 c=Config(config_path);root=c.path('panel_root');out=c.path('output_dir')
 inputs=[c.input('audit_contract'),c.input('original_audit_script'),c.input('support_definition',SUPPORT_SHA),c.input('iqtree_binary',IQ_SHA)]
 for role,rel,h in [('chain_exit','control/CHAIN_EXIT.json',None),('inference_exit','inference/EXIT_RECEIPT.json',None),('inference_state','inference_results/STATE.json',None),('jobs','generation/INFERENCE_JOBS.json',JOBS_SHA),('generation_seal','generation/GENERATION_OUTPUT_SEAL.json',GEN_SEAL_SHA),('panel_freeze','generation/PANEL_INPUT_FREEZE.json',INPUT_FREEZE_SHA),('program_freeze','inference/PROGRAM_FREEZE.json',PROGRAM_FREEZE_SHA),('output_manifest','inference_results/OUTPUT_MANIFEST.json',None),('panel_results','inference_results/PANEL_RESULTS.json',None),('input_binding','inference_results/INPUT_BINDING.json',None)]:
  inputs.append({'role':role,'path':str(root/rel),'sha256':h})
 if 'scope_receipt' in c.data.get('inputs',{}):inputs.append(c.input('scope_receipt'))
 return {'script':'audit_panel','inputs':inputs,'outputs':[str(out/n) for n in ['RESULT.json','MODEL_AND_EDGE_DETAILS.json','WARNINGS.json','PER_GENE.tsv','OUTPUT_SEAL.json']],
   'required_threads':1,'tools':{},'python_packages':['tskit'],'transitive_inputs':'The generation/program/output seals bind all per-gene FASTA, tree sequences, IQ-TREE logs/reports/bootstrap trees and receipts',
   'scientific_contract':'Completed 173-gene / 173000-bootstrap-tree panel only; fixed seed and hash bindings; audit 2 GiB / no swap / one configured CPU / 600 seconds; no jobs launched'}

def main():
 global ROOT,PORTABLE,AUDIT_CPU,IQ_BINARY
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True,type=Path);p.add_argument('--execute',action='store_true');a=p.parse_args()
 PORTABLE=Config(a.config);spec=plan(a.config)
 if not a.execute:return report(spec)
 require_preflight(spec);ROOT=PORTABLE.path('panel_root');a.out=PORTABLE.path('output_dir');a.scope_receipt=PORTABLE.input('scope_receipt')['path']
 AUDIT_CPU=PORTABLE.data.get('audit_cpu')
 if not isinstance(AUDIT_CPU,int) or isinstance(AUDIT_CPU,bool) or AUDIT_CPU<0:raise ValueError('audit_cpu must be an explicit nonnegative CPU ID')
 IQ_BINARY=PORTABLE.input('iqtree_binary',IQ_SHA)['path'];execute(a);return 0
if __name__=='__main__':raise SystemExit(main())
