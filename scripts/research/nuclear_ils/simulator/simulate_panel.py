#!/usr/bin/env python3
"""Full12 S0 MSC-to-AA input generation only. Defaults to a source-checking plan.

No msprime import until guarded --execute. Re-inference is emitted as jobs only.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,re,signal,subprocess,sys,time,traceback
AA_ORDER='ARNDCQEGHILKMFPSTWYV';AA=set(AA_ORDER)
STAGES=('genealogy','AA','inference');SEED_MAX=2147483647

def require(ok,why):
    if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p=Path(p);temp=p.with_suffix(p.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(p)
def fresh_save(p,obj):
    require(not Path(p).exists(),'Refuse overwrite: '+str(p));save(p,obj)
def finite_positive(value,name):
    value=float(value);require(math.isfinite(value) and value>0,name+' must be positive finite');return value
def fasta(path):
    rows={};label=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            label=line[1:].split()[0];require(label not in rows,'Duplicate FASTA label');rows[label]=''
        elif line.strip():require(label is not None,'Sequence before FASTA header');rows[label]+=line.strip()
    require(bool(rows) and len({len(x) for x in rows.values()})==1,'Empty/ragged FASTA');return rows
def writefa(path,rows):
    with Path(path).open('x') as out:
        for label,seq in rows.items():out.write('>'+label+'\n'+seq+'\n')
def mask_digest(rows):
    return hashlib.sha256(json.dumps([(t,''.join('1' if c in AA else c for c in s)) for t,s in rows.items()],separators=(',',':')).encode()).hexdigest()
def bound_source(entry,base):
    require(set(entry)>={'path','sha256'},'Source requires path and sha256');p=Path(entry.get('resolved_path') or entry['path']);p=p if p.is_absolute() else base/p
    require(p.is_file() and sha(p)==entry['sha256'],'Source SHA mismatch: '+str(p));return p.resolve()

def collapse_population_nodes(nodes,taxa):
    """Contract exact same-age internal edges; do not create epsilon intervals."""
    require(len(taxa)==12 and len(set(taxa))==12,'Exactly twelve distinct labels')
    require(all(re.fullmatch(r'[A-Za-z0-9_.-]+',x) for x in taxa),'Unsupported Newick label')
    by={str(n['id']):dict(n) for n in nodes};require(len(by)==len(nodes),'Duplicate population id')
    for key,n in by.items():
        n['id']=key;n['parent_id']=None if n['parent_id'] is None else str(n['parent_id']);n['children']=[str(x) for x in n['children']]
        n['time_CU']=float(n['time_CU']);require(math.isfinite(n['time_CU']) and n['time_CU']>=0,'Bad population time')
        require(len(n['children'])==len(set(n['children'])),'Duplicate child')
    roots=[i for i,n in by.items() if n['parent_id'] is None];require(len(roots)==1,'One rooted population tree required');root=roots[0];visited=set();active=set()
    def walk(i):
        require(i in by and i not in active and i not in visited,'Cycle/shared/missing child');active.add(i);n=by[i]
        if not n['children']:
            require(len(n['descendant_taxa'])==1 and n['descendant_taxa'][0] in taxa and n['time_CU']==0,'Extant leaf must have one taxon at time0');desc=set(n['descendant_taxa'])
        else:
            require(len(n['children'])>=2,'Unary source population node');desc=set()
            for child in n['children']:
                require(child in by and by[child]['parent_id']==i,'Parent/child disagreement')
                require(n['time_CU']>=by[child]['time_CU'],'Negative population interval')
                if not by[child]['children'] or i==root:require(n['time_CU']>by[child]['time_CU'],'Terminal and root gap must be strictly positive')
                child_desc=walk(child);require(not desc&child_desc,'Repeated descendant taxon');desc|=child_desc
        require(desc==set(n['descendant_taxa']) and len(desc)==len(n['descendant_taxa']),'Descendant ledger mismatch')
        active.remove(i);visited.add(i);return desc
    require(walk(root)==set(taxa) and visited==set(by),'Incomplete/disconnected population tree')
    # Union only directly related internal populations with exactly equal times.
    parent={i:i for i in by}
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    collapsed=[]
    for i,n in by.items():
        p=n['parent_id']
        if p is not None and n['children'] and by[p]['children'] and n['time_CU']==by[p]['time_CU']:
            parent[find(i)]=find(p);collapsed.append([p,i])
    groups={}
    for i in by:groups.setdefault(find(i),[]).append(i)
    result=[]
    for representative,members in groups.items():
        # Choose the ancestral-most member of this zero-edge component.
        top=next(i for i in members if by[i]['parent_id'] not in members);n=by[top]
        kids=sorted({find(k) for i in members for k in by[i]['children'] if find(k)!=representative})
        result.append({'id':representative,'parent_id':None if n['parent_id'] is None else find(n['parent_id']),'children':kids,'descendant_taxa':list(n['descendant_taxa']),'time_CU':n['time_CU'],'source_node_ids':sorted(members)})
    result.sort(key=lambda n:(n['time_CU'],n['id']));lookup={n['id']:n for n in result}
    for n in result:
        for kid in n['children']:require(n['time_CU']>lookup[kid]['time_CU'],'Zero interval survived contraction')
    return result,collapsed

def deterministic_seeds(task,scenario,panel,genes):
    ledger={};used={}
    for gene in genes:
        ledger[gene]={}
        for stage in STAGES:
            payload=json.dumps([task,scenario,panel,gene,stage],ensure_ascii=True,separators=(',',':')).encode();digest=hashlib.sha256(payload).hexdigest();seed=int(digest,16)%SEED_MAX+1
            require(seed not in used,'Seed collision; stop rather than changing seeds');used[seed]=[gene,stage]
            ledger[gene][stage]={'seed':seed,'sha256':digest}
    return ledger

def model_string(gene):
    require(gene['model_matrix'] in ['JTT','WAG','LG'],'Unexpected source AA matrix')
    require(gene.get('frequency_order',AA_ORDER)==AA_ORDER,'AA frequency state order differs')
    freq=gene['ordered_normalized_frequencies'];require(len(freq)==20 and all(math.isfinite(float(x)) and float(x)>=0 for x in freq),'Twenty nonnegative finite AA frequencies required')
    require(abs(math.fsum(freq)-1)<1e-12 and sum(x>0 for x in freq)>=2,'Frequencies not normalized; no silent normalization')
    alpha=finite_positive(gene['alpha'],'alpha')
    return gene['model_matrix']+'+F{'+','.join(format(float(x),'.17g') for x in freq)+'}+G4{'+format(alpha,'.17g')+'}'

def original_command(gene):
    cmd=list(gene['original_inference_argv']);require(bool(cmd),'Missing original inference command')
    wanted={'-m':'MFP','-mset':'JTT,WAG,LG','-mfreq':'F','-mrate':'G4','-B':'1000','-alrt':'1000','-T':'1'}
    for key,val in wanted.items():require(cmd.count(key)==1 and cmd[cmd.index(key)+1]==val,'Original inference setting differs: '+key)
    for key in ['-bnni','-safe','--boot-trees']:require(cmd.count(key)==1,'Original inference flag missing '+key)
    for key in ['-s','-seed','--prefix']:require(cmd.count(key)==1,'Inference replacement field ambiguous')
    return cmd

def load_fitted(path,expected_sha=None):
    path=Path(path).resolve()
    if expected_sha:require(sha(path)==expected_sha,'Fitted record SHA differs')
    fitted=json.loads(path.read_text());require(fitted['schema_version']=='dioscorea.conditional_s0.v1' and fitted['scenario']=='S0' and fitted['status']=='FITTED_CONDITIONAL_SCENARIO','Requires a genuinely fitted S0 record in the agreed interface')
    taxa=fitted['taxa'];require(len(taxa)==12,'full12 required')
    nodes,collapsed=collapse_population_nodes(fitted['rooted_population_nodes'],taxa)
    kappa=finite_positive(fitted['kappa'],'kappa');genes=fitted['gene_entries'];order=[g['gene'] for g in genes]
    require(len(order)==173 and len(set(order))==173 and '2661at3193' not in order,'Original173 gene set required')
    require(sum(g['length'] for g in genes)==68796,'Original 68796 AA columns required')
    rates=fitted['gene_relative_rates'];require(set(rates)==set(order),'Missing/extra gene rates')
    for value in rates.values():finite_positive(value,'relative gene rate')
    require(abs(math.fsum(math.log(float(rates[g])) for g in order)/173)<1e-8,'Gene relative-rate geometric mean must be1')
    bindings=[];masks={};originals={};next_column=1
    for g in genes:
        require(isinstance(g['length'],int) and g['length']>0,'Invalid gene length');model_string(g);original_command(g)
        require(re.fullmatch(r'[A-Za-z0-9_.-]+',g['gene']) is not None and g['gene'] not in ['.','..'],'Unsafe gene identifier')
        sources={k:g[k] for k in ['source_locus','source_model','source_report','source_mask','source_receipt']}
        mapped={kind:bound_source(v,path.parent) for kind,v in sources.items()};matrix=fasta(mapped['source_mask']);mask=g['source_mask'];start=int(mask['matrix_start_1based']);end=int(mask['matrix_end_1based']);require(1<=start<=end<=68796 and end-start+1==g['length'],'Invalid original matrix slice');require({len(s) for s in matrix.values()}=={68796},'source_mask must be original full12 complete matrix, not ambiguous pre-slice')
        require(start==next_column,'Gene order or original mask has gap/overlap');next_column=end+1
        full={t:seq[start-1:end] for t,seq in matrix.items()};original=fasta(mapped['source_locus'])
        require(list(full)==taxa and {len(s) for s in full.values()}=={g['length']},'Full12 source mask order/dimensions differ')
        kept={t:s for t,s in full.items() if set(s)&AA}
        require(list(kept)==list(original) and kept==original,'Original omitted-row policy/source sequences differ')
        masks[g['gene']]=full;originals[g['gene']]=original;bindings.append({'gene':g['gene'],'source_bindings':sources,'resolved_paths':{k:str(v) for k,v in mapped.items()},'mask_digest':mask_digest(full)})
    require(next_column==68797,'Mask does not exhaust original matrix')
    return fitted,{'taxa12':taxa,'nodes':nodes,'collapsed_zero_edges':collapsed,'kappa':kappa,'gene_order':order,'bindings':bindings,'masks':masks,'originals':originals}

def build_demography(msprime,nodes,taxa):
    demography=msprime.Demography();names={n['id']:'P'+str(i) for i,n in enumerate(nodes)}
    tips={n['descendant_taxa'][0]:n['id'] for n in nodes if not n['children']}
    for n in nodes:demography.add_population(name=names[n['id']],initial_size=0.5,growth_rate=0)
    for n in sorted((n for n in nodes if n['children']),key=lambda n:(n['time_CU'],n['id'])):
        demography.add_population_split(time=n['time_CU'],derived=[names[k] for k in n['children']],ancestral=names[n['id']])
    demography.sort_events();demography.validate()
    require(all(float(x)==0 for row in demography.migration_matrix for x in row),'Unexpected migration')
    samples=[msprime.SampleSet(num_samples=1,population=names[tips[t]],time=0,ploidy=1) for t in taxa]
    population_taxa={names[tips[t]]:t for t in taxa}
    return demography,samples,population_taxa

def genealogy_newicks(ts,demography,population_taxa,taxa,kappa,rate):
    require(ts.num_samples==12 and ts.num_trees==1 and ts.sequence_length==1,'MSC did not yield one full12 genealogy')
    require(ts.time_units=='generations','Unexpected software time unit')
    tree=ts.first();require(tree.num_roots==1,'Genealogy not fully coalesced')
    labels={}
    for sample in ts.samples():
        node=ts.node(sample);name=demography.populations[node.population].name
        require(name in population_taxa and node.time==0,'Unexpected sample population/time');labels[int(sample)]=population_taxa[name]
    require(set(labels.values())==set(taxa) and len(set(labels.values()))==12,'Sample labels/count differ')
    factor=finite_positive(kappa*rate,'gene CU-to-AA factor');edges=[]
    def encode(u,scale):
        children=list(tree.children(u))
        if not children:require(u in labels,'Non-sample terminal node');word=labels[u]
        else:word='('+','.join(sorted(encode(v,scale) for v in children))+')'
        if u!=tree.root:
            parent=tree.parent(u);delta=float(tree.time(parent)-tree.time(u));finite_positive(delta,'genealogy interval');length=finite_positive(delta*scale,'serialized branch length');word+=':'+format(length,'.17g')
        return word
    for u in tree.nodes():
        if u!=tree.root:
            p=tree.parent(u);delta=float(tree.time(p)-tree.time(u));finite_positive(delta,'genealogy interval');edges.append({'parent':int(p),'child':int(u),'CU':delta,'AA':finite_positive(delta*factor,'AA branch')})
    return encode(tree.root,1)+';',encode(tree.root,factor)+';',{'sample_node_labels':labels,'software_time_units':ts.time_units,'interpretation':'Software generations are dimensionless CU because global ploidy2 and every population initial_size0.5 give two-lineage mean waiting time1. No physical Ne or generation time inferred.','global_ploidy':2,'sample_ploidy':1,'population_initial_size':0.5,'CU_to_AA_factor':factor,'per_gene_total_length_rescaling':False,'edges':edges}

def apply_mask(generated,source,original):
    require(set(generated)==set(source) and {len(s) for s in generated.values()}=={len(next(iter(source.values())))},'Generated full12 labels/length mismatch')
    require(all(set(s)<=AA for s in generated.values()),'Noncanonical simulated AA')
    masked={t:''.join(generated[t][i] if c in AA else c for i,c in enumerate(seq)) for t,seq in source.items()}
    require(mask_digest(masked)==mask_digest(source),'Full-cell missing mask changed');kept={t:s for t,s in masked.items() if set(s)&AA}
    require(list(kept)==list(original),'All-missing-row deletion differs');return masked,kept

def execution_guard(contract_path,receipt_path,fitted_path,out,inference_root):
    require(sys.platform=='linux','Linux required');contract=json.loads(Path(contract_path).read_text());receipt=json.loads(Path(receipt_path).read_text())
    require(contract['status']=='AUTHORIZED_LINUX_EXECUTION','No authorized execution contract')
    require(sha(fitted_path)==contract['fitted_sha256'] and sha(__file__)==contract['script_sha256'],'Runtime source binding differs')
    root=Path(contract['linux_root']).resolve();require(str(root).startswith('/home/') and out.is_relative_to(root) and inference_root.is_relative_to(root),'Declared Linux execution root mismatch')
    require(out!=inference_root and not out.is_relative_to(inference_root) and not inference_root.is_relative_to(out),'Generation and inference output roots overlap')
    require(receipt['status']=='VERIFIED_BEFORE_EXEC','Invalid scope receipt')
    rel=next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'));cg=Path('/sys/fs/cgroup')/rel.lstrip('/')
    require(cg.resolve()==Path(receipt['cgroup']).resolve(),'Stale/different cgroup receipt')
    require((cg/'memory.max').read_text().strip()=='4294967296' and (cg/'memory.swap.max').read_text().strip()=='0','Require effective4GiB/swap0')
    require(receipt['memory_max_bytes']==4294967296 and receipt['swap_max_bytes']==0,'Receipt memory binding mismatch')
    require(sorted(os.sched_getaffinity(0))==receipt['logical_cpu_ids']==contract['logical_cpu_ids'] and len(receipt['logical_cpu_ids'])==1,'Require exact declared one-CPU binding')
    expected=[str(Path(sys.executable).resolve()),str(Path(sys.argv[0]).resolve()),*sys.argv[1:]];actual=list(receipt['command']);actual[:2]=[str(Path(x).resolve()) for x in actual[:2]]
    require(actual==expected,'Scope was not created for this exact command')
    require(sha(contract['iqtree'])==contract['iqtree_sha256'],'IQTREE binary changed')
    require(0<contract['generation_timeout_seconds']<=600 and 0<contract['total_timeout_seconds']<=3600,'Undeclared/excess execution timeout')
    return contract,{'contract_sha256':sha(contract_path),'scope_receipt_sha256':sha(receipt_path),'current_cgroup':str(cg),'scope':receipt}

def timed_alisim(command,dest,seconds):
    record={'argv':command,'started_epoch':time.time(),'timeout_seconds':seconds}
    with (dest/'alisim.stdout').open('x') as stdout,(dest/'alisim.stderr').open('x') as stderr:
        p=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'});record['pid']=p.pid
        try:
            p.wait(timeout=seconds);require(p.returncode==0,'AliSim returned nonzero')
        except BaseException:
            if p.poll() is None:
                os.killpg(p.pid,signal.SIGTERM)
                try:p.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
            record['error']=traceback.format_exc();raise
        finally:record.update(exit_code=p.returncode,elapsed_seconds=time.time()-record['started_epoch']);save(dest/'ALISIM_RECEIPT.json',record)

def execute(args,fitted,checked,seeds):
    out=args.out.resolve();inference_root=args.inference_root.resolve();require(not out.exists(),'Refuse existing generation directory')
    require(not inference_root.exists() or not any(inference_root.iterdir()),'Refuse existing inference outputs')
    require(args.fitted_sha256 and args.execution_contract and args.scope_receipt,'Execute requires fitted SHA, execution contract and live scope receipt')
    contract,guard=execution_guard(args.execution_contract,args.scope_receipt,args.fitted,out,inference_root)
    # Actual dependency/API validation is intentionally deferred until authorized Linux execution.
    import msprime
    demography,samples,population_taxa=build_demography(msprime,checked['nodes'],checked['taxa12'])
    out.mkdir(parents=True);start=time.monotonic();deadline=start+contract['total_timeout_seconds'];state={'status':'GENERATING_INPUTS_NOT_A_REINFERRED_PANEL','scenario_id':'S0','panel_id':args.panel_id,'completed_genes':0,'inference_completed':0,'formal_P_value':None,'driver_pid':os.getpid(),'msprime_version':msprime.__version__,'guard':guard};save(out/'STATE.json',state)
    plan={'task_id':args.task_id,'fitted_path':str(args.fitted.resolve()),'fitted_sha256':sha(args.fitted),'script_sha256':sha(__file__),'scenario_id':'S0','panel_id':args.panel_id,'taxa12':checked['taxa12'],'gene_order':checked['gene_order'],'gene_sources':checked['bindings'],'seeds':seeds,'collapsed_nodes':checked['nodes'],'collapsed_zero_edges':checked['collapsed_zero_edges'],'demography_recipe':{'populations':checked['nodes'],'migration_matrix':'all zero; no migration events','constant_population_initial_size':0.5},'initial_size':0.5,'global_ploidy':2,'SampleSet_ploidy':1,'migration_rate':0,'recombination_rate':0,'sequence_length':1,'model':'hudson'}
    save(out/'PANEL_INPUT_FREEZE.json',plan);jobs=[]
    try:
        for entry in fitted['gene_entries']:
            require(time.monotonic()<deadline,'Panel generation deadline reached');g=entry['gene'];dest=out/'genes'/g;dest.mkdir(parents=True)
            ts=msprime.sim_ancestry(samples=samples,demography=demography,ploidy=2,model='hudson',recombination_rate=0,sequence_length=1,random_seed=seeds[g]['genealogy']['seed'],record_migrations=False)
            require(time.monotonic()<deadline,'Deadline exceeded during genealogy');cu,aa,receipt=genealogy_newicks(ts,demography,population_taxa,checked['taxa12'],checked['kappa'],fitted['gene_relative_rates'][g]);ts.dump(str(dest/'true_genealogy.trees'))
            (dest/'true_genealogy.CU.nwk').write_text(cu+'\n');(dest/'generation.AA.nwk').write_text(aa+'\n');save(dest/'GENEALOGY_RECEIPT.json',receipt)
            command=[contract['iqtree'],'--alisim',str(dest/'generated'),'-t',str(dest/'generation.AA.nwk'),'--length',str(entry['length']),'--seqtype','AA','-m',model_string(entry),'--seed',str(seeds[g]['AA']['seed']),'-T','1','-af','fasta']
            remaining=deadline-time.monotonic();require(remaining>0,'Deadline before AliSim');timed_alisim(command,dest,min(contract['generation_timeout_seconds'],remaining))
            generated_path=dest/'generated.fa';generated=fasta(generated_path);masked,kept=apply_mask(generated,checked['masks'][g],checked['originals'][g]);full=dest/'masked.full12.faa';infer=dest/'masked.inference.faa';writefa(full,masked);writefa(infer,kept)
            cmd=original_command(entry);prefix=inference_root/g/'inferred';cmd[0]=contract['iqtree'];cmd[cmd.index('-s')+1]=str(infer);cmd[cmd.index('-seed')+1]=str(seeds[g]['inference']['seed']);cmd[cmd.index('--prefix')+1]=str(prefix)
            job={'gene':g,'input_fasta':str(infer),'input_sha256':sha(infer),'expected_inference_taxa':list(kept),'original_full_length':entry['length'],'argv':cmd,'seed':seeds[g]['inference']['seed'],'output_prefix':str(prefix),'original_argv':list(entry['original_inference_argv']),'true_tree_path':str(dest/'generation.AA.nwk'),'true_tree_sha256':sha(dest/'generation.AA.nwk'),'true_tree_CU':{'path':str(dest/'true_genealogy.CU.nwk'),'sha256':sha(dest/'true_genealogy.CU.nwk')},'source_bindings':[{'role':k,'path':str(bound_source(entry[k],args.fitted.resolve().parent)),'sha256':entry[k]['sha256'],'source_path':entry[k]['path']} for k in ['source_locus','source_model','source_report','source_mask','source_receipt']],'generated_full12':{'path':str(generated_path),'sha256':sha(generated_path)},'masked_full12':{'path':str(full),'sha256':sha(full)},'mask_digest':mask_digest(masked),'true_genealogy':{'path':str(dest/'true_genealogy.trees'),'sha256':sha(dest/'true_genealogy.trees')}}
            jobs.append(job);save(dest/'INPUT_RECEIPT.json',job);state['completed_genes']=len(jobs);save(out/'STATE.json',state)
        fresh_save(out/'GENERATION_OUTPUT_SEAL.json',{'scope':'Immutable generated gene inputs and truth, before inference job table; excludes mutable state and final manifests to avoid circular hashes.','fitted_sha256':sha(args.fitted),'panel_input_freeze_sha256':sha(out/'PANEL_INPUT_FREEZE.json'),'files':{str(p.relative_to(out)):sha(p) for p in sorted((out/'genes').rglob('*')) if p.is_file()}})
        for job in jobs:
            job['source_bindings'].append({'role':'generation_output_seal','path':str(out/'GENERATION_OUTPUT_SEAL.json'),'sha256':sha(out/'GENERATION_OUTPUT_SEAL.json')})
            job['source_bindings'].append({'role':'panel_input_freeze','path':str(out/'PANEL_INPUT_FREEZE.json'),'sha256':sha(out/'PANEL_INPUT_FREEZE.json')})
            job['source_bindings'].append({'role':'execution_contract','path':str(args.execution_contract.resolve()),'sha256':sha(args.execution_contract)})
        fresh_save(out/'INFERENCE_JOBS.json',{'status':'INPUTS_GENERATED_REINFERENCE_NOT_STARTED','scenario_id':'S0','panel_id':args.panel_id,'taxa12':checked['taxa12'],'gene_order':checked['gene_order'],'jobs':jobs,'initial_reinference_order':['209041at3193','140075at3193','49273at3193'],'remaining_order_policy':'Original gene_order excluding initial three; no outcome-based selection','inference_root':str(inference_root),'prefix_relocation_policy':'Root may create separately hashed execution table replacing only the unique --prefix value; preserve this source table.','source_bindings':[{'role':'fitted','path':str(args.fitted.resolve()),'sha256':sha(args.fitted)},{'role':'panel_input_freeze','path':str(out/'PANEL_INPUT_FREEZE.json'),'sha256':sha(out/'PANEL_INPUT_FREEZE.json')},{'role':'generation_output_seal','path':str(out/'GENERATION_OUTPUT_SEAL.json'),'sha256':sha(out/'GENERATION_OUTPUT_SEAL.json')},{'role':'execution_contract','path':str(args.execution_contract.resolve()),'sha256':sha(args.execution_contract)}],'not_complete_panel_replicate':True,'formal_P_value':None})
        state['status']='ALL173_INPUTS_GENERATED_REINFERENCE_NOT_STARTED'
    except BaseException:
        state.update(status='FAILED_OR_PARTIAL_NO_RETRY',error=traceback.format_exc());save(out/'PARTIAL_INFERENCE_JOBS.json',{'completed_jobs':jobs,'not_a_complete_panel':True});raise
    finally:
        state['elapsed_seconds']=time.monotonic()-start;save(out/'STATE.json',state);save(out/'OUTPUT_MANIFEST.json',{'files':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},'status':state['status']})
    return state

def main():
    p=argparse.ArgumentParser();p.add_argument('--fitted',type=Path,required=True);p.add_argument('--fitted-sha256');p.add_argument('--panel-id',default='S0_panel_0001');p.add_argument('--task-id',default='dioscorea_full12_MSC_S0_01');p.add_argument('--out',type=Path);p.add_argument('--inference-root',type=Path);p.add_argument('--execution-contract',type=Path);p.add_argument('--scope-receipt',type=Path);p.add_argument('--execute',action='store_true');args=p.parse_args()
    fitted,checked=load_fitted(args.fitted,args.fitted_sha256);seeds=deterministic_seeds(args.task_id,'S0',args.panel_id,checked['gene_order'])
    if args.execute:
        require(args.out and args.inference_root,'Execute requires separate generation/inference paths');result=execute(args,fitted,checked,seeds)
    else:result={'status':'PLAN_AND_SOURCE_CHECK_ONLY_NO_MSPRIME_OR_IQTREE','fitted_sha256':sha(args.fitted),'taxa':12,'genes':173,'columns':68796,'zero_edges_collapsed':checked['collapsed_zero_edges'],'seed_count':len(seeds)*3,'seed_collisions':0,'msprime_installed_or_API_run':'NOT VERIFIED','generation_plan':'One independent full12 Hudson genealogy per gene; no migration/recombination; scale each CU edge by kappa*r_gene; original AA model and full12 cell mask; emit original inference command jobs only.','formal_P_value':None}
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
