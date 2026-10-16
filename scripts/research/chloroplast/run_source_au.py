#!/usr/bin/env python3
"""Sequential paired-source searches -> frozen topology union -> paired AU.

No inference at import or by default. Existing directories are never resumed.
"""
from pathlib import Path
import argparse,json,os,shutil,signal,subprocess,time,traceback
from source_au_common import *

def verify_inputs(p):
    m=json.loads((p/'PREPARED.json').read_text());c=json.loads((p/'CONFIG.json').read_text())
    for name,r in m['files'].items():assert sha(p/name)==r['sha256'],name
    assert sha(c['iqtree'])==c['iqtree_sha256'],'IQ-TREE binary changed'
    assert m['columns']==122234 and m['taxa_per_variant']==10
    assert c['maximum_rescues_per_variant']==2 and c['rescue_seeds']==[20261010,20261011]
    assert c['candidate_max_before_dedup']==18 and c['nuclear_unchanged']
    return m,c

def base_args(variant,inputs,prefix,seed,c):
    # Same model classes and numerical search flags as old shared10 full01.
    return [c['iqtree'],'-s',str(inputs/(variant+'.fa')),'-p',str(inputs/'joint.partitions.nex'),'-T',str(c['threads']),'-seed',str(seed),'-keep-ident','--prefix',str(prefix)]

def stop(child):
    if child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
        try:child.wait(timeout=10)
        except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()

def run_command(cmd,prefix,timeout=None):
    prefix.parent.mkdir(parents=True,exist_ok=True)
    assert not list(prefix.parent.glob(prefix.name+'.*')),'Existing prefix; never retry or overwrite'
    r={'status':'RUNNING','argv':cmd,'started_epoch':time.time(),'timeout_seconds':timeout};save(prefix.with_suffix('.runtime.json'),r)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='4',MKL_NUM_THREADS='1');child=None
    try:
        with prefix.with_suffix('.stdout').open('x') as out,prefix.with_suffix('.stderr').open('x') as err:
            child=subprocess.Popen(cmd,stdout=out,stderr=err,env=env,start_new_session=True);r['pid']=child.pid;save(prefix.with_suffix('.runtime.json'),r)
            try:code=child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:stop(child);code=child.returncode;r['status']='PILOT_TIME_LIMIT_NO_INFERENCE_CLAIM'
        r['exit_code']=code
        if r['status']=='RUNNING':
            if code!=0:raise RuntimeError('Command exit '+str(code))
            r['status']='COMPLETE'
    except BaseException:
        if child is not None:stop(child)
        r.update(status='FAILED_NO_AUTO_RESTART',error=traceback.format_exc());raise
    finally:
        r.update(finished_epoch=time.time(),elapsed_seconds=time.time()-r['started_epoch'])
        r['files']={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(prefix.parent.glob(prefix.name+'.*')) if p.is_file() and not p.name.endswith('.runtime.json')};save(prefix.with_suffix('.runtime.json'),r)
    return r

def search_result(prefix,variant,hypothesis,seed,c):
    text=prefix.with_suffix('.iqtree').read_text();match=re.search(r'Log-likelihood of the tree:\s*([-+0-9.eE]+)',text);assert match
    newick=prefix.with_suffix('.treefile').read_text().strip();tips,_=topology(newick);assert tips==frozenset(c['tips'])
    if hypothesis in c['constraints']:assert has_split(newick,c['constraints'][hypothesis],c['tips'])
    return {'variant':variant,'hypothesis':hypothesis,'seed':seed,'prefix':str(prefix),'logL_search':float(match[1]),'newick':newick}

def rescue_needed(rows,c):
    best=max(rows,key=lambda x:x['logL_search']);best_u=max((x for x in rows if x['hypothesis']=='unconstrained'),key=lambda x:x['logL_search'])
    return best['logL_search']>best_u['logL_search']+c['rescue_trigger_logL_margin'],best,best_u

def audit_au(prefix,entries,columns):
    table=user_tree_table(prefix.with_suffix('.iqtree'));ns,rows,bounds=sitelh(prefix.with_suffix('.sitelh'))
    assert ns==columns and len(rows)==len(entries) and set(table)==set(range(1,len(entries)+1))
    checks=[]
    for i,((label,values),bound) in enumerate(zip(rows,bounds),1):
        total=math.fsum(values);report=table[i]['logL'];allow=bound+float(sig_radius(report,10))+1e-7
        assert abs(total-report)<=allow,'Site sum inconsistent with six-significant-digit output bound'
        checks.append({'tree_index_1based':i,'topology_sha256':entries[i-1]['topology_sha256'],'label':label,'site_sum':total,'reported_logL':report,'site_rounding_radius':bound,'report_rounding_radius':float(sig_radius(report,10)),'floating_sum_slack':1e-7,'absolute_difference':abs(total-report)})
    # IQ-TREE writes evaluated trees in supplied candidate order.
    lines=[line.strip() for line in prefix.with_suffix('.trees').read_text().splitlines() if line.strip()]
    assert len(lines)==len(entries)
    for line,e in zip(lines,entries):assert topology_id(line)==e['topology_sha256'],'Evaluated candidate order/topology changed'
    return {'status':'PASS_HASH_ORDER_TOPOLOGY_AND_ROUNDING','rows':checks,'table':table}

def execute(inputs,out,mode,m,c):
    cg=cgroup_check(c);assert not out.exists(),'Existing output; do not auto-resume'
    assert shutil.disk_usage(inputs).free>=c['minimum_free_disk_bytes']
    out.mkdir(parents=True,exist_ok=False)
    state={'status':'RUNNING','mode':mode,'started_epoch':time.time(),'input_receipt_sha256':sha(inputs/'PREPARED.json'),'driver_sha256':sha(__file__),'cgroup':cg,'jobs':[],'rescue_decisions':[]};save(out/'RUN_STATUS.json',state)
    try:
        if mode=='pilot':
            for variant in c['variants']:
                prefix=out/variant/'pilot';cmd=base_args(variant,inputs,prefix,c['search_seeds'][0],c)+['-g',str(inputs/'HT_cherry.constraint.nwk'),'-n',str(c['pilot_iterations'])]
                state['jobs'].append(run_command(cmd,prefix,c['pilot_timeout_per_variant']));save(out/'RUN_STATUS.json',state)
            state['status']='PILOT_TIMING_ONLY_NO_AUTO_FULL'
        else:
            all_records=[];best_unconstrained={}
            for variant in c['variants']:
                rows=[]
                for hypothesis in ['unconstrained',*c['constraints']]:
                    seeds=c['search_seeds'] if hypothesis=='unconstrained' else c['constrained_search_seeds']
                    for seed in seeds:
                        prefix=out/variant/(hypothesis+'_'+str(seed));cmd=base_args(variant,inputs,prefix,seed,c)
                        if hypothesis!='unconstrained':cmd+=['-g',str(inputs/(hypothesis+'.constraint.nwk'))]
                        rec=run_command(cmd,prefix);row=search_result(prefix,variant,hypothesis,seed,c);rows.append(row);state['jobs'].append(rec);save(out/variant/'SEARCH_RESULTS.json',rows);save(out/'RUN_STATUS.json',state)
                for k,seed in enumerate(c['rescue_seeds']):
                    needed,best,best_u=rescue_needed(rows,c)
                    decision={'variant':variant,'rescue_number':k+1,'seed':seed,'triggered':needed,'best_logL':best['logL_search'],'best_unconstrained_logL':best_u['logL_search'],'margin':c['rescue_trigger_logL_margin'],'criterion':'search likelihood only; no AU available yet'}
                    state['rescue_decisions'].append(decision);save(out/'RUN_STATUS.json',state)
                    if not needed:break
                    prefix=out/variant/('unconstrained_rescue_'+str(k+1));cmd=base_args(variant,inputs,prefix,seed,c)+['-t',best['prefix']+'.treefile']
                    rec=run_command(cmd,prefix);rows.append(search_result(prefix,variant,'unconstrained',seed,c));state['jobs'].append(rec);save(out/variant/'SEARCH_RESULTS.json',rows);save(out/'RUN_STATUS.json',state)
                needed,best,best_u=rescue_needed(rows,c)
                assert not needed,'SEARCH_INADEQUATE_AFTER_BOUNDED_RESCUE; stop interpretation'
                best_unconstrained[variant]=best_u;all_records.extend(rows)
            old_manifest=json.loads((inputs/'historical.CANDIDATE_MANIFEST.json').read_text())
            for e in old_manifest['candidates']:
                all_records.append({'variant':'historical_OQ_124448','hypothesis':'historical_candidate','old_tree_index_1based':e['tree_index_1based'],'newick':e['newick']})
            entries=freeze_union(all_records,c['tips'],c['constraints'],c['candidate_max_before_dedup'])
            candidate_file=out/'UNION.candidates.nwk';candidate_file.write_text('\n'.join(e['newick'] for e in entries)+'\n')
            manifest={'status':'FROZEN_BEFORE_ANY_SOURCE_AU','candidates':entries,'candidate_count':len(entries),'pre_dedup_count':len(all_records),'candidate_file_sha256':sha(candidate_file),'best_unconstrained':best_unconstrained,'no_selection_using_AU':True}
            save(out/'UNION_CANDIDATE_MANIFEST.json',manifest);state.update(status='SEARCH_COMPLETE_CANDIDATE_UNION_FROZEN',union_candidate_count=len(entries),candidate_file_sha256=sha(candidate_file));save(out/'RUN_STATUS.json',state)
            bindings={}
            for variant in c['variants']:
                assert sha(candidate_file)==manifest['candidate_file_sha256']
                prefix=out/variant/'au';cmd=base_args(variant,inputs,prefix,c['AU_seed'],c)+['-te',best_unconstrained[variant]['prefix']+'.treefile','-z',str(candidate_file),'-zb',str(c['AU_RELL_replicates']),'-au','-wsl']
                rec=run_command(cmd,prefix);state['jobs'].append(rec);save(out/'RUN_STATUS.json',state)
                audit=audit_au(prefix,entries,m['columns']);save(out/variant/'AU_NUMERIC_ORDER_AUDIT.json',audit)
                binding={'variant':variant,'union_candidate_file_sha256':sha(candidate_file),'union_manifest_sha256':sha(out/'UNION_CANDIDATE_MANIFEST.json'),'sitelh_sha256':sha(prefix.with_suffix('.sitelh')),'iqtree_sha256':sha(prefix.with_suffix('.iqtree')),'candidate_count':len(entries),'columns':m['columns'],'audit_sha256':sha(out/variant/'AU_NUMERIC_ORDER_AUDIT.json'),'model_policy':'Each variant estimates model parameters on its best unconstrained topology, shares those estimates across candidate evaluations, and optimizes candidate branch lengths.'}
                save(out/variant/'EVALUATION_BINDING.json',binding);bindings[variant]=binding
            assert len({b['union_candidate_file_sha256'] for b in bindings.values()})==1
            state['status']='PAIRED_UNION_AU_COMPLETE_SUMMARY_AND_SCIENTIFIC_REVIEW_PENDING'
    except BaseException:state.update(status='FAILED_NO_AUTO_RESTART',error=traceback.format_exc());raise
    finally:
        state.update(finished_epoch=time.time(),elapsed_seconds=time.time()-state['started_epoch']);cgroot=Path(cg['path'])
        cg['memory_events_after']=(cgroot/'memory.events').read_text()
        if (cgroot/'memory.peak').exists():cg['memory_peak_bytes']=(cgroot/'memory.peak').read_text().strip()
        save(out/'RUN_STATUS.json',state)
    return state

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--mode',choices=['pilot','full'],default='full');ap.add_argument('--execute',action='store_true');a=ap.parse_args();a.inputs=a.inputs.resolve();a.out=a.out.resolve();m,c=verify_inputs(a.inputs)
    if not a.execute:
        print(json.dumps({'status':'DRY_RUN_NO_INFERENCE','mode':a.mode,'variants':c['variants'],'columns':m['columns'],'taxa':10,'regular_search_count':12,'maximum_total_rescues':4,'maximum_candidates_before_dedup':18,'AU_jobs':2,'threads':4,'memory_bytes':c['cgroup_memory_max_bytes'],'sequential':True,'no_automatic_follow_on':True},indent=2));return
    def interrupted(signum,frame):raise InterruptedError('Signal '+str(signum))
    signal.signal(signal.SIGTERM,interrupted)
    print(json.dumps(execute(a.inputs,a.out,a.mode,m,c),indent=2))
if __name__=='__main__':main()
