#!/usr/bin/env python3
"""Prepare the frozen chloroplast ILS inputs; default is read-only preflight."""
from pathlib import Path
import argparse, hashlib, json, re, shutil
try:
    from .portable_inputs import Config, report, require_preflight
except ImportError:
    from portable_inputs import Config, report, require_preflight
expected={'source_alignment':'72a7839d0a41d124a4d02e119941b0f7464c92f43ce205914a257bd9a7dd8852','cp98_tree':'4b0352705d1dfdecf5e166995658610b7d3924344465cd4247a14374c78c25dc','cp97_tree':'8504215ab28c68ad6e0ce58e51d64c9f52107e556571e8f7518a795ae6a078dd'}

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def plan(config_path):
    c=Config(config_path);out=c.path('output_dir')/'inputs'
    roles=['source_alignment','cp98_tree','cp97_tree','nuclear_S0']
    return {'script':'prepare_inputs','inputs':[c.input(k,expected.get(k)) for k in roles],
            'outputs':[str(out/n) for n in ['hime2_hachi_joint_acgt.fa','joint_acgt_columns_1based.txt','cp98_tree.nwk','cp97_tree.nwk','guide.json','observed.json','SOURCE_BINDINGS.json']],
            'required_threads':1,'tools':{},'python_packages':[], 'scientific_contract':'Frozen 98 sequences x 124448 columns; hash-bound tree sources; no simulation'}

def execute(config_path):
    c=Config(config_path);spec=plan(config_path);require_preflight(spec)
    R=c.path('output_dir');destination=R/'inputs'
    if destination.exists():raise ValueError('Input preparation requires a fresh output inputs directory')
    paths={row['role']:Path(row['path']) for row in spec['inputs']}
    expected={row['role']:row['sha256'] for row in spec['inputs'] if row['sha256']}
    destination.mkdir(parents=True)
    def save(name,obj): (destination/name).write_text(json.dumps(obj,indent=2)+'\n')
    for k,v in expected.items():assert sha(paths[k])==v,(k,sha(paths[k]))
    seq={}
    for line in paths['source_alignment'].read_text().splitlines():
     if line.startswith('>'):name=line[1:].split()[0];assert name not in seq;seq[name]=''
     else:seq[name]+=line.strip().upper()
    assert len(seq)==98 and set(map(len,seq.values()))=={124448}
    labels={'H':['DT_local_path2','PV612491.1','PZ211775.1'],'T':['OQ525999.1','Dtokoro_Kita1_2022_author_declared_public_CP'],'G':['OQ525995.1'],'Z':['KP899622.1','OQ526000.1']}
    a,b=seq['DT_local_path2'],seq['Edo_CP_read_validated_partial_coordinate_projection']
    mask=[i for i,(x,y) in enumerate(zip(a,b)) if x in 'ACGT' and y in 'ACGT']
    diffs=[i+1 for i in mask if a[i]!=b[i]]
    (R/'inputs'/'hime2_hachi_joint_acgt.fa').write_text('>hime2\n'+''.join(a[i] for i in mask)+'\n>Hachi-1\n'+''.join(b[i] for i in mask)+'\n')
    (R/'inputs'/'joint_acgt_columns_1based.txt').write_text('\n'.join(str(i+1) for i in mask)+'\n')
    # Minimal Newick parser: immutable source trees; internal labels ignored, lengths retained.
    def treeparse(s):
     tok=re.findall(r"\(|\)|,|:|;|[^\s(),:;]+",s); idx=0; nodes=[]
     def node(parent=None):
      nonlocal idx
      n=len(nodes);nodes.append({'parent':parent,'name':None,'length':0.,'children':[]})
      if tok[idx]=='(':
       idx+=1
       while True:
        nodes[n]['children'].append(node(n))
        if tok[idx]!=',':break
        idx+=1
       assert tok[idx]==')';idx+=1
       if tok[idx] not in [':',',',')',';']:idx+=1
      else:nodes[n]['name']=tok[idx];idx+=1
      if tok[idx]==':':idx+=1;nodes[n]['length']=float(tok[idx]);idx+=1
      return n
     node();assert tok[idx]==';';return nodes
    observations=[]
    for tree_id in ['cp98_tree','cp97_tree']:
     nodes=treeparse(paths[tree_id].read_text());leaf={n['name']:i for i,n in enumerate(nodes) if n['name']}
     def distance(a,b):
      anc={};d=0.;u=leaf[a]
      while u is not None:anc[u]=d;d+=nodes[u]['length'];u=nodes[u]['parent']
      d=0.;u=leaf[b]
      while u not in anc:d+=nodes[u]['length'];u=nodes[u]['parent']
      return d+anc[u]
     for h in labels['H']:
      for t in labels['T']:
       for z in labels['Z']:
        v={'H':h,'T':t,'G':labels['G'][0],'Z':z}
        sums={x:distance(v['H'],v[x[1]])+distance(v[other[0]],v[other[1]]) for x,other in [('HT','GZ'),('HG','TZ'),('HZ','TG')]}
        order=sorted(sums,key=sums.get);gap=sums[order[1]]-sums[order[0]]
        assert gap>1e-12,(v,sums)
        observations.append({'tree':tree_id,'labels':v,'unrooted_quartet':order[0],'pair_sums_substitutions_per_site':sums,'split_gap':gap,'same_individual_binding_verified':False,'source_scenario_not_independent_replicate':True})
     shutil.copyfile(paths[tree_id],R/'inputs'/(tree_id+'.nwk'))
    s0=json.loads(paths['nuclear_S0'].read_text());allnodes=s0['rooted_population_nodes'];idnode={n['id']:n for n in allnodes}
    ntaxa={'H':'Dioscorea_tenuipes_historical','T':'Dioscorea_tokoro','G':'Dioscorea_gracillima','Z':'Dioscorea_zingiberensis'}
    ntimes={k:idnode[n]['time_CU'] for k,n in [('HT','pop_003'),('GZ','pop_004'),('ROOT','pop_002')]}
    x=2*ntimes['ROOT']-ntimes['HT']-ntimes['GZ']
    guide={'nuclear_taxa':ntaxa,'split_times_nuclear_CU':ntimes,'unrooted_internal_length_nuclear_CU':x,'coalescent_mean_nuclear_CU':1,'cp_rate_multipliers':[1,2,4],'multiplier_one_role':'nuclear_like_sensitivity_control','primary_source_Z':'KP899622.1','physical_time_Ne_DNA_mutation_rate':'NOT_ESTIMATED','aa_kappa_not_used':s0['kappa'],'fixed_nuclear_topology':'((H,T),(G,Z))','limitations':s0['uncalibrated_limitations']}
    save('guide.json',guide)
    save('observed.json',{'source_alignment_length':124448,'joint_acgt_length':len(mask),'differences':len(diffs),'difference_source_columns_1based':diffs,'p_distance':len(diffs)/len(mask),'excluded_unknown_or_gap':124448-len(mask),'IR_policy':'source frozen alignment retains one IR; LSC80600 IR25351 SSC18497','distance_pair':['hime2','Hachi-1'],'quartets':observations})
    save('SOURCE_BINDINGS.json',{'sources':{k:{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for k,p in paths.items()},'assumptions':['H historical nuclear to DT_local_path2 CP is a conditional label mapping; sequence binding not specimen authentication.','Public PV/PZ and both tokoro/zingiberensis sequences are source substitutions, not matched nuclear individuals.','Hachi and panthaica are not nuclear guide tips.'],'derived':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in (R/'inputs').iterdir() if p.is_file() and p.name!='SOURCE_BINDINGS.json'}})
    print(json.dumps({'L':len(mask),'K':len(diffs),'quartet_scenarios':len(observations),'quartet_summary':[(o['tree'],o['labels']['H'],o['labels']['T'],o['labels']['Z'],o['unrooted_quartet']) for o in observations],'x':x}))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    if args.execute:execute(args.config);return 0
    return report(plan(args.config))

if __name__=='__main__':raise SystemExit(main())
