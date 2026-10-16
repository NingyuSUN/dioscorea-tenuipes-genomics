#!/usr/bin/env python3
"""Conditional plastid genealogies and JC69 distance sensitivity, NOT an empirical ILS P test."""
import argparse,csv,gzip,hashlib,json,math,os,platform,time
from pathlib import Path
import numpy as np
import msprime
P=argparse.ArgumentParser();P.add_argument('--root',required=True);P.add_argument('--phase',choices=['pilot','full'],required=True);args=P.parse_args()
R=Path(args.root);O=R/('results_'+args.phase);O.mkdir(exist_ok=False)
start=time.monotonic();guide=json.loads((R/'inputs/guide.json').read_text());obs=json.loads((R/'inputs/observed.json').read_text())
seed=2026100301+(0 if args.phase=='full' else 100000)
ntrees=20000 if args.phase=='full' else 500
ndist=100000 if args.phase=='full' else 2000
nseq=10 if args.phase=='full' else 1
L=obs['joint_acgt_length'];Kobs=obs['differences'];times=guide['split_times_nuclear_CU']
def save(n,d): (O/n).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def table(n,rows):
 with (O/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def wilson(k,n):
 z=1.959963984540054;q=k/n;den=1+z*z/n;c=(q+z*z/(2*n))/den
 h=z*math.sqrt(q*(1-q)/n+z*z/(4*n*n))/den
 return max(0,c-h),min(1,c+h)
def jc(d):return -.75*np.expm1(-4*np.asarray(d)/3)
# One nonrecombining haploid sample per population, with mean pairwise coalescence 1/m nuclear CU.
counts=[];topology_draws=[]
for i,m in enumerate(guide['cp_rate_multipliers']):
 d=msprime.Demography()
 for n in ['H','T','G','Z','HT','GZ','ROOT']:d.add_population(name=n,initial_size=1/m)
 d.add_population_split(time=times['HT'],derived=['H','T'],ancestral='HT')
 d.add_population_split(time=times['GZ'],derived=['G','Z'],ancestral='GZ')
 d.add_population_split(time=times['ROOT'],derived=['HT','GZ'],ancestral='ROOT');d.sort_events();d.validate()
 c={s:0 for s in ['HT','HG','HZ']};label={0:'H',1:'T',2:'G',3:'Z'}
 with gzip.open(O/f'genealogies_m{m}.nwk.gz','wt') as nw:
  for j,ts in enumerate(msprime.sim_ancestry(samples={s:1 for s in ['H','T','G','Z']},demography=d,ploidy=1,sequence_length=1,recombination_rate=0,num_replicates=ntrees,random_seed=seed+i)):
   tr=ts.first();assert ts.num_trees==1 and ts.num_samples==4
   assert [ts.population(ts.node(x).population).metadata['name'] for x in ts.samples()]==['H','T','G','Z']
   s={'HT':2*tr.tmrca(0,1)+2*tr.tmrca(2,3),'HG':2*tr.tmrca(0,2)+2*tr.tmrca(1,3),'HZ':2*tr.tmrca(0,3)+2*tr.tmrca(1,2)}
   order=sorted(s,key=s.get);assert s[order[1]]>s[order[0]]+1e-12
   c[order[0]]+=1
   topology_draws.append({'multiplier':m,'replicate':j+1,'topology':order[0]})
   nw.write(tr.as_newick(node_labels=label,precision=12)+'\n')
 x=guide['unrooted_internal_length_nuclear_CU']*m
 for top,k in c.items():
  p=1-2/3*math.exp(-x) if top=='HT' else math.exp(-x)/3
  lo,hi=wilson(k,ntrees);z=(k-ntrees*p)/math.sqrt(ntrees*p*(1-p))
  counts.append({'multiplier':m,'topology':top,'count':k,'n':ntrees,'frequency':k/ntrees,'wilson95_lower':lo,'wilson95_upper':hi,'analytic_probability':p,'analytic_z':z})
 print(json.dumps({'stage':'genealogies_complete','m':m,'counts':c,'elapsed_s':time.monotonic()-start}),flush=True)
table('topology_counts.csv',counts);table('topology_draws.csv',topology_draws)
comparisons=[]
for ob in obs['quartets']:
 for c in counts:
  if c['topology']==ob['unrooted_quartet']:comparisons.append({'source_tree':ob['tree'],**ob['labels'],**c,'empirical_ILS_P':None})
table('source_quartet_comparisons.csv',comparisons)
# a=2*mu*tau; b=2*mu*c. These are hypothetical DNA mutation scales, never fitted to the target pair.
a_grid=[0.,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3]
b_grid=[0.,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2]
rows=[];fullrows=[];alphabet=np.frombuffer(b'ACGT',dtype='S1');rng=np.random.default_rng(seed+100)
with gzip.open(O/'example_simulated_pairs.fa.gz','wt') as fasta:
 for ai,a in enumerate(a_grid):
  for bi,b in enumerate(b_grid):
   D=a+(rng.exponential(b,ndist) if b else np.zeros(ndist));p=jc(D)
   K=rng.binomial(L,p);k=int(np.count_nonzero(K<=Kobs));lo,hi=wilson(k,ndist)
   rows.append({'a':a,'b':b,'joint_sites':L,'threshold_observed_K':Kobs,'draws':ndist,'draws_K_le_observed':k,'conditional_tail_frequency':k/ndist,'wilson95_lower':lo,'wilson95_upper':hi,'mean_K':float(K.mean()),'q025_K':float(np.quantile(K,.025)),'median_K':float(np.median(K)),'q975_K':float(np.quantile(K,.975)),'empirical_ILS_P':None})
   for r in range(nseq):
    dd=a+(rng.exponential(b) if b else 0.);pr=float(jc(dd/2))
    root=rng.integers(0,4,L,dtype=np.int8)
    # Two independent JC69 branches D/2 from one ancestral sequence.
    pair=[]
    for branch in range(2):
     change=rng.random(L)<pr;offset=rng.integers(1,4,L,dtype=np.int8)
     pair.append((root+change*offset)%4)
    kk=int(np.count_nonzero(pair[0]!=pair[1]));pp=float(jc(dd))
    fullrows.append({'a_index':ai,'b_index':bi,'replicate':r+1,'a':a,'b':b,'realized_D':dd,'L':L,'K':kk,'conditional_expected_K':L*pp,'conditional_variance_K':L*pp*(1-pp),'saved_fasta':r==0})
    if r==0:
     for j,seq in enumerate(pair):
      fasta.write(f'>a{ai}_b{bi}_rep1_tip{j+1}\n'+alphabet[seq].tobytes().decode('ascii')+'\n')
table('distance_grid.csv',rows);table('full_sequence_checks.csv',fullrows)
save('RUN_RECEIPT.json',{'status':'SIMULATION_COMPLETE_VALIDATION_PENDING','phase':args.phase,'seed':seed,'independent_genealogies_per_scenario':ntrees,'genealogy_scenarios':len(guide['cp_rate_multipliers']),'distance_grid_cells':len(rows),'distance_draws_per_cell':ndist,'actual_full_sequence_pairs':len(fullrows),'saved_full_sequence_pairs':len(rows),'one_plastome_one_linked_genealogy':True,'same_retained_site_count_as_observed':True,'original_mask_coordinates_saved_but_not_reapplied_to_simulated_compressed_sequences':True,'sequence_model':'JC69 homogeneous independent sites conditional on one shared genealogy per plastome pair','posterior_predictive_JML':False,'empirical_ILS_P':None,'parameters_estimated_from_target_pair':False,'elapsed_s':time.monotonic()-start,'python':platform.python_version(),'numpy':np.__version__,'msprime':msprime.__version__,'cpu_affinity':sorted(os.sched_getaffinity(0)),'inputs':{str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (R/'inputs').iterdir() if p.is_file()}})
print(json.dumps({'stage':'simulation_complete','phase':args.phase,'elapsed_s':time.monotonic()-start}),flush=True)
