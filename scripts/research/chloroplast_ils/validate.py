#!/usr/bin/env python3
"""Independent process: SciPy analytic integration; Bio.Phylo source tree verification."""
import argparse,csv,gzip,json,math,hashlib,time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.stats import binom
from scipy.optimize import brentq
from Bio import Phylo
P=argparse.ArgumentParser();P.add_argument('--root',required=True);P.add_argument('--phase',required=True);a=P.parse_args();R=Path(a.root);O=R/('results_'+a.phase)
def read(n):return list(csv.DictReader((O/n).open()))
obs=json.loads((R/'inputs/observed.json').read_text());guide=json.loads((R/'inputs/guide.json').read_text());checks=[]
def check(name,ok,detail=None):checks.append({'name':name,'pass':bool(ok),'detail':detail})
L=obs['joint_acgt_length'];K=obs['differences']
expected_n=20000 if a.phase=='full' else 500
expected_d=100000 if a.phase=='full' else 2000
expected_s=10 if a.phase=='full' else 1
counts=read('topology_counts.csv');draws_all=read('topology_draws.csv')
check('topology_row_count',len(counts)==9 and len({(x['multiplier'],x['topology']) for x in counts})==9)
check('topology_draw_count',len(draws_all)==expected_n*3)
for m in [1,2,4]:
 rows=[x for x in counts if int(x['multiplier'])==m]
 check(f'fixed_replicates_m{m}',sum(int(x['count']) for x in rows)==expected_n and all(int(x['n'])==expected_n for x in rows))
 for row in rows:
  check(f'draw_summary_parity_m{m}_{row["topology"]}',sum(int(x['multiplier'])==m and x['topology']==row['topology'] for x in draws_all)==int(row['count']))
grid=read('distance_grid.csv')
check('distance_grid_count',len(grid)==80 and len({(x['a'],x['b']) for x in grid})==80)
check('fixed_distance_replicates',all(int(x['draws'])==expected_d for x in grid))
f=(R/'inputs/hime2_hachi_joint_acgt.fa').read_text().splitlines();check('observed_pair_recount',len(f[1])==len(f[3])==L and sum(x!=y for x,y in zip(f[1],f[3]))==K)
for name in ['cp98_tree','cp97_tree']:
 tr=Phylo.read(str(R/'inputs'/(name+'.nwk')),'newick')
 for i,ob in enumerate(o for o in obs['quartets'] if o['tree']==name):
  v=ob['labels'];s={x:tr.distance(v['H'],v[x[1]])+tr.distance(v[o[0]],v[o[1]]) for x,o in [('HT','GZ'),('HG','TZ'),('HZ','TG')]}
  check(f'{name}_independent_projection_{i}',min(s,key=s.get)==ob['unrooted_quartet'])
for row in read('topology_counts.csv'):
 m=int(row['multiplier']);p=1-2/3*math.exp(-m*guide['unrooted_internal_length_nuclear_CU']) if row['topology']=='HT' else math.exp(-m*guide['unrooted_internal_length_nuclear_CU'])/3
 n=int(row['n']);k=int(row['count']);check(f'topology_m{m}_{row["topology"]}_analytic',abs(k/n-p)<=6*math.sqrt(p*(1-p)/n)+6/n,{'expected':p,'observed':k/n})
# Recount genealogy files with Bio.Phylo, including projection of a deterministic subset.
for m in guide['cp_rate_multipliers']:
 path=O/f'genealogies_m{m}.nwk.gz'
 with gzip.open(path,'rt') as f:n=sum(1 for _ in f)
 check(f'genealogy_file_count_m{m}',n==int(next(x['n'] for x in read('topology_counts.csv') if int(x['multiplier'])==m)))
 from io import StringIO
 draws=[x for x in read('topology_draws.csv') if int(x['multiplier'])==m]
 with gzip.open(path,'rt') as f:
  for i,line in enumerate(f):
   if i>=50:break
   tr=Phylo.read(StringIO(line),'newick');s={x:tr.distance('H',x[1])+tr.distance(o[0],o[1]) for x,o in [('HT','GZ'),('HG','TZ'),('HZ','TG')]}
   check(f'genealogy_independent_projection_m{m}_{i}',min(s,key=s.get)==draws[i]['topology'])
def pjc(d):return -.75*math.expm1(-4*d/3)
def tail(aa,bb):
 if bb==0:return float(binom.cdf(K,L,pjc(aa)))
 fn=lambda u:math.exp(-u)*float(binom.cdf(K,L,pjc(aa+bb*u)))
 cut=min(50.,50./(L*bb))
 v,e=quad(fn,0,cut,epsabs=1e-12,epsrel=1e-8,limit=200)
 v2,e2=quad(fn,cut,np.inf,epsabs=1e-12,epsrel=1e-8,limit=200)
 assert e+e2<1e-8
 return v+v2
analytic=[]
for row in read('distance_grid.csv'):
 aa=float(row['a']);bb=float(row['b']);q=tail(aa,bb);n=int(row['draws']);k=int(row['draws_K_le_observed'])
 check(f'distance_analytic_{aa}_{bb}',abs(k/n-q)<=6*math.sqrt(q*(1-q)/n)+6/n,{'analytic':q,'observed':k/n})
 analytic.append({**row,'conditional_tail_probability_analytic':q})
with (O/'distance_grid_validated.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(analytic[0]));w.writeheader();w.writerows(analytic)
full=read('full_sequence_checks.csv')
check('fixed_full_sequence_pair_count',len(full)==80*expected_s)
res=0.;var=0.
for aa,bb in sorted({(x['a'],x['b']) for x in full}):
 rows=[x for x in full if x['a']==aa and x['b']==bb];rr=0.;vv=0.
 check(f'full_sequence_cell_count_{aa}_{bb}',len(rows)==expected_s)
 for row in rows:
  pp=pjc(float(row['realized_D']));ex=L*pp;va=L*pp*(1-pp)
  check(f'full_sequence_expected_recomputed_{aa}_{bb}_{row["replicate"]}',abs(ex-float(row['conditional_expected_K']))<1e-7 and abs(va-float(row['conditional_variance_K']))<1e-7)
  rr+=int(row['K'])-ex;vv+=va
 check(f'full_sequence_cell_distribution_{aa}_{bb}',abs(rr)<=6*math.sqrt(vv)+6)
 res+=rr;var+=vv
check('actual_full_sequences_vs_JC69_distance_distribution',abs(res)<=6*math.sqrt(var)+6,{'residual':res,'variance':var,'z':res/math.sqrt(var)})
# Every saved full pair is independently read and matches the logged site distance.
with gzip.open(O/'example_simulated_pairs.fa.gz','rt') as f:lines=[x.strip() for x in f]
saved=[x for x in full if x['saved_fasta']=='True'];check('saved_fasta_pair_count',len(lines)==4*len(saved))
for i,row in enumerate(saved):
 s,t=lines[4*i+1],lines[4*i+3];check(f'saved_pair_{i}',len(s)==len(t)==L and sum(x!=y for x,y in zip(s,t))==int(row['K']))
threshold=brentq(lambda v:float(binom.cdf(K,L,pjc(v)))-.05,0,.01,xtol=1e-15)
result={'status':'VALIDATED_CONDITIONAL_SIMULATIONS' if all(c['pass'] for c in checks) else 'VALIDATION_FAILED','checks':checks,'pass_count':sum(c['pass'] for c in checks),'fail_count':sum(not c['pass'] for c in checks),'distance_upper_bound_b0_a_at_tail_0_05':threshold,'threshold_is_model_sensitivity_not_species_divergence_estimate':True,'empirical_ILS_P':None,'JML_reproduction':False,'independent_engines':'msprime/numpy simulation; scipy integration and Bio.Phylo validation'}
(O/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='checks'}),flush=True)
if result['fail_count']:raise SystemExit(1)
