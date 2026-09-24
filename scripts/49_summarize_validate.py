#!/usr/bin/env python3
"""Validate archived analysis results and generate concise tables/figures."""
import argparse,json,hashlib,itertools,sys
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from relationship_core import read_bed,ibs_distance,checked_neighbor_net
p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();r=a.root;res=r/'results';checks={}
d=pd.read_csv(res/'panel_stability.tsv',sep='\t');half=pd.read_csv(res/'half_contig_stability.tsv',sep='\t')
assert len(d)==47 and len(half)==100
checks['panels']=len(d);checks['half_contig_replicates']=len(half)
checks['geographic_NJ_groups_preserved_all_panels']=bool(d.filter(regex='_own_split_NJ$').all().all() and half.filter(regex='_own_split_NJ$').all().all())
checks['all_NNet_group_splits_positive']=bool((d.filter(regex='_own_split_NNet_weight$')>0).all().all() and (half.filter(regex='_own_split_NNet_weight$')>0).all().all())
ind=pd.read_csv(r/'independent_R/cross_implementation.tsv',sep='\t');assert len(ind)==147 and ind.RF_ape_Biopython.eq(0).all();checks['NJ_cross_implementation_trees']=len(ind);checks['NJ_max_patristic_error']=float(ind.max_patristic_difference.max())
checks['network_cross_implementation']=json.loads((r/'network_cross_validation.json').read_text());assert checks['network_cross_implementation']['maximum_weight_difference']<1e-6
errs=[]
for folder in sorted((r/'data').iterdir()):
 if not (folder/'pruned.bed').exists():continue
 names=[x.split()[1] for x in (folder/'pruned.fam').read_text().splitlines()];bim=pd.read_csv(folder/'pruned.bim',sep=r'\s+',header=None,dtype=str);assert bim[1].is_unique and not (bim[1]=='.').any();g=read_bed(folder/'pruned.bed',len(names),len(bim));dist,joint=ibs_distance(g)
 flat=folder/'flat.mdist' if folder.name!='supplement17' else r/'data/primary17_flat.mdist'
 if flat.exists():
  error=float(abs(dist-np.loadtxt(flat)).max());assert error<5.1e-7;errs.append({'panel':folder.name,'maximum_difference':error})
 if (folder/'panel.json').exists():
  meta=json.loads((folder/'panel.json').read_text());assert len(g)==meta['pruned'];called=(g>=0).sum(axis=1);ac=np.where(g<0,0,g).sum(axis=1);assert (np.minimum(ac,2*called-ac)>=2).all();assert (2*called>=meta['minimum_AN']).all()
 if folder.name in ['supplement17','strict17','complete17']:
  direct=pd.read_csv(r/'relatedness'/folder.name/'first50_genotypes.tsv',sep='\t').fillna(-1).to_numpy();assert np.array_equal(g[:50],direct)
  king=pd.read_csv(r/'relatedness'/folder.name/'king_matrix.tsv',sep='\t',index_col=0).loc[names,names].to_numpy()
  for i,j in itertools.combinations(range(len(names)),2):
   ok=(g[:,i]>=0)&(g[:,j]>=0);x,y=g[ok,i],g[ok,j];v=.5-np.sum((x-y)**2)/(4*min((x==1).sum(),(y==1).sum()));assert abs(v-king[i,j])<1e-12
checks['PLINK_distance_cross_checks']=errs
for f in res.glob('*/network_solver.json'):assert json.loads(f.read_text())['kkt_violation']<1e-9
checks['NNLS_optimality_checks']=len(list(res.glob('*/network_solver.json')));assert checks['NNLS_optimality_checks']==147
m=pd.read_csv(res/'supplement17/distance.tsv',sep='\t',index_col=0);_,s,_=checked_neighbor_net(list(m.index),m.values);_,s2,_=checked_neighbor_net(list(m.index),m.values);assert [(x.part1(),x.weight) for x in s]==[(x.part1(),x.weight) for x in s2]
checks['NNLS_repeatability']=True
# Collect every successful PC-Relate sensitivity, with explicit training composition.
rows=[];training=[]
for panel in ['supplement17','strict17','complete17']:
 for th in ['022','044','088']:
  folder=r/'relatedness'/panel/('threshold_'+th);train=(folder/'training_set.txt').read_text().splitlines();training.append({'panel':panel,'threshold':float(th)/1000,'training_n':len(train),**{gr:sum(s.startswith(gr) for s in train) for gr in ['AIC','HCO','IZU','DRR']}})
  for k in [1,2]:
   x=pd.read_csv(folder/f'PC{k}_pairs.tsv',sep='\t');assert len(x)==136 and np.isfinite(x.kin).all();x['panel']=panel;x['threshold']=float(th)/1000;x['PCs']=k;rows.append(x)
allpc=pd.concat(rows,ignore_index=True);allpc.to_csv(res/'pcrelate_all_sensitivities.tsv',sep='\t',index=False);pd.DataFrame(training).to_csv(res/'pcrelate_training_composition.tsv',sep='\t',index=False)
checks['PCRelate_successful_configurations']=len(rows)
# Priority pairs based on convergent within-locality rankings, not pedigree calls.
qc=pd.read_csv(res/'pairwise_qc.tsv',sep='\t');prior=[]
for i,j in [('HCO-F9-S7','HCO-M6-S5'),('IZU-F5-S2','IZU-F8-S4')]:
 sel=allpc[(allpc.ID1.eq(i)&allpc.ID2.eq(j))|(allpc.ID1.eq(j)&allpc.ID2.eq(i))];q=qc[(qc.sample1.eq(i)&qc.sample2.eq(j))|(qc.sample1.eq(j)&qc.sample2.eq(i))].iloc[0]
 vals=[]
 for panel in ['supplement17','strict17','complete17']:
  km=pd.read_csv(r/'relatedness'/panel/'king_matrix.tsv',sep='\t',index_col=0);vals.append(km.loc[i,j])
 prior.append({'sample1':i,'sample2':j,'KING_min':min(vals),'KING_max':max(vals),'PCRelate_min':sel.kin.min(),'PCRelate_max':sel.kin.max(),'primary_IBS0':q.IBS0_fraction,'interpretation':'priority relatedness candidate; pedigree degree unresolved'})
pd.DataFrame(prior).to_csv(res/'priority_relatedness_pairs.tsv',sep='\t',index=False)
# Two compact figures used with network and NJ.
fig,ax=plt.subplots(figsize=(10,9));im=ax.imshow(m,cmap='viridis',vmin=0);ax.set_xticks(range(17),m.index,rotation=90,fontsize=8);ax.set_yticks(range(17),m.index,fontsize=8);fig.colorbar(im,ax=ax,label='Pairwise-called diploid IBS distance');ax.set_title('17 DT individuals: shared distance matrix');fig.tight_layout();fig.savefig(res/'distance_heatmap.png',dpi=200);fig.savefig(res/'distance_heatmap.pdf');plt.close(fig)
bs=pd.read_csv(res/'baseline_split_subsampling.tsv',sep='\t');fig,ax=plt.subplots(figsize=(10,8));labs=bs.side.str.replace(',', ' + ',regex=False);ax.barh(range(len(bs)),bs.NJ_recoveries,color='#377c9b');ax.set_yticks(range(len(bs)),labs,fontsize=7);ax.set_xlim(0,105);ax.set_xlabel('Recovery count in 100 random half-contig subsets (not bootstrap confidence)');ax.set_title('Broad groups stable; some within-locality splits variable');
for i,x in enumerate(bs.NJ_recoveries):ax.text(x+.5,i,str(x),va='center',fontsize=8)
fig.tight_layout();fig.savefig(res/'split_stability.png',dpi=200,bbox_inches='tight');fig.savefig(res/'split_stability.pdf',bbox_inches='tight');plt.close(fig)
checks['source_hash_bound']=json.loads((r/'data/panels.json').read_text())['source_sha256'];checks['status']='PASS_ENGINEERING_CHECKS_NOT_BIOLOGICAL_VALIDATION'
(r/'VALIDATION.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
