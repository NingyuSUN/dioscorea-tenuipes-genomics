#!/usr/bin/env python3
"""Reproducible descriptive DT relationships. No root or pedigree inference."""
import argparse,copy,hashlib,itertools,json,sys
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
from splitspy.splits.basic_split import split_dist
from splitspy.splits.splits_io import print_splits_nexus
from splitspy.outlines.outline_algo import compute
from relationship_core import read_bed,ibs_distance,nj_tree,tree_splits,canonical_split,residual_fit,checked_neighbor_net
COLORS={'AIC':'#c65a29','HCO':'#287f9a','IZU':'#679b39','public':'#794ca1'}
PUBLIC={'DRR351754':'Mie1','DRR351755':'Yuga','DRR351756':'Utsu'}
def group(s):return s[:3] if s[:3] in COLORS else 'public'
def label(s):return s+(' / '+PUBLIC[s] if s in PUBLIC else '')
def load_panel(p):
    fam=pd.read_csv(p/'pruned.fam',sep=r'\s+',header=None);bim=pd.read_csv(p/'pruned.bim',sep=r'\s+',header=None,dtype=str)
    assert bim[1].is_unique and not (bim[1]=='.').any()
    names=fam[1].tolist();g=read_bed(p/'pruned.bed',len(names),len(bim));return names,g,bim

def network_plot(names,cycle,splits,out,title):
    graph,_=compute(names,cycle,splits,rooted=False)
    fig,ax=plt.subplots(figsize=(12,9))
    for e in graph.edges():
        a,b=e.src().pos,e.tar().pos;ax.plot([a[0],b[0]],[a[1],b[1]],color='#777777',lw=.7,alpha=.75)
    coords=np.array([v.pos for v in graph.nodes() if v.label]);centre=coords.mean(axis=0)
    for v in graph.nodes():
        if v.label:
            x,y=v.pos;c=COLORS[group(v.label)];ax.scatter(x,y,s=26,color=c,zorder=3)
            ax.annotate(label(v.label),(x,y),xytext=(5 if x>=centre[0] else -5,3),textcoords='offset points',ha='left' if x>=centre[0] else 'right',color=c,fontsize=9)
    xmin,xmax,ymin,ymax=graph.bbox();ax.plot([xmin,xmin+.05],[ymin-.035,ymin-.035],color='black',lw=1);ax.text(xmin+.025,ymin-.045,'0.05 IBS distance',ha='center',fontsize=8)
    ax.set_aspect('equal');ax.margins(.30);ax.axis('off');ax.set_title(title+'\nNeighbor-Net split system, phylogenetic outline (unrooted)',fontsize=13)
    fig.text(.05,.025,'Pairwise-called diploid IBS distance; boxes are conflicting splits, not a test of hybridization.',fontsize=9)
    for ext in ['png','pdf','svg']:fig.savefig(out/f'network.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

def projected_splits(splits,names,own):
    result={}
    for s in splits:
        key=canonical_split({names[i-1] for i in s.part1()}&own,own)
        if key:result[key]=result.get(key,0.)+float(s.weight)
    return result

def analyze(name,names,g,out,reference=None,plot=False,distance=None):
    d,joint=ibs_distance(g) if distance is None else distance
    tree=nj_tree(names,d);cycle,splits,solver=checked_neighbor_net(names,d)
    fit=residual_fit(d,np.asarray(split_dist(len(names),splits))[1:,1:])
    own={s for s in names if group(s)!='public'};t14=copy.deepcopy(tree)
    for s in set(names)-own:t14.prune(s)
    ts=tree_splits(t14);ns=projected_splits(splits,names,own)
    row={'panel':name,'n':len(names),'snps':len(g),'network_fit_pct':fit,'min_joint_sites':int(joint[np.triu_indices(len(names),1)].min()),'negative_nj_branches':sum(c.branch_length is not None and c.branch_length< -1e-12 for c in tree.find_clades())}
    for gr in ['AIC','HCO','IZU']:
        members={s for s in own if group(s)==gr};key=canonical_split(members,own)
        row[f'{gr}_own_split_NJ']=key in ts if key else None
        row[f'{gr}_own_split_NNet_weight']=ns.get(key,0.) if key else None
        conflict=sum(w for sk,w in ns.items() if all([members&set(sk),members-set(sk), (own-members)&set(sk),(own-members)-set(sk)]))
        row[f'{gr}_conflicting_NNet_weight']=conflict
        within=[d[i,j] for i,j in itertools.combinations(range(len(names)),2) if names[i] in members and names[j] in members]
        between=[d[i,j] for i in range(len(names)) for j in range(len(names)) if names[i] in members and names[j] in own-members]
        row[f'{gr}_separation_margin']=float(min(between)-max(within)) if within and between else None
    for p in PUBLIC:
        if p in names:
            means={gr:float(np.mean([d[names.index(p),j] for j,s in enumerate(names) if group(s)==gr])) for gr in ['AIC','HCO','IZU']}
            row[p+'_nearest_locality']=min(means,key=means.get)
            for gr,v in means.items():row[p+'_'+gr+'_mean_distance']=v
    if reference is not None:
        rt=copy.deepcopy(reference)
        for s in {x.name for x in rt.get_terminals()}-set(names):rt.prune(s)
        a,b=tree_splits(rt),tree_splits(tree)
        row['RF_normalized']=len(a^b)/(len(a)+len(b)) if a or b else 0.
        row['baseline_split_recovery']=len(a&b)/len(a) if a else 1.
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(d,index=names,columns=names).to_csv(out/'distance.tsv',sep='\t')
        pd.DataFrame(joint,index=names,columns=names).to_csv(out/'joint_sites.tsv',sep='\t')
        Phylo.write(tree,out/'nj.nwk','newick',format_branch_length='%1.12f')
        records=[{'side1':','.join(names[i-1] for i in sorted(s.part1())),'side2':','.join(names[i-1] for i in sorted(s.part2())),'weight':float(s.weight)} for s in splits]
        pd.DataFrame(records).to_csv(out/'network_splits.tsv',sep='\t',index=False)
        print_splits_nexus(names,splits,cycle,fit,filename=str(out/'network.nex'))
        (out/'summary.json').write_text(json.dumps(row,indent=2,default=lambda x:x.item())+'\n')
        (out/'network_solver.json').write_text(json.dumps(solver,indent=2)+'\n')
        if plot:network_plot(names,cycle,splits,out,name)
    return row,tree,splits

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--replicates',type=int,default=100);a=p.parse_args();out=a.root/'results';out.mkdir(exist_ok=True)
    names,g,bim=load_panel(a.root/'data/supplement17');rows=[];manifest=[]
    row,ref,sp=analyze('supplement17',names,g,out/'supplement17',plot=True);rows.append(row)
    pairs=[]
    for i,j in itertools.combinations(range(len(names)),2):
        ok=(g[:,i]>=0)&(g[:,j]>=0);x,y=g[ok,i],g[ok,j]
        pairs.append({'sample1':names[i],'sample2':names[j],'group1':group(names[i]),'group2':group(names[j]),'joint_sites':int(ok.sum()),'IBS_distance':float(np.mean(np.abs(x-y))/2),'IBS0_fraction':float(np.mean(np.abs(x-y)==2)),'both_het_fraction':float(np.mean((x==1)&(y==1))),'het1_pair':int((x==1).sum()),'het2_pair':int((y==1).sum())})
    pd.DataFrame(pairs).to_csv(out/'pairwise_qc.tsv',sep='\t',index=False)
    pd.DataFrame([{'sample':s,'group':group(s),'missing_fraction':float(np.mean(g[:,i]<0)),'heterozygosity_called':float(np.mean(g[g[:,i]>=0,i]==1))} for i,s in enumerate(names)]).to_csv(out/'sample_qc.tsv',sep='\t',index=False)
    for folder in sorted((a.root/'data').iterdir()):
        if not (folder/'pruned.bed').exists():continue
        for suffix in ['bed','bim','fam']:
            f=folder/f'pruned.{suffix}';manifest.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
        if folder.name=='supplement17':continue
        ns,gg,bb=load_panel(folder);row,_,_=analyze(folder.name,ns,gg,out/folder.name,ref,plot=not folder.name.startswith('loo_'));rows.append(row)
    # Complete calls at exactly the primary selected loci: isolates missingness from refiltering.
    complete=np.all(g>=0,axis=1);row,_,_=analyze('primary17_complete_loci',names,g[complete],out/'primary17_complete_loci',ref);rows.append(row)
    best={}
    for i,r in bim.iterrows():
        h=hashlib.sha256(('20260914:'+r[1]).encode()).hexdigest()
        if r[0] not in best or h<best[r[0]][0]:best[r[0]]=(h,i)
    idx=sorted(v[1] for v in best.values());row,_,_=analyze('one_per_contig17',names,g[idx],out/'one_per_contig17',ref,True);rows.append(row)
    own_idx=[i for i,s in enumerate(names) if group(s)!='public'];row,_,_=analyze('primary17_projected14',[names[i] for i in own_idx],g[:,own_idx],out/'primary17_projected14',ref);rows.append(row)
    # Explicit exploratory pair-representative controls; not pedigree-based exclusions.
    for h in ['HCO-F9-S7','HCO-M6-S5']:
        for z in ['IZU-F5-S2','IZU-F8-S4']:
            ids=[j for j,s in enumerate(names) if s not in {h,z}];nm='pair_control_'+h+'_'+z
            row,_,_=analyze(nm,[names[j] for j in ids],g[:,ids],out/nm,ref);row['omitted']=h+','+z;rows.append(row)
    for i,omit in enumerate(names):
        ids=[j for j in range(len(names)) if j!=i];row,_,_=analyze(f'fixed_loo_{i:02d}',[names[j] for j in ids],g[:,ids],out/f'fixed_loo_{i:02d}',ref);row['omitted']=omit;rows.append(row)
    rng=np.random.default_rng(20260915);contigs=np.array(sorted(best));bs=[];base_splits=tree_splits(ref);support={k:0 for k in base_splits}
    for i in range(a.replicates):
        selected=rng.choice(contigs,len(contigs)//2,replace=False);keep=bim[0].isin(selected).to_numpy()
        row,t,s=analyze(f'half_contigs_{i:03d}',names,g[keep],out/f'half_contigs_{i:03d}',ref)
        row['selected_contigs']=len(selected);row['contig_selection_sha256']=hashlib.sha256('\n'.join(sorted(selected)).encode()).hexdigest();bs.append(row)
        for k in tree_splits(t)&base_splits:support[k]+=1
    pd.DataFrame(rows).to_csv(out/'panel_stability.tsv',sep='\t',index=False);pd.DataFrame(bs).to_csv(out/'half_contig_stability.tsv',sep='\t',index=False)
    pd.DataFrame([{'side':','.join(k),'NJ_recoveries':v,'replicates':a.replicates} for k,v in sorted(support.items())]).to_csv(out/'baseline_split_subsampling.tsv',sep='\t',index=False)
    (out/'input_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(pd.DataFrame(rows).to_string(index=False));print('Completed',len(bs),'half-contig subsets')
if __name__=='__main__':main()
