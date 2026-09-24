#!/usr/bin/env python3
"""Create figures/tables with explicitly descriptive inference boundaries."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS={'AIC':'#d97706','HCO':'#2563eb','IZU':'#059669','PUBLIC_INDIVIDUAL':'#9333ea'}
PANELS=('main14','strict14','complete14','supplement17')
SEEDS=(20260914,20260915,20260916)


def read_bed(prefix):
    prefix=Path(prefix)
    with prefix.with_suffix('.fam').open() as f:n=sum(1 for _ in f)
    with prefix.with_suffix('.bim').open() as f:m=sum(1 for _ in f)
    raw=np.fromfile(prefix.with_suffix('.bed'),dtype=np.uint8)
    assert list(raw[:3])==[108,27,1]
    assert len(raw)==3+m*((n+3)//4)
    packed=raw[3:].reshape(m,(n+3)//4)
    codes=((packed[:,:,None] >> (2*np.arange(4,dtype=np.uint8))) & 3).reshape(m,-1)[:,:n]
    return np.array([2.,np.nan,1.,0.])[codes]


def align_q(q,reference):
    perm=min(itertools.permutations(range(q.shape[1])),key=lambda p:float(np.square(q[:,p]-reference).sum()))
    return q[:,perm]


def correlation(a,b):
    if np.std(a)==0 or np.std(b)==0:return None
    return float(np.corrcoef(a,b)[0,1])


def read_pca(directory,stem='pca'):
    df=pd.read_csv(directory/f'{stem}.eigenvec',sep=r'\s+')
    df=df.rename(columns={'IID':'sample'})
    ev=np.loadtxt(directory/f'{stem}.eigenval')
    return df,ev/ev[ev>0].sum()*100


def plot_pca(ax,df,ev,metadata,title):
    for g,color in COLORS.items():
        rows=df[df['sample'].map(metadata['group'])==g]
        if rows.empty:continue
        ax.scatter(rows.PC1,rows.PC2,c=color,label=g.replace('_INDIVIDUAL',''),s=48,edgecolor='white',linewidth=.6)
        if g=='PUBLIC_INDIVIDUAL':
            for _,r in rows.iterrows():
                label={'DRR351754':'Mie1','DRR351755':'Yuga','DRR351756':'Utsu'}[r['sample']]
                ax.annotate(label,(r.PC1,r.PC2),xytext=(6,5),textcoords='offset points',fontsize=9)
        else:
            ax.annotate(f'{g} (n={len(rows)})',(rows.PC1.mean(),rows.PC2.mean()),xytext=(10,12),textcoords='offset points',fontsize=10,fontweight='bold',color=color)
    # Pairwise missing-data relationship matrices can have negative eigenvalues.
    # Do not label positive-spectrum fractions as literal explained variance.
    ax.set(xlabel='PC1',ylabel='PC2',title=title)
    ax.margins(.18);ax.spines[['top','right']].set_visible(False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    assert (a.analysis/'STRUCTURE_COMPLETE').exists()
    metadata=pd.read_csv(a.data/'supplement17.manifest.tsv',sep='\t').set_index('sample')
    original_qc=pd.read_csv(a.audit/'main.sample_qc.bcftools.tsv',sep='\t').set_index('sample')
    counts=json.loads((a.analysis/'pca_counts.json').read_text())
    filter_summary=json.loads((a.data/'filter_summary.json').read_text())
    summary={'filter':filter_summary,'counts':counts,'pca':{},'sensitivity':{},'admixture':{},'fst':[]}
    fig,axs=plt.subplots(2,2,figsize=(12,10),layout='constrained')
    pcas={};sampleqcs=[]
    titles={'main14':'Primary: 14 study samples','strict14':'Stricter genotype filters: 14 samples',
        'complete14':'Complete genotypes: 14 samples','supplement17':'Supplement: + 3 public individuals'}
    for name,ax in zip(PANELS,axs.flat):
        d=a.analysis/name
        df,ev=read_pca(d);pcas[name]=df.set_index('sample')
        fam=pd.read_csv(d/'pruned.fam',sep=r'\s+',header=None)
        assert list(df['sample'])==list(fam[1])
        plot_pca(ax,df,ev,metadata,titles[name])
        df.to_csv(a.out/f'{name}.pca.tsv',sep='\t',index=False)
        imiss=pd.read_csv(d/'pca.imiss',sep=r'\s+').set_index('IID').loc[df['sample']]
        allmiss=pd.read_csv(d/'qc.imiss',sep=r'\s+').set_index('IID').loc[df['sample']]
        mat=read_bed(d/'pruned')
        missing=np.isnan(mat).astype(float);centered=missing-missing.mean(axis=1,keepdims=True)
        vals,vecs=np.linalg.eigh(centered.T@centered)
        metrics={f'PC{i}':{'missing_fraction_r':correlation(df[f'PC{i}'],imiss.F_MISS),
            'original_SNP_mean_depth_r':correlation(df[f'PC{i}'],original_qc.loc[df['sample'],'mean_depth_bcftools']),
            'missing_pattern_PC1_r':correlation(df[f'PC{i}'],vecs[:,-1]) if vals[-1]>1e-12 else None} for i in (1,2)}
        for j,sample in enumerate(df['sample']):
            sampleqcs.append({'panel':name,'sample':sample,'group':metadata.loc[sample,'group'],
                'missing_after_site_filter_MAC2':float(allmiss.iloc[j].F_MISS),
                'missing_on_PCA_panel':float(imiss.iloc[j].F_MISS),
                'original_SNP_mean_depth':float(original_qc.loc[sample,'mean_depth_bcftools']),
                'heterozygosity_on_PCA_SNPs':float(np.sum(mat[:,j]==1)/np.sum(~np.isnan(mat[:,j])))})
        one,oneev=read_pca(d,'one_contig_pca')
        one=one.set_index('sample').loc[df['sample']]
        summary['pca'][name]={'PC1_positive_spectrum_percent':float(ev[0]),'PC2_positive_spectrum_percent':float(ev[1]),
            'raw_relationship_eigenvalues':np.loadtxt(d/'pca.eigenval').tolist(),
            'qc_correlations':metrics,'one_per_contig_distance_correlation':correlation(pdist(df[['PC1','PC2']]),pdist(one[['PC1','PC2']]))}
    axs[0,0].legend(fontsize=8,frameon=False)
    fig.suptitle('Dioscorea tenuipes | Population structure sensitivity checks',fontsize=15)
    for ext in ('png','pdf','svg'):fig.savefig(a.out/f'pca_comparison.{ext}',dpi=200)
    plt.close(fig)
    pd.DataFrame(sampleqcs).to_csv(a.out/'sample_qc.tsv',sep='\t',index=False)
    main=pcas['main14']
    for name in ('strict14','complete14','supplement17'):
        other=pcas[name].loc[main.index]
        summary['sensitivity'][name]={'PC12_pairwise_distance_r':correlation(pdist(main[['PC1','PC2']]),pdist(other[['PC1','PC2']]))}
    # Report descriptive between-group and within-group PCA distances.
    groups=metadata.loc[main.index,'group'].to_numpy();coords=main[['PC1','PC2']].to_numpy()
    distances=[]
    for i in range(len(main)):
        for j in range(i+1,len(main)):
            distances.append({'sample1':main.index[i],'sample2':main.index[j],
                'group_pair':'-'.join(sorted([groups[i],groups[j]])),
                'PC12_distance':float(np.linalg.norm(coords[i]-coords[j]))})
    pd.DataFrame(distances).to_csv(a.out/'pca_pairwise_distances.tsv',sep='\t',index=False)

    cvrows=[]
    qbest={}
    for name in ('main14','strict14','supplement17'):
        byk={}
        for k in range(1,5):
            runs=[]
            for seed in SEEDS:
                d=a.analysis/name/'admixture'/f'K{k}_seed{seed}';log=(d/'run.log').read_text()
                cv=float(re.search(r'CV error \(K=\d+\):\s*([\d.eE+-]+)',log).group(1))
                ll=float(re.findall(r'Loglikelihood:\s*([\d.eE+-]+)',log)[-1])
                q=np.loadtxt(d/f'panel.{k}.Q',ndmin=2)
                assert q.shape==(counts[name]['samples'],k)
                assert np.allclose(q.sum(axis=1),1,atol=1e-5)
                runs.append({'seed':seed,'cv':cv,'loglikelihood':ll,'q':q})
                cvrows.append({'panel':name,'K':k,'seed':seed,'CV_error':cv,'loglikelihood':ll})
            best=max(runs,key=lambda r:r['loglikelihood'])
            fam=pd.read_csv(a.analysis/name/'pruned.fam',sep=r'\s+',header=None)
            if name!='main14':
                refids=pd.read_csv(a.analysis/'main14'/'pruned.fam',sep=r'\s+',header=None)[1].tolist()
                common=[fam[1].tolist().index(s) for s in refids]
                perm=min(itertools.permutations(range(k)),key=lambda x:float(np.square(best['q'][common][:,x]-qbest['main14',k]).sum()))
                best['q']=best['q'][:,perm]
            qbest[name,k]=best['q']
            byk[str(k)]={'CV_mean':float(np.mean([r['cv'] for r in runs])),
                'CV_min':min(r['cv'] for r in runs),'CV_max':max(r['cv'] for r in runs),
                'best_full_fit_seed':best['seed'],
                'largest_aligned_Q_difference':max(float(np.max(np.abs(align_q(r['q'],best['q'])-best['q']))) for r in runs)}
            qdf=pd.DataFrame(best['q'],columns=[f'component{i+1}' for i in range(k)])
            qdf.insert(0,'sample',fam[1]);qdf.to_csv(a.out/f'{name}.K{k}.bestQ.tsv',sep='\t',index=False)
        summary['admixture'][name]={'K':byk,'lowest_mean_CV_K':min(range(1,5),key=lambda k:byk[str(k)]['CV_mean'])}
    pd.DataFrame(cvrows).to_csv(a.out/'admixture_runs.tsv',sep='\t',index=False)
    fig,axs=plt.subplots(4,2,figsize=(13,12),layout='constrained')
    for col,name in enumerate(('main14','supplement17')):
        fam=pd.read_csv(a.analysis/name/'pruned.fam',sep=r'\s+',header=None)
        for row,k in enumerate((1,2,3,4)):
            ax=axs[row,col];q=qbest[name,k];bottom=np.zeros(len(q))
            for component in range(k):
                ax.bar(np.arange(len(q)),q[:,component],bottom=bottom,width=.9,label=f'C{component+1}')
                bottom+=q[:,component]
            ax.set(ylim=(0,1),ylabel=f'K={k} | component fraction',title=titles[name] if row==0 else '')
            ax.set_xticks(np.arange(len(q)),fam[1],rotation=70,ha='right',fontsize=7)
            ax.spines[['top','right']].set_visible(False)
    fig.suptitle('ADMIXTURE | CV favors K=1; K=2-4 show conditional fits, not proven ancestry',fontsize=12)
    for ext in ('png','pdf','svg'):fig.savefig(a.out/f'admixture_comparison.{ext}',dpi=200)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for name in ('main14','strict14','supplement17'):
        ks=summary['admixture'][name]['K'];mean=[ks[str(k)]['CV_mean'] for k in range(1,5)]
        ax.plot(range(1,5),mean,'o-',label=name)
        ax.fill_between(range(1,5),[ks[str(k)]['CV_min'] for k in range(1,5)],[ks[str(k)]['CV_max'] for k in range(1,5)],alpha=.15)
    ax.set(xlabel='K',ylabel='5-fold CV error',xticks=range(1,5));ax.legend();fig.savefig(a.out/'admixture_cv.png',dpi=200);plt.close(fig)

    for name in ('main14','strict14'):
        for pair in ('AIC_HCO','AIC_IZU','HCO_IZU'):
            log=(a.analysis/name/'fst'/f'{pair}.stdout.log').read_text()
            weighted=re.search(r'weighted Fst estimate:\s*([\d.eE+-]+)',log)
            mean=re.search(r'mean Fst estimate:\s*([\d.eE+-]+)',log)
            assert weighted and mean
            summary['fst'].append({'panel':name,'pair':pair,'weighted_WC84_FST':float(weighted.group(1)),
                'mean_site_FST_not_primary':float(mean.group(1)),
                'scope':'filtered polymorphic SNPs; descriptive; no callable-genome or demographic inference'})
    pd.DataFrame(summary['fst']).to_csv(a.out/'pairwise_fst.tsv',sep='\t',index=False)
    (a.out/'results_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (a.out/'SUMMARY_COMPLETE').write_text('All counts, sample-order/Q checks, figures and summary tables complete.\n')


if __name__=='__main__':main()
