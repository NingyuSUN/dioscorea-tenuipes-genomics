#!/usr/bin/env python3
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
root=Path(sys.argv[1]);p=root/'published_comparison';s=json.loads((p/'comparison_summary.json').read_text())
fig,axes=plt.subplots(3,1,figsize=(10,7),layout='constrained',sharex=True,gridspec_kw={'height_ratios':[1.2,1,1.5]})
colors=['#88b9da','#d8c494','#94c5a0','#d8c494'];names=['LSC','IR','SSC','IR']
for y,who,label in [(1,'paper','Published: PV612491.1'),(0,'own','Recovered: Izu candidate')]:
 r=s['region_lengths'][who];lengths=[r['LSC'],r['IR_each'],r['SSC'],r['IR_each']];pos=0
 for length,c,n in zip(lengths,colors,names):
  axes[0].barh(y,length/1000,left=pos/1000,height=.5,color=c,edgecolor='white');axes[0].text((pos+length/2)/1000,y,n,ha='center',va='center',fontsize=9);pos+=length
 axes[0].text(pos/1000+1,y,f'{pos:,} bp',va='center',fontsize=9)
axes[0].set_yticks([0,1],['Izu candidate','PV612491.1']);axes[0].set_title('Plastome comparison: structure, sequence differences and retained-read support',loc='left',fontsize=12)
d=pd.read_csv(p/'substitutions_for_read_validation.tsv',sep='\t');bins=np.arange(0,155001,1000);counts,_=np.histogram(d['own_path2_position'],bins=bins)
axes[1].bar(bins[:-1]/1000,counts,width=1,align='edge',color='#ad5865');axes[1].set_ylabel('Substitutions\nper 1 kb');axes[1].set_title(f"{s['difference_counts']['substitutions']} aligned substitutions; both IR copies included",loc='left',fontsize=10)
z=pd.read_csv(root/'results/plastome.depth.tsv',sep='\t')
mask=(z['position_1based']>109153)&(z['position_1based']<=128097)
for col in ['depth_BQ20','depth_BQ20_MAPQ20']:z.loc[mask,col]=z.loc[mask,col].to_numpy()[::-1]
g=z.groupby((z['position_1based']-1)//250).mean(numeric_only=True)
axes[2].plot(g.position_1based/1000,g.depth_BQ20,label='BQ >= 20, all mapping qualities',color='#357f99',lw=1)
axes[2].plot(g.position_1based/1000,g.depth_BQ20_MAPQ20,label='BQ >= 20, MAPQ >= 20',color='#bd8643',lw=1)
axes[2].set_ylabel('Read depth');axes[2].set_xlabel('Position (kb; panels B/C use own plastome path 2)');axes[2].legend(fontsize=8,loc='lower left');axes[2].set_title('Retained assembly reads; 250-bp mean. Low unique depth in repeated regions.',loc='left',fontsize=10)
for ax in axes:
 ax.spines[['top','right']].set_visible(False);ax.set_xlim(0,172)
fig.text(.01,-.015,'SSC orientation standardized for comparison. Read validation is not independent of assembly; no gene annotation shown.',fontsize=8)
for suffix in ['png','pdf','svg']:fig.savefig(p/f'plastome_comparison.{suffix}',dpi=180,bbox_inches='tight')
