from pathlib import Path
import csv,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
r=Path(__file__).resolve().parents[1];o=r/'results/full';f=r/'results/figures';f.mkdir(exist_ok=True)
rows=list(csv.DictReader((o/'distance_grid_validated.csv').open()));aa=sorted({float(x['a']) for x in rows});bb=sorted({float(x['b']) for x in rows});z=np.zeros((len(aa),len(bb)))
for x in rows:z[aa.index(float(x['a'])),bb.index(float(x['b']))]=float(x['conditional_tail_probability_analytic'])
plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(np.maximum(z,1e-6),norm=LogNorm(vmin=1e-6,vmax=1),cmap='viridis',aspect='auto')
for i in range(len(aa)):
 for j in range(len(bb)):
  val=z[i,j];label='<1e-6' if val<1e-6 else f'{val:.3g}'
  ax.text(j,i,label,ha='center',va='center',fontsize=8,color='black' if val>.05 else 'white')
fmt=lambda v:'0 (limit)' if v==0 else f'{v:.0e}'
ax.set_xticks(range(len(bb)),[fmt(v) for v in bb],rotation=45,ha='right');ax.set_yticks(range(len(aa)),[fmt(v) for v in aa]);ax.set_xlabel('b = mean additional ancestral coalescent mutation length');ax.set_ylabel('a = minimum divergence mutation length');ax.set_title('Conditional probability of <=1 difference in 123,256 sites\nJC69 sensitivity grid; no empirical ILS P value')
fig.colorbar(im,ax=ax,label='Conditional tail probability (analytic)');fig.text(.02,.01,'a and b are hypothetical substitutions/site parameters. Zero counts are not zero probabilities.\nThis does not estimate Hachi-1 history, mutation rate, or divergence time.',fontsize=9);fig.tight_layout(rect=(0,.065,1,1));fig.savefig(f/'distance_compatibility.png',dpi=180);fig.savefig(f/'distance_compatibility.pdf');plt.close(fig)
rows=list(csv.DictReader((o/'topology_counts.csv').open()));fig,axes=plt.subplots(1,3,figsize=(10,4),sharey=True)
for ax,m in zip(axes,[1,2,4]):
 rr=[x for x in rows if int(x['multiplier'])==m];vals=[float(x['frequency']) for x in rr];ax.bar([x['topology'] for x in rr],vals,color=['#3c78a8','#d1843a','#65a573'])
 for i,v in enumerate(vals):ax.text(i,v+.015,(f'{v:.3%}' if v<.01 else f'{v:.2%}'),ha='center')
 ax.set_ylim(0,1.09);ax.set_title(f'm={m}; n=20,000');ax.set_xlabel('Unrooted quartet')
axes[0].set_ylabel('Fraction of simulated plastid genealogies');fig.suptitle('Conditional fixed nuclear history; H=hime2 label, T=tokoro, G=gracillima, Z=zingiberensis');fig.text(.02,.015,'m=1 is a nuclear-like control; m=2/4 are organellar assumptions, not estimated values.\nObserved H relationship is HG with tokoro OQ525999.1, HT with Kita1; source identity remains conditional.',fontsize=9);fig.tight_layout(rect=(0,.13,1,.93));fig.savefig(f/'genealogy_frequencies.png',dpi=180);fig.savefig(f/'genealogy_frequencies.pdf');plt.close(fig)
print(str(f))
