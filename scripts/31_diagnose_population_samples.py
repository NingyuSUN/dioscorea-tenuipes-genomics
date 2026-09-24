#!/usr/bin/env python3
"""Deterministic, geographically dispersed genomic windows for VCF QC diagnostics.
Not an LD-pruned population panel; thresholds are diagnostic, not final filtering.
"""
import argparse,bisect,collections,csv,json,random
from pathlib import Path
import pysam,numpy as np
p=argparse.ArgumentParser();p.add_argument('--vcf',required=True);p.add_argument('--fai',required=True);p.add_argument('--out',required=True);a=p.parse_args()
o=Path(a.out);o.mkdir(parents=True,exist_ok=False)
lengths=[(x[0],int(x[1])) for l in open(a.fai) if (x:=l.split())];cum=[];total=0
for _,n in lengths:total+=n;cum.append(total)
rng=random.Random(20260914);regions=collections.defaultdict(list)
for pos in rng.sample(range(total),1000):
 i=bisect.bisect_right(cum,pos);start=max(0,pos-(cum[i-1] if i else 0)-500);end=min(lengths[i][1],start+1000);regions[lengths[i][0]].append((start,end))
merged=[]
for chrom,_ in lengths:
 for start,end in sorted(regions[chrom]):
  if merged and merged[-1][0]==chrom and start<=merged[-1][2]:merged[-1]=(chrom,merged[-1][1],max(end,merged[-1][2]))
  else:merged.append((chrom,start,end))
(o/'windows.bed').write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in merged))
v=pysam.VariantFile(a.vcf);names=list(v.header.samples);qc=[collections.Counter() for _ in names];balance=[[] for _ in names];dosages=[];n=0
for chrom,start,end in merged:
 for r in v.fetch(chrom,start,end):
  if r.start<start or r.start>=end:continue
  n+=1;row=[]
  for i,s in enumerate(r.samples.values()):
   gt=s.get('GT');dp=s.get('DP');gq=s.get('GQ');ad=s.get('AD');q=qc[i];q['records']+=1
   called=gt is not None and len(gt)==2 and None not in gt
   good=called and dp is not None and dp>=5 and gq is not None and gq>=20
   q['called']+=called;q['good_dp5_gq20']+=good;row.append(sum(gt) if good else -1)
   q['homref_dp0_or_missing']+=called and gt==(0,0) and (dp is None or dp==0)
   if good and sum(gt)==1 and dp>=10 and ad is not None and len(ad)==2 and all(x is not None and x>=0 for x in ad) and sum(ad)>=10:
    balance[i].append(ad[1]/sum(ad))
  dosages.append(row)
data=np.array(dosages,dtype=np.int8);pairs=[]
# Pairwise similarity is only a duplicate/outlier screen, not a kinship estimator.
for i in range(len(names)):
 for j in range(i):
  ok=(data[:,i]>=0)&(data[:,j]>=0);nr=ok&((data[:,i]>0)|(data[:,j]>0))
  pairs.append(dict(sample1=names[j],sample2=names[i],shared_good=int(ok.sum()),discordant=int((data[ok,i]!=data[ok,j]).sum()),shared_with_any_nonref=int(nr.sum()),discordant_nonref=int((data[nr,i]!=data[nr,j]).sum())))
rows=[]
for name,q,ab in zip(names,qc,balance):
 rows.append(dict(sample=name,**dict(q),allele_balance_hets=len(ab),allele_balance_median=float(np.median(ab)) if ab else None,ab_fraction_035_065=float(np.mean((np.array(ab)>=.35)&(np.array(ab)<=.65))) if ab else None))
for filename,records in [('sample_diagnostics.tsv',rows),('pairwise_screen.tsv',pairs)]:
 with (o/filename).open('w') as f:
  w=csv.DictWriter(f,fieldnames=records[0],delimiter='\t');w.writeheader();w.writerows(records)
with (o/'allele_balance_histogram.tsv').open('w') as f:
 f.write('sample\tbin_start\tbin_end\theterozygous_calls\n')
 for name,ab in zip(names,balance):
  hist,edges=np.histogram(ab,bins=np.linspace(0,1,21))
  for k,count in enumerate(hist):f.write(f'{name}\t{edges[k]:.2f}\t{edges[k+1]:.2f}\t{count}\n')
(o/'summary.json').write_text(json.dumps(dict(complete=True,vcf=a.vcf,seed=20260914,windows=len(merged),covered_reference_bases=sum(e-s for _,s,e in merged),records=n,filter='GT diploid/nonmissing, DP>=5, GQ>=20; AB: heterozygous plus DP>=10 and AD sum>=10',scope='Genomic-window diagnostic only, not LD-pruned, not a ploidy/contamination or kinship test',samples=rows),indent=2)+'\n')
print('complete',n,'records',len(merged),'windows')
