#!/usr/bin/env python3
"""Build the authorized38 PCA panels and six full VCFtools FST checks.

All genotypes originate from the unchanged masked17 source. PLINK filters are
cross-checked by reproducing the original full14 PCA marker set and genotypes.
"""
import argparse,concurrent.futures,hashlib,itertools,json,shlex,shutil,subprocess
from pathlib import Path
import numpy as np

def run(cmd,log):
 with log.open('w') as f:
  f.write(shlex.join(list(map(str,cmd)))+'\n');f.flush()
  subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True)

def nlines(p):
 with p.open() as f:return sum(1 for _ in f)

def marker_hash(p):
 h=hashlib.sha256()
 with p.open() as f:
  for line in f:
   x=line.split();h.update(('\t'.join([x[0],x[3],*sorted(x[4:6])])+'\n').encode())
 return h.hexdigest()

def make_scenarios(samples):
 def entry(name,kind,mode,keep):
  return {'name':name,'kind':kind,'mode':mode,'keep':keep,'samples':keep,'omitted':[s for s in samples if s not in keep],'group_counts':{g:sum(s.startswith(g+'-') for s in keep) for g in ['AIC','HCO','IZU']}}
 panels=[entry('baseline_main','baseline','fixed',samples),entry('baseline_complete','baseline','complete',samples)]
 for i,omit in enumerate(samples):
  for mode in ['fixed','refilter']:panels.append(entry(f'loo_{i:02d}_{mode}','loo',mode,[s for s in samples if s!=omit]))
 hco=['HCO-F9-S7','HCO-M6-S5'];izu=['IZU-F5-S2','IZU-F8-S4']
 assert set(hco+izu)<=set(samples)
 for i,(h,z) in enumerate(itertools.product(hco,izu)):
  keep=[s for s in samples if s not in set(hco+izu)-{h,z}]
  for mode in ['fixed','refilter']:panels.append(entry(f'pair_{i:02d}_{mode}','pair',mode,keep))
 assert len(panels)==38 and len({p['name'] for p in panels})==38
 return panels

def main():
 p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--tools',type=Path,required=True);a=p.parse_args()
 d=a.data.resolve();o=a.out.resolve();o.mkdir(parents=True,exist_ok=False);(o/'panels').mkdir();(o/'vcftools').mkdir();(o/'baseline_original').mkdir()
 samples=(d/'main14.samples.txt').read_text().splitlines();assert len(samples)==14 and len(set(samples))==14;assert all(s.startswith(('AIC-','HCO-','IZU-')) for s in samples)
 assert {'AIC':sum(s.startswith('AIC-') for s in samples),'HCO':sum(s.startswith('HCO-') for s in samples),'IZU':sum(s.startswith('IZU-') for s in samples)}=={'AIC':4,'HCO':4,'IZU':6}
 (o/'own14.keep').write_text(''.join(f'{s}\t{s}\n' for s in samples))
 panels=make_scenarios(samples);spec={'sample_order':samples,'populations':{g:[s for s in samples if s.startswith(g+'-')] for g in ['AIC','HCO','IZU']},'panels':panels,'expected_baseline_sites':5564311,'expected_complete_sites':3367199,'scope':'14LOO +4paircombinations, fixed/refilter; full14 complete FST; no new rawreads'}
 (o/'panels.json').write_text(json.dumps(spec,indent=2)+'\n')
 for suffix in ['eigenvec','eigenval','log']:
  shutil.copy2(d/f'structure/main14/pca.{suffix}',o/f'baseline_original/pca.{suffix}')
 inputs=[d/'masked17.bcf',d/'main14.bcf',d/'complete14.bcf',d/'main14.samples.txt']
 inputs +=[d/f'structure/{panel}/pruned.{ext}' for panel in ['main14','complete14'] for ext in ['bed','bim','fam']]
 hashes={str(x):hashlib.file_digest(x.open('rb'),'sha256').hexdigest() for x in inputs}
 (o/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
 base=[a.tools/'plink','--allow-extra-chr','--keep-allele-order','--threads','2','--memory','3500']
 def source_bed():
  run(base+['--bcf',d/'masked17.bcf','--double-id','--keep',o/'own14.keep','--set-missing-var-ids','@:#:$1:$2','--make-bed','--out',o/'source_own14'],o/'source_import.log')
  observed=[x.split()[1] for x in (o/'source_own14.fam').read_text().splitlines()];assert observed==samples
  assert (o/'source_own14.bed').stat().st_size==3+nlines(o/'source_own14.bim')*((len(samples)+3)//4)
  (o/'SOURCE_BED_COMPLETE').write_text('Unfiltered SNP-major own14 genotypes imported from masked17; sample order and dimensions checked.\n');print('SOURCE_BED_COMPLETE',flush=True)
 def fst_check(name,pair):
  sub=o/'vcftools'/name;sub.mkdir(exist_ok=True);src=d/('main14.bcf' if name=='baseline_main' else 'complete14.bcf');prefix=sub/pair
  g1,g2=pair.split('_');cmd1=['bcftools','view',str(src),'-Ov'];cmd2=[str(a.tools/'vcftools'),'--vcf','-','--weir-fst-pop',str(d/f'{g1}.samples.txt'),'--weir-fst-pop',str(d/f'{g2}.samples.txt'),'--fst-window-size','100000','--fst-window-step','100000','--out',str(prefix)]
  with prefix.with_suffix('.stdout.log').open('w') as log:
   log.write(shlex.join(cmd1)+' | '+shlex.join(cmd2)+'\n');log.flush();p1=subprocess.Popen(cmd1,stdout=subprocess.PIPE,stderr=log);p2=subprocess.Popen(cmd2,stdin=p1.stdout,stdout=log,stderr=log);p1.stdout.close();rc2=p2.wait();rc1=p1.wait();assert rc1==rc2==0
  print('VCFTOOLS_COMPLETE',name,pair,flush=True)
 def pca_panel(info):
  name=info['name'];out=o/'panels'/name;out.mkdir();keep=out/'keep.txt';keep.write_text(''.join(f'{s}\t{s}\n' for s in info['keep']))
  if info['mode']=='refilter' or name=='baseline_main':
   run(base+['--bfile',o/'source_own14','--keep',keep,'--geno','0.1','--mac','2','--make-bed','--out',out/'all_mac2'],out/'filter.log')
   run(base+['--bfile',out/'all_mac2','--bp-space','1000','--make-bed','--out',out/'spaced'],out/'space.log')
   run(base+['--bfile',out/'spaced','--indep-pairwise','50kb','1','0.2','--out',out/'ld'],out/'ld.stdout.log')
   run(base+['--bfile',out/'spaced','--extract',out/'ld.prune.in','--make-bed','--out',out/'pruned'],out/'extract.log')
   assert nlines(out/'pruned.bim')==nlines(out/'ld.prune.in')
  else:
   old='complete14' if info['mode']=='complete' else 'main14'
   run(base+['--bfile',d/f'structure/{old}/pruned','--keep',keep,'--make-bed','--out',out/'pruned'],out/'extract.log')
  found=[x.split()[1] for x in (out/'pruned.fam').read_text().splitlines()];assert found==info['keep']
  run(base+['--bfile',out/'pruned','--pca',str(len(found)),'header','tabs','--missing','--freqx','--out',out/'pca'],out/'pca.stdout.log')
  q=dict(info,sample_count=len(found),pca_snps=nlines(out/'pruned.bim'),normalized_marker_sha256=marker_hash(out/'pruned.bim'))
  (out/'panel.json').write_text(json.dumps(q,indent=2)+'\n');print('PCA_COMPLETE',name,q['sample_count'],q['pca_snps'],flush=True)
  return q
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  source=pool.submit(source_bed)
  fullfst=[pool.submit(fst_check,name,pair) for name in ['baseline_main','baseline_complete'] for pair in ['AIC_HCO','AIC_IZU','HCO_IZU']]
  source.result()
  baseline=pca_panel(panels[0]);old=d/'structure/main14'
  assert baseline['pca_snps']==37566
  assert baseline['normalized_marker_sha256']==marker_hash(old/'pruned.bim')
  assert (o/'panels/baseline_main/pruned.bed').read_bytes()==(old/'pruned.bed').read_bytes()
  x=np.loadtxt(o/'panels/baseline_main/pca.eigenvec',skiprows=1,usecols=(2,3));y=np.loadtxt(old/'pca.eigenvec',skiprows=1,usecols=(2,3));dx=np.sqrt(((x[:,None]-x[None,:])**2).sum(2));dy=np.sqrt(((y[:,None]-y[None,:])**2).sum(2));assert np.max(np.abs(dx-dy))<1e-5
  check={'sample_order_matches':True,'normalized_markers_match':True,'BED_bytes_match':True,'PCA_distance_max_abs_difference':float(np.max(np.abs(dx-dy))),'pca_snps':37566}
  (o/'baseline_pca_validation.json').write_text(json.dumps(check,indent=2)+'\n');print('BASELINE_PCA_VERIFIED',check,flush=True)
  [f.result() for f in fullfst[:3]]
  newjobs=[pool.submit(pca_panel,info) for info in panels[1:]]
  results=[baseline]+[f.result() for f in newjobs]
  [f.result() for f in fullfst]
  assert len(results)==38
 assert {str(x):hashlib.file_digest(x.open('rb'),'sha256').hexdigest() for x in inputs}==hashes
 (o/'pca_panels_summary.json').write_text(json.dumps(results,indent=2)+'\n')
 (o/'PCA_COMPLETE').write_text('All38PCA panels and6fullVCFtools checks complete; original inputs unchanged.\n');print('PCA_ALL_COMPLETE',flush=True)
if __name__=='__main__':main()
