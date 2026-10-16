#!/usr/bin/env python3
"""PCA/FST masks with frozen original methods and explicit matched-site controls.
Requires a FULL prepare_ab_beds.py result, NumPy, PLINK1.9 and archived FST core.
No ADMIXTURE/demography/selection analysis; no original file overwritten.
"""
import argparse,collections,csv,hashlib,importlib.util,itertools,json,math,subprocess,time,sys
from pathlib import Path
import numpy as np

NAMES=['baseline','ad_control','ab20_80','ab30_70']
EXPECTED_FST_CORE_SHA='58d033294bd33d0cfb4616eba7481e8c6dfbc685c6df75326248a91dd00661d8'
EXPECTED_FST={'AIC_HCO':(.3312913257919642,4538382),'AIC_IZU':(.332169705556553,4776008),'HCO_IZU':(.195661227449393,4300810)}


def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()


def decode(data,n):
 width=(n+3)//4
 packed=np.frombuffer(data,dtype=np.uint8).reshape(-1,width)
 return np.array([0,-1,1,2],dtype=np.int8)[(packed[:,np.arange(n)//4]>>(2*(np.arange(n)%4)))&3]


def marker_key(line):
 fields=line.split()
 return (fields[0],fields[3],*sorted(fields[4:6]))


def validate_baseline_pca_inputs(new_prefix,old_prefix,samples):
 """IDs may differ, but positions, sample order and physical alleles may not."""
 paths=lambda prefix:{e:Path(str(prefix)+'.'+e) for e in ['bim','bed','fam']}
 new,old=paths(new_prefix),paths(old_prefix)
 newlines=new['bim'].read_text().splitlines();oldlines=old['bim'].read_text().splitlines()
 checks={
  'canonical_locus_order_matches':[marker_key(x) for x in newlines]==[marker_key(x) for x in oldlines],
  'ordered_A1_A2_match':[x.split()[4:6] for x in newlines]==[x.split()[4:6] for x in oldlines],
  'new_sample_order_matches':[x.split()[1] for x in new['fam'].read_text().splitlines()]==samples,
  'old_sample_order_matches':[x.split()[1] for x in old['fam'].read_text().splitlines()]==samples,
  'BED_bytes_match':sha(new['bed'])==sha(old['bed'])}
 if not all(checks.values()):raise ValueError('Baseline PCA input identity failure: '+','.join(k for k,v in checks.items() if not v))
 return {'checks':checks,'canonical_loci':len(newlines),
  'literal_variant_ids_identical':[x.split()[1] for x in newlines]==[x.split()[1] for x in oldlines],
  'new_sha256':{k:sha(p) for k,p in new.items()},'old_sha256':{k:sha(p) for k,p in old.items()}}


def preflight(args):
 plink=args.tools/'plink'
 required=[args.fst_core,args.masks/'MASK_RECEIPT.json',plink]
 required += [args.original_main/('pruned.'+e) for e in ['bim','bed','fam']]
 required += [args.original_main/'pca.eigenvec']
 required += [args.masks/name/('source.'+e) for name in NAMES for e in ['bim','bed','fam']]
 missing=[str(p) for p in required if not p.is_file()]
 if missing:raise FileNotFoundError('Missing analysis input/dependency: '+', '.join(missing))
 version=subprocess.check_output([str(plink),'--version'],text=True).strip()
 if 'v1.90b6.21' not in version:raise ValueError('Expected historical PLINK1.90b6.21')
 return {'required_paths':[str(p) for p in required],'plink_version':version,
  'original_main_sha256':{p.name:sha(p) for p in required if p.parent==args.original_main}}


def validate_nested_masks(genotypes):
 base=genotypes['baseline']
 for name,g in genotypes.items():
  if g.shape!=base.shape or np.any((g>=0)&(g!=base)):
   raise ValueError('Changed genotype or missing-call imputation: '+name)
  if np.any((base>=0)&(base!=1)&(g<0)):
   raise ValueError('Only baseline heterozygotes may be additionally masked: '+name)
 for wider,narrower in zip(NAMES,NAMES[1:]):
  if np.any((genotypes[narrower]>=0)&(genotypes[wider]<0)):
   raise ValueError('Non-nested masks: '+wider+' -> '+narrower)


def loss_summary(genotypes,selection=None):
 base=genotypes['baseline'] if selection is None else genotypes['baseline'][selection]
 base_called=(base>=0).sum(0);base_hets=(base==1).sum(0)
 result={}
 for name,g in genotypes.items():
  g=g if selection is None else g[selection]
  result[name]=np.array([base_called,base_hets,(g>=0).sum(0),((base>=0)&(g<0)).sum(0)])
 return len(base),result


def require_estimable_fst(accumulator):
 if accumulator.defined_sites<=0 or not math.isfinite(accumulator.denominator) or accumulator.denominator==0 or not math.isfinite(accumulator.weighted_fst):
  raise ValueError('NOT_ESTIMABLE: no defined finite weighted FST with nonzero denominator')


def metrics(current,reference):
 if current.shape!=reference.shape:raise ValueError('PCA coordinate shapes differ')
 if not np.isfinite(current).all() or not np.isfinite(reference).all():raise ValueError('NOT_ESTIMABLE: nonfinite PCA coordinates')
 d=lambda x:np.sqrt(((x[:,None]-x[None,:])**2).sum(2))[np.triu_indices(len(x),1)]
 dc,dr=d(current),d(reference)
 if len(dc)<2 or any(np.ptp(v)<=10*np.finfo(float).eps*max(1,float(np.max(np.abs(v)))) for v in [dc,dr]):raise ValueError('NOT_ESTIMABLE: constant pairwise distances')
 x=current-current.mean(0);y=reference-reference.mean(0)
 if np.linalg.norm(x)==0 or np.linalg.norm(y)==0:raise ValueError('Degenerate PCA configuration')
 u,s,vt=np.linalg.svd(x.T@y);rotation=u@vt;scale=s.sum()/np.sum(x*x)
 result={'pairwise_distance_r':float(np.corrcoef(dc,dr)[0,1]),'relative_Procrustes_residual':float(np.linalg.norm(x@rotation*scale-y)/np.linalg.norm(y)),'maximum_pairwise_distance_difference':float(np.max(np.abs(dc-dr)))}
 if not all(math.isfinite(v) for v in result.values()):raise ValueError('NOT_ESTIMABLE: nonfinite PCA summary')
 return result


def coordinates(path,samples):
 lines=Path(path).read_text().splitlines();rows=[x.split() for x in lines[1:]]
 if [r[1] for r in rows]!=samples:raise ValueError('PCA sample order differs')
 return np.array([[float(x[2]),float(x[3])] for x in rows])


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--masks',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
 p.add_argument('--fst-core',type=Path,required=True)
 p.add_argument('--tools',type=Path,default=Path('workspace/wgs/popgen_tools_20260914/bin'))
 p.add_argument('--original-main',type=Path,default=Path('workspace/wgs/population_analysis_20260914/structure/main14'))
 p.add_argument('--chunk-sites',type=int,default=200000)
 a=p.parse_args();o=a.out.resolve();o.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 receipt={'status':'RUNNING','commands':[],'script_sha256':sha(__file__),'numpy':np.__version__}
 try:
  receipt['preflight']=preflight(a)
  receipt['commands'].append([str(a.tools/'plink'),'--version'])
  receipt['fst_core_sha256']=sha(a.fst_core)
  receipt['source_mask_receipt_sha256']=sha(a.masks/'MASK_RECEIPT.json')
  if receipt['fst_core_sha256']!=EXPECTED_FST_CORE_SHA:raise ValueError('FST core differs from the archived VCFtools-validated implementation')
  cfg=json.loads((a.masks/'MASK_RECEIPT.json').read_text())
  if cfg['status']!='MASK_COMPLETE' or cfg['scope']!='FULL_SELECTED_SNP_SOURCE' or cfg['sites']!=20625692:raise ValueError('Full validated mask result required; pilot is not scientific input')
  samples=cfg['samples'];n=len(samples);width=(n+3)//4
  if n!=14:raise ValueError('Expected14samples')
  if not 1000<=a.chunk_sites<=1000000:raise ValueError('Invalid chunk size')
  source={name:{e:a.masks/name/('source.'+e) for e in ['bed','bim','fam']} for name in NAMES}
  for name in NAMES:
   if sha(source[name]['bed'])!=cfg['output_BED_sha256'][name]:raise ValueError('Prepared BED hash changed')
   if sha(source[name]['bim'])!=cfg['source_sha256']['bim'] or sha(source[name]['fam'])!=cfg['source_sha256']['fam']:raise ValueError('Shared BIM/FAM changed')
  spec=importlib.util.spec_from_file_location('ab_existing_fst_core',a.fst_core);core=importlib.util.module_from_spec(spec);sys.modules[spec.name]=core;spec.loader.exec_module(core)
  populations={g:[i for i,s in enumerate(samples) if s.startswith(g+'-')] for g in ['AIC','HCO','IZU']};pairs=list(itertools.combinations(populations,2))
  if {g:len(v) for g,v in populations.items()}!={'AIC':4,'HCO':4,'IZU':6}:raise ValueError('Primary group mismatch')
  modes=['fixed_baseline_positions','refiltered','common_positions','common_defined_for_pair']
  accum={(name,mode,pair):core.WeightedFstAccumulator() for name in NAMES for mode in modes for pair in pairs}
  old_bim=a.original_main/'pruned.bim';original_lines=old_bim.read_text().splitlines()
  original_keys=[marker_key(line) for line in original_lines];oldids=set(original_keys)
  if len(original_keys)!=len(oldids) or len(original_keys)!=37566:raise ValueError('Expected 37,566 unique original PCA marker identities')
  original_id_seen=set();common_ids=[];sitecounts=collections.Counter()
  loss_stats=collections.defaultdict(lambda:np.zeros((4,n),dtype=np.int64))
  loss_sites=collections.Counter()
  def add_loss_scope(genotypes,scope,selection=None):
   count,values=loss_summary(genotypes,selection)
   loss_sites[scope]+=count
   for name,v in values.items():loss_stats[(name,scope)]+=v
  hs={name:source[name]['bed'].open('rb') for name in NAMES}
  for h in hs.values():
   if h.read(3)!=b'\x6c\x1b\x01':raise ValueError('Expected SNP-major BED')
  processed=0
  with source['baseline']['bim'].open() as bim:
   while True:
    b=hs['baseline'].read(width*a.chunk_sites)
    if not b:break
    allg={'baseline':decode(b,n)}
    for name in NAMES[1:]:
     other=hs[name].read(len(b))
     if len(other)!=len(b):raise ValueError('BED lengths differ')
     allg[name]=decode(other,n)
    size=len(allg['baseline']);lines=[bim.readline() for _ in range(size)]
    if any(not line for line in lines):raise ValueError('BIM shorter than BED')
    processed+=size
    base=allg['baseline'];rawstats=core.population_counts(base,range(n));eligible=core.mask_from_counts(rawstats.called,rawstats.allele_count,n,.1,1)
    sitecounts['baseline_MAC1_missing10']+=int(eligible.sum())
    for name,g in allg.items():
     core.validate_genotypes(g)
    validate_nested_masks(allg)
    add_loss_scope(allg,'all_source_positions')
    # No genotype change can rescue baseline-failing MAC/missingness sites.
    # Restrict FST work to existing baseline MAC1 set; PCA shared MAC2 is a subset.
    clipped={name:g[eligible] for name,g in allg.items()}
    add_loss_scope(clipped,'baseline_MAC1_positions')
    counts={name:{g:core.population_counts(x,cols) for g,cols in populations.items()} for name,x in clipped.items()}
    masks={};pca_common=np.ones(int(eligible.sum()),dtype=bool)
    for name,c in counts.items():
     called=sum(v.called for v in c.values());ac=sum(v.allele_count for v in c.values())
     masks[name]=core.mask_from_counts(called,ac,n,.1,1)
     pca_common &= core.mask_from_counts(called,ac,n,.1,2)
     sitecounts[name+'_refilter_MAC1']+=int(masks[name].sum())
    common=np.logical_and.reduce(list(masks.values()));sitecounts['common_MAC1_positions']+=int(common.sum())
    add_loss_scope(clipped,'common_MAC1_positions',common)
    common_full=np.zeros(size,dtype=bool);common_full[np.flatnonzero(eligible)]=pca_common
    for i,line in enumerate(lines):
     identifier=line.split()[1];key=marker_key(line)
     if key in oldids:
      if key in original_id_seen:raise ValueError('Duplicated original marker identity')
      original_id_seen.add(key)
      if common_full[i]:common_ids.append(identifier)
    for pair in pairs:
     comp={name:core.wc84_pair(c[pair[0]],c[pair[1]]) for name,c in counts.items()}
     common_defined=common&np.logical_and.reduce([np.isfinite(c.fst) for c in comp.values()])
     for name in NAMES:
      for mode,mask in [('fixed_baseline_positions',np.ones_like(common)),('refiltered',masks[name]),('common_positions',common),('common_defined_for_pair',common_defined)]:accum[(name,mode,pair)].update(comp[name],mask)
    if processed%1000000==0:print('FST processed',processed,'elapsed_s',round(time.monotonic()-start,1),flush=True)
   if bim.readline():raise ValueError('BIM longer than BED')
  for h in hs.values():
   if h.read(1):raise ValueError('Extra BED data')
   h.close()
  if processed!=20625692 or sitecounts['baseline_MAC1_missing10']!=5564311 or original_id_seen!=oldids:raise ValueError('Frozen baseline site validation failed')
  baseline_checks={}
  for pair in pairs:
   v=accum[('baseline','fixed_baseline_positions',pair)];expected,defined=EXPECTED_FST['_'.join(pair)]
   baseline_checks['_'.join(pair)]={'value':v.weighted_fst,'defined_sites':v.defined_sites,'matches':abs(v.weighted_fst-expected)<1e-10 and v.defined_sites==defined}
  if not all(x['matches'] for x in baseline_checks.values()):raise ValueError('Baseline FST no longer matches archived actual-VCFtools-validated result')
  rows=[]
  for (name,mode,pair),v in accum.items():
   if v.infinite_sites:raise ValueError('Infinite FST sites; inspect before interpretation')
   require_estimable_fst(v)
   rows.append({'mask':name,'mode':mode,'pair':'_'.join(pair),**v.as_dict()})
  with (o/'fst_summary.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
  with (o/'sample_mask_loss.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['mask','scope','sample','selected_positions','baseline_called','baseline_heterozygotes','retained_called','calls_set_missing','loss_fraction_of_baseline_called','loss_fraction_of_baseline_heterozygotes'])
   for (name,scope),values in sorted(loss_stats.items()):
    for i,s in enumerate(samples):
     called,hets,retained,lost=map(int,values[:,i])
     if retained+lost!=called:raise ValueError('Sample loss denominator mismatch')
     w.writerow([name,scope,s,loss_sites[scope],called,hets,retained,lost,lost/called if called else '',lost/hets if hets else ''])
  if len(common_ids)<14:raise ValueError('Too few common original PCA markers')
  commonfile=o/'common_original_PCA_markers.txt';commonfile.write_text(''.join(s+'\n' for s in common_ids))
  receipt['FST_baseline_validation']=baseline_checks;receipt['sitecounts']=dict(sitecounts);receipt['common_original_PCA_markers']=len(common_ids)
  receipt['nested_masks_verified_all_sites']=True
  plink=a.tools/'plink';receipt['plink_version']=receipt['preflight']['plink_version']
  basecmd=[str(plink),'--allow-extra-chr','--keep-allele-order','--threads','2','--memory','3500']
  def run(args,log):
   command=basecmd+list(map(str,args));receipt['commands'].append(command)
   with log.open('w') as outlog:subprocess.run(command,stdout=outlog,stderr=subprocess.STDOUT,check=True)
  pca=[]
  for mode in ['refiltered','common_original_markers']:
   for name in NAMES:
    d=o/'pca'/mode/name;d.mkdir(parents=True)
    prefix=a.masks/name/'source'
    if mode=='refiltered':
     run(['--bfile',prefix,'--geno','0.1','--mac','2','--make-bed','--out',d/'all_mac2'],d/'01_filter.log')
     run(['--bfile',d/'all_mac2','--bp-space','1000','--make-bed','--out',d/'spaced'],d/'02_space.log')
     run(['--bfile',d/'spaced','--indep-pairwise','50kb','1','0.2','--out',d/'ld'],d/'03_ld.log')
     run(['--bfile',d/'spaced','--extract',d/'ld.prune.in','--make-bed','--out',d/'pruned'],d/'04_extract.log')
    else:run(['--bfile',prefix,'--extract',commonfile,'--make-bed','--out',d/'pruned'],d/'04_extract.log')
    found=[x.split()[1] for x in (d/'pruned.fam').read_text().splitlines()]
    if found!=samples:raise ValueError('PCA input sample order drift')
    markerids=[line.split()[1] for line in (d/'pruned.bim').read_text().splitlines()]
    if mode=='common_original_markers' and markerids!=common_ids:raise ValueError('Common PCA marker set/order mismatch')
    if name=='baseline' and mode=='refiltered':
     receipt['baseline_PCA_input_identity']=validate_baseline_pca_inputs(d/'pruned',a.original_main/'pruned',samples)
    run(['--bfile',d/'pruned','--pca','14','header','tabs','--missing','--freqx','--out',d/'pca'],d/'05_pca.log')
    pca.append({'mask':name,'mode':mode,'markers':len(markerids),'BED_sha256':sha(d/'pruned.bed'),'BIM_sha256':sha(d/'pruned.bim'),'directory':str(d)})
  refs={mode:coordinates(o/'pca'/mode/'baseline/pca.eigenvec',samples) for mode in ['refiltered','common_original_markers']}
  ad_refs={mode:coordinates(o/'pca'/mode/'ad_control/pca.eigenvec',samples) for mode in ['refiltered','common_original_markers']}
  original=coordinates(a.original_main/'pca.eigenvec',samples)
  if metrics(refs['refiltered'],original)['maximum_pairwise_distance_difference']>=1e-5:raise ValueError('Archived PCA not reproduced')
  for row in pca:
   current=coordinates(Path(row['directory'])/'pca.eigenvec',samples)
   row['vs_same_mode_baseline']=metrics(current,refs[row['mode']])
   row['vs_same_mode_AD_control']=metrics(current,ad_refs[row['mode']])
  receipt['baseline_common_vs_refiltered']=metrics(refs['common_original_markers'],refs['refiltered'])
  (o/'pca_comparisons.json').write_text(json.dumps(pca,indent=2)+'\n')
  receipt['prepared_input_checks_after']={name:{'bed':sha(source[name]['bed'])==cfg['output_BED_sha256'][name],
   'bim':sha(source[name]['bim'])==cfg['source_sha256']['bim'],
   'fam':sha(source[name]['fam'])==cfg['source_sha256']['fam']} for name in NAMES}
  receipt['prepared_inputs_unchanged']=all(all(c.values()) for c in receipt['prepared_input_checks_after'].values())
  if not receipt['prepared_inputs_unchanged']:raise ValueError('Prepared input mutated')
  receipt['original_main_unchanged']=all(sha(a.original_main/name)==h for name,h in receipt['preflight']['original_main_sha256'].items())
  if not receipt['original_main_unchanged']:raise ValueError('Original PCA reference mutated during analysis')
  receipt.update(status='SENSITIVITY_COMPLETE',FST_rows=len(rows),PCA_panels=len(pca),limitations=['Masking-only test conditional on existing discovered SNPs and the diploid caller.','Fixed FST positions can include sites failing post-mask MAC/missingness; refiltered and common-position results are distinct.','Common-defined FST additionally requires each pair-specific site estimate defined in every mask.','Common PCA markers are a shared eligible subset of the ORIGINAL pruned set; LD is not recalculated in that comparison.','Sensitivity ranges are not confidence intervals; no biological ploidy or causal reference-bias inference.'])
 except Exception as e:receipt.update(status='FAILED',error=repr(e));raise
 finally:
  receipt['elapsed_seconds']=time.monotonic()-start;receipt['script_exit_code']=0 if receipt['status']=='SENSITIVITY_COMPLETE' else 1
  (o/'RUN_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'status':receipt['status'],'FST_rows':receipt['FST_rows'],'PCA_panels':receipt['PCA_panels'],'elapsed_seconds':receipt['elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
