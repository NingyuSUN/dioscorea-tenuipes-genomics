"""Descriptive same-material comparison of published genotype STATES and CP SNPs."""
from pathlib import Path
import collections,datetime,gzip,hashlib,json,resource,time
import numpy as np
R=Path(__file__).resolve().parents[1];I=R/'inputs/published';O=R/'published_comparison';assert not O.exists();O.mkdir()
start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
names=['dt_kita1','dt_waka1','dt_mzawa','ed_hachi1','tn_yuga','tn_utsu','tn_mie1','dq_mie'];N=1054866
path=I/'Fig_4_nighbor-net__vcf_pp_kita_waka_mzawa_edo_yuga_utsu_mie_quinq.no_missing.m2M2vsnps.vcf.gz.gt.gz';A=np.empty((N,8),dtype=np.int8);chrom=collections.Counter();coord=hashlib.sha256();last={}
with gzip.open(path,'rt') as f:
 for i,line in enumerate(f):
  a=line.split();assert len(a)==10 and i<N;g=list(map(int,a[2:]));assert all(x in [0,1,2] for x in g);A[i]=g;chrom[a[0]]+=1;position=int(a[1]);assert position>last.get(a[0],0);last[a[0]]=position;coord.update(('\t'.join(a[:2])+'\n').encode())
assert i+1==N
assert not any(c.upper() in ['CP','MT','CHLOROPLAST','MITOCHONDRIA'] for c in chrom),chrom
pair=[]
for i in range(8):
 for j in range(i+1,8):
  d=np.abs(A[:,i]-A[:,j]);pair.append({'a':names[i],'b':names[j],'columns':N,'identical_encoded_genotypes':int((d==0).sum()),'different_encoded_genotypes':int((d>0).sum()),'encoded_state_mismatch_fraction':float((d>0).mean()),'encoded_halfstep_distance':float(d.mean()/2)})
tok=np.all(A[:,:3]==A[:,[0]],axis=1)&np.isin(A[:,0],[0,2]);ten=np.all(A[:,4:7]==A[:,[4]],axis=1)&np.isin(A[:,4],[0,2]);diag=tok&ten&(A[:,0]!=A[:,4]);h=A[diag,3];t=A[diag,0];u=A[diag,4]
D={'columns':int(diag.sum()),'hachi_same_tokoro_homozygous_code':int((h==t).sum()),'hachi_heterozygous_code':int((h==1).sum()),'hachi_same_tenuipes_homozygous_code':int((h==u).sum()),'criterion':'all3tokoro and all3tenuipes each share opposite homozygous encoded calls; selection based on sixcomparators, notHachi; no recomputedDP/AD/GQ filter','scope':'descriptive reanalysis of published SNP-selected calls, not new independent hybrid-origin validation or ancestry-fraction estimate'}
assert sum(D[k] for k in ['hachi_same_tokoro_homozygous_code','hachi_heterozygous_code','hachi_same_tenuipes_homozygous_code'])==D['columns']
cp={};key=None
for line in (I/'Fig_6_organellar_trees__62_fasta_files__SNPs_CP.fasta').read_text().splitlines():
 if line.startswith('>'):key=names[int(line[1:].strip())];cp[key]=''
 else:cp[key]+=line.strip()
assert list(cp)==names and all(len(x)==829 and set(x)<=set('ACGT') for x in cp.values())
CP=[{'a':names[i],'b':names[j],'selected_SNP_columns':829,'mismatches':sum(x!=y for x,y in zip(cp[names[i]],cp[names[j]]))} for i in range(8) for j in range(i+1,8)]
# Independent aggregate check via integer contingency tables, not vectorized subtraction.
checks=[]
for row in pair:
 i=names.index(row['a']);j=names.index(row['b']);ct=np.bincount(A[:,i]*3+A[:,j],minlength=9).reshape(3,3);mismatch=int(ct.sum()-np.trace(ct));step=sum(abs(x-y)*int(ct[x,y]) for x in range(3) for y in range(3))/(2*N);assert mismatch==row['different_encoded_genotypes'] and abs(step-row['encoded_halfstep_distance'])<1e-12;checks.append({'a':row['a'],'b':row['b'],'contingency_verified':True})
O.joinpath('NUCLEAR_PAIRWISE.json').write_text(json.dumps(pair,indent=2)+'\n');O.joinpath('CP_PAIRWISE.json').write_text(json.dumps(CP,indent=2)+'\n');O.joinpath('DIAGNOSTIC_STATE_COUNTS.json').write_text(json.dumps(D,indent=2)+'\n')
fields=list(pair[0]);O.joinpath('NUCLEAR_PAIRWISE.tsv').write_text('\t'.join(fields)+'\n'+''.join('\t'.join(str(r[k]) for k in fields)+'\n' for r in pair))
fields=list(CP[0]);O.joinpath('CP_PAIRWISE.tsv').write_text('\t'.join(fields)+'\n'+''.join('\t'.join(str(r[k]) for k in fields)+'\n' for r in CP))
receipt={'status':'COMPLETE_DESCRIPTIVE_REANALYSIS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'nuclear_columns':N,'cp_columns':829,'samples':names,'chromosomes':dict(chrom),'coordinate_sha256':coord.hexdigest(),'GT_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'independent_contingency_checks':checks,'elapsed_seconds':time.monotonic()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limitations':['Original source calls not revalidated; no nuclearDP/AD/GQ','IUPAC and0/1/2 reflect author encoding, not reliable genotype dosage for everyploidy','Distance is descriptive in SNP-selected columns, not substitutions pergenome site or admixturefraction','CP829 locus is linked; nuclear/CP distances have different units and ascertainment','hime2 is absent from published8panel; results do not locate hime2 nuclear ancestry','No tree or mechanistic likelihood test performed']};O.joinpath('RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'hachi_nuclear':[r for r in pair if 'ed_hachi1' in [r['a'],r['b']]],'diagnostic':D,'elapsed_seconds':receipt['elapsed_seconds']},indent=2))
