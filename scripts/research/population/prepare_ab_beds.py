#!/usr/bin/env python3
"""Prepare three additional masks on the frozen DP/GQ-masked own14 BED.
Read original AD once. Verify every incoming locus and baseline genotype.
A pilot times a prefix only and MUST NOT be used for biological analysis.
"""
import argparse,collections,csv,hashlib,json,re,subprocess,time,sys
from pathlib import Path

POLICIES={'ad_control':None,'ab20_80':(1,4,5),'ab30_70':(3,7,10)}
HISTORICAL_VCF_SHA='43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9'


def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()


def fields(raw):
 gt,dp,gq,ad=raw.split(':')
 if '.' in gt or dp=='.' or gq=='.':return -1,ad
 alleles=re.split(r'[/|]',gt)
 if len(alleles)!=2 or any(v not in ('0','1') for v in alleles):raise ValueError('Unexpected non-diploid/biallelic GT')
 if not (5<=int(dp)<=100 and int(gq)>=20):return -1,ad
 return sum(map(int,alleles)),ad


def bed_code(alt_dosage,a1,a2,ref,alt):
 if {a1,a2}!={ref,alt}:raise ValueError('BED and VCF alleles differ')
 if alt_dosage==-1:return 1
 if alt_dosage==1:return 2
 allele=ref if alt_dosage==0 else alt
 return 0 if allele==a1 else 3


def additional_masks(ad):
 values=ad.split(',')
 if len(values)!=2 or any(v=='.' for v in values):return dict.fromkeys(POLICIES,'AD_missing_or_wrong_count')
 try:ref,alt=map(int,values)
 except ValueError:return dict.fromkeys(POLICIES,'AD_invalid')
 if ref<0 or alt<0:return dict.fromkeys(POLICIES,'AD_invalid')
 total=ref+alt
 if total<10:return dict.fromkeys(POLICIES,'ADsum_below_10')
 result={'ad_control':None}
 for name,(low,high,scale) in list(POLICIES.items())[1:]:
  result[name]=None if low*total<=scale*alt<=high*total else 'outside_symmetric_AB_interval'
 return result


def mask_record(packed, calls, alleles):
 ref,alt,a1,a2=alleles
 result={k:bytearray(packed) for k in POLICIES};events=[]
 for i,raw in enumerate(calls):
  dosage,ad=fields(raw);code=(packed[i//4]>>(2*(i%4)))&3
  if code!=bed_code(dosage,a1,a2,ref,alt):raise ValueError('Baseline GT differs from raw DP5-100/GQ20 rule, sample index '+str(i))
  if code!=2:continue
  for policy,reason in additional_masks(ad).items():
   if reason:
    shift=2*(i%4);result[policy][i//4]=(result[policy][i//4]&~(3<<shift))|(1<<shift)
   events.append((policy,i,reason))
 return result,events


def index_bim(path,width):
 """Only one entry per contiguous contig, not a dictionary of all variants."""
 index={};count=0;current=None;previous=0;offset=0
 with Path(path).open('rb') as f:
  while True:
   line=f.readline()
   if not line:break
   fields=line.decode().split()
   if len(fields)!=6:raise ValueError('Malformed BIM')
   chrom=fields[0];pos=int(fields[3])
   if chrom!=current:
    if chrom in index:raise ValueError('BIM contig is not contiguous: '+chrom)
    index[chrom]={'bim_byte_offset':offset,'bed_byte_offset':3+count*width,'variant_offset':count,'count':0}
    current=chrom;previous=0
   if pos<=previous:raise ValueError('BIM positions not strictly increasing: '+chrom)
   previous=pos;index[chrom]['count']+=1;count+=1;offset+=len(line)
 return index,count


class ContigReader:
 """Bind each queried locus to its original BIM/BED segment despite contig order."""
 def __init__(self,bed,bim,index,width):
  self.bed,self.bim,self.index,self.width=bed,bim,index,width
  self.current=None;self.used=0;self.seen=set()
 def check_segment(self):
  if self.current is not None and self.used!=self.index[self.current]['count']:
   raise ValueError('Query did not cover the complete original contig segment: '+self.current)
 def read(self,chrom,pos):
  changed=chrom!=self.current
  if changed:
   self.check_segment()
   if chrom not in self.index or chrom in self.seen:raise ValueError('Unknown or repeated query contig: '+chrom)
   self.current=chrom;self.used=0;self.seen.add(chrom)
   self.bed.seek(self.index[chrom]['bed_byte_offset']);self.bim.seek(self.index[chrom]['bim_byte_offset'])
  line=self.bim.readline().decode();fields=line.split();packed=self.bed.read(self.width)
  if len(fields)!=6 or len(packed)!=self.width or (fields[0],fields[3])!=(chrom,pos):
   raise ValueError('BIM/VCF locus differs within bound contig: '+chrom+':'+pos)
  self.used+=1
  if self.used>self.index[chrom]['count']:raise ValueError('Query exceeds BIM contig length')
  return line,fields,packed,changed
 def finish(self):
  self.check_segment()
  if self.seen!=set(self.index):raise ValueError('Query omitted some BIM contigs')


def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--vcf',type=Path,default=Path('workspace/wgs/population_audit_20260914/inputs/Dtenuipes_population.vcf.gz'))
 p.add_argument('--baseline-prefix',type=Path,default=Path('workspace/wgs/population_influence_20260915/source_own14'))
 p.add_argument('--out',type=Path,required=True);p.add_argument('--bcftools',default='bcftools')
 p.add_argument('--pilot-sites',type=int,default=0,help='First N PASS source sites, timing/engineering only; 0 means full analysis')
 a=p.parse_args();o=a.out.resolve();o.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 source={e:Path(str(a.baseline_prefix)+'.'+e) for e in ['bed','bim','fam']}
 receipt={'status':'RUNNING','scope':'PILOT_ONLY' if a.pilot_sites else 'FULL_SELECTED_SNP_SOURCE','source_paths':{k:str(v) for k,v in source.items()},'source_VCF':str(a.vcf),'historical_source_VCF_sha256':HISTORICAL_VCF_SHA,'VCF_full_hash_recomputed':False,'script_sha256':sha(__file__),'policies':POLICIES,'commands':[]}
 proc=None;handles={};n=0
 try:
  if a.pilot_sites<0:raise ValueError('Negative pilot size')
  families=[line.split() for line in source['fam'].read_text().splitlines()];samples=[x[1] for x in families]
  if len(samples)!=14 or len(set(samples))!=14 or collections.Counter(s[:3] for s in samples)!={'AIC':4,'HCO':4,'IZU':6}:raise ValueError('Expected primary14 FAM')
  if {'AIC-F20-S13','AIC-M6-S11'}&set(samples):raise ValueError('Excluded sample in baseline')
  samplefile=o/'samples.txt';samplefile.write_text(''.join(s+'\n' for s in samples));receipt['samples']=samples
  receipt['source_sha256']={k:sha(v) for k,v in source.items()}
  width=(len(samples)+3)//4
  if source['bed'].stat().st_size!=3+width*20625692:raise ValueError('Expected frozen 20,625,692-site BED')
  contig_index,total=index_bim(source['bim'],width)
  if total!=20625692:raise ValueError('Expected frozen BIM variant count')
  (o/'contig_offsets.json').write_text(json.dumps(contig_index,indent=2)+'\n')
  receipt['contig_offset_sha256']=sha(o/'contig_offsets.json')
  receipt['contig_segments']=len(contig_index)
  header_cmd=[a.bcftools,'view','--no-version','-h',str(a.vcf)];receipt['commands'].append(header_cmd)
  header=subprocess.check_output(header_cmd);receipt['header_sha256']=hashlib.sha256(header).hexdigest()
  if receipt['header_sha256']!='f9a99bc5d89c3b69463fd459f015dafa94d3d54284a592e166552651998566f9':raise ValueError('Source header drift')
  for name in ['baseline',*POLICIES]:
   d=o/name;d.mkdir()
   (d/'source.fam').symlink_to(source['fam'].resolve())
   if not a.pilot_sites:(d/'source.bim').symlink_to(source['bim'].resolve())
   if name=='baseline' and not a.pilot_sites:(d/'source.bed').symlink_to(source['bed'].resolve())
   else:
    handles[name]=(d/'source.bed').open('w+b');handles[name].write(b'\x6c\x1b\x01')
    if not a.pilot_sites:handles[name].truncate(3+width*total)
  pilot_bim=(o/'pilot.bim').open('w') if a.pilot_sites else None
  counts=collections.Counter()
  command=[a.bcftools,'query','-S',str(samplefile),'-i','FILTER="PASS"','-f','%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP:%GQ:%AD]\n',str(a.vcf)]
  receipt['commands'].append(command)
  receipt['input_setup_seconds']=time.monotonic()-start
  stream_start=time.monotonic()
  with source['bed'].open('rb') as bed,source['bim'].open('rb') as bim,(o/'query.stderr').open('w') as err:
   if bed.read(3)!=b'\x6c\x1b\x01':raise ValueError('Expected SNP-major BED')
   reader=ContigReader(bed,bim,contig_index,width)
   proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=err,text=True)
   for line in proc.stdout:
    columns=line.rstrip('\n').split('\t')
    if len(columns)!=4+len(samples):raise ValueError('Input query shape mismatch')
    chrom,pos,ref,alt=columns[:4]
    locus,site,packed,new_contig=reader.read(chrom,pos)
    if len(ref)!=1 or len(alt)!=1:raise ValueError('Source VCF is not biallelic SNP')
    if new_contig and not a.pilot_sites:
     for handle in handles.values():handle.seek(contig_index[chrom]['bed_byte_offset'])
    try:masked,events=mask_record(packed,columns[4:],(ref,alt,site[4],site[5]))
    except ValueError as e:raise ValueError(f'{chrom}:{pos}: {e}') from e
    for name,data in masked.items():handles[name].write(data)
    for policy,i,reason in events:
     counts[(policy,samples[i],'baseline_eligible_heterozygotes')]+=1
     counts[(policy,samples[i],'retained' if reason is None else reason)]+=1
    n+=1
    if a.pilot_sites:
     handles['baseline'].write(packed);pilot_bim.write(locus)
    if n%100000==0:print(f'validated_sites={n}; elapsed_s={time.monotonic()-start:.1f}',flush=True)
    if a.pilot_sites and n>=a.pilot_sites:
     proc.stdout.close();proc.terminate();break
   receipt['query_exit_code']=proc.wait()
   if not a.pilot_sites:
    if receipt['query_exit_code']:raise RuntimeError('VCF query failed')
    reader.finish()
  receipt['stream_seconds']=time.monotonic()-stream_start
  for h in handles.values():h.close()
  if pilot_bim:
   pilot_bim.close()
   for name in ['baseline',*POLICIES]:(o/name/'source.bim').symlink_to((o/'pilot.bim').resolve())
  if not a.pilot_sites and n!=20625692:raise ValueError('Expected all 20,625,692 source sites')
  with (o/'mask_counts.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['policy','sample','category','calls']);w.writerows([(*key,value) for key,value in sorted(counts.items())])
  receipt['source_unchanged']={k:sha(v)==receipt['source_sha256'][k] for k,v in source.items()}
  if not all(receipt['source_unchanged'].values()):raise RuntimeError('Baseline input changed')
  receipt.update(status='PILOT_COMPLETE_NOT_BIOLOGICAL' if a.pilot_sites else 'MASK_COMPLETE',sites=n,baseline_calls_verified=n*len(samples),output_BED_sha256={name:sha(o/name/'source.bed') for name in ['baseline',*POLICIES]})
 except Exception as e:
  receipt.update(status='FAILED',error=repr(e),sites=n)
  if proc is not None and proc.poll() is None:proc.terminate();proc.wait()
  raise
 finally:
  for h in handles.values():
   if not h.closed:h.close()
  receipt['elapsed_seconds']=time.monotonic()-start
  receipt['script_exit_code']=0 if receipt['status'] in ('MASK_COMPLETE','PILOT_COMPLETE_NOT_BIOLOGICAL') else 1
  (o/'MASK_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
