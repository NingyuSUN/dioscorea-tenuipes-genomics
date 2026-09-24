#!/usr/bin/env python3
"""Actual-read prefix pilot; new outputs only, not full-run validation."""
import argparse,concurrent.futures,csv,datetime,gzip,hashlib,json,shlex,subprocess,urllib.request
from pathlib import Path
import numpy as np
from wgs_probe_core import LimitedReader,read_record,read_id

def run(cmd,log):
 with log.open('w') as f:
  f.write(shlex.join(list(map(str,cmd)))+'\n');f.flush();subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True)
def fetch(url,out,n,limit,expected):
 req=urllib.request.Request(url,headers={'Range':f'bytes=0-{limit-1}','User-Agent':'DT-research-read-probe/1.0'})
 hist=np.zeros(94,dtype=np.int64);base={chr(k):0 for k in b'ACGTN'};lengths={};namehash=hashlib.sha256();batch=[];first=last=None
 with urllib.request.urlopen(req,timeout=90) as response:
  headers={k:response.headers.get(k) for k in ['Content-Length','Content-Range','ETag','Last-Modified']};counter=LimitedReader(response,limit)
  with gzip.GzipFile(fileobj=counter) as source,gzip.open(out,'wb',compresslevel=1) as dest:
   for i in range(n):
    rec=read_record(source)
    if rec is None:raise ValueError('Too few records')
    h,s,p,q=rec;rid=read_id(h);namehash.update(rid+b'\n');first=first or rid.decode();last=rid.decode()
    dest.write(b'\n'.join(rec)+b'\n');lengths[len(s)]=lengths.get(len(s),0)+1
    for b in base:base[b]+=s.count(b.encode())
    batch.append(q)
    if len(batch)>=5000:
     hist+=np.bincount(np.frombuffer(b''.join(batch),dtype=np.uint8)-33,minlength=94);batch=[]
   if batch:hist+=np.bincount(np.frombuffer(b''.join(batch),dtype=np.uint8)-33,minlength=94)
 total=int(hist.sum());summary={'url':url,'records':n,'bases':total,'length_counts':lengths,'base_counts':base,'q20_fraction':float(hist[20:].sum()/total),'q30_fraction':float(hist[30:].sum()/total),'mean_phred':float(np.dot(np.arange(94),hist)/total),'phred_histogram':hist.tolist(),'read_id_sha256':namehash.hexdigest(),'first_id':first,'last_id':last,'compressed_bytes_read':counter.count,'compressed_prefix_sha256':counter.sha.hexdigest(),'response_headers':headers,'full_file_expected':expected,'full_file_MD5_verified':False,'prefix_fastq_sha256':hashlib.file_digest(out.open('rb'),'sha256').hexdigest()}
 out.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n');return summary

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--ena',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--plastome',type=Path,required=True);p.add_argument('--pairs',type=int,default=1000000);p.add_argument('--byte-limit-per-file',type=int,default=180000000);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 meta=list(csv.DictReader(a.ena.open(),delimiter='\t'));assert len(meta)==1;meta=meta[0];assert meta['run_accession']=='SRR33602354';urls=['https://'+u for u in meta['fastq_ftp'].split(';')];assert len(urls)==2
 inputhash={str(x):hashlib.file_digest(x.open('rb'),'sha256').hexdigest() for x in [a.reference,a.plastome]}
 refs=a.out/'reference';refs.mkdir();combined=refs/'combined.fa'
 with combined.open('wb') as f:
  with a.reference.open('rb') as src:
   for chunk in iter(lambda:src.read(1024*1024),b''):f.write(chunk)
  seq=''.join(x.strip() for x in a.plastome.read_text().splitlines() if not x.startswith('>'));assert len(seq)==154018
  f.write(b'\n>decoy_PV612491.1\n'+seq.encode()+b'\n')
 run(['samtools','faidx',combined],refs/'faidx.log')
 fastqs=[a.out/f'prefix_R{i}.fastq.gz' for i in [1,2]]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  index=pool.submit(run,['bwa','index','-p',refs/'combined',combined],refs/'bwa_index.log')
  jobs=[pool.submit(fetch,url,fq,a.pairs,a.byte_limit_per_file,{'md5':meta['fastq_md5'].split(';')[i],'bytes':int(meta['fastq_bytes'].split(';')[i])}) for i,(url,fq) in enumerate(zip(urls,fastqs))]
  results=[j.result() for j in jobs];assert results[0]['read_id_sha256']==results[1]['read_id_sha256']
  (a.out/'raw_prefix_qc.json').write_text(json.dumps({'selection':'first records; deterministic prefix, NOT random subsample','pairs':a.pairs,'pair_ids_match':True,'mates':results,'source_input_hashes':inputhash},indent=2)+'\n');print('PREFIX_QC_COMPLETE',[(x['records'],x['q30_fraction']) for x in results],flush=True);index.result()
 bam=a.out/'prefix.bam';cmd=['bwa','mem','-t','8','-R','@RG\\tID:SRR33602354_prefix\\tSM:SRR33602354\\tPL:DNBSEQ',refs/'combined',*fastqs]
 with (a.out/'bwa_mem.log').open('w') as log:
  log.write(shlex.join(list(map(str,cmd)))+'\n');log.flush();align=subprocess.Popen(list(map(str,cmd)),stdout=subprocess.PIPE,stderr=log)
  sort=subprocess.Popen(['samtools','sort','-@','2','-m','1G','-o',str(bam),'-'],stdin=align.stdout,stdout=log,stderr=log);align.stdout.close();rc2=sort.wait();rc1=align.wait();assert rc1==rc2==0,(rc1,rc2)
 run(['samtools','quickcheck','-v',bam],a.out/'quickcheck.log');run(['samtools','index',bam],a.out/'index.log')
 for args,name in [(['flagstat',bam],'flagstat.txt'),(['stats',bam],'samtools_stats.txt'),(['idxstats',bam],'idxstats.tsv'),(['coverage','--ff','0x904','-q','20','-Q','20',bam],'coverage_Q20_BQ20.tsv')]:run(['samtools',*args],a.out/name)
 counts={}
 for name,flags in [('primary',['-F','0x900']),('mapped_primary',['-F','0x904']),('mapped_Q20_primary',['-F','0x904','-q','20']),('proper_pair_primary',['-F','0x904','-f','2']),('unmapped_primary',['-F','0x900','-f','4']),('supplementary',['-f','0x800'])]:counts[name]=int(subprocess.check_output(['samtools','view','-c',*flags,str(bam)],text=True))
 assert counts['primary']==2*a.pairs
 (a.out/'alignment_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
 after={str(x):hashlib.file_digest(x.open('rb'),'sha256').hexdigest() for x in [a.reference,a.plastome]};assert after==inputhash
 (a.out/'COMPLETE').write_text('Real prefix reads, read IDs and BAM structure checked; source references unchanged. Not a full-run or ploidy validation.\n');print('COMPLETE',counts,flush=True)
if __name__=='__main__':main()
