#!/usr/bin/env python3
"""Paired, matched-ID diagnostic of raw vs trimmed primary alignments."""
import argparse,gzip,json,re,subprocess,shlex,hashlib
from pathlib import Path
from wgs_probe_core import read_record,read_id
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();d=a.out;t=d/'trim_diagnostic';t.mkdir(exist_ok=False)
names=set();fq=[t/f'clean_R{i}.fastq.gz' for i in [1,2]]
with gzip.open(d/'clean_prefix_R1.fastq.gz','rb') as r1,gzip.open(d/'clean_prefix_R2.fastq.gz','rb') as r2,gzip.open(fq[0],'wb') as w1,gzip.open(fq[1],'wb') as w2:
 for _ in range(100000):
  x,y=read_record(r1),read_record(r2);assert x and y;assert read_id(x[0])==read_id(y[0]);names.add(read_id(x[0]).decode())
  for z,w in [(x,w1),(y,w2)]:w.write(b'\n'.join(z)+b'\n')
assert len(names)==100000
(t/'read_names.txt').write_text('\n'.join(sorted(names))+'\n')
cmd=['bwa','mem','-t','8',str(d/'reference/combined'),*map(str,fq)]
with (t/'bwa_mem.log').open('w') as log:
 log.write(shlex.join(cmd)+'\n');log.flush();pr=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log)
 cv=subprocess.run(['samtools','view','-b','-o',str(t/'clean.bam'),'-'],stdin=pr.stdout,stderr=log);pr.stdout.close();assert cv.returncode==0 and pr.wait()==0
subprocess.run(['samtools','quickcheck',str(t/'clean.bam')],check=True)
def summary(bam,subset):
 r={'primary_reads':0,'mapped_reads':0,'nuclear_reads':0,'nuclear_Q20_reads':0,'nuclear_query_bases':0,'nuclear_primary_softclip_bases':0,'nuclear_aligned_query_bases':0,'nuclear_NM':0}
 pr=subprocess.Popen(['samtools','view','-F','0x900',str(bam)],stdout=subprocess.PIPE,text=True)
 for line in pr.stdout:
  f=line.split('\t')
  if subset and f[0] not in names:continue
  r['primary_reads']+=1
  if int(f[1])&4:continue
  r['mapped_reads']+=1
  if f[2]=='decoy_PV612491.1':continue
  r['nuclear_reads']+=1;r['nuclear_Q20_reads']+=int(int(f[4])>=20);r['nuclear_query_bases']+=len(f[9])
  for size,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
   if op=='S':r['nuclear_primary_softclip_bases']+=int(size)
   elif op in 'MI=X':r['nuclear_aligned_query_bases']+=int(size)
  nm=[v for v in f[11:] if v.startswith('NM:i:')];assert len(nm)==1;r['nuclear_NM']+=int(nm[0][5:])
 assert pr.wait()==0 and r['primary_reads']==200000
 r['nuclear_primary_softclip_fraction']=r['nuclear_primary_softclip_bases']/r['nuclear_query_bases']
 return r
res={'scope':'first100000 retained pairs; SAME names compared, NOT random','selected_ids_sha256':hashlib.sha256((t/'read_names.txt').read_bytes()).hexdigest(),'raw':summary(d/'prefix.bam',True),'trimmed':summary(t/'clean.bam',False),'limitations':['primary soft clips may also have supplementary alignments','cannot identify cause of mismatch or sample taxonomy','not a full-run or joint-call validation']}
(t/'trim_sensitivity.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
