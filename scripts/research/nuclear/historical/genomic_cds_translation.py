#!/usr/bin/env python3
"""One-pass genomic translation at pre-existing miniprot GFF CDS coordinates.
No search, alignment, rescue, annotation editing, or genome download.
"""
from pathlib import Path
import json,hashlib,collections,resource,time,datetime,os,traceback
D=Path(__file__).resolve().parent;I=D.parent
GENOME=Path('workspace/wgs/genome.fa')
EXPECTED='36fecc3e8b3054728d9779d22d61957062fd914bf24cd04f09e27d363d97c87b'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(n,v):(D/n).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def rc(b):return b.translate(bytes.maketrans(b'ACGTNacgtn',b'TGCANtgcan'))[::-1]
def main():
 assert not (D/'STATUS.json').exists(),'Existing attempt; never overwrite/retry'
 started=time.monotonic();status={'status':'RUNNING_COORDINATE_ONLY_CDS_TRANSLATION','started_utc':now(),'pid':os.getpid(),'single_thread_python':True}
 save('STATUS.json',status)
 try:
  grp=next(x.split(':',2)[2] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
  gp=Path('/sys/fs/cgroup')/grp.lstrip('/');limits={k:((gp/k).read_text().strip() if (gp/k).exists() else None) for k in ['memory.max','memory.swap.max','cpu.max','cpuset.cpus.effective']}
  assert int(limits['memory.max'])<=1073741824 and limits['memory.swap.max']=='0'
  os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
  assert len(os.sched_getaffinity(0))==1
  if limits['cpu.max'] is not None:
   quota,period=map(int,limits['cpu.max'].split());assert quota<=period
  save('CGROUP_OBSERVED.json',{'path':str(gp),'limits':limits,'process_affinity':sorted(os.sched_getaffinity(0))})
  gff=I/'alignments.gff';queryfile=I/'historical_173.faa';status0=json.loads((I/'STATUS.json').read_text());manifest=json.loads((I/'QUERY_MANIFEST.json').read_text())
  assert status0['exit_code']==0 and sha(gff)==status0['output_sha256'] and sha(queryfile)==manifest['query_sha256']
  queries={}
  for block in queryfile.read_text().split('>')[1:]:
   a=block.splitlines();queries[a[0].split()[0]]=''.join(a[1:])
  assert len(queries)==173 and sum(map(len,queries.values()))==122636
  alignments=[];cur=None
  for line in gff.read_text().splitlines():
   if line.startswith('##PAF\t'):
    if cur:alignments.append(cur)
    a=line.split('\t')[1:];cur={'query':a[0],'query_length':int(a[1]),'query_start0':int(a[2]),'query_end0':int(a[3]),'strand':a[4],'contig':a[5],'contig_length':int(a[6]),'features':[]}
   elif line.startswith('##STA\t'):cur['STA']=line.split('\t',1)[1]
   elif line and not line.startswith('#'):
    a=line.split('\t');at=dict(x.split('=',1) for x in a[8].split(';') if x)
    f={'type':a[2],'start':int(a[3]),'end':int(a[4]),'strand':a[6],'phase':a[7],'attributes':at};cur['features'].append(f)
    if f['type']=='mRNA':cur.update(alignment_id=at['ID'],rank=int(at['Rank']))
  if cur:alignments.append(cur)
  requests=collections.defaultdict(list);pieces={};coords={}
  def request(contig,start,end,key):
   assert 0<=start<end;assert key not in pieces
   requests[contig].append((start,end,key));pieces[key]=[];coords[key]=(start,end)
  for a in alignments:
   a['CDS']=sorted([f for f in a['features'] if f['type']=='CDS'],key=lambda f:f['start'],reverse=a['strand']=='-')
   for i,f in enumerate(a['CDS']):request(a['contig'],f['start']-1,f['end'],a['alignment_id']+':cds:'+str(i))
   a['junctions']=[]
   for i,(x,y) in enumerate(zip(a['CDS'],a['CDS'][1:])):
    start,end=(x['end'],y['start']-1) if a['strand']=='+' else (y['end'],x['start']-1)
    assert end-start>=4
    request(a['contig'],start,start+2,a['alignment_id']+':spliceL:'+str(i));request(a['contig'],end-2,end,a['alignment_id']+':spliceR:'+str(i))
    a['junctions'].append({'genomic_start0':start,'genomic_end0':end,'length_nt':end-start})
  for contig in requests:requests[contig].sort()
  h=hashlib.sha256();contig=None;pos=0;todo=[];active=[];next_req=0;seen=set();contig_lengths={};nt=0
  with GENOME.open('rb') as handle:
   for raw in handle:
    h.update(raw)
    if raw.startswith(b'>'):
     if contig is not None:contig_lengths[contig]=pos
     contig=raw[1:].split()[0].decode();assert contig not in seen;seen.add(contig);pos=0;todo=requests.get(contig,[]);active=[];next_req=0
     continue
    line=raw.strip();end=pos+len(line);nt+=len(line)
    if todo:
     while next_req<len(todo) and todo[next_req][0]<end:active.append(todo[next_req]);next_req+=1
     active=[z for z in active if z[1]>pos]
     for start,stop,key in active:
      b=line[max(start,pos)-pos:min(stop,end)-pos]
      if b:pieces[key].append(b)
    pos=end
   if contig is not None:contig_lengths[contig]=pos
  assert h.hexdigest()==EXPECTED
  assert set(requests)<=seen
  pieces={k:b''.join(v).upper() for k,v in pieces.items()}
  for key,b in pieces.items():assert len(b)==coords[key][1]-coords[key][0]
  bases='TCAG';aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
  codons=[a+b+c for a in bases for b in bases for c in bases];code=dict(zip(codons,aas));assert len(code)==64
  results=[]
  for a in alignments:
   aid=a['alignment_id'];oriented=[];phases=[];cum=0
   first_phase=int(a['CDS'][0]['phase']);assert first_phase in [0,1,2]
   for i,f in enumerate(a['CDS']):
    b=pieces[aid+':cds:'+str(i)];b=rc(b) if a['strand']=='-' else b
    expected_phase=first_phase if i==0 else (3-cum%3)%3
    phases.append({'CDS_index_transcription_order':i+1,'observed_phase':int(f['phase']),'expected_phase':expected_phase,'matches':int(f['phase'])==expected_phase})
    oriented.append(b);cum+=len(b)-(first_phase if i==0 else 0)
   dna=b''.join(oriented)[first_phase:];prot=''.join(code.get(dna[i:i+3].decode(),'X') for i in range(0,len(dna)-2,3));core=prot[:-1] if prot.endswith('*') else prot
   q=queries[a['query']];full=a['query_start0']==0 and a['query_end0']==len(q);sta=a['STA'];stacore=sta[:-1] if sta.endswith('*') else sta
   junctions=[]
   for i,j in enumerate(a['junctions']):
    left=pieces[aid+':spliceL:'+str(i)];right=pieces[aid+':spliceR:'+str(i)]
    donor,acceptor=(left,right) if a['strand']=='+' else (rc(right),rc(left))
    junctions.append({**j,'donor':donor.decode(),'acceptor':acceptor.decode()})
   normal=not(set(dna)-set(b'ACGT'));phaseok=all(x['matches'] for x in phases)
   exact=full and core==q and len(dna)%3==0 and normal and phaseok
   results.append({'alignment_id':aid,'query':a['query'],'rank':a['rank'],'contig':a['contig'],'strand':a['strand'],'CDS_segments':len(a['CDS']),'CDS_nt_count_after_first_phase':len(dna),'first_phase_bases_removed_once':first_phase,'trailing_nt_not_translated':len(dna)%3,'phase_consistency_all':phaseok,'phase_checks':phases,'non_ACGT_CDS_nt_count':sum(x not in b'ACGT' for x in dna),'full_query_span':full,'translated_length_including_terminal_stop':len(prot),'translated_terminal_stop':prot.endswith('*'),'translated_internal_stop_count':core.count('*'),'translated_X_count':prot.count('X'),'translated_protein':prot,'translation_core_sha256':hashlib.sha256(core.encode()).hexdigest(),'CDS_nt_sha256':hashlib.sha256(dna).hexdigest(),'query_sha256':hashlib.sha256(q.encode()).hexdigest(),'translation_core_equals_full_query':core==q,'unambiguous_full_length_exact':exact,'translation_core_equals_STA_core':core==stacore,'splice_junctions':junctions,'non_GT_AG_junctions':sum(x['donor']+x['acceptor']!='GTAG' for x in junctions)})
  primary=[r for r in results if r['rank']==1]
  summary={'status':'COMPLETE_COORDINATE_ONLY_GENOMIC_TRANSLATION','genome_sha256':h.hexdigest(),'genome_total_nt':nt,'genome_contigs':len(seen),'alignments':len(results),'primary_alignments':len(primary),'primary_unambiguous_full_length_exact':sum(r['unambiguous_full_length_exact'] for r in primary),'all_unambiguous_full_length_exact':sum(r['unambiguous_full_length_exact'] for r in results),'primary_translation_equals_STA_core':sum(r['translation_core_equals_STA_core'] for r in primary),'all_translation_equals_STA_core':sum(r['translation_core_equals_STA_core'] for r in results),'primary_phase_inconsistent':[r['query'] for r in primary if not r['phase_consistency_all']],'primary_trailing_nt_queries':[r['query'] for r in primary if r['trailing_nt_not_translated']],'primary_internal_stop_queries':[r['query'] for r in primary if r['translated_internal_stop_count']],'primary_non_ACGT_queries':[r['query'] for r in primary if r['non_ACGT_CDS_nt_count']],'source_sha256':{str(gff):sha(gff),str(queryfile):sha(queryfile),str(D/'genomic_cds_translation.py'):sha(D/'genomic_cds_translation.py')},'limits':'Genomic extraction and standard-code translation are independent of miniprot STA. Coordinates/splice structure remain selected by miniprot, so this is not an independent alignment or an authentication of sample identity.'}
  save('CDS_TRANSLATION.json',{'summary':summary,'records':results});status.update(status='COMPLETE_EXIT0',summary=summary,output_sha256=sha(D/'CDS_TRANSLATION.json'),output_bytes=(D/'CDS_TRANSLATION.json').stat().st_size)
 except BaseException:
  status.update(status='FAILED_NO_AUTOMATIC_RETRY',error=traceback.format_exc());raise
 finally:
  status.update(finished_utc=now(),elapsed_seconds=time.monotonic()-started,peak_self_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);save('STATUS.json',status)
 print(json.dumps({k:v for k,v in status.get('summary',{}).items() if k!='source_sha256'}))
if __name__=='__main__':main()
