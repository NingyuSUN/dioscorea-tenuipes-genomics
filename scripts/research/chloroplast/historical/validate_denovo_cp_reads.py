#!/usr/bin/env python3
"""Read-back validation of CP-like de novo contigs, before any tree placement.
Coordinates remain each original contig; no reference sequence supplies a base.
"""
from pathlib import Path
import datetime,hashlib,json,os,resource,statistics,subprocess,time
import recover_hachi_cp_callable as call
R=Path('workspace/wgs/dt_research_p57_20260930');O=R/'denovo_readback'
CONTIG=R/'denovo_1m/assembly/final.contigs.fa'
M='workspace/miniconda3/envs/ragtag/bin/minimap2';S='workspace/miniconda3/envs/busco/bin/samtools'
def sha(p):return call.sha256_file(p)
def save(n,x): (O/n).write_text(json.dumps(x,indent=2)+'\n')
def fa(p):
 d={};k=None
 for l in Path(p).read_text().splitlines():
  if l.startswith('>'):k=l[1:].split()[0];d[k]=''
  else:d[k]+=l.strip()
 return d
def limit():resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
def run():
 assert not O.exists();O.mkdir()
 start=time.monotonic();state={'status':'RUNNING','stage':'PREFLIGHT','started_utc':call.utcnow(),'pid':os.getpid(),'proc_start_ticks':Path('/proc/self/stat').read_text().split()[21]};save('STATUS.json',state)
 try:
  assembly=json.loads((R/'denovo_1m/STATUS.json').read_text());assert assembly['status']=='COMPLETE_REVIEW_PENDING' and sha(CONTIG)==assembly['contigs_sha256']
  screen=json.loads((R/'denovo_cp_audit/SCREEN.json').read_text());assert screen['source_contigs_sha256']==sha(CONTIG)
  ids=list(dict.fromkeys(x['contig'] for x in screen['alignments'] if x['qualifying_cp_like']))
  refs=fa(CONTIG);assert len(ids)==3 and all(x in refs for x in ids)
  reads=[R/f'pilot_1m_qualityonly/clean_R{i}.fastq.gz' for i in [1,2]]
  inputs={str(p):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [CONTIG,*reads]}
  save('FROZEN_PARAMETERS.json',{'created_before_readback':True,'source_input_hashes':inputs,'CP_like_contigs':ids,'all_assembly_contigs_in_mapping':len(refs),'base_calling_thresholds':call.THRESHOLDS,'additional_filters':{'contig_end_mask_each_end':200,'min_fraction_of_contig_median_depth':0.25,'max_fraction_of_contig_median_depth':4.0,'median_scope':'all internal positions with ACGT template depth>0, after BQ/MAPQ/template filters','assembly_disagreement':'mask N, retain observed raw dominant and count'}})
  bam=O/'all_contigs.bam';cmd=[M,'-ax','sr','-t','4',str(CONTIG),*map(str,reads)];sort=[S,'sort','-@','2','-m','256M','-o',str(bam),'-']
  save('MAPPING_COMMANDS.json',{'minimap2':cmd,'sort':sort,'limit_per_process_bytes':4*1024**3})
  state['stage']='MAP_ALL_DE_NOVO_CONTIGS';save('STATUS.json',state)
  with (O/'minimap.stderr').open('w') as err1,(O/'sort.stderr').open('w') as err2:
   p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err1,preexec_fn=limit)
   q=subprocess.Popen(sort,stdin=p.stdout,stderr=err2,preexec_fn=limit);p.stdout.close();codeq=q.wait();codep=p.wait();assert codep==codeq==0,(codep,codeq)
  subprocess.run([S,'index','-@','2',str(bam)],check=True,preexec_fn=limit)
  import pysam
  pysam.quickcheck(str(bam))
  resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
  reports=[]
  for name in ids:
   state['stage']='CALL_CONTIG';state['contig']=name;state['completed_contigs']=len(reports);save('STATUS.json',state)
   seq=refs[name];length=len(seq)
   with pysam.AlignmentFile(str(bam),'rb',require_index=True,threads=1) as b:evidence,filtered=call.collect_evidence(b,name,length)
   masks,triggers=call.make_masks(evidence,length)
   depths=[call.column_depth(evidence,i) for i in range(length)]
   observed=[depths[i] for i in range(200,length-200) if depths[i]>0];median=statistics.median(observed) if observed else 0
   result=[];counts={};disagreements=[]
   cols=['position_1based','assembly_base','observed_dominant','validated_base','DP','dominant_fraction','dominant_fwd','dominant_rev','failure_reasons']
   with (O/(name+'.per_base.tsv')).open('w') as h:
    h.write('\t'.join(cols)+'\n')
    for i in range(length):
     base,reasons,dp,denom,dominant,support=call.evaluate_position(evidence,i,length,masks)
     if not median or dp<median*.25:reasons.append('below0.25x_contig_median')
     if median and dp>median*4:reasons.append('above4x_contig_median')
     if base!='N' and base!=seq[i]:reasons.append('assembly_disagreement');disagreements.append({'position1':i+1,'assembly':seq[i],'dominant':base,'DP':dp,'support':support})
     if reasons:base='N'
     result.append(base)
     for reason in reasons:counts[reason]=counts.get(reason,0)+1
     h.write('\t'.join(map(str,[i+1,seq[i],dominant,base,dp,support/denom if denom else 'NA',evidence[dominant+'_fwd'][i],evidence[dominant+'_rev'][i],';'.join(reasons) or 'PASS']))+'\n')
   output=''.join(result)
   (O/(name+'.validated.fa')).write_text('>'+name+'_read_validated_partial\n'+output+'\n')
   item={'contig':name,'length':length,'median_nonzero_internal_template_DP':median,'callable_bases':length-output.count('N'),'unknown_N_bases':output.count('N'),'failure_counts_overlapping':counts,'assembly_disagreements':disagreements,'filtered_alignment_observations':filtered,'sequence_sha256':hashlib.sha256(output.encode()).hexdigest(),'indel_risk_trigger_count':len(triggers)}
   reports.append(item);save(name+'.metrics.json',item)
  state.update(status='COMPLETE_READBACK_PARTIAL_CONTIGS',stage='COMPLETE',completed_utc=call.utcnow(),elapsed_seconds=time.monotonic()-start,max_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,child_peak_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024,contigs=reports,BAM_sha256=sha(bam),limitations=['All contigs from same nonrandom first1M paired reads; not independent biological validation.','CP-like origin inferred from similarity/depth; nuclear or mitochondrial CP-like material not mathematically excluded.','No physical circular structure, IR copy count or heteroplasmy validation.','Unmarked PCR duplicates remain; overlapping mates collapsed only.','Contig ends and indel neighborhoods masked. No reference base filling.'])
  save('STATUS.json',state);print(json.dumps(state,indent=2))
 except BaseException as e:state.update(status='FAILED_PRESERVED',error=repr(e),completed_utc=call.utcnow());save('STATUS.json',state);raise
if __name__=='__main__':run()
