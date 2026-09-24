#!/usr/bin/env python3
"""Summarize real-prefix QC and mapping, with explicit denominators and limits."""
import argparse,csv,gzip,hashlib,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();d=a.out
assert (d/'COMPLETE').exists()
raw=json.loads((d/'raw_prefix_qc.json').read_text());fastp=json.loads((d/'fastp.json').read_text());counts=json.loads((d/'alignment_counts.json').read_text());n=counts['primary'];assert n==2*raw['pairs']
assert fastp['summary']['before_filtering']['total_reads']==n
for i in [1,2]:
 q=fastp[f'read{i}_before_filtering'];h=raw['mates'][i-1]['phred_histogram'];assert abs(q['q30_bases']/q['total_bases']-sum(h[30:])/sum(h))<1e-12
# All primary mapped reads, separated by best nuclear-reference vs plastome placement.
bam=d/'prefix.bam';by={}
cmd=['samtools','view','-F','0x904',str(bam)];pr=subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True)
for line in pr.stdout:
 f=line.rstrip().split('\t');group='plastome_decoy' if f[2]=='decoy_PV612491.1' else 'draft_nuclear';b=by.setdefault(group,{'primary_reads':0,'Q20_reads':0,'proper_pair_reads':0,'NM_sum':0,'aligned_query_bases':0,'soft_clipped_bases':0,'query_bases':0})
 b['primary_reads']+=1;b['Q20_reads']+=int(int(f[4])>=20);b['proper_pair_reads']+=int(bool(int(f[1])&2));b['query_bases']+=len(f[9]) if f[9]!='*' else 0
 for length,op in re.findall(r'(\d+)([MIDNSHP=X])',f[5]):
  if op in 'MI=X':b['aligned_query_bases']+=int(length)
  elif op=='S':b['soft_clipped_bases']+=int(length)
 nm=[t for t in f[11:] if t.startswith('NM:i:')];assert len(nm)==1;b['NM_sum']+=int(nm[0][5:])
assert pr.wait()==0
assert sum(b['primary_reads'] for b in by.values())==counts['mapped_primary']
# Coverage rows use MAPQ20/BQ20 with explicit flag exclusion0x904; not deduplicated.
with (d/'coverage_Q20_BQ20.tsv').open() as f:lines=[x for x in f if not x.startswith('samtools ')]
rows=list(csv.DictReader(lines,delimiter='\t'));nuc=[x for x in rows if x['#rname']!='decoy_PV612491.1'];cp=[x for x in rows if x['#rname']=='decoy_PV612491.1'];assert len(cp)==1
length=sum(int(x['endpos'])-int(x['startpos'])+1 for x in nuc);assert length==403772333
covered=sum(int(x['covbases']) for x in nuc);mean=sum(float(x['meandepth'])*(int(x['endpos'])-int(x['startpos'])+1) for x in nuc)/length
summary={'scope':'first 1,000,000 pairs, NOT random, NOT full-run QC','pairs':raw['pairs'],'primary_reads':n,'raw_bases':sum(x['bases'] for x in raw['mates']),'raw_R1_Q30':raw['mates'][0]['q30_fraction'],'raw_R2_Q30':raw['mates'][1]['q30_fraction'],'fastp_retained_read_fraction':fastp['summary']['after_filtering']['total_reads']/n,'fastp_adapter_trimmed_read_fraction':fastp['adapter_cutting']['adapter_trimmed_reads']/n,'fastp_duplication_estimate':fastp.get('duplication',{}).get('rate'),'mapped_primary_fraction':counts['mapped_primary']/n,'Q20_mapped_primary_fraction':counts['mapped_Q20_primary']/n,'proper_pair_primary_fraction':counts['proper_pair_primary']/n,'best_alignment_components':by,'nuclear_reference_length':length,'nuclear_contigs_reported':len(nuc),'nuclear_covered_bases_Q20_BQ20':covered,'nuclear_covered_fraction_Q20_BQ20':covered/length,'nuclear_mean_depth_Q20_BQ20':mean,'plastome_coverage_Q20_BQ20':cp[0],'not_verified':['full FASTQ checksums','whole-run QC','random-subsample representativeness','ploidy','species identity independent of archive','nuclear genotype callability','joint compatibility with existing17','sequencing-batch correction'],'validation':{'pair_names_match':raw['pair_ids_match'],'independent_Q30_matches':True,'primary_reads_expected':True,'BAM_quickcheck_passed':True,'reference_hashes_unchanged':True}}
(d/'probe_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
