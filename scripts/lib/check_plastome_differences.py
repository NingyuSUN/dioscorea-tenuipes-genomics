#!/usr/bin/env python3
"""Read support at alignment differences; counts are descriptive, not independent replication."""
import csv,json,sys
from collections import Counter
from pathlib import Path
import pysam
root=Path(sys.argv[1]);p=root/'published_comparison';rows=list(csv.DictReader((p/'substitutions_for_read_validation.tsv').open(),delimiter='\t'))
out=[]
bam=pysam.AlignmentFile(str(root/'validation/native.all.bam'),'rb')
for row in rows:
    pos=int(row['native_path1_position'])-1
    # Use the rotated reference near the artificial native ends.
    rotated=pos<1000 or pos>=bam.lengths[0]-1000
    b=pysam.AlignmentFile(str(root/'validation/rotated.all.bam'),'rb') if rotated else bam
    querypos=(pos-bam.lengths[0]//2)%bam.lengths[0] if rotated else pos
    counts=Counter();unique=Counter()
    for col in b.pileup(b.references[0],querypos,querypos+1,truncate=True,stepper='nofilter',min_base_quality=0,max_depth=100000,ignore_overlaps=False,ignore_orphans=False):
        for pr in col.pileups:
            r=pr.alignment;i=pr.query_position
            if i is None or r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail:continue
            if r.query_qualities is None or r.query_qualities[i]<20:continue
            base=r.query_sequence[i].upper();counts[base]+=1
            if r.mapping_quality>=20:unique[base]+=1
    expected=row['native_path1_base'];n=sum(counts.values());u=sum(unique.values())
    row.update(depth_BQ20=n,own_base_support=counts[expected],own_base_fraction=counts[expected]/n if n else None,unique_depth=u,unique_own_base_fraction=unique[expected]/u if u else None,counts=json.dumps(dict(counts),sort_keys=True))
    out.append(row)
    if rotated:b.close()
bam.close()
with (p/'substitution_read_support.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
s={'substitution_sites':len(out),'sites_depth_ge20_own_fraction_ge0_9':sum(r['depth_BQ20']>=20 and r['own_base_fraction']>=.9 for r in out),'sites_unique_depth_ge20_own_fraction_ge0_9':sum(r['unique_depth']>=20 and r['unique_own_base_fraction']>=.9 for r in out),'min_depth':min(r['depth_BQ20'] for r in out),'min_own_base_fraction':min(r['own_base_fraction'] for r in out if r['own_base_fraction'] is not None),'limitations':['Assembly-recruited reads, no independent library validation.','IR multi-mapping retained in total counts; unique mapping count reported separately.','Counts include overlapping mates and do not represent independent molecules.','Indels and gene annotation not validated by this test.']}
(p/'substitution_read_support.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s))
