#!/usr/bin/env python3
from pathlib import Path
from collections import Counter, defaultdict
import argparse, hashlib, json
from datetime import datetime, timezone

parser=argparse.ArgumentParser(description='Read-only audit of preserved D. tenuipes BUSCO, AGP and hints outputs; does not launch pipelines.')
parser.add_argument('--wgs', type=Path, required=True, help='Existing server wgs directory containing both historical BUSCO runs and RagTag outputs')
root=parser.parse_args().wgs.resolve()
result={'checked_at_utc':datetime.now(timezone.utc).isoformat()}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def agp(rel):
    placements={}; gaps=0
    for line in (root/rel).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t')
        if f[4] in ['N','U']:gaps+=int(f[5]);continue
        assert f[5] not in placements, f[5]
        placements[f[5]]=(f[0],int(f[1]),int(f[2]),int(f[6]),int(f[7]),f[8])
    return placements,gaps
kita,kgap=agp('ragtag_out_Kita1/ragtag.scaffold.agp')
waka,wgap=agp('ragtag_out_Waka1/ragtag.scaffold.agp')
assert set(kita)==set(waka)
matrix=defaultdict(Counter)
for line in (root/'kita_vs_waka.paf').read_text().splitlines():
    f=line.split('\t');matrix[f[0]][f[5]]+=int(f[9])
mapping={q:c.most_common(1)[0][0] for q,c in matrix.items()}
agreement=Counter(); examples=[]; bp=Counter()
for contig,k in kita.items():
    w=waka[contig];ka=k[0].removesuffix('_RagTag');wa=w[0].removesuffix('_RagTag')
    if ka in mapping and wa in mapping.values():
        cat='both_anchored_same_chr' if mapping[ka]==wa else 'both_anchored_different_chr'
        if cat.endswith('different_chr') and len(examples)<10:examples.append([contig,ka,wa])
    elif ka in mapping:cat='Kita_only_anchored'
    elif wa in mapping.values():cat='Waka_only_anchored'
    else:cat='both_unplaced'
    agreement[cat]+=1;bp[cat]+=k[4]-k[3]+1
result['scaffold_comparison']={'Kita_to_Waka_from_PAF':mapping,'contigs':dict(agreement),'bp':dict(bp),
    'discordance_examples':examples,'Kita_gap_bp':kgap,'Waka_gap_bp':wgap,
    'boundary':'Chromosome membership comparison only. Neither sequence order nor orientation nor biological joins independently validated.',
    'PAF_best_to_second_ratios':{q: c.most_common()[0][1]/c.most_common()[1][1] for q,c in matrix.items() if len(c)>1}}

result['lift_validation']={}
for basename in ['hints_fixed_introns','hints_fixed_exonparts']:
    expected=Counter();actual=Counter();orient=Counter()
    for line in (root/(basename+'.gff')).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t'); seq,start,end,cstart,cend,strand=kita[f[0]]
        s,e=int(f[3]),int(f[4]);assert cstart<=s<=e<=cend
        if strand=='+':ns,ne=start+s-cstart,start+e-cstart
        elif strand=='-':
            ns,ne=start+cend-e,start+cend-s
            f[6]={'+':'-','-':'+','.':'.','?':'?'}[f[6]]
        else:raise ValueError(strand)
        orient[strand]+=1; f[0]=seq;f[3]=str(ns);f[4]=str(ne)
        expected[tuple(f)]+=1
    for line in (root/(basename+'.scaffolded.gff')).read_text().splitlines():
        if line and not line.startswith('#'):actual[tuple(line.split('\t'))]+=1
    result['lift_validation'][basename]={'input_records':sum(expected.values()),'output_records':sum(actual.values()),
        'missing_expected':sum((expected-actual).values()),'unexpected_output':sum((actual-expected).values()),
        'component_orientation_counts':dict(orient),'input_sha256':sha(root/(basename+'.gff')),
        'output_sha256':sha(root/(basename+'.scaffolded.gff'))}

result['busco']={}; calls={}
for key,name,fasta in [('current','Dtenuipes_proteins','braker.aa'),('old','Dtenuipes_proteins_ALT','D.tenuipes_protein.fa')]:
    d=root/'busco_runs'/name
    summary=d/f'short_summary.specific.embryophyta_odb10.{name}.txt'
    table=d/'run_embryophyta_odb10/full_table.tsv'
    byid={}
    for line in table.read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t');old=byid.setdefault(f[0],f[1]);assert old==f[1]
    assert len(byid)==1614
    calls[key]=byid
    result['busco'][key]={'summary':summary.read_text(),'summary_sha256':sha(summary),
         'table_sha256':sha(table),'input_sha256':sha(root/fasta),'status_counts':dict(Counter(byid.values()))}
assert set(calls['old'])==set(calls['current'])
transitions=Counter((calls['old'][b],calls['current'][b]) for b in calls['old'])
result['busco']['old_to_current_transitions']=[{'old':a,'current':b,'count':n} for (a,b),n in sorted(transitions.items())]
complete={'Complete','Duplicated'}
result['busco']['old_complete_current_not_complete_ids']=[b for b in sorted(calls['old']) if calls['old'][b] in complete and calls['current'][b] not in complete]
result['busco']['current_complete_old_not_complete_ids']=[b for b in sorted(calls['old']) if calls['current'][b] in complete and calls['old'][b] not in complete]
result['file_hashes']={str(p.relative_to(root)):sha(p) for p in [root/'ragtag_out_Kita1/ragtag.scaffold.agp',root/'ragtag_out_Waka1/ragtag.scaffold.agp',root/'kita_vs_waka.paf',root/'fix_hints.py',root/'run_repeat.sh',root/'busco_alt_proteins.sh']}
print(json.dumps(result,ensure_ascii=False,indent=2))
