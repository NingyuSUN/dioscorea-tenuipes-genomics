"""Post-terminal descriptive gCF sensitivity, excluding every locus with any logged warning.

Does not infer new gene trees or change the frozen primary concordance analysis.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,hashlib,json

CLASSES=('concordant','alternative_1','alternative_2','other_or_unresolved','not_decisive')
PAIRS=CLASSES[:3]


def tally(rows):
    c=Counter();seen=set()
    for r in rows:
        assert r['gene'] not in seen;seen.add(r['gene'])
        k=r['classification'];assert k in CLASSES;c[k]+=1
        v=r['resolving_edge_dual_high_support'];assert v in ('True','False')
        if v=='True':assert k in PAIRS;c['high_'+k]+=1
    den=len(rows)-c['not_decisive']
    out={'total_loci':len(rows),'decisive_loci':den,**{k:c[k] for k in CLASSES},
         **{'high_'+k:c['high_'+k] for k in PAIRS},
         'gCF_percent':100*c['concordant']/den if den else None,
         'decisive_without_dual_high_pairing':den-sum(c['high_'+k] for k in PAIRS)}
    assert den==sum(c[k] for k in CLASSES[:4])
    assert out['decisive_without_dual_high_pairing']>=0
    return out


def run(source,out):
    assert not out.exists(),'Refuse overwrite'
    status=json.loads((source/'FULL_CONCORDANCE_STATUS.json').read_text())
    assert status['status']=='COMPLETE_SCIENTIFIC_REVIEW_PENDING','Full run must finish first'
    audit=json.loads((source/'validation/ALL_LOCUS_AUDIT.json').read_text())
    cross=json.loads((source/'validation/ALL_GCF_CROSSCHECK.json').read_text())
    assert cross['status']=='PASS_INDEPENDENT_IQTREE_CROSSCHECK'
    assert audit['status']=='PASS_CORPUS_WITH_WARNINGS_RETAINED' and audit['loci']==173
    warnings={r['gene']:r['warnings'] for r in audit['records'] if r['warnings']}
    manifest=json.loads((source/'metadata/LOCUS_MANIFEST.json').read_text())
    genes={r['gene'] for r in manifest['loci']};assert len(genes)==173 and '2661at3193' not in genes
    assert {r['gene'] for r in audit['records']}==genes
    rows=list(csv.DictReader((source/'concordance_independent/LOCUS_BRANCH_CLASSIFICATIONS.tsv').open(),delimiter='\t'))
    by_branch=defaultdict(list)
    for r in rows:by_branch[r['reference_split']].append(r)
    original={r['reference_split']:r for r in csv.DictReader((source/'concordance_independent/BRANCH_CONCORDANCE.tsv').open(),delimiter='\t')}
    assert set(by_branch)==set(original) and len(by_branch)==9
    result=[]
    for branch,records in sorted(by_branch.items()):
        assert {r['gene'] for r in records}==genes and len(records)==173
        full=tally(records)
        for k,v in full.items():
            old=original[branch][k]
            if v is None:assert old==''
            elif k=='gCF_percent':assert abs(float(old)-v)<1e-9
            else:assert int(old)==v
        filtered=tally([r for r in records if r['gene'] not in warnings])
        delta=filtered['gCF_percent']-full['gCF_percent'] if filtered['gCF_percent'] is not None and full['gCF_percent'] is not None else None
        result.append({'reference_split':json.loads(branch),'all_loci':full,'without_warning_loci':filtered,'gCF_change_percentage_points':delta})
    out.mkdir(parents=True)
    packet={'status':'COMPLETE_DESCRIPTIVE_SENSITIVITY_REVIEW_PENDING','criterion':'Exclude all loci with any retained log warning; no topology or result-based selection',
        'posthoc':True,'primary_result_unchanged':True,'total_loci':173,'excluded_loci':warnings,'retained_loci':173-len(warnings),'branches':result,
        'source_sha256':{str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source/'validation/ALL_LOCUS_AUDIT.json',source/'validation/ALL_GCF_CROSSCHECK.json',source/'concordance_independent/LOCUS_BRANCH_CLASSIFICATIONS.tsv',source/'concordance_independent/BRANCH_CONCORDANCE.tsv']},
        'limits':['Post hoc QC sensitivity, not a replacement main analysis','No new concatenated tree or gene-tree inference','Warning-free outputs do not prove model adequacy or biological correctness','Raw gCF retains low-support resolutions and is not a probability of hybridization']}
    (out/'WARNING_SENSITIVITY.json').write_text(json.dumps(packet,indent=2)+'\n')
    print(json.dumps({'excluded':len(warnings),'retained':173-len(warnings),'branches':len(result)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.source,a.out)
