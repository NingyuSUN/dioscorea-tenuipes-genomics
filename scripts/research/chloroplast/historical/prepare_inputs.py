#!/usr/bin/env python3
"""Copy already completed CP matrices and bind exact fixed-column/model controls."""
from pathlib import Path
import argparse, collections, hashlib, json, re, shutil, time

ROOT=Path('workspace/wgs/dt_research_p57_20260930/cp_trees')
SOURCE=Path('workspace/historical_research/historical_Kita2022_prepare03_circular_both')
P19=Path('workspace/wgs/dt_research_p19_20260928/whole_plastome')
IQ=Path('workspace/wgs/dt_nuclear_expanded_20260915/tools/bin/iqtree2')
EXPECTED={
 'plus_Kita94.fixed_columns.fa':'59437ccd800c687f94ff39bf84dd99191930688fe50c08b6aca710d6ee51cdd8',
 'plus_Kita_Edo95.fixed_columns.fa':'45c23615b75b8a641b5036cbc40faafee39107c15f270653d966fc34241c0cf5',
 'fixed_regions.nex':'37336bcd8f6c1fede5806ec893b2baabbe20bbcc91d938649b8bf623a70b1f89',
 'PREPARATION_RECEIPT.json':'0fd1a4e26ebef9b111363ce78f7e63bc4ec50fdca0afc1fd8a7fdb3e209fd9a2'}
P19_EXPECTED={'whole_plastome_singleIR.fa':'82ad8eccc0a7d5494b025a14f1c2d66a3f77f8b86fb192edfbaa7eb003e7752a',
 'main.best_scheme.nex':'c6ccbe85a37c76bead22cd651334294d49b03989a6998408345e5c8d6fd3959b',
 'main.best_model.nex':'d4aef60d16ca2f80fee5ad2f33f3a909ead4e0a0c3447793367c2cc4a307ffdf',
 'main.log':'03c601a583058428d005e81938d0ba9f0fe1eb87179c1247d4a6b9b76fdeafb7',
 'main.treefile':'54c0fc60fa3df340240026368101622f8de6256a291536414e83934af806ed27'}
KITA='Dtokoro_Kita1_2022_author_declared_public_CP'
HACHI='Edo_CP_read_validated_partial_coordinate_projection'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')
def fasta(p):
    rows={};name=None
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            name=line[1:].strip();assert name and name not in rows;rows[name]=''
        elif line.strip():
            assert name;rows[name]+=line.strip().upper()
    assert rows and all(set(s)<=set('ACGTRYSWKMBDHVN-?') for s in rows.values())
    return rows
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'inputs');ap.add_argument('--execute-copy',action='store_true');a=ap.parse_args()
    for n,h in EXPECTED.items():assert sha(SOURCE/n)==h,n
    for n,h in P19_EXPECTED.items():assert sha(P19/n)==h,n
    receipt=json.loads((SOURCE/'PREPARATION_RECEIPT.json').read_text());assert receipt['status']=='HISTORICAL_2022_PREPARED_SOURCE_SENSITIVITY_NO_TREE'
    for n in ['plus_Kita94.fixed_columns.fa','plus_Kita_Edo95.fixed_columns.fa','fixed_regions.nex']:
        assert receipt['files'][n]['sha256']==EXPECTED[n]
    original=fasta(P19/'whole_plastome_singleIR.fa');p94=fasta(SOURCE/'plus_Kita94.fixed_columns.fa');p95=fasta(SOURCE/'plus_Kita_Edo95.fixed_columns.fa')
    assert len(original)==93 and len(p94)==94 and len(p95)==95
    assert set(p94)-set(original)=={KITA} and set(p95)-set(p94)=={HACHI}
    assert all(p94[k]==v==p95[k] for k,v in original.items())
    assert all(p95[k]==v for k,v in p94.items())
    assert all(len(s)==124448 for panel in [original,p94,p95] for s in panel.values())
    assert 'DT_local_path2' in original and 'OQ525999.1' in original
    assert all(x in original for x in ['KX171420.1','NC_036658.1'])
    model=(P19/'main.best_scheme.nex').read_text()
    regions=re.findall(r'charset\s+(\w+)\s*=\s*(\d+)-(\d+)',model,re.I)
    assert regions==[('LSC','1','80600'),('IR','80601','105951'),('SSC','105952','124448')]
    assert re.findall(r'(TVM\+F\+G4):\s*(LSC|IR|SSC)',model)==[('TVM+F+G4','LSC'),('TVM+F+G4','IR'),('TVM+F+G4','SSC')]
    command=next(l for l in (P19/'main.log').read_text().splitlines() if l.startswith('Command: '))
    assert ' -p ' in command and ' -m MFP+MERGE ' in command
    panels={}
    for name,rows in [('control94',p94),('main95',p95)]:
        groups=collections.defaultdict(list)
        for k,s in rows.items():groups[hashlib.sha256(s.encode()).hexdigest()].append(k)
        panels[name]={'taxa':len(rows),'columns':124448,'tips':list(rows),'identical_sequence_groups':[v for v in groups.values() if len(v)>1],
                      'new_tip_callable':{k:{'ACGT':sum(b in 'ACGT' for b in rows[k]),'N':rows[k].count('N'),'gap':rows[k].count('-')} for k in [KITA,HACHI] if k in rows}}
    summary={'status':'VERIFIED_PREPARE_ONLY_NO_TREE','original93_unchanged_in_both_panels':True,'all94_rows_identical_in_95':True,
             'same_columns':124448,'column_occupancy_not_recomputed':True,'panels':panels,
             'nuclear_inputs':'unchanged; 2661 excluded in existing nuclear baseline',
             'P19_original_command':command,'model_policy':'P19 selected model classes fixed: three TVM+F+G4 partitions, -p edge-linked proportional; numerical model parameters and branch lengths re-estimated independently per panel.',
             'source_types':{KITA:'Natsume2022 author-declared Kita1 public CP assembly; retained original provenance label',HACHI:'DRR351751/SAMD00444457 Hachi1 (Edo-dokoro) first-1M-pair de novo/readback partial coordinate projection'},
             'no_source_replacement':'OQ525999.1 and all original93 remain; historical public Kita and Hachi are additional individuals, not replicate species observations.',
             'prepared_epoch':time.time(),'prepare_script_sha256':sha(__file__)}
    if not a.execute_copy:print(json.dumps(summary,indent=2));return
    a.out.mkdir(parents=True,exist_ok=False)
    copies={'control93.fixed_columns.fa':P19/'whole_plastome_singleIR.fa',
            'control94.fixed_columns.fa':SOURCE/'plus_Kita94.fixed_columns.fa','main95.fixed_columns.fa':SOURCE/'plus_Kita_Edo95.fixed_columns.fa',
            'fixed_regions.nex':SOURCE/'fixed_regions.nex','fixed_models.nex':P19/'main.best_scheme.nex',
            'SOURCE_PREPARATION_RECEIPT.json':SOURCE/'PREPARATION_RECEIPT.json',
            'P19.main.log':P19/'main.log','P19.main.best_model.nex':P19/'main.best_model.nex','P19.main.treefile':P19/'main.treefile'}
    for name,src in copies.items():shutil.copy2(src,a.out/name);assert sha(a.out/name)==sha(src)
    config={'iqtree':str(IQ),'iqtree_sha256':sha(IQ),'threads':4,'cgroup_memory_max_bytes':6*1024**3,'cgroup_memory_swap_max_bytes':0,
            'cpu_affinity_max_cpus':4,'sequential_panels':['main95','control94'],'outgroups':['KX171420.1','NC_036658.1'],
            'seed':20260930,'pilot_iterations':10,'pilot_total_budget_seconds':300,'full_bootstrap':1000,'full_alrt':1000,'full_bnni':True,
            'full_stage_does_not_follow_pilot_automatically':True,'model_file':'fixed_models.nex','minimum_free_disk_bytes':2*1024**3}
    save(a.out/'CONFIG.json',config)
    summary['files']={p.name:{'sha256':sha(p),'bytes':p.stat().st_size,'source':str(copies.get(p.name,'generated CONFIG'))} for p in a.out.iterdir() if p.is_file()}
    save(a.out/'PREPARED.json',summary)
    print(json.dumps({'status':summary['status'],'output':str(a.out),'receipt_sha256':sha(a.out/'PREPARED.json'),'panel_counts':{k:v['taxa'] for k,v in panels.items()}},indent=2))
if __name__=='__main__':main()
