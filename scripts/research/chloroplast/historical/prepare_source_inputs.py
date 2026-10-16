#!/usr/bin/env python3
"""Create strict paired shared10 matrices; never start inference."""
from pathlib import Path
import argparse,csv,json,shutil
import portable_cli as portable
def configure(config_path):
    return portable.configure(config_path, ('design_dir','source_dir','historical_dir','iqtree','out'))

def plan(config_path):
    c=configure(config_path);w=c['design_dir'];src=c['source_dir'];hist=c['historical_dir']
    inputs={n:w/n for n in ['INPUT_METADATA.json','DESIGN.json','joint.partitions.nex',
                            'HT_cherry.constraint.nwk','CP_TGZ_split.constraint.nwk']}
    inputs.update({n:src/n for n in ['Kita_OQ_shared.control.fa','Kita_OQ_shared.replacement.fa','Kita_OQ_shared.columns.tsv']})
    if (w/'INPUT_METADATA.json').is_file():
        for n in json.loads((w/'INPUT_METADATA.json').read_text())['source_file_hashes']:
            if Path(n).is_absolute() or '..' in Path(n).parts: raise ValueError('Unsafe source manifest member')
            inputs['source/'+n]=src/n
    if (w/'DESIGN.json').is_file():
        for n in json.loads((w/'DESIGN.json').read_text()).get('constraints',{}):
            if Path(n).name!=n or n in ('.','..'): raise ValueError('Unsafe constraint filename')
            inputs[n+'.constraint.nwk']=w/(n+'.constraint.nwk')
    inputs.update({n:hist/n for n in ['candidates.nwk','CANDIDATE_MANIFEST.json']})
    inputs.update(script=Path(__file__),helper=Path(__file__).with_name('source_au_common.py'))
    names=['OQ_joint.fa','Kita_joint.fa','joint.partitions.nex','HT_cherry.constraint.nwk',
           'CP_TGZ_split.constraint.nwk','INPUT_METADATA.json','DESIGN.json','joint.columns.tsv',
           'historical.candidates.nwk','historical.CANDIDATE_MANIFEST.json','joint.blocks.json','CONFIG.json','PREPARED.json']
    return portable.contract(c['_config'],inputs,{n:c['out']/n for n in names},
                             {'iqtree':c['iqtree']},c['out'],taxa=10,columns=122234)

def main():
    a=portable.parser(__doc__).parse_args();checked=portable.preflight(plan(a.config))
    if not a.execute: portable.emit(checked); return
    from source_au_common import sha, fasta, topology_id, has_split, save
    import hashlib
    c=configure(a.config);w=c['design_dir'];src=c['source_dir'];HIST=c['historical_dir'];IQ=c['iqtree'];a.out=c['out']
    meta=json.loads((w/'INPUT_METADATA.json').read_text());design=json.loads((w/'DESIGN.json').read_text())
    for n,h in meta['source_file_hashes'].items():assert sha(src/n)==h
    old=fasta(src/'Kita_OQ_shared.control.fa');new=fasta(src/'Kita_OQ_shared.replacement.fa')
    assert len(old)==len(new)==93 and all(old[k]==new[k] for k in set(old)&set(new))
    variants={}
    for label,rows in [('OQ_joint',old),('Kita_joint',new)]:
        newtip=meta['prospective_matrix_hashes_and_rows'][label]['source_tokoro_tip']
        seqs={alias:rows[newtip if alias=='tokoro' else original] for alias,original in meta['representatives'].items()}
        assert len(seqs)==10 and all(len(s)==122234 for s in seqs.values())
        text=''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items());assert hashlib.sha256(text.encode()).hexdigest()==meta['prospective_matrix_hashes_and_rows'][label]['serialized_single_line_fasta_sha256_if_prepared']
        variants[label]=text
    assert sha(HIST/'candidates.nwk')==design['historical_candidate_sha256']
    old_entries=json.loads((HIST/'CANDIDATE_MANIFEST.json').read_text())['candidates'];assert len(old_entries)==2
    historical_lines=(HIST/'candidates.nwk').read_text().splitlines();assert len(historical_lines)==2
    for line,e in zip(historical_lines,old_entries):assert topology_id(line)==topology_id(e['newick'])
    for n in design['constraints']:
        content=(w/(n+'.constraint.nwk')).read_text();assert has_split(content,design['constraints'][n],meta['shared_labels'])
    config={'status':'FROZEN_EXECUTABLE_INPUTS_NOT_LAUNCHED','variants':['OQ_joint','Kita_joint'],'tips':meta['shared_labels'],'constraints':design['constraints'],
            'search_seeds':design['search_seeds'],'constrained_search_seeds':design['constrained_search_seeds'],'AU_seed':design['AU_seed'],'AU_RELL_replicates':10000,
            'rescue_seeds':[20261010,20261011],'rescue_trigger_logL_margin':0.01,'maximum_rescues_per_variant':2,'candidate_max_before_dedup':18,
            'threads':4,'cgroup_memory_max_bytes':4*1024**3,'cpu_affinity_max_cpus':4,'iqtree':str(IQ),'iqtree_sha256':sha(IQ),'pilot_iterations':10,'pilot_timeout_per_variant':150,
            'nuclear_unchanged':True,'nuclear_excluded_loci':['2661at3193'],'minimum_free_disk_bytes':2*1024**3}
    receipt={'status':'PREPARED_PAIRED_SOURCE_INPUTS_NO_INFERENCE','columns':122234,'taxa_per_variant':10,'all_other_nine_sequences_identical':True,
             'source_metadata_sha256':sha(w/'INPUT_METADATA.json'),'design_sha256':sha(w/'DESIGN.json'),'prepare_script_sha256':sha(__file__),
             'variants':meta['prospective_matrix_hashes_and_rows'],'source_files':meta['source_file_hashes'],'historical_candidate_count':2,'nuclear_input_or_results_changed':False}
    a.out.mkdir(parents=True,exist_ok=False)
    for label,text in variants.items():(a.out/(label+'.fa')).write_text(text)
    for n in ['joint.partitions.nex','HT_cherry.constraint.nwk','CP_TGZ_split.constraint.nwk','INPUT_METADATA.json','DESIGN.json']:shutil.copy2(w/n,a.out/n)
    shutil.copy2(src/'Kita_OQ_shared.columns.tsv',a.out/'joint.columns.tsv');shutil.copy2(HIST/'candidates.nwk',a.out/'historical.candidates.nwk');shutil.copy2(HIST/'CANDIDATE_MANIFEST.json',a.out/'historical.CANDIDATE_MANIFEST.json')
    blocks=[{'name':name,'start_1based':v[0],'end_1based':v[1],'sites':v[1]-v[0]+1} for name,v in meta['region_boundaries_1based'].items()];save(a.out/'joint.blocks.json',blocks);save(a.out/'CONFIG.json',config)
    receipt['files']={p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(a.out.iterdir()) if p.is_file()};save(a.out/'PREPARED.json',receipt)
    print(json.dumps({'status':receipt['status'],'path':str(a.out),'receipt_sha256':sha(a.out/'PREPARED.json'),'matrix_sha256':{v:sha(a.out/(v+'.fa')) for v in variants}},indent=2))
if __name__=='__main__':main()
