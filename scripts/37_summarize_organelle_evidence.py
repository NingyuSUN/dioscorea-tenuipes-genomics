#!/usr/bin/env python3
"""Summarize plastome validation and extract unjoined mitochondrial homology candidates."""
import argparse, importlib.util, json, statistics, hashlib
from pathlib import Path
from collections import defaultdict


def union(intervals):
    out=[]
    for start,end in sorted(intervals):
        if out and start<=out[-1][1]:out[-1][1]=max(out[-1][1],end)
        else:out.append([start,end])
    return out


def span(intervals):return sum(b-a for a,b in union(intervals))


def readpaf(path):
    rows=[]
    for line in Path(path).open():
        f=line.split()
        if len(f)<12:continue
        row=dict(query=f[0],qlen=int(f[1]),qs=int(f[2]),qe=int(f[3]),strand=f[4],target=f[5],tlen=int(f[6]),ts=int(f[7]),te=int(f[8]),matches=int(f[9]),alnlen=int(f[10]),mapq=int(f[11]))
        if row['alnlen']>=1000 and row['matches']/row['alnlen']>=.9:rows.append(row)
    return rows


def depth(path,shift=0):
    values=[]
    for l in Path(path).open():values.append(int(l.split()[2]))
    if shift:values=values[-shift:]+values[:-shift]
    return values


def dstat(v):
    return dict(bases=len(v),min=min(v),median=statistics.median(v),mean=statistics.mean(v),zero=sum(x==0 for x in v),below10=sum(x<10 for x in v),below20=sum(x<20 for x in v))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('--genome',required=True,type=Path);a=ap.parse_args()
    root=a.root;v=root/'validation';out=root/'results';out.mkdir(exist_ok=False)
    assert (v/'ALIGNMENT_COMPLETE').exists(),'Read validation not complete'
    spec=importlib.util.spec_from_file_location('recover',root/'35_recover_organelle_sequences.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    summary={'recovery':json.loads((root/'recovered/recovery.json').read_text()),'paf_filter':{'min_alignment_bp':1000,'min_identity':.90,'purpose':'Exploratory homology screening; not organelle provenance certification'}}
    n=summary['recovery']['chloroplast_paths'][0]['length'];shift=n//2
    native=depth(v/'native.depth.tsv');rot=depth(v/'rotated.depth.tsv',shift)
    assert len(native)==len(rot)==n
    # Prefer native mapping away from its artificial ends; rotated mapping near them.
    combined=[rot[i] if i<1000 or i>=n-1000 else native[i] for i in range(n)]
    un=depth(v/'native.unique_depth.tsv');ur=depth(v/'rotated.unique_depth.tsv',shift)
    uc=[ur[i] if i<1000 or i>=n-1000 else un[i] for i in range(n)]
    summary['read_depth']={'native':dstat(native),'rotated':dstat(rot),'circular_boundary_aware':dstat(combined),'unique_mapping_circular_boundary_aware':dstat(uc),'base_quality_min':20,'all_mapping_mapq_min':0,'unique_mapping_mapq_min':20,'rotation_boundary_window_bp':1000,'caveat':'Same assembly-recruited reads; duplicated/inverted repeat ambiguity and recruitment bias remain. No duplicate removal.'}
    sp=[x.split('\t')[0] for x in (v/'circular_junction_support.tsv').read_text().splitlines()[1:]]
    summary['circular_junction']={'spanning_alignments':len(sp),'unique_read_names':len(set(sp)),'min_flank_bp':30,'min_mapq':20,'caveat':'Sequence-junction support does not prove a single physical circular molecule.'}
    with (out/'plastome.depth.tsv').open('w') as f:
        f.write('position_1based\tdepth_BQ20\tdepth_BQ20_MAPQ20\n')
        for i,(d,u) in enumerate(zip(combined,uc),1):f.write(f'{i}\t{d}\t{u}\n')
    by=defaultdict(lambda:defaultdict(list));target=defaultdict(list);lengths={}
    for row in readpaf(v/'genome_vs_organelle.paf'):
        by[row['query']][row['target']].append((row['qs'],row['qe']));target[row['target']].append((row['ts'],row['te']));lengths[row['query']]=row['qlen']
    rows=[];selected=set()
    for q,d in by.items():
        cp=span(d.get('D.tokoro_chro',[]));mt=span(d.get('D.tokoro_mt',[]));ln=lengths[q]
        # This is a retained-candidate filter, NOT a mitochondrial genome call.
        keep=mt>=5000 and mt/ln>=.5 and cp/ln<.2
        if keep:selected.add(q)
        rows.append((q,ln,mt,mt/ln,cp,cp/ln,int(keep)))
    rows.sort(key=lambda x:x[2],reverse=True)
    with (out/'genome_organelle_homology.tsv').open('w') as f:
        f.write('contig\tlength\tmt_homology_bp\tmt_fraction\tcp_homology_bp\tcp_fraction\tretained_candidate\n')
        for row in rows:f.write('\t'.join(map(str,row))+'\n')
    # Streaming extraction preserves complete source contigs; no reference sequence is copied.
    count=0;bp=0
    with a.genome.open() as fin,(out/'mitochondrial_homology_candidates.UNVALIDATED.fasta').open('w') as fout:
        keep=False
        for line in fin:
            if line.startswith('>'):
                keep=line[1:].split()[0] in selected
                if keep:count+=1
            if keep:
                fout.write(line)
                if not line.startswith('>'):bp+=len(line.strip())
    summary['mitochondria']={'status':'NO_VALIDATED_MITOGENOME','candidate_contigs':count,'candidate_bp':bp,'retention_filter':'mt homology >=5000 bp and >=50% of contig; cp homology <20%','reference_homology_covered_bp':{k:span(z) for k,z in target.items()},'limitations':['Candidate contigs may be nuclear mitochondrial insertions (NUMTs) or mixed/chimeric sequences.','Homology coverage of D. tokoro is not D. tenuipes mitogenome completeness.','No whole-genome read set or mitochondrial assembly graph recovered; do not join candidates or use as an authenticated mitochondrial genome.','Historical mitochondrial log reports incomplete/complicated graph, not a completed circular mitogenome.']}
    h=hashlib.sha256()
    with a.genome.open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):h.update(chunk)
    summary['genome_sha256']=h.hexdigest()
    (out/'organelle_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='recovery'}))

if __name__=='__main__':main()
