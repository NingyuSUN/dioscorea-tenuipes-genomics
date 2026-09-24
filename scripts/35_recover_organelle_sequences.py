#!/usr/bin/env python3
"""Recover graph-supported plastome paths; audit retained FASTQ, never fill a reference."""
import argparse, gzip, hashlib, json, re
from collections import Counter
from pathlib import Path


def fasta(path):
    records = {}; name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                name = line[1:].strip(); assert name not in records
                records[name] = ''
            elif line.strip():
                assert name is not None
                records[name] += line.strip().upper()
    assert records and all(records.values())
    return records


def rc(seq):
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def write_fasta(path, name, seq):
    with open(path, 'w') as f:
        f.write('>' + name + '\n')
        for i in range(0, len(seq), 80): f.write(seq[i:i+80] + '\n')


def graph(path):
    segments = {}; edges = set(); ks = set()
    for line in Path(path).read_text().splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            assert f[1] not in segments and f[2] != '*'
            segments[f[1]] = f[2].upper()
        elif f[0] == 'L':
            assert re.fullmatch(r'\d+M', f[5]), 'Unsupported overlap CIGAR'
            k = int(f[5][:-1]); ks.add(k)
            edges.add((f[1]+f[2], f[3]+f[4]))
            flip = lambda x: '+' if x == '-' else '-'
            edges.add((f[3]+flip(f[4]), f[1]+flip(f[2])))
    assert len(ks) == 1
    return segments, edges, ks.pop()


def recover(header, seq, segments, edges, k):
    assert header.endswith('(circular)')
    path = header.removesuffix('(circular)').split(',')
    pieces = [segments[x[:-1]] if x[-1] == '+' else rc(segments[x[:-1]]) for x in path]
    for i in range(len(path)):
        j = (i+1) % len(path)
        assert (path[i], path[j]) in edges, 'Path edge absent'
        assert pieces[i][-k:] == pieces[j][:k], 'Graph overlap mismatch'
    joined = ''.join(p[:-k] for p in pieces)
    assert len(joined) == len(seq)
    offset = (joined+joined).find(seq)
    assert 0 <= offset < len(seq), 'FASTA differs from graph-derived path'
    return dict(path=path, overlap=k, rotation_offset=offset, length=len(seq), graph_matches=True)


def reads(path):
    with gzip.open(path, 'rt') as h:
        while True:
            name = h.readline()
            if not name: break
            seq, plus, qual = h.readline().strip(), h.readline(), h.readline().rstrip('\r\n')
            assert name.startswith('@') and plus.startswith('+') and len(seq) == len(qual) and seq
            yield name.split()[0].removesuffix('/1').removesuffix('/2'), len(seq)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('root', type=Path); args=ap.parse_args()
    root=args.root; inp=root/'inputs'; out=root/'recovered'; out.mkdir(exist_ok=False)
    manifest=json.loads((root/'source_manifest.json').read_text())
    for row in manifest:
        p=inp/row['name']
        assert p.stat().st_size==row['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256']
    segments, edges, k=graph(inp/'embplant_pt.K105.complete.graph1.selected_graph.gfa')
    summary={'status':'RECOVERED_PLASTOME_WITH_GRAPH_SUPPORT_NOT_FINAL_PUBLICATION_VALIDATION',
             'provenance':'Historical GetOrganelle v1.7.6.1 output from Dt_Izu_clean paired reads; accession and current nuclear assembly binding unresolved',
             'input_hashes_verified':len(manifest),'chloroplast_paths':[]}
    for i, name in enumerate(['D.tenuipes_chloroplast.fasta','embplant_pt.K105.complete.graph1.2.path_sequence.fasta'],1):
        fs=fasta(inp/name); assert len(fs)==1; header, seq=next(iter(fs.items()))
        result=recover(header,seq,segments,edges,k)
        result.update(source=name,n_bases=seq.count('N'),gc_percent=100*(seq.count('G')+seq.count('C'))/len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest())
        assert not set(seq)-set('ACGT')
        output=f'Dtenuipes_Izu_reference_candidate.plastome.path{i}.fasta'
        write_fasta(out/output,f'Dtenuipes_plastome_path{i}',seq); result['output']=output
        if i==1:
            shift=len(seq)//2
            write_fasta(out/'plastome.rotated_validation.fasta','plastome_rotated',seq[shift:]+seq[:shift])
            write_fasta(out/'circular_junction.fasta','circular_junction',seq[-1000:]+seq[:1000])
            result['rotation_validation_offset']=shift
        summary['chloroplast_paths'].append(result)
    summary['structural_difference']={'same_graph':True,'paths_differ_only_in_SSC_vertex_orientation': summary['chloroplast_paths'][0]['path'][2][:-1]==summary['chloroplast_paths'][1]['path'][2][:-1] and all(summary['chloroplast_paths'][0]['path'][j]==summary['chloroplast_paths'][1]['path'][j] for j in (0,1,3))}
    r1=reads(inp/'extended_1_paired.fq.gz');r2=reads(inp/'extended_2_paired.fq.gz');counts=Counter();n=0;bases=0
    while True:
        a=next(r1,None);b=next(r2,None)
        assert (a is None)==(b is None),'Unequal paired file counts'
        if a is None:break
        assert a[0]==b[0], 'Paired names mismatch'
        n+=1; bases+=a[1]+b[1]; counts.update([a[1],b[1]])
    summary['retained_reads']={'paired_fragments':n,'paired_bases':bases,'read_length_counts':dict(counts),'unpaired':{}}
    for i in (1,2):
        rows=list(reads(inp/f'extended_{i}_unpaired.fq.gz'))
        summary['retained_reads']['unpaired'][str(i)]={'reads':len(rows),'bases':sum(x[1] for x in rows)}
    summary['limitations']=['Reads were recruited for plastome assembly; not independent validation or unbiased mitochondrial data.', 'Two flip-flop paths are structural alternatives, not two individuals or evidence of heteroplasmy.', 'No mitochondrial genome reconstructed by this script.']
    (out/'recovery.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__=='__main__':main()
