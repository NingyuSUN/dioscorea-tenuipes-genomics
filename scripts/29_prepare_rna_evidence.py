#!/usr/bin/env python3
"""Prepare traceable structural RNA evidence from complete StringTie GTF.

Python 3.9+ standard library. Output GFF3 is structural evidence, not a drop-in
BRAKER-ETP RNA input. No read-support multiplicity or transcript abundance is
inferred. A malformed transcript is quarantined in full, never exon by exon.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import quote

COMPLEMENT=str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn','TGCAYRKMVHDBNtgcayrkmvhdbn')
CANONICAL={'GT-AG','GC-AG','AT-AC'}

def read_lines(path):
    with path.open() as stream:
        yield from stream

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def fasta(path):
    seqs={};name=None;parts=[]
    def save():
        if name is not None:
            if name in seqs:raise ValueError('duplicate FASTA ID: '+name)
            seqs[name]=''.join(parts).upper()
            if not seqs[name]:raise ValueError('empty FASTA record: '+name)
    for line in read_lines(path):
        if line.startswith('>'):
            save();name=line[1:].split()[0];parts=[]
        elif line.strip():
            if name is None:raise ValueError('sequence before FASTA header')
            parts.append(line.strip())
    save()
    if not seqs:raise ValueError('empty FASTA')
    return seqs

def parse_gtf(path,genome):
    models=defaultdict(lambda:{'rows':[],'exons':[],'transcripts':[],'issues':set()})
    for n,line in enumerate(read_lines(path),1):
        if not line.strip() or line.startswith('#'):continue
        f=line.rstrip('\r\n').split('\t')
        if len(f)!=9:raise ValueError(f'GTF line {n}: expected 9 fields; cannot safely identify affected transcript')
        match=re.findall(r'(?:^|;)\s*transcript_id\s+"([^"\t\n]+)"',f[8])
        if len(match)!=1:raise ValueError(f'GTF line {n}: expected one transcript_id')
        m=models[match[0]];m['rows'].append(f)
        if f[2] not in {'transcript','exon'}:raise ValueError(f'GTF line {n}: unsupported feature {f[2]}')
        try:s,e=int(f[3]),int(f[4])
        except ValueError:
            m['issues'].add('invalid_coordinate');continue
        if f[0] not in genome or not 1<=s<=e<=len(genome[f[0]]):m['issues'].add('coordinate_out_of_bounds')
        if f[6] not in '+-.':m['issues'].add('invalid_strand')
        item=(f[0],s,e,f[6])
        m['exons' if f[2]=='exon' else 'transcripts'].append(item)
    valid={};issues={}
    for tid,m in models.items():
        exons=sorted(set(m['exons']))
        if len(m['transcripts'])!=1:m['issues'].add('missing_or_multiple_transcript_record')
        if not exons:m['issues'].add('no_exons')
        if len({(x[0],x[3]) for x in m['transcripts']+exons})!=1:m['issues'].add('mixed_contig_or_strand')
        if exons:
            if any(a[2]>=b[1]-1 for a,b in zip(exons,exons[1:])):m['issues'].add('overlapping_or_adjacent_exons')
            if len(m['transcripts'])==1 and (m['transcripts'][0][1],m['transcripts'][0][2])!=(exons[0][1],exons[-1][2]):m['issues'].add('transcript_bounds_disagree')
            if len(exons)>1 and exons[0][3]=='.':m['issues'].add('unstranded_multi_exon')
        if m['issues']:issues[tid]=sorted(m['issues'])
        else:
            m['duplicate_exon_rows']=len(m['exons'])-len(exons)
            m['exons']=exons;valid[tid]=m
    return valid,issues

def derive_introns(valid):
    introns=defaultdict(set)
    for tid,m in valid.items():
        for a,b in zip(m['exons'],m['exons'][1:]):
            introns[(a[0],a[2]+1,b[1]-1,a[3])].add(tid)
    return introns

def load_agp(path,genome,scaffold):
    placements={};end_by_object={};part_by_object={}
    for n,line in enumerate(read_lines(path),1):
        if not line.strip() or line.startswith('#'):continue
        f=line.rstrip().split('\t')
        if len(f)!=9:raise ValueError(f'AGP line {n}: expected 9 fields')
        obj,start,end,part=f[0],int(f[1]),int(f[2]),int(f[3])
        if obj not in scaffold or not 1<=start<=end<=len(scaffold[obj]):raise ValueError('AGP object out of bounds')
        if start!=end_by_object.get(obj,0)+1 or part!=part_by_object.get(obj,0)+1:raise ValueError('AGP object discontinuity')
        end_by_object[obj]=end;part_by_object[obj]=part
        if f[4] in {'N','U'}:
            if int(f[5])!=end-start+1 or set(scaffold[obj][start-1:end])!={'N'}:raise ValueError('AGP gap mismatch')
            continue
        name,cs,ce,strand=f[5],int(f[6]),int(f[7]),f[8]
        if name in placements:raise ValueError('repeated AGP component')
        if name not in genome or not 1<=cs<=ce<=len(genome[name]) or ce-cs!=end-start:raise ValueError('AGP component out of bounds')
        if strand not in {'+','-'}:raise ValueError('unsupported AGP orientation')
        seq=genome[name][cs-1:ce]
        if strand=='-':seq=seq.translate(COMPLEMENT)[::-1]
        if seq!=scaffold[obj][start-1:end]:raise ValueError('AGP sequence mismatch: '+name)
        placements[name]=(obj,start,end,cs,ce,strand)
    if set(placements)!=set(genome):raise ValueError('AGP does not account for every source contig')
    if end_by_object!={k:len(v) for k,v in scaffold.items()}:raise ValueError('AGP does not cover complete scaffold FASTA')
    return placements

def lift(key,placements):
    contig,s,e,strand=key
    obj,start,end,cs,ce,direction=placements[contig]
    if not cs<=s<=e<=ce:raise ValueError('feature outside AGP component')
    if direction=='+':return obj,start+s-cs,start+e-cs,strand
    return obj,start+ce-e,start+ce-s,{'+':'-','-':'+','.':'.'}[strand]

def classify(key,genome,min_intron):
    contig,s,e,strand=key
    seq=genome[contig][s-1:e]
    if strand=='-':seq=seq.translate(COMPLEMENT)[::-1]
    motif=seq[:2]+'-'+seq[-2:]
    if len(seq)<min_intron:return motif,'short_intron'
    if set(seq)-set('ACGT'):return motif,'ambiguous_sequence'
    if strand not in {'+','-'}:return motif,'unknown_strand'
    return motif,'canonical' if motif in CANONICAL else 'noncanonical'

def gff_row(identifier,key,motif,model_count):
    c,s,e,strand=key
    return f'{c}\tStringTie_structural\tintron\t{s}\t{e}\t.\t{strand}\t.\tID={identifier};motif={motif};transcript_models={model_count};evidence=assembled_transcript\n'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for arg in ['gtf','genome','agp','scaffold','out']:ap.add_argument('--'+arg,type=Path,required=True)
    ap.add_argument('--old-introns',type=Path,help='Optional previous invalid reconstruction, for coordinate comparison only')
    ap.add_argument('--min-intron',type=int,default=20,help='Conservative candidate cutoff; shorter records remain in review TSV (default 20 bp)')
    args=ap.parse_args()
    if args.min_intron<4:ap.error('--min-intron must be >=4')
    if args.out.exists():ap.error('output directory already exists; use a new versioned path')
    sources={k:{'path':str(v.resolve()),'sha256':sha256(v),'bytes':v.stat().st_size} for k,v in vars(args).items() if isinstance(v,Path) and k!='out' and v is not None}
    genome=fasta(args.genome);scaffold=fasta(args.scaffold)
    placements=load_agp(args.agp,genome,scaffold)
    valid,issues=parse_gtf(args.gtf,genome)
    introns=derive_introns(valid)
    if not valid or not introns:raise ValueError('no usable transcript/intron structures')
    # All input parsing and sequence checks precede output creation.
    args.out.mkdir(parents=True)
    stats=Counter();motifs=Counter();canonical=set();reasons={}
    with (args.out/'introns.all.tsv').open('w') as table, (args.out/'introns.canonical.contigs.gff3').open('w') as gc, (args.out/'introns.canonical.Kita1.gff3').open('w') as gs:
        table.write('evidence_id\tcontig\tstart\tend\tstrand\tlength\tmotif\tcategory\ttranscript_model_count\ttranscript_ids\tKita1_seqid\tKita1_start\tKita1_end\tKita1_strand\n')
        for f in (gc,gs):f.write('##gff-version 3\n# Structural transcript evidence; not read-support counts or ready-made BRAKER ETP input.\n')
        for n,(key,tids) in enumerate(sorted(introns.items()),1):
            identifier=f'RNAI{n:07d}';new=lift(key,placements)
            motif,category=classify(key,genome,args.min_intron)
            if classify(new,scaffold,args.min_intron)!=(motif,category):raise ValueError('lift changed transcript-oriented motif/category')
            c,s,e,strand=key;stats[category]+=1;motifs[motif]+=1;reasons[key]=category
            table.write('\t'.join(map(str,[identifier,c,s,e,strand,e-s+1,motif,category,len(tids),','.join(quote(t,safe='._-') for t in sorted(tids)),*new]))+'\n')
            if category=='canonical':
                canonical.add(key);gc.write(gff_row(identifier,key,motif,len(tids)));gs.write(gff_row(identifier,new,motif,len(tids)))
    with (args.out/'transcripts.structural.contigs.gtf').open('w') as original, (args.out/'transcripts.structural.Kita1.gtf').open('w') as lifted:
        for tid,m in sorted(valid.items()):
            seen=set()
            for f in m['rows']:
                if tuple(f) in seen:continue
                seen.add(tuple(f));original.write('\t'.join(f)+'\n')
                new=lift((f[0],int(f[3]),int(f[4]),f[6]),placements)
                a=f.copy();a[0],a[3],a[4],a[6]=new[0],str(new[1]),str(new[2]),new[3]
                lifted.write('\t'.join(a)+'\n')
    with (args.out/'transcripts.quarantine.tsv').open('w') as f:
        f.write('transcript_id\treasons\n')
        for tid,rs in sorted(issues.items()):f.write(tid+'\t'+','.join(rs)+'\n')
    old_stats=Counter()
    if args.old_introns:
        with (args.out/'previous_introns.audit.tsv').open('w') as out:
            out.write('old_line\tcontig\tstart\tend\tstrand\tcomparison\n')
            for n,line in enumerate(read_lines(args.old_introns),1):
                if not line.strip() or line.startswith('#'):continue
                a=line.rstrip().split('\t')
                if len(a)!=9 or a[2]!='intron':raise ValueError('unexpected old intron record')
                key=(a[0],int(a[3]),int(a[4]),a[6])
                label='canonical_supported_structure' if key in canonical else 'supported_structure_review' if key in introns else 'absent_from_valid_full_transcripts'
                old_stats[label]+=1;out.write('\t'.join(map(str,[n,*key,label]))+'\n')
    summary={'generated_at_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha256(Path(__file__)),
        'inputs':sources,'min_intron_bp':args.min_intron,
        'transcripts_total':len(valid)+len(issues),'transcripts_structurally_valid':len(valid),'transcripts_quarantined':len(issues),
        'single_exon_transcripts':sum(len(m['exons'])==1 for m in valid.values()),
        'unstranded_single_exon_transcripts':sum(m['exons'][0][3]=='.' for m in valid.values()),
        'exon_occurrences':sum(len(m['exons']) for m in valid.values()),
        'duplicate_exon_rows':sum(m['duplicate_exon_rows'] for m in valid.values()),
        'intron_occurrences':sum(len(m['exons'])-1 for m in valid.values()),
        'unique_introns':len(introns),'intron_categories':dict(stats),'motifs':dict(motifs),'old_intron_comparison':dict(old_stats),
        'validation':{'all_AGP_component_sequences_match':True,'AGP_component_count':len(placements),'all_intron_motifs_preserved_after_lift':True},
        'limitations':['Not raw-read validated. No read counts, expression abundance or coding status inferred.',
        'Canonical sequence and valid transcript structure do not establish biological truth.',
        'Kita1 coordinates follow the existing AGP; they do not validate biological scaffold joins.',
        'GFF3 is structural evidence, not a drop-in BRAKER ETP input; choose pipeline-specific integration separately.']}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.out/'README.md').write_text('''# RNA transcript structural evidence

Primary files: introns.canonical.contigs.gff3 and introns.canonical.Kita1.gff3.
These contain deduplicated introns between consecutive exons in complete, structurally valid StringTie transcript models, with GT-AG/GC-AG/AT-AC motifs, unambiguous intron sequence, and the minimum length recorded in summary.json.

introns.all.tsv retains ALL derived introns, the review category, source transcript IDs and both coordinate systems. Noncanonical or excluded records are review candidates, not automatically false. transcripts.structural.*.gtf preserves structurally valid full transcript models, including single-exon models and models with noncanonical introns; this is not a coding-gene annotation. Any malformed transcript is quarantined whole, without reconnecting surviving exons.

transcript_models counts assembled models, NOT independent reads or samples. Source GTF scores are preserved in the GTF but are not interpreted as abundance or read support. No mult attribute is generated. No BAM/FASTQ was available for this step. GFF3 files must not be passed directly as complete RNA evidence to BRAKER ETP; pipeline-specific integration and weighting remain to be designed.

Every AGP component was checked against both source and scaffold FASTAs, including reverse complement and N gaps. Every intron has matching transcript-oriented sequence category/motif after conversion. No intron joins different input contigs. This checks coordinate consistency, not the biological correctness of the reference-guided assembly.

The success marker is summary.json plus a valid SHA256SUMS.txt. A failure after creating the output directory can leave partial files: never use an incomplete directory; rerun only into a new path. See summary.json for input hashes, exact counts and limitations. Input files are not copied into this output directory.
''')
    files=sorted(p for p in args.out.iterdir() if p.is_file())
    (args.out/'SHA256SUMS.txt').write_text(''.join(sha256(p)+'  '+p.name+'\n' for p in files))
    print(json.dumps({k:summary[k] for k in ['transcripts_total','transcripts_structurally_valid','transcripts_quarantined','unique_introns','intron_categories','old_intron_comparison']},indent=2))

if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as exc:sys.exit('ERROR: '+str(exc))
