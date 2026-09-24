#!/usr/bin/env python3
"""Audit small terminal boundary shifts and re-examine repeat indels with 101-mers.

Leaves all prior files unchanged. Refined annotations remain provisional.
"""
import argparse
import collections
import csv
import gzip
import importlib.util
import json
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord


def load(name, filename):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).parent/filename)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    ann=load('annotation','40_annotate_plastome.py')
    ind=load('indels','41_validate_plastome_indels.py')
    record=SeqIO.read(args.out/'Dtenuipes.path2.provisional.gb','genbank')
    ref=SeqIO.read(args.root/'published_comparison/PV612491.gb','genbank')
    refinements=[]
    patterns={}
    for f in record.features:
        if f.type!='CDS':
            continue
        old,_,_=ann.cds_audit(f,record.seq)
        if old['orf_status']=='intact_ORF':
            continue
        assert len(f.location.parts)==1 and old['internal_stop_count']==0
        s,e,strand=int(f.location.start),int(f.location.end),f.location.strand
        trial=FeatureLocation(s-30 if strand==-1 else s,e+30 if strand==1 else e,strand=strand)
        extended=trial.extract(record.seq)
        peptide=str(extended[:len(extended)//3*3].translate(table=11))
        stop=peptide.find('*')
        assert stop>=0
        ntlen=(stop+1)*3
        delta=ntlen-len(f)
        assert 0 < delta <= 30
        new=FeatureLocation(e-ntlen,e,strand=strand) if strand==-1 else FeatureLocation(s,s+ntlen,strand=strand)
        gene=f.qualifiers['gene'][0]
        fid=f.qualifiers['label'][0]
        stop_start=int(new.start) if strand==-1 else int(new.end)-3
        pattern=ind.circular(str(record.seq),stop_start-14,stop_start+17)
        patterns[fid]=(pattern,ind.rc(pattern))
        row=dict(feature_id=fid,gene=gene,old_location=str(f.location),new_location=str(new),
                 extension_bp=delta,old_nt_length=len(f),new_nt_length=ntlen,
                 stop_codon=str(new.extract(record.seq)[-3:]),stop_genomic_start0=stop_start,
                 stop_context_31mer=pattern,stop_context_read_names=0,
                 status='pending_read_evidence')
        refinements.append(row)
        for g in record.features:
            if g.type=='gene' and g.qualifiers.get('gene')==f.qualifiers['gene'] and g.location==f.location:
                g.location=new
                g.qualifiers.setdefault('note',[]).append('3-prime boundary refined with associated CDS; see terminal_refinement.tsv.')
        f.location=new
        check,_,prot=ann.cds_audit(f,record.seq)
        assert check['orf_status']=='intact_ORF'
        f.qualifiers['translation']=[prot]
        f.qualifiers['note']=[x for x in f.qualifiers['note'] if not x.startswith('ORF_REVIEW:')]
        f.qualifiers['note'].append('3-prime boundary extended by %d bp to first in-frame stop; candidate subject to read support; source projection retained separately.' % delta)
    support=collections.defaultdict(set)
    for path in sorted((args.root/'inputs').glob('extended_*fq.gz')):
        with gzip.open(path,'rt') as fh:
            while True:
                name=fh.readline().strip()
                if not name:
                    break
                sequence,plus,quality=fh.readline().strip().upper(),fh.readline(),fh.readline().strip()
                assert len(sequence)==len(quality)
                name=name[1:].split()[0]
                if name.endswith('/1') or name.endswith('/2'):
                    name=name[:-2]
                for fid,pair in patterns.items():
                    for pattern in pair:
                        start=sequence.find(pattern)
                        while start>=0:
                            if min(quality[start:start+len(pattern)])>='5':
                                support[fid].add(name)
                            start=sequence.find(pattern,start+1)
    for row in refinements:
        n=len(support[row['feature_id']])
        row['stop_context_read_names']=n
        assert n>=5, 'Boundary refinement lacks read support'
        row['status']='ORF_boundary_refined_with_read_supported_stop_context'
    ann.write_tsv(args.out/'terminal_refinement.tsv',refinements)
    features=[f for f in record.features if f.type=='CDS']
    reference_features=[f for f in ref.features if f.type=='CDS']
    audit=[]
    cds=[]
    proteins=[]
    for f,rf in zip(features,reference_features):
        assert f.qualifiers['gene']==rf.qualifiers['gene']
        info,nt,p=ann.cds_audit(f,record.seq)
        ri,rnt,rp=ann.cds_audit(rf,ref.seq)
        assert info['orf_status']=='intact_ORF'
        fid=f.qualifiers['label'][0]
        gene=f.qualifiers['gene'][0]
        alignment=pairwise2.align.globalms(rp,p,2,-1,-5,-0.5,one_alignment_only=True)[0]
        aa_equal=sum(x==y for x,y in zip(alignment.seqA,alignment.seqB))
        audit.append(dict(feature_id=fid,gene=gene,own_location=str(f.location),
                          own_nt_length=len(nt),reference_nt_length=len(rnt),
                          own_protein_length=len(p),reference_protein_length=len(rp),
                          identical_protein=p==rp,amino_acid_identity_including_gaps=aa_equal/len(alignment.seqA),
                          own_aligned_protein=alignment.seqB,reference_aligned_protein=alignment.seqA,
                          orf_status=info['orf_status']))
        cds.append(SeqRecord(Seq(nt),id=fid,description=gene+' intact_ORF provisional_annotation'))
        proteins.append(SeqRecord(Seq(p),id=fid,description=gene+' intact_ORF provisional_annotation'))
    SeqIO.write(record,args.out/'Dtenuipes.path2.refined_provisional.gb','genbank')
    SeqIO.write(cds,args.out/'CDS.refined_provisional.fna','fasta')
    SeqIO.write(proteins,args.out/'proteins.refined_provisional.faa','fasta')
    ann.write_tsv(args.out/'CDS.refined_audit.tsv',audit)
    with open(args.out/'Dtenuipes.path2.refined_provisional.gff3','w') as fh:
        fh.write('##gff-version 3\n##sequence-region %s 1 %d\n'%(record.id,len(record)))
        fh.write('# Reference projection with audited terminal refinement; compound parts retain biological order.\n')
        for f in record.features:
            if f.type=='source':
                continue
            fid=f.qualifiers['label'][0]
            consumed=0
            offset=int(f.qualifiers.get('codon_start',['1'])[0])-1
            for i,part in enumerate(f.location.parts):
                phase=str(offset if i==0 else (offset-consumed)%3) if f.type=='CDS' else '.'
                attrs='ID='+fid
                for key in ['gene','product']:
                    if key in f.qualifiers:
                        attrs+=';'+key+'='+quote(f.qualifiers[key][0],safe='')
                if 'trans_splicing' in f.qualifiers:
                    attrs+=';exception=trans-splicing'
                fh.write('\t'.join(map(str,[record.id,'audited_reference_projection',f.type,int(part.start)+1,int(part.end),'.','+' if part.strand==1 else '-',phase,attrs]))+'\n')
                consumed+=len(part)
    reread=SeqIO.read(args.out/'Dtenuipes.path2.refined_provisional.gb','genbank')
    assert reread.seq==record.seq
    for f,g in zip(record.features,reread.features):
        assert f.extract(record.seq)==g.extract(reread.seq)
    # Retain initial assay; longer markers are a separate sensitivity analysis.
    initial=list(csv.DictReader(open(args.out/'indel_validation.tsv'),delimiter='\t'))
    unresolved=[dict(e) for e in initial if e['validation_status']!='candidate_supported']
    lookup=ind.build_assays(str(record.seq),unresolved,101)
    supports,metadata=ind.scan_fastq(sorted((args.root/'inputs').glob('extended_*fq.gz')),lookup,unresolved,101)
    for e,s in zip(unresolved,supports):
        ind.classify(e,s)
        e['assay_kmer_length']=101
    if unresolved:
        ann.write_tsv(args.out/'indel_repeat_recheck_101mer.tsv',unresolved)
    final=[]
    rechecked={e['event_id']:e for e in unresolved}
    for e in initial:
        e['initial_31mer_status']=e['validation_status']
        e['final_status']=rechecked.get(e['event_id'],e)['validation_status']
        e['final_assay_kmer_length']=101 if e['event_id'] in rechecked else 31
        best=rechecked.get(e['event_id'],e)
        e['final_own_read_names']=best['own_total_read_names']
        e['final_public_read_names']=best['public_total_read_names']
        e['final_conflicting_read_names']=best['conflicting_read_names']
        e['final_repeat_copy_resolution']=best['repeat_copy_resolution']
        e['final_own_both_junctions_read_names']=best['own_both_junctions_read_names']
        # Impact considers BOTH published and refined annotations, including deleted CDS ends.
        hits=[f for f in record.features if f.type in {'CDS','tRNA','rRNA'} and ann.overlaps_event(f.location,int(e['own_start0']),int(e['own_end0']))]
        e['refined_overlapping_features']=';'.join(f.type+':'+f.qualifiers.get('gene',f.qualifiers.get('product',['']))[0] for f in hits)
        final.append(e)
    ann.write_tsv(args.out/'indel_validation.final.tsv',final)
    summary=dict(intact_CDS_after_boundary_refinement=len(features),reference_CDS_count=len(reference_features),
                 terminal_refinements=refinements,
                 unique_functional_gene_names=113,  # 79 CDS + 30 tRNA + 4 rRNA; source audit.
                 final_indel_status_counts=dict(collections.Counter(e['final_status'] for e in final)),
                 indels_with_equivalent_gap_placements=sum(int(e['equivalent_gap_shift_min_bp'])!=0 or int(e['equivalent_gap_shift_max_bp'])!=0 for e in final),
                 indels_copy_unresolved=sum(e['final_repeat_copy_resolution']=='copy_unresolved' for e in final),
                 maximum_indel_length=max(int(e['length']) for e in final),
                 proteins_identical_to_reference=sum(r['identical_protein'] for r in audit),
                 limitations=['All annotations remain reference-guided and provisional; no independent RNA structure/gene-discovery validation.',
                              'Intact ORFs and read-supported stop contexts do not demonstrate expression or biological function.',
                              'Read evidence is from assembly-recruited FASTQs; IR copies and PCR duplicates can share evidence.',
                              'Equivalent indel gap positions are not separately confirmed events.'])
    (args.out/'refined_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.out/'repeat_recheck_run.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
