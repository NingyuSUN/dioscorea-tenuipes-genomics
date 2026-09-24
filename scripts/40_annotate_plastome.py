#!/usr/bin/env python3
"""Project PV612491.1 annotations through a verified alignment; do not repair sequence.

Biopython required. Coordinates in audit TSVs are zero-based, half-open, unless
explicitly named otherwise. GenBank/GFF follow their own one-based conventions.
"""
import argparse
import collections
import copy
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


def write_tsv(path, rows):
    with open(path, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def alignment_map(a, b):
    assert len(a) == len(b)
    mapping, events = {}, []
    rp = qp = col = 0
    while col < len(a):
        x, y = a[col], b[col]
        assert x != '-' or y != '-'
        if x == '-' or y == '-':
            kind = 'insertion' if x == '-' else 'deletion'
            start, r0, q0 = col, rp, qp
            while col < len(a) and ((a[col] == '-') if kind == 'insertion' else (b[col] == '-')):
                if a[col] != '-':
                    mapping[rp] = None
                    rp += 1
                if b[col] != '-':
                    qp += 1
                col += 1
            events.append(dict(event_id='INDEL_%03d' % (len(events)+1), type=kind,
                               ref_start0=r0, ref_end0=rp, own_start0=q0, own_end0=qp,
                               published_allele=a[start:col].replace('-', ''),
                               own_allele=b[start:col].replace('-', ''), length=col-start))
        else:
            mapping[rp] = qp
            rp += 1
            qp += 1
            col += 1
    return mapping, events


def lift_location(loc, mapping):
    """Retain biological part order/strand, including trans-spliced rps12.

    Internal insertions are included, while insertions exactly outside an exon
    boundary are excluded. Deleted boundary bases are reported for review.
    """
    parts, issues = [], []
    for p in loc.parts:
        s, e = int(p.start), int(p.end)
        mapped = [mapping[i] for i in range(s, e) if mapping[i] is not None]
        if not mapped:
            raise ValueError('entire feature part deleted')
        if mapping[s] is None or mapping[e-1] is None:
            issues.append('deleted_feature_boundary')
        parts.append(FeatureLocation(mapped[0], mapped[-1]+1, strand=p.strand))
    return (parts[0] if len(parts) == 1 else CompoundLocation(parts, operator=loc.operator)), issues


def cds_audit(feature, sequence):
    nt = feature.extract(sequence)[int(feature.qualifiers.get('codon_start', ['1'])[0])-1:]
    table = int(feature.qualifiers.get('transl_table', ['11'])[0])
    protein = str(nt[:len(nt)//3*3].translate(table=table))
    try:
        complete = str(feature.translate(sequence))
        status, reason = 'intact_ORF', ''
    except Exception as exc:
        complete = protein.rstrip('*')
        status, reason = 'needs_review', str(exc)
    return dict(cds_nt_length=len(nt), frame_remainder=len(nt) % 3,
                internal_stop_count=protein[:-1].count('*'), start_codon=str(nt[:3]),
                terminal_codon=str(nt[-3:]), orf_status=status, orf_issue=reason), str(nt), complete


def overlaps_event(location, s, e):
    if s == e:
        return any(int(p.start) < s < int(p.end) for p in location.parts)
    return any(max(int(p.start), s) < min(int(p.end), e) for p in location.parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    refpath = args.root/'published_comparison/PV612491.gb'
    ref = SeqIO.read(refpath, 'genbank')
    aln = list(SeqIO.parse(args.root/'published_comparison/full_alignment.fasta', 'fasta'))
    assert len(aln) == 2
    a, b = map(lambda r: str(r.seq).upper(), aln)
    ownseq = b.replace('-', '')
    assert a.replace('-', '') == str(ref.seq).upper()
    ownpath = args.root/'recovered/Dtenuipes_Izu_reference_candidate.plastome.path2.fasta'
    if not ownpath.exists():
        ownpath = args.root/ownpath.name
    assert str(SeqIO.read(ownpath, 'fasta').seq).upper() == ownseq
    mapping, events = alignment_map(a, b)
    assert len(events) == 171
    record = SeqRecord(Seq(ownseq), id='Dtenuipes_plastome_path2', name='Dtenuipes_cp',
                       description='Dioscorea tenuipes chloroplast; provisional reference-guided annotation; path2 coordinates')
    record.annotations = dict(molecule_type='DNA', topology='circular', organism='Dioscorea tenuipes')
    record.features = [SeqFeature(FeatureLocation(0, len(record)), type='source', qualifiers={
        'organism': ['Dioscorea tenuipes'], 'organelle': ['plastid:chloroplast'],
        'mol_type': ['genomic DNA'], 'note': ['Recovered assembly; sample provenance pending final binding; not the PV612491.1 specimen.']})]
    audit, cds, proteins, rnas, refcounts = [], [], [], [], collections.Counter()
    features_for_overlap = []
    for index, feature in enumerate(ref.features):
        if feature.type == 'source':
            continue
        fid = 'DTCp_%04d' % index
        gene = feature.qualifiers.get('gene', [''])[0]
        loc, issues = lift_location(feature.location, mapping)
        q = {k: copy.deepcopy(v) for k, v in feature.qualifiers.items()
             if k in {'gene', 'product', 'codon_start', 'transl_table', 'trans_splicing', 'pseudo'}}
        q['note'] = ['Provisional projection from PV612491.1; not independently curated. Source feature index %d.' % index]
        q['inference'] = ['similar to DNA sequence:INSD:PV612491.1']
        q['label'] = [fid]
        if issues:
            q['note'].extend(issues)
        f = SeqFeature(loc, type=feature.type, qualifiers=q)
        row = dict(feature_id=fid, feature_type=f.type, gene=gene,
                   reference_location=str(feature.location), own_path2_location=str(loc),
                   reference_nt_length=len(feature), own_nt_length=len(f),
                   length_change=len(f)-len(feature), boundary_issue=';'.join(issues),
                   cds_nt_length='', frame_remainder='', internal_stop_count='', start_codon='',
                   terminal_codon='', orf_status='', orf_issue='', reference_orf_status='',
                   source_translation_matches='', identical_protein_to_reference='')
        if feature.type == 'CDS':
            info, nt, prot = cds_audit(f, record.seq)
            refinfo, refnt, refprot = cds_audit(feature, ref.seq)
            row.update(info)
            row['reference_orf_status'] = refinfo['orf_status']
            row['source_translation_matches'] = refprot == feature.qualifiers.get('translation', [''])[0]
            row['identical_protein_to_reference'] = prot == refprot
            cds.append(SeqRecord(Seq(nt), id=fid, description=gene+' '+info['orf_status']))
            proteins.append(SeqRecord(Seq(prot), id=fid, description=gene+' '+info['orf_status']))
            if info['orf_status'] == 'intact_ORF':
                q['translation'] = [prot]
            else:
                q['note'].append('ORF_REVIEW: '+info['orf_issue'])
        elif feature.type in {'tRNA', 'rRNA'}:
            rnas.append(SeqRecord(f.extract(record.seq), id=fid, description=gene))
        if feature.type in {'CDS', 'tRNA', 'rRNA'}:
            refcounts[feature.type] += 1
        audit.append(row)
        record.features.append(f)
        features_for_overlap.append((fid, f))
    SeqIO.write(record, args.out/'Dtenuipes.path2.provisional.gb', 'genbank')
    SeqIO.write(record, args.out/'Dtenuipes.path2.fasta', 'fasta')
    SeqIO.write(cds, args.out/'CDS.all_provisional.fna', 'fasta')
    SeqIO.write(proteins, args.out/'proteins.all_provisional.faa', 'fasta')
    SeqIO.write([p for p in proteins if p.description.endswith('intact_ORF')], args.out/'proteins.intact_ORF.faa', 'fasta')
    SeqIO.write(rnas, args.out/'RNAs.provisional.fna', 'fasta')
    write_tsv(args.out/'annotation_audit.tsv', audit)
    with open(args.out/'Dtenuipes.path2.provisional.gff3', 'w') as fh:
        fh.write('##gff-version 3\n##sequence-region %s 1 %d\n' % (record.id, len(record)))
        fh.write('# Provisional reference projection; CDS/RNA features have independent IDs; compound parts share an ID.\n')
        for fid, f in features_for_overlap:
            consumed = 0
            offset = int(f.qualifiers.get('codon_start', ['1'])[0])-1
            for part_index, p in enumerate(f.location.parts):
                phase = str(offset if part_index == 0 else (offset-consumed) % 3) if f.type == 'CDS' else '.'
                attrs = 'ID='+fid
                for k in ('gene', 'product'):
                    if k in f.qualifiers:
                        attrs += ';'+k+'='+quote(f.qualifiers[k][0], safe='')
                if 'trans_splicing' in f.qualifiers:
                    attrs += ';exception=trans-splicing'
                fh.write('\t'.join(map(str, [record.id, 'reference_projection', f.type, int(p.start)+1, int(p.end), '.', '+' if p.strand == 1 else '-', phase, attrs]))+'\n')
                consumed += len(p)
    for event in events:
        s, e = event['own_start0'], event['own_end0']
        hit = [(fid, f) for fid, f in features_for_overlap if overlaps_event(f.location, s, e)]
        event['overlapping_features'] = ';'.join(fid+':'+f.type+':'+f.qualifiers.get('gene', [''])[0] for fid, f in hit)
        refhit = [f for f in ref.features if f.type != 'source' and overlaps_event(f.location, event['ref_start0'], event['ref_end0'])]
        event['reference_overlapping_features'] = ';'.join(f.type+':'+f.qualifiers.get('gene', f.qualifiers.get('product',['']))[0] for f in refhit)
        event['coding_genes'] = ';'.join(sorted({f.qualifiers.get('gene', [''])[0] for f in [x[1] for x in hit]+refhit if f.type == 'CDS'}))
        event['context'] = 'CDS' if event['coding_genes'] else ('RNA' if any(f.type in {'tRNA','rRNA'} for _,f in hit) else ('gene_intron_or_span' if any(f.type=='gene' for _,f in hit) else 'intergenic'))
        event['region'] = next(label for x,y,label in [(0,83600,'LSC'),(83600,109153,'IR1'),(109153,128097,'SSC'),(128097,153650,'IR2')] if x <= s < y)
        # A gap at an exon boundary may have no interval overlap; retain an explicit flag.
        event['near_feature_boundary_2bp'] = ';'.join(sorted({fid for fid,f in features_for_overlap if f.type in {'CDS','tRNA','rRNA'} and any(min(abs(s-int(p.start)),abs(s-int(p.end)),abs(e-int(p.start)),abs(e-int(p.end))) <= 2 for p in f.location.parts)}))
    write_tsv(args.out/'indel_events.annotated.tsv', events)
    unique = {t: len({f.qualifiers.get('gene', f.qualifiers.get('product',['']))[0] for f in ref.features if f.type==t}) for t in refcounts}
    summary = dict(sequence_sha256=hashlib.sha256(ownseq.encode()).hexdigest(), coordinate_basis='path2; TSV zero-based half-open; GenBank/GFF one-based',
                   annotation_source='PV612491.1', reference_file_sha256=hashlib.sha256(refpath.read_bytes()).hexdigest(),
                   transferred_feature_counts=dict(collections.Counter(f.type for f in record.features)),
                   functional_feature_counts=dict(refcounts), unique_gene_names_by_type=unique,
                   intact_CDS=sum(r['orf_status']=='intact_ORF' for r in audit),
                   review_CDS=[r for r in audit if r['orf_status']=='needs_review'],
                   reference_translation_mismatches=[r['feature_id'] for r in audit if r['source_translation_matches'] is False],
                   boundary_issues=[r['feature_id'] for r in audit if r['boundary_issue']],
                   indel_context_counts=dict(collections.Counter(e['context'] for e in events)),
                   limitation='Reference-guided provisional annotation, no de novo gene discovery or independent RNA evidence. Missing reference features cannot be discovered by transfer.')
    (args.out/'annotation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    # GenBank round-trip must preserve each feature's extracted sequence and part order.
    reread = SeqIO.read(args.out/'Dtenuipes.path2.provisional.gb', 'genbank')
    assert len(reread.features) == len(record.features)
    for old,new in zip(record.features, reread.features):
        assert str(old.extract(record.seq)) == str(new.extract(reread.seq))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
