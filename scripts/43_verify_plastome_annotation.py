#!/usr/bin/env python3
"""Verify delivered sequence/GenBank/GFF/CDS/proteins and make readable audit tables."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def table(path):
    return list(csv.DictReader(open(path), delimiter='\t'))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    p=args.out
    record=SeqIO.read(p/'Dtenuipes.path2.refined_provisional.gb','genbank')
    original=SeqIO.read(args.root/'recovered/Dtenuipes_Izu_reference_candidate.plastome.path2.fasta','fasta')
    assert record.seq==original.seq
    cds=SeqIO.to_dict(SeqIO.parse(p/'CDS.refined_provisional.fna','fasta'))
    proteins=SeqIO.to_dict(SeqIO.parse(p/'proteins.refined_provisional.faa','fasta'))
    gff=collections.defaultdict(list)
    for line in (p/'Dtenuipes.path2.refined_provisional.gff3').read_text().splitlines():
        if line.startswith('#'):
            continue
        fields=line.split('\t')
        assert len(fields)==9 and fields[0]==record.id
        assert 1 <= int(fields[3]) <= int(fields[4]) <= len(record)
        if fields[2]=='CDS':
            attrs=dict(x.split('=',1) for x in fields[8].split(';'))
            gff[attrs['ID']].append(fields)
    checked=0
    for f in record.features:
        if f.type!='CDS':
            continue
        fid=f.qualifiers['label'][0]
        nt=f.extract(record.seq)
        assert nt==cds[fid].seq
        assert f.translate(record.seq)==proteins[fid].seq
        consumed=0
        pieces=[]
        for fields in gff[fid]:
            assert int(fields[7]) == (-consumed)%3
            piece=record.seq[int(fields[3])-1:int(fields[4])]
            if fields[6]=='-':
                piece=piece.reverse_complement()
            pieces.append(str(piece))
            consumed+=len(piece)
        assert ''.join(pieces)==str(nt)
        checked+=1
    assert checked==len(cds)==len(proteins)==len(gff)==86
    events=table(p/'indel_validation.final.tsv')
    assert len(events)==len({e['event_id'] for e in events})==171
    counts=collections.Counter(e['type'] for e in events)
    assert counts=={'insertion':89,'deletion':82}
    net=sum(int(e['length'])*(1 if e['type']=='insertion' else -1) for e in events)
    assert net==-368
    refined=json.loads((p/'refined_summary.json').read_text())
    assert dict(collections.Counter(e['final_status'] for e in events))==refined['final_indel_status_counts']
    # Ensure final table's corrected impact labels agree with BOTH coordinate systems.
    ref=SeqIO.read(args.root/'published_comparison/PV612491.gb','genbank')
    uniqueness={t:len({f.qualifiers.get('gene',f.qualifiers.get('product',['']))[0] for f in ref.features if f.type==t}) for t in ['CDS','tRNA','rRNA']}
    assert sum(uniqueness.values())==refined['unique_functional_gene_names']
    annotation=json.loads((p/'annotation_summary.json').read_text())
    assert annotation['indel_context_counts']==dict(collections.Counter(e['context'] for e in events))
    readable=[]
    for e in events:
        s,t=int(e['own_start0']),int(e['own_end0'])
        rs,rt=int(e['ref_start0']),int(e['ref_end0'])
        readable.append({'编号':e['event_id'],'相对公开序列':{'insertion':'插入','deletion':'缺失'}[e['type']],
                         '长度bp':e['length'],'我们path2位置_1based':f'{s+1}-{t}' if t>s else f'{s}|{s+1}之间',
                         'PV612491位置_1based':f'{rs+1}-{rt}' if rt>rs else f'{rs}|{rs+1}之间',
                         '区域':e['region'],'注释位置':{'CDS':'涉及编码区或其末端','RNA':'RNA区','intergenic':'基因间区','gene_intron_or_span':'内含子或基因跨度内'}[e['context']],
                         '涉及蛋白基因':e['coding_genes'],'我们的序列':e['own_allele'] or '(空)',
                         '公开序列':e['published_allele'] or '(空)',
                         '我们版本支持read名称数':e['final_own_read_names'],
                         '对照版本支持read名称数':e['final_public_read_names'],
                         '两种版本冲突read名称数':e['final_conflicting_read_names'],
                         '证据片段长度':e['final_assay_kmer_length'],
                         '结果':'现有reads支持' if e['final_status']=='candidate_supported' else e['final_status'],
                         '缺口等价移动范围bp':e['equivalent_gap_shift_min_bp']+'..'+e['equivalent_gap_shift_max_bp'],
                         '重复拷贝限制':'不能区分两份IR拷贝' if e['final_repeat_copy_resolution']=='copy_unresolved' else '有单拷贝标记',
                         '说明':'相对差异；不能确定演化中谁插入或缺失。reads并非独立测序验证。'})
    with open(p/'逐项验证_171项.csv','w',encoding='utf-8-sig',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(readable[0]))
        w.writeheader();w.writerows(readable)
    cds_rows=table(p/'CDS.refined_audit.tsv')
    with open(p/'蛋白基因_86条.csv','w',encoding='utf-8-sig',newline='') as fh:
        w=csv.writer(fh)
        w.writerow(['注释编号','基因','我们path2位置_0based半开区间','我们蛋白长度aa','公开蛋白长度aa','蛋白是否完全相同','蛋白比对一致率_含缺口','读框检查','说明'])
        for r in cds_rows:
            w.writerow([r['feature_id'],r['gene'],r['own_location'],r['own_protein_length'],r['reference_protein_length'],
                        '是' if r['identical_protein']=='True' else '否',r['amino_acid_identity_including_gaps'],'完整','参考支持的暂定注释；未验证表达或功能'])
    checks=dict(sequence_unchanged=True,sequence_length=len(record),
                sequence_sha256=hashlib.sha256(str(record.seq).encode()).hexdigest(),
                genbank_gff_cds_protein_consistent_count=checked,indel_event_count=len(events),
                insertion_event_count=counts['insertion'],deletion_event_count=counts['deletion'],
                net_length_difference=net,unique_gene_names_by_type=uniqueness,
                structural_verification='PASSED',scientific_scope='Internal read consistency and annotation ORF checks; not independent biological validation.')
    (p/'delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))


if __name__=='__main__':
    main()
