#!/usr/bin/env python3
"""Validate every alignment column, complete short end alignment, export differences."""
import argparse, csv, hashlib, json, re
from pathlib import Path
from collections import Counter


def fasta(p):return ''.join(l.strip() for l in p.read_text().splitlines() if not l.startswith('>')).upper()


def short_align(t,q):
    assert max(len(t),len(q))<1000
    dp=[[0]*(len(q)+1) for _ in range(len(t)+1)]
    for i in range(len(t)+1):dp[i][0]=-3*i
    for j in range(len(q)+1):dp[0][j]=-3*j
    for i in range(1,len(t)+1):
        for j in range(1,len(q)+1):dp[i][j]=max(dp[i-1][j-1]+(2 if t[i-1]==q[j-1] else -2),dp[i-1][j]-3,dp[i][j-1]-3)
    i,j=len(t),len(q);at=[];aq=[]
    while i or j:
        if i and j and dp[i][j]==dp[i-1][j-1]+(2 if t[i-1]==q[j-1] else -2):at.append(t[i-1]);aq.append(q[j-1]);i-=1;j-=1
        elif i and dp[i][j]==dp[i-1][j]-3:at.append(t[i-1]);aq.append('-');i-=1
        else:at.append('-');aq.append(q[j-1]);j-=1
    return ''.join(at[::-1]),''.join(aq[::-1])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);a=ap.parse_args();root=a.root;p=root/'published_comparison'
    own=root if (root/'Dtenuipes_Izu_reference_candidate.plastome.path1.fasta').exists() else root/'recovered'
    t=fasta(p/'PV612491.fasta');q=fasta(own/'Dtenuipes_Izu_reference_candidate.plastome.path2.fasta');q1=fasta(own/'Dtenuipes_Izu_reference_candidate.plastome.path1.fasta')
    lines=(p/'path2.paf').read_text().splitlines();assert len(lines)==1
    f=lines[0].split('\t');assert f[4]=='+'
    ts,te,qs,qe=map(int,[f[7],f[8],f[2],f[3]])
    cs=next(z[5:] for z in f if z.startswith('cs:Z:'));tok=re.findall(r'=[A-Za-z]+|:[0-9]+|\*[A-Za-z]{2}|\+[A-Za-z]+|-[A-Za-z]+',cs);assert ''.join(tok)==cs
    at=[];aq=[];ti,qi=ts,qs
    for z in tok:
        op=z[0];v=z[1:].upper()
        if op in '=:':
            n=len(v) if op=='=' else int(v);ta=t[ti:ti+n];qa=q[qi:qi+n];assert ta==qa
            if op=='=':assert ta==v
            at.append(ta);aq.append(qa);ti+=n;qi+=n
        elif op=='*':
            assert t[ti]==v[0] and q[qi]==v[1];at.append(v[0]);aq.append(v[1]);ti+=1;qi+=1
        elif op=='+':assert q[qi:qi+len(v)]==v;at.append('-'*len(v));aq.append(v);qi+=len(v)
        elif op=='-':assert t[ti:ti+len(v)]==v;at.append(v);aq.append('-'*len(v));ti+=len(v)
    assert (ti,qi)==(te,qe)
    pt,pq=short_align(t[:ts],q[:qs]);st,sq=short_align(t[te:],q[qe:])
    at=pt+''.join(at)+st;aq=pq+''.join(aq)+sq
    assert at.replace('-','')==t and aq.replace('-','')==q
    counts=Counter();ti=qi=0;events=[];snps=[];last=None
    rc=lambda z:z.translate(str.maketrans('ACGT','TGCA'))[::-1]
    ss,ee=109153,128097
    assert q1[:ss]==q[:ss] and q1[ee:]==q[ee:] and q1[ss:ee]==rc(q[ss:ee])
    for col,(x,y) in enumerate(zip(at,aq)):
        if x==y:counts['matches']+=1;typ='match'
        elif x=='-':counts['insertion_bp']+=1;typ='insertion'
        elif y=='-':counts['deletion_bp']+=1;typ='deletion'
        else:
            counts['substitutions']+=1;typ='substitution'
            pos=ss+ee-1-qi if ss<=qi<ee else qi
            alt=rc(y) if ss<=qi<ee else y
            assert q1[pos]==alt
            snps.append([ti+1,qi+1,x,y,pos+1,alt])
        if typ not in ['match','substitution'] and typ!=last:counts[typ+'_events']+=1
        if typ!='match':events.append([col+1,ti+1 if x!='-' else '',qi+1 if y!='-' else '',typ,x,y])
        if x!='-':ti+=1
        if y!='-':qi+=1
        last=typ
    counts['aligned_columns']=len(at)
    assert len(q)-len(t)==counts['insertion_bp']-counts['deletion_bp']
    with (p/'full_alignment.fasta').open('w') as h:
        for name,seq in [('PV612491.1',at),('Dtenuipes_Izu_candidate_path2',aq)]:
            h.write('>'+name+'\n');h.writelines(seq[i:i+80]+'\n' for i in range(0,len(seq),80))
    for name,headers,rows in [('alignment_differences.tsv',['alignment_col','published_position','own_path2_position','type','published_base','own_base'],events),('substitutions_for_read_validation.tsv',['published_position','own_path2_position','published_base','own_base','native_path1_position','native_path1_base'],snps)]:
        with (p/name).open('w') as h:w=csv.writer(h,delimiter='\t');w.writerow(headers);w.writerows(rows)
    summary={'accession':'PV612491.1','published_length':len(t),'own_length':len(q),'length_difference_own_minus_published':len(q)-len(t),'own_gc_percent':100*(q.count('G')+q.count('C'))/len(q),'published_gc_percent':100*(t.count('G')+t.count('C'))/len(t),'difference_counts':dict(counts),'identity_including_gap_columns_percent':100*counts['matches']/len(at),'identity_paired_bases_only_percent':100*counts['matches']/(counts['matches']+counts['substitutions']),'unaligned_bases':0,'full_alignment_recovers_both_inputs':True,'SSC_flip_only_between_own_paths':True,'method':'minimap2 asm5 cs checked base-by-base; short terminal sequences globally aligned (match +2 mismatch -2 gap -3), no change to input sequences','limitations':['Alignment-derived differences, not yet all read-validated variants.','Indel event counts depend on alignment and gap placement.','Both IR copies included; repeated differences are not independent observations.','Different samples; sequence discrepancy alone does not determine misidentification, hybridization, or assembly error.'],'region_lengths':{'own':{'LSC':83600,'SSC':18944,'IR_each':25553},'paper':{'LSC':83966,'SSC':18910,'IR_each':25571}}}
    (p/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
if __name__=='__main__':main()
