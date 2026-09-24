#!/usr/bin/env python3
"""Read-only, streaming exact-allele VCF comparison and sample QC; requires pysam.
Reference is read into memory (~400 MB here); input files and indices are never modified.
Exact matching intentionally does not infer normalized equivalence of multiallelic sites.
--limit restricts coordinate groups, explicitly marking the result incomplete.
"""
import argparse, collections, csv, itertools, json, re, sys, time
from pathlib import Path
import pysam

def alias(s):
    return s[4:-5].replace('_','-') if s.startswith('Dte_') and s.endswith('_L001') else s

def missing(gt): return not gt or any(x is None for x in gt)
def canonical(gt): return tuple(sorted(-1 if x is None else x for x in (gt or ())))
def homref(gt): return bool(gt) and all(x == 0 for x in gt)

def reference(path):
    seqs={}; name=None; parts=[]
    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                if name is not None: seqs[name]=''.join(parts).upper()
                name=l[1:].split()[0]; parts=[]
            else: parts.append(l.strip())
    if name is not None: seqs[name]=''.join(parts).upper()
    return seqs

class Audit:
    def __init__(self, path, ref):
        self.path=path; self.v=pysam.VariantFile(path,threads=2); self.ref=ref
        self.names=list(self.v.header.samples); self.aliases=list(map(alias,self.names))
        if len(set(self.aliases)) != len(self.aliases): raise ValueError('Sample alias collision')
        self.ranks={c:i for i,c in enumerate(self.v.header.contigs)}
        self.d=collections.Counter(); self.filters=collections.Counter(); self.ploidy=collections.Counter()
        self.stats=[collections.Counter() for _ in self.names]; self.previous=None; self.seen_keys=set(); self.lastpos=None
        self.dict_errors=[c for c,x in self.v.header.contigs.items() if c not in ref or x.length != len(ref[c])]
    def key(self,r): return (self.ranks[r.contig], r.pos)
    def consume(self,r):
        self.d['records']+=1; key=self.key(r)
        if self.previous is not None and key<self.previous: raise ValueError(f'Unsorted VCF {self.path} {r.contig}:{r.pos}')
        self.previous=key
        alleles=(r.ref,r.alts)
        if self.lastpos!=key:self.seen_keys=set();self.lastpos=key
        if alleles in self.seen_keys: raise ValueError(f'Duplicate exact variant {self.path} {r.contig}:{r.pos}')
        self.seen_keys.add(alleles)
        self.filters[';'.join(r.filter) or '.']+=1
        if not r.alts:self.d['reference_only']+=1
        elif len(r.alts)==1 and len(r.ref)==1 and len(r.alts[0])==1 and r.alts[0] in 'ACGT':self.d['biallelic_snp']+=1
        if self.ref.get(r.contig,'')[r.start:r.start+len(r.ref)]!=r.ref.upper():self.d['ref_mismatch']+=1
        called_alleles=0
        for s, st in zip(r.samples.values(),self.stats):
            gt=s.get('GT'); dp=s.get('DP');gq=s.get('GQ'); st['records']+=1
            st['missing']+=missing(gt);st['fully_missing']+=not gt or all(x is None for x in gt)
            st['non_diploid']+=bool(gt) and len(gt)!=2
            if gt: called_alleles+=sum(x is not None for x in gt)
            if dp is None: st['dp_missing']+=1
            else:st['dp_observed']+=1;st['dp_sum']+=dp;st['dp_zero']+=dp==0;st['dp_lt5']+=dp<5;st['dp_gt100']+=dp>100
            if gq is None:st['gq_missing']+=1
            else:st['gq_observed']+=1;st['gq_sum']+=gq;st['gq_lt20']+=gq<20
            if not missing(gt):
                st['called']+=1
                hr=homref(gt);het=len(set(gt))>1
                st['homref']+=hr;st['het']+=het;st['homalt']+=not hr and not het
                st['homref_dp_zero_or_missing']+=hr and (dp is None or dp==0)
                st['called_dp_lt5_or_missing']+=dp is None or dp<5
                st['called_gq_lt20_or_missing']+=gq is None or gq<20
                st['called_dp5_gq20']+=dp is not None and dp>=5 and gq is not None and gq>=20
        if 'AN' in r.info:
            self.d['info_an_present']+=1;self.d['info_an_disagrees_gt']+=r.info['AN']!=called_alleles
    def report(self):
        samples=[]
        for name,al,st in zip(self.names,self.aliases,self.stats):
            row={'sample':name,'alias':al,**dict(st)}
            row['missing_fraction']=st['missing']/st['records'] if st['records'] else None
            row['heterozygosity_called']=st['het']/st['called'] if st['called'] else None
            row['mean_dp_observed']=st['dp_sum']/st['dp_observed'] if st['dp_observed'] else None
            samples.append(row)
        return {'path':self.path,'records':self.d['records'],'ref_mismatch':self.d['ref_mismatch'],**dict(self.d),'reference_dictionary_errors':self.dict_errors,'header_contigs':len(self.ranks),'filters':dict(self.filters),'samples':samples}

def groups(a):
    for key, group in itertools.groupby(a.v,key=a.key):
        records=list(group)
        for r in records:a.consume(r)
        yield key,records

def run(left,right,refpath,out,limit=0):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'audit.json').exists():raise ValueError('Output already exists; choose a fresh output directory')
    ref=reference(refpath);a=Audit(left,ref);b=Audit(right,ref)
    if list(a.ranks)!=list(b.ranks):raise ValueError('Different contig orders; comparison requires identical dictionaries')
    if set(a.aliases)!=set(b.aliases):raise ValueError('Different sample sets')
    namepairs=[(x,b.names[b.aliases.index(alias(x))]) for x in a.names]
    cnt=collections.Counter(); per={alias(n):collections.Counter() for n in a.names};examples=[]
    ia=groups(a);ib=groups(b);ra=next(ia,None);rb=next(ib,None);n=0;start=time.time()
    while ra is not None or rb is not None:
        if limit and n>=limit:break
        n+=1
        if rb is None or (ra is not None and ra[0]<rb[0]):cnt['left_only']+=len(ra[1]);ra=next(ia,None)
        elif ra is None or rb[0]<ra[0]:cnt['right_only']+=len(rb[1]);rb=next(ib,None)
        else:
            da={(r.ref,r.alts):r for r in ra[1]};db={(r.ref,r.alts):r for r in rb[1]}
            cnt['left_only']+=len(da.keys()-db.keys());cnt['right_only']+=len(db.keys()-da.keys())
            for k in da.keys() & db.keys():
                x=da[k];y=db[k];cnt['matched_records']+=1
                for ln,rn in namepairs:
                    u=x.samples[ln];v=y.samples[rn];g=u.get('GT');h=v.get('GT');al=alias(ln)
                    if canonical(g)!=canonical(h):
                        cnt['changed_gt']+=1;per[al]['changed_gt']+=1
                        label='missing_to_homref' if missing(g) and homref(h) else 'other_gt_change'
                        cnt[label]+=1;per[al][label]+=1
                        if u.get('DP')==v.get('DP') and u.get('GQ')==v.get('GQ'):cnt[label+'_same_dp_gq']+=1
                        if len(examples)<100:examples.append({'chrom':x.contig,'pos':x.pos,'ref':x.ref,'alt':x.alts,'sample':al,'left_GT':g,'right_GT':h,'left_DP':u.get('DP'),'right_DP':v.get('DP'),'left_GQ':u.get('GQ'),'right_GQ':v.get('GQ'),'change':label})
            ra=next(ia,None);rb=next(ib,None)
        if n%100000==0:print(json.dumps({'coordinate_groups':n,'left_records':a.d['records'],'right_records':b.d['records'],'matched':cnt['matched_records'],'changed':cnt['changed_gt'],'elapsed_seconds':round(time.time()-start)}),file=sys.stderr,flush=True)
    result={'complete':ra is None and rb is None,'limit_coordinate_groups':limit,'scope_note':'Limited runs include one lookahead coordinate group in per-file QC; comparison counts only processed groups.','pysam':pysam.__version__,'reference':refpath,'left':a.report(),'right':b.report(),'comparison':dict(cnt),'changes_by_sample':{k:dict(v) for k,v in per.items()},'examples':examples,'elapsed_seconds':time.time()-start}
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    for label in ['left','right']:
        rows=result[label]['samples']; fields=list(dict.fromkeys(k for row in rows for k in row))
        with (out/f'{label}.samples.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    a.v.close();b.v.close();return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--left',required=True);p.add_argument('--right',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',required=True);p.add_argument('--limit',type=int,default=0)
    args=p.parse_args();r=run(args.left,args.right,args.reference,args.out,args.limit);print(json.dumps({'complete':r['complete'],'comparison':r['comparison']}))
