#!/usr/bin/env python3
"""New strict/complete/leave-one-out panels; original completed panels read-only."""
import argparse, concurrent.futures, hashlib, json, math, shlex, subprocess
from pathlib import Path

def run(cmd, log):
    with log.open('w') as f:
        f.write(shlex.join(map(str,cmd))+'\n');f.flush()
        subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--plink',required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    source=a.data/'masked17.bcf'
    samples=subprocess.check_output(['bcftools','query','-l',source],text=True).splitlines()
    assert len(samples)==17 and len(set(samples))==17
    assert not {'AIC-F20-S13','AIC-M6-S11'}&set(samples)
    source_hash=hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
    tasks=[('strict17',samples,True,False),('complete17',samples,False,True)]
    tasks += [(f'loo_{i:02d}',[s for s in samples if s!=omit],False,False) for i,omit in enumerate(samples)]
    def prepare(task):
        name,keep,strict,complete=task;d=a.out/name;d.mkdir();(d/'keep.txt').write_text('\n'.join(keep)+'\n')
        min_an=2*(len(keep) if complete else math.ceil(.9*len(keep)))
        commands=[['bcftools','view','-S',str(d/'keep.txt'),str(source),'-Ou']]
        if strict:commands.append(['bcftools','+setGT','-Ou','--','-t','q','-n','.','-i','FMT/DP="." | FMT/GQ="." | FMT/DP<8 | FMT/DP>60 | FMT/GQ<30'])
        commands += [['bcftools','+fill-tags','-Ou','--','-t','AC,AN'],['bcftools','view','-i',f'AN>={min_an} && AC>=2 && (AN-AC)>=2','-Ob','-o',str(d/'candidates.bcf')]]
        with (d/'filter.log').open('w') as f:
            f.write(' | '.join(shlex.join(c) for c in commands)+'\n');f.flush()
            processes=[];prev=None
            for j,c in enumerate(commands):
                pr=subprocess.Popen(c,stdin=prev,stdout=subprocess.PIPE if j<len(commands)-1 else f,stderr=f)
                if prev:prev.close()
                prev=pr.stdout;processes.append(pr)
            codes=[pr.wait() for pr in processes]
            if any(codes):raise RuntimeError((name,codes))
        base=[a.plink,'--allow-extra-chr','--keep-allele-order','--threads','2','--memory','4000']
        run(base+['--bcf',d/'candidates.bcf','--double-id','--set-missing-var-ids','@:#:$1:$2','--make-bed','--out',d/'all_mac2'],d/'import.log')
        run(base+['--bfile',d/'all_mac2','--bp-space','1000','--make-bed','--out',d/'spaced'],d/'space.log')
        run(base+['--bfile',d/'spaced','--indep-pairwise','50kb','1','0.2','--out',d/'ld'],d/'ld.stdout.log')
        run(base+['--bfile',d/'spaced','--extract',d/'ld.prune.in','--make-bed','--out',d/'pruned'],d/'extract.log')
        run(base+['--bfile',d/'pruned','--distance','square','1-ibs','flat-missing','--out',d/'flat'],d/'distance.log')
        ids=[x.split()[1] for x in (d/'pruned.bim').read_text().splitlines()]
        assert len(ids)==len(set(ids)) and '.' not in ids
        assert len(ids)==sum(1 for _ in (d/'ld.prune.in').open())
        actual=[x.split()[1] for x in (d/'pruned.fam').read_text().splitlines()]
        assert actual==keep,(name,actual,keep)
        info={'panel':name,'samples':keep,'omitted':sorted(set(samples)-set(keep)),'strict':strict,'complete':complete,'minimum_AN':min_an,'candidates':sum(1 for _ in (d/'all_mac2.bim').open()),'pruned':sum(1 for _ in (d/'pruned.bim').open())}
        (d/'panel.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info),flush=True)
        return info
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:result=list(pool.map(prepare,tasks))
    after=hashlib.file_digest(source.open('rb'),'sha256').hexdigest();assert after==source_hash
    (a.out/'panels.json').write_text(json.dumps({'source':str(source),'source_sha256':source_hash,'panels':result},indent=2)+'\n')
    (a.out/'COMPLETE').write_text('All 19 new panels completed; source unchanged.\n')
if __name__=='__main__':main()
