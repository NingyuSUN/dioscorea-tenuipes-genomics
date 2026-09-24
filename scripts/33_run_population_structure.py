#!/usr/bin/env python3
"""PCA, repeated ADMIXTURE, and descriptive pairwise WC84 FST.

Requires outputs of script 32 and PLINK 1.90b6.21, ADMIXTURE 1.3.0,
VCFtools 0.1.16. Physical coordinates stay on original contigs. Original
PLINK files are never chromosome-renamed. ADMIXTURE gets a disposable
numeric-contig copy AFTER LD pruning, solely for its input parser.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import os

PANELS=('main14','strict14','complete14','supplement17')
SEEDS=(20260914,20260915,20260916)


def run(cmd,log,cwd=None):
    with Path(log).open('w') as f:
        f.write(shlex.join([str(x) for x in cmd])+'\n');f.flush()
        subprocess.run([str(x) for x in cmd],stdout=f,stderr=subprocess.STDOUT,cwd=cwd,check=True)


def nlines(path):
    with Path(path).open() as f:
        return sum(1 for _ in f)


def prepare_admixture(prefix,out):
    out.mkdir(exist_ok=False)
    for suffix in ('.bed','.fam'):
        (out/('panel'+suffix)).symlink_to(prefix.with_suffix(suffix).resolve())
    mapping={}
    with prefix.with_suffix('.bim').open() as src,(out/'panel.bim').open('w') as dst:
        for line in src:
            fields=line.split()
            code=mapping.setdefault(fields[0],len(mapping)+1)
            fields[0]=str(code)
            dst.write('\t'.join(fields)+'\n')
    (out/'contig_code_map.tsv').write_text('original_contig\tadmixture_parser_code\n'+''.join(f'{c}\t{v}\n' for c,v in mapping.items()))
    return out/'panel.bed'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--tools',required=True,type=Path)
    a=parser.parse_args()
    a.data=a.data.resolve();a.out=a.out.resolve()
    assert (a.data/'FILTER_COMPLETE').exists()
    a.out.mkdir(parents=True,exist_ok=False)
    plink=a.tools/'plink';admixture=a.tools/'admixture';vcftools=a.tools/'vcftools'
    base=[plink,'--allow-extra-chr','--keep-allele-order','--threads','2','--memory','5000']

    def panel(name):
        d=a.out/name;d.mkdir()
        # Keep MAC1 in BCF for FST. Only PCA/ADMIXTURE use MAC2.
        run(base+['--bcf',a.data/f'{name}.bcf','--double-id','--mac','2','--make-bed','--out',d/'all_mac2'],d/'import.stdout.log')
        run(base+['--bfile',d/'all_mac2','--missing','--het','--out',d/'qc'],d/'qc.stdout.log')
        run(base+['--bfile',d/'all_mac2','--bp-space','1000','--make-bed','--out',d/'spaced'],d/'space.stdout.log')
        run(base+['--bfile',d/'spaced','--indep-pairwise','50kb','1','0.2','--out',d/'ld'],d/'ld.stdout.log')
        run(base+['--bfile',d/'spaced','--extract',d/'ld.prune.in','--make-bed','--out',d/'pruned'],d/'extract.stdout.log')
        n=nlines(d/'pruned.fam')
        run(base+['--bfile',d/'pruned','--pca',str(n),'header','tabs','--distance','square','ibs','--missing','--out',d/'pca'],d/'pca.stdout.log')
        # Sensitivity: one deterministic hash-selected SNP per represented contig.
        by_contig={}
        for line in (d/'pruned.bim').read_text().splitlines():
            fields=line.split();key=hashlib.sha256(('20260914:'+fields[1]).encode()).hexdigest()
            if fields[0] not in by_contig or key<by_contig[fields[0]][0]:
                by_contig[fields[0]]=(key,fields[1])
        (d/'one_per_contig.txt').write_text(''.join(v[1]+'\n' for v in by_contig.values()))
        run(base+['--bfile',d/'pruned','--extract',d/'one_per_contig.txt','--pca',str(n),'header','tabs','--out',d/'one_contig_pca'],d/'one_contig_pca.stdout.log')
        counts={'samples':n,'MAC2_snps':nlines(d/'all_mac2.bim'),'spaced_snps':nlines(d/'spaced.bim'),'LD_pruned_snps':nlines(d/'pruned.bim'),'one_per_contig_snps':len(by_contig)}
        assert counts['LD_pruned_snps']>0
        print('PCA complete',name,counts,flush=True)
        return name,counts

    counts={}
    with ThreadPoolExecutor(max_workers=3) as pool:
        for name,c in pool.map(panel,PANELS):counts[name]=c
    (a.out/'pca_counts.json').write_text(json.dumps(counts,indent=2)+'\n')
    (a.out/'PCA_COMPLETE').write_text('All four panels and one-per-contig sensitivities completed.\n')

    jobs=[]
    for name in ('main14','strict14','supplement17'):
        adm=a.out/name/'admixture';adm.mkdir()
        for k in range(1,5):
            for seed in SEEDS:
                d=adm/f'K{k}_seed{seed}'
                bed=prepare_admixture(a.out/name/'pruned',d)
                jobs.append((name,k,seed,d,bed))

    def admix(job):
        name,k,seed,d,bed=job
        run([admixture,'--cv=5','-s',str(seed),'-j2',bed,str(k)],d/'run.log',d)
        assert (d/f'panel.{k}.Q').exists()
        print('ADMIXTURE complete',name,k,seed,flush=True)

    def fst(job):
        name,g1,g2=job;d=a.out/name/'fst';d.mkdir(exist_ok=True)
        prefix=d/f'{g1}_{g2}'
        commands=[['bcftools','view',str(a.data/f'{name}.bcf'),'-Ov'],
            [str(vcftools),'--vcf','-','--weir-fst-pop',str(a.data/f'{g1}.samples.txt'),
             '--weir-fst-pop',str(a.data/f'{g2}.samples.txt'),'--fst-window-size','100000',
             '--fst-window-step','100000','--out',str(prefix)]]
        with prefix.with_suffix('.stdout.log').open('w') as log:
            log.write(' | '.join(shlex.join(c) for c in commands)+'\n');log.flush()
            p=subprocess.Popen(commands[0],stdout=subprocess.PIPE,stderr=log)
            q=subprocess.Popen(commands[1],stdin=p.stdout,stdout=log,stderr=log)
            p.stdout.close()
            codes=(p.wait(),q.wait())
            if any(codes):raise RuntimeError(f'FST failed {job}: {codes}')
        print('FST complete',name,g1,g2,flush=True)

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures=[pool.submit(admix,j) for j in jobs]
        futures += [pool.submit(fst,(name,g1,g2)) for name in ('main14','strict14') for g1,g2 in [('AIC','HCO'),('AIC','IZU'),('HCO','IZU')]]
        for f in futures:f.result()
    (a.out/'STRUCTURE_COMPLETE').write_text('PCA, all 36 ADMIXTURE runs, and six pairwise FST estimates completed.\n')


if __name__=='__main__':main()
