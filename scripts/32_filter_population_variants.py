#!/usr/bin/env python3
"""Create an explicitly excluded, genotype-masked population analysis version.

Original VCF and original contig coordinates are preserved. Requires bcftools
1.13 with setGT/fill-tags. No missing-to-reference imputation is performed.
"""
import argparse
import csv
import json
from pathlib import Path
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor

EXCLUDED = {'AIC-M6-S11', 'AIC-F20-S13'}


def run_pipe(commands, output, log):
    """Run a checked pipe; publish output only if every process succeeds."""
    output, log = Path(output), Path(log)
    partial = output.with_name(output.name + '.partial')
    if output.exists():
        raise FileExistsError(output)
    processes = []
    with log.open('w') as err, partial.open('wb') as dest:
        err.write(' | '.join(shlex.join(c) for c in commands) + '\n')
        err.flush()
        previous = None
        for i, cmd in enumerate(commands):
            proc = subprocess.Popen(cmd, stdin=previous, stdout=dest if i == len(commands)-1 else subprocess.PIPE, stderr=err)
            if previous is not None:
                previous.close()
            previous = proc.stdout
            processes.append(proc)
        codes = [p.wait() for p in processes]
    if any(codes):
        raise RuntimeError(f'Pipeline failed {codes}; see {log}')
    partial.rename(output)


def mask_commands(source, dp, gq, maxdp):
    # Single | operates per sample; || could select samples across a whole site.
    expr = f'FMT/DP="." | FMT/GQ="." | FMT/DP<{dp} | FMT/DP>{maxdp} | FMT/GQ<{gq}'
    return [
        ['bcftools','+setGT',str(source),'-Ou','--','-t','.','-n','.'],
        ['bcftools','+setGT','-Ou','--','-t','q','-n','.','-i',expr],
        ['bcftools','annotate','-x','INFO,^FORMAT/GT,FORMAT/DP,FORMAT/GQ','-Ob','--threads','4'],
    ]


def panel_commands(source, samples, missing):
    return [
        ['bcftools','view','-S',str(samples),str(source),'-Ou'],
        ['bcftools','+fill-tags','-Ou','--','-t','AC,AN,AF,NS,F_MISSING'],
        ['bcftools','view','-i',f'F_MISSING<={missing}','-c','1:minor','-Ou'],
        ['bcftools','annotate','--set-id','%CHROM:%POS:%REF:%FIRST_ALT','-Ob','--threads','4'],
    ]


def index(path):
    subprocess.run(['bcftools','index',str(path)],check=True)
    return int(subprocess.check_output(['bcftools','index','-n',str(path)],text=True))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',required=True,type=Path)
    p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
    ids=[r['sample'] for r in rows]
    observed=subprocess.check_output(['bcftools','query','-l',str(a.vcf)],text=True).splitlines()
    assert ids==observed and len(set(ids))==19
    assert EXCLUDED <= set(ids)
    own=[r for r in rows if r['sample'] not in EXCLUDED and r['group'] in ('AIC','HCO','IZU')]
    all17=[r for r in rows if r['sample'] not in EXCLUDED]
    assert len(own)==14 and len(all17)==17
    assert {g:sum(r['group']==g for r in own) for g in ('AIC','HCO','IZU')}=={'AIC':4,'HCO':4,'IZU':6}
    for name, selected in [('main14',own),('supplement17',all17)]:
        (a.out/f'{name}.samples.txt').write_text(''.join(r['sample']+'\n' for r in selected))
        with (a.out/f'{name}.manifest.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
    for g in ('AIC','HCO','IZU'):
        (a.out/f'{g}.samples.txt').write_text(''.join(r['sample']+'\n' for r in own if r['group']==g))
    (a.out/'excluded.tsv').write_text('sample\treason\nAIC-M6-S11\tvery_low_supported_genotype_fraction;user_authorized_exclusion\nAIC-F20-S13\textreme_heterozygosity;user_authorized_exclusion\n')
    masked=a.out/'masked17.bcf'
    cmds=[['bcftools','view','-S',str(a.out/'supplement17.samples.txt'),'-f','PASS','-m2','-M2','-v','snps',str(a.vcf),'-Ou']]
    cmds+=mask_commands('-',5,20,100)
    print('Masking 17 samples: DP 5-100, GQ >=20',flush=True)
    run_pipe(cmds,masked,a.out/'mask.log')
    counts={'masked17':index(masked)}
    print('Masking complete',counts,flush=True)

    def panel(name,sample_set,missing,strict=False):
        output=a.out/f'{name}.bcf'
        cmds=panel_commands(masked,a.out/f'{sample_set}.samples.txt',missing)
        if strict:
            cmds=mask_commands(masked,8,30,60)[:-1]+panel_commands('-',a.out/f'{sample_set}.samples.txt',missing)
        run_pipe(cmds,output,a.out/f'{name}.filter.log')
        n=index(output)
        print(name,n,flush=True)
        return name,n

    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(panel,*x) for x in [
            ('main14','main14',.1),('strict14','main14',.1,True),
            ('complete14','main14',0),('supplement17','supplement17',.1)]]
        for future in futures:
            name,n=future.result();counts[name]=n
    config={'source':str(a.vcf),'source_sha256_from_completed_audit':'43472edd28b406493d5745d71a1d7d403d2086d810e59d4f50194858a176adb9',
        'excluded':sorted(EXCLUDED),'coordinates':'original contigs; no relabeling or liftover',
        'primary':{'min_DP':5,'max_DP':100,'min_GQ':20,'max_missing':.1,'min_MAC':1},
        'strict':{'min_DP':8,'max_DP':60,'min_GQ':30,'max_missing':.1,'min_MAC':1},
        'complete14':'primary genotype thresholds; zero missing genotypes',
        'PCA_additional_filter':'MAC>=2; within-contig LD pruning applied by script 33',
        'counts':counts,'bcftools':subprocess.check_output(['bcftools','--version'],text=True)}
    (a.out/'filter_summary.json').write_text(json.dumps(config,indent=2)+'\n')
    (a.out/'FILTER_COMPLETE').write_text('All pipelines and indexes completed successfully.\n')


if __name__=='__main__':
    main()
