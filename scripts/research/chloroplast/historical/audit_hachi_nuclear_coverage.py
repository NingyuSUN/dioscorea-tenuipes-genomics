#!/usr/bin/env python3
"""Hachi first-1M-pair nuclear coding-region coverage; no genotypes or consensus."""
import argparse
from collections import Counter
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import time

import portable_cli as portable
READ_HASHES = ['c7a1b8d568f4b46144f334cb30ba35e272f2d3b9cb2b760fe8d3c96d9a8c63ce',
               '812047cd5b64c383a9a63b1ae6dc08c78567fdf11234f7451104438454194e76']


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''): h.update(b)
    return h.hexdigest()


def collapse(observations):
    """Count each eligible template once; disagreeing mates are not ACGT depth."""
    by_template = {}
    for key, base in observations:
        by_template.setdefault(key, set()).add(base)
    return sum(len(v) == 1 and next(iter(v)) in 'ACGT' for v in by_template.values())


def configure(config_path):
    c=portable.configure(config_path, ('binding','genome','cp','minimap2','samtools','out'),('reads',))
    if len(c['reads'])!=2: raise ValueError('Exactly two ordered read mates are required')
    return c

def plan(config_path):
    c=configure(config_path)
    inputs={k:c[k] for k in ('binding','genome','cp')};inputs.update(read1=c['reads'][0],read2=c['reads'][1],script=Path(__file__))
    names=['COVERAGE_RECEIPT.json','STATUS.json','CDS_TEMPLATE_DEPTH.tsv','Kita2022_nuclear_plus_CP.fa',
           'Kita2022_nuclear_plus_CP.fa.fai','Hachi1.first1M.Kita2022.bam','Hachi1.first1M.Kita2022.bam.bai']
    return portable.contract(c['_config'],inputs,{n:c['out']/n for n in names},
                             {k:c[k] for k in ('minimap2','samtools')},c['out'],markers=143,required_threads=8)

def main():
    a=portable.parser(__doc__).parse_args();checked=portable.preflight(plan(a.config))
    if not a.execute: portable.emit(checked); return
    c=configure(a.config);OUT=c['out'];ROOT=OUT.parent;BINDING=c['binding'];CP=c['cp'];READS=c['reads'];MM=str(c['minimap2']);SAM=str(c['samtools'])
    if c.get('required_hostname') and socket.gethostname()!=c['required_hostname']:
        raise ValueError('Execution hostname differs from the configured resource contract')
    assert not OUT.exists(), 'Inspect existing run; no overwrite or automatic restart'
    assert shutil.disk_usage(portable.existing_parent(ROOT)).free > 8 * 1024**3
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3,) * 2)
    OUT.mkdir(parents=True,exist_ok=False)
    start = time.monotonic()
    state = dict(status='RUNNING', stage='VERIFY_INPUTS', pid=os.getpid(),
                 proc_start_ticks=Path('/proc/self/stat').read_text().split()[21],
                 started_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                 script_sha256=sha(__file__), binding_sha256=sha(BINDING),
                 source_scope='Nonrandom first1000000 paired reads of DRR351751;997147 retained pairs',
                 genotype_calling=False, per_process_address_space_limit_bytes=8*1024**3,
                 commands=[])

    def save(**kw):
        state.update(kw)
        state['updated_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
        t = OUT / 'STATUS.tmp'
        t.write_text(json.dumps(state, indent=2) + '\n')
        t.replace(OUT / 'STATUS.json')

    def run(cmd, label):
        state['commands'].append(cmd)
        save(stage=label)
        with (OUT / (label + '.stdout')).open('w') as o, (OUT / (label + '.stderr')).open('w') as e:
            p = subprocess.Popen(cmd, stdout=o, stderr=e)
            save(child_pid=p.pid)
            code = p.wait()
        save(child_pid=None)
        if code: raise RuntimeError(label + ' exit ' + str(code))

    save()
    try:
        d = json.loads(BINDING.read_text())
        genome = c['genome']  # Relocation preserves the independently checked source genome digest.
        assert sha(genome) == d['sources']['source_genome_sha256']
        for p, expected in zip(READS, READ_HASHES): assert sha(p) == expected
        assert sha(CP) == '1870fe09da36eee2c3048cd52d0aa28f50b788016c151c15e969eb852e454485'
        targets = [x for x in d['target_results'] if (x.get('selected_translation') or {}).get('exact_selected_protein_match')]
        assert len(targets) == 143
        assert all(x['marker'] != '2661at3193' for x in targets)
        save(stage='PREPARE_WHOLE_GENOME_REFERENCE')
        ref = OUT / 'Kita2022_nuclear_plus_CP.fa'
        with ref.open('wb') as f, gzip.open(genome, 'rb') as g:
            shutil.copyfileobj(g, f)
            with gzip.open(CP, 'rb') as c:
                seq = b''.join(x.strip() for x in c if not x.startswith(b'>'))
            f.write(b'\n>CP_Kita2022\n' + seq + b'\n')
        state['reference'] = {'path': str(ref), 'sha256': sha(ref), 'source_genome_sha256': sha(genome),
                              'CP_source_sha256': sha(CP), 'scope': 'Legacy protein-matched nuclear reference plus CP competitor; not latest2026 primary CP comparison'}
        run([SAM, 'faidx', str(ref)], 'FAIDX')
        bam = OUT / 'Hachi1.first1M.Kita2022.bam'
        mc = [MM, '-ax', 'sr', '--secondary=yes', '-t', '4', str(ref), *map(str, READS)]
        sc = [SAM, 'sort', '-@', '2', '-m', '256M', '-o', str(bam), '-']
        state['commands'].extend([mc, sc]); save(stage='MAP_WHOLE_REFERENCE')
        with (OUT / 'minimap.stderr').open('w') as me, (OUT / 'sort.stderr').open('w') as se:
            p = subprocess.Popen(mc, stdout=subprocess.PIPE, stderr=me)
            q = subprocess.Popen(sc, stdin=p.stdout, stderr=se)
            p.stdout.close(); save(mapping_pid=p.pid, sorting_pid=q.pid)
            qc, pc = q.wait(), p.wait()
        assert qc == pc == 0, (pc, qc)
        save(mapping_pid=None, sorting_pid=None)
        run([SAM, 'quickcheck', '-v', str(bam)], 'QUICKCHECK')
        run([SAM, 'index', '-@', '2', str(bam)], 'INDEX')
        import pysam
        rows = []
        unique = set()
        with pysam.AlignmentFile(bam) as b, (OUT/'CDS_TEMPLATE_DEPTH.tsv').open('w') as f:
            f.write('marker\tchrom\tposition_1based\tACGT_template_depth\n')
            for i, target in enumerate(targets):
                t = target['selected_translation']; chrom = t['seqid']
                positions = set()
                for first, last, phase in t['segments']: positions.update(range(first-1,last))
                depths = dict.fromkeys(positions, 0)
                for first, last, phase in t['segments']:
                    for pile in b.pileup(chrom,first-1,last,truncate=True,stepper='nofilter',
                                         max_depth=100000,min_base_quality=0,ignore_overlaps=False,
                                         ignore_orphans=False,compute_baq=False):
                        obs=[]
                        for pr in pile.pileups:
                            rd=pr.alignment
                            if rd.flag & 3844 or rd.mapping_quality < 20 or rd.mapping_quality == 255: continue
                            if pr.is_del or pr.is_refskip or pr.query_position is None: continue
                            j=pr.query_position
                            if rd.query_qualities[j] < 20: continue
                            key=(rd.get_tag('RG') if rd.has_tag('RG') else '',rd.query_name)
                            obs.append((key,rd.query_sequence[j].upper()))
                        depths[pile.reference_pos]=collapse(obs)
                hist=Counter(depths.values()); n=len(depths)
                row=dict(marker=target['marker'],transcript=target['selected_transcript'],CDS_bases=n,
                         mean_ACGT_template_depth=sum(depths.values())/n,
                         depth_histogram=dict(sorted(hist.items())),
                         bases_ge={str(k):sum(v>=k for v in depths.values()) for k in [1,2,5,10,20]})
                rows.append(row)
                for pos,val in sorted(depths.items()):
                    f.write('%s\t%s\t%d\t%d\n'%(target['marker'],chrom,pos+1,val)); unique.add((chrom,pos))
                save(stage='CDS_COVERAGE',completed_markers=i+1,total_markers=143)
        totals={str(k):sum(x['bases_ge'][str(k)] for x in rows) for k in [1,2,5,10,20]}
        result=dict(status='COMPLETE_COVERAGE_ONLY',marker_count=143,original_panel_loci=173,
                    marker_base_observations=sum(x['CDS_bases'] for x in rows),unique_genomic_positions=len(unique),
                    bases_ge=totals,markers_with_at_least_90percent_CDS_depth10=sum(x['bases_ge']['10']>=.9*x['CDS_bases'] for x in rows),
                    loci=rows,filter='BQ20 MAPQ20; unknownMAPQ255/excludedflags3844 removed; RG+QNAME collapse, discordant eligible bases excluded',
                    limitations=['Coverage is not validated diploid genotyping or orthology.','Reference allele bias is not eliminated by whole-genome mapping.',
                                 'Unmarked PCR duplicates may remain; templates are not independent molecules.',
                                 'Nonrandom first1M pairs do not characterize the full public run.','143 coordinate-bound loci; 30 of original173 unassessed.',
                                 'No genotype, consensus, haplotype, or phylogenetic tip generated.'])
        (OUT/'COVERAGE_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
        save(status='COMPLETE_COVERAGE_ONLY',stage='COMPLETE',elapsed_seconds=time.monotonic()-start,
             BAM_sha256=sha(bam),receipt_sha256=sha(OUT/'COVERAGE_RECEIPT.json'),
             max_child_peak_RSS_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024)
    except BaseException as ex:
        save(status='FAILED',error=repr(ex),elapsed_seconds=time.monotonic()-start)
        raise


if __name__ == '__main__': main()
