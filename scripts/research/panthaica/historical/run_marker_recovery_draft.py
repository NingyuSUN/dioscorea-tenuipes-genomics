#!/usr/bin/env python3
"""Prepared, not production-run: BUSCO 5.8.3 transcriptome + fixed 173-locus candidates.
Default dry-run does not read production FASTA, write files, or start BUSCO.
No assembly, alignment, tree, candidate collapsing, or topology-driven selection.
"""
import argparse
from collections import Counter, defaultdict
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT=None
ENV=None
PORTABLE=None
TOOLS={}
EVIDENCE_SHA='de212e8f6f5502f5d75e68f8d509433c7c9e9cc6638cae304b6ce16ef456d660'
AA=set('ACDEFGHIKLMNPQRSTVWY')


try:
    from .portable_inputs import Config, report, require_preflight, run_directory
except ImportError:
    from portable_inputs import Config, report, require_preflight, run_directory

def now():return dt.datetime.now(dt.timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
    return h.hexdigest()

def save(path,data):
    temp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    with temp.open('w') as f:json.dump(data,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temp,path)

def evidence():
    row=PORTABLE.input('marker_evidence',EVIDENCE_SHA)
    p=Path(row['path'])
    if sha(p)!=row['sha256']:raise ValueError('Readiness evidence hash mismatch')
    e=PORTABLE.relocated(json.loads(p.read_text()))
    assert e['baseline']['genes']==173 and len(e['baseline']['markers'])==173
    assert e['baseline']['main_sha256']=='982addd02b6a6d0a188aa0de9b8e688b95bec926703094ab7f9238a8a91cc01e'
    return e

def check_lineage(root,e):
    expected={x['relative_path'] for x in e['lineage']['files']}
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if actual!=expected:raise ValueError('Lineage file set changed')
    for row in e['lineage']['files']:
        p=root/row['relative_path']
        if p.is_symlink() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:
            raise ValueError('Frozen lineage mismatch: '+str(p))

def fasta(path):
    records={};ident=None;seq=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if ident is not None:records[ident]=''.join(seq).upper()
                ident=line[1:].split()[0]
                if ident in records:raise ValueError('Duplicate FASTA ID '+ident)
                seq=[]
            elif line.strip():
                if ident is None:raise ValueError('Sequence before FASTA header')
                seq.append(line.strip())
    if ident is not None:records[ident]=''.join(seq).upper()
    return records

def strip_one_terminal_stop(sequence):
    return sequence[:-1] if sequence.endswith('*') else sequence

def quality_reason(sequence,median_length):
    s=strip_one_terminal_stop(sequence)
    if '*' in s:return 'internal_stop'
    if len(s)<50:return 'shorter_than_50_aa'
    if sum(c in AA for c in s)/len(s)<0.95:return 'canonical_aa_fraction_below_0.95'
    if not 0.70<=len(s)/median_length<=1.50:return 'length_ratio_outside_0.70_to_1.50'
    return None

def extract(run,e,out):
    # The server BUSCO Python has Biopython; no new dependency is installed.
    from Bio.Seq import Seq
    table=run/'run_embryophyta_odb10/full_table.tsv'
    entries=defaultdict(list)
    for line in table.read_text().splitlines():
        if line and not line.startswith('#'):
            fields=line.split('\t');entries[fields[0]].append(fields)
    if len(entries)!=1614:raise ValueError('Expected exactly 1614 unique BUSCO marker statuses')
    states={g:rows[0][1] for g,rows in entries.items()}
    for g,rows in entries.items():
        if len(rows)>1 and any(r[1]!='Duplicated' for r in rows):raise ValueError('Unexpected repeated marker rows: '+g)
    folder=run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences'
    expected={g for g,s in states.items() if s=='Complete'}
    if {p.stem for p in folder.glob('*.faa')}!=expected:raise ValueError('Single-copy FASTA/table set mismatch')
    accepted=out/'candidates';accepted.mkdir()
    audit=[]
    for ref in e['baseline']['markers']:
        g=ref['marker'];rows=entries[g]
        row={'marker':g,'busco_status':states[g],'source_rows':rows,'baseline_raw_median_aa':ref['baseline_raw_median_aa'],'reason':None}
        if states[g]!='Complete':row['reason']='not_complete_single_copy'
        else:
            faa=folder/(g+'.faa');fna=folder/(g+'.fna')
            aa=fasta(faa);nt=fasta(fna) if fna.exists() else {}
            if len(aa)!=1 or len(nt)!=1:row['reason']='not_exactly_one_bound_protein_and_cds'
            else:
                ident,seq=next(iter(aa.items()));dna=nt.get(ident)
                row.update(id=ident,aa_length=len(strip_one_terminal_stop(seq)),faa_path=str(faa),faa_sha256=sha(faa),fna_path=str(fna),fna_sha256=sha(fna),length_ratio=len(strip_one_terminal_stop(seq))/ref['baseline_raw_median_aa'])
                if ident!=rows[0][2] or dna is None:row['reason']='table_protein_cds_id_mismatch'
                elif len(dna)%3:row['reason']='cds_not_divisible_by_three'
                else:
                    translated=str(Seq(dna).translate(table=1))
                    if strip_one_terminal_stop(translated)!=strip_one_terminal_stop(seq):row['reason']='cds_translation_not_exactly_bound_to_protein'
                    elif '*' in translated[:-1]:row['reason']='internal_stop_in_cds_translation'
                    else:row['reason']=quality_reason(seq,ref['baseline_raw_median_aa'])
                if row['reason'] is None:
                    # Sequence identity is retained in the audit; one fixed taxon label in FASTA.
                    clean=''.join(c if c in AA else 'X' for c in strip_one_terminal_stop(seq))
                    (accepted/(g+'.faa')).write_text('>Dioscorea_panthaica_DIPA1_transcriptome\n'+clean+'\n')
        audit.append(row)
    report={'status':'CANDIDATES_EXTRACTED_ALIGNMENT_AND_ORTHOLOGY_REVIEW_PENDING','fixed_loci_denominator':173,'accepted_candidates':sum(x['reason'] is None for x in audit),'status_counts_all_1614':dict(Counter(states.values())),'full_table_sha256':sha(table),'candidate_audit':audit,'limits':['BUSCO single-copy is a call within this leaf transcript assembly, not proven one genomic copy.','All BUSCO raw/intermediate/multicopy/fragmented outputs remain in the run directory.','No alignment, contamination test, read support check, paralog resolution or phylogeny has been performed.']}
    save(out/'CANDIDATE_RECEIPT.json',report)
    return report

def validate_assembly_scope(assembly):
    if assembly.get('input_scope')!='FULL_CRR515167' or assembly.get('source_raw_pairs')!=28665951 or assembly.get('formal_marker_eligible') is not True:
        raise ValueError('Pilot/prefix/unknown assembly cannot enter formal marker recovery')

def execute(a):
    if ROOT is None or ROOT.resolve()!=ROOT or not ROOT.is_dir():raise RuntimeError('Configured task root is not an existing absolute directory')
    e=evidence()
    transcript=Path(a.transcripts).resolve();assembly_receipt=Path(a.assembly_receipt).resolve()
    if sha(transcript)!=a.transcripts_sha256 or sha(assembly_receipt)!=a.assembly_receipt_sha256:raise ValueError('Assembly input/provenance hash mismatch')
    assembly=PORTABLE.relocated(json.loads(assembly_receipt.read_text()))
    if assembly.get('status')!='COMPLETE':raise ValueError('Normalized assembly receipt must have status COMPLETE')
    validate_assembly_scope(assembly)
    binding=assembly.get('transcripts',{})
    if binding.get('sha256')!=a.transcripts_sha256 or Path(binding.get('path','')).resolve()!=transcript:
        raise ValueError('Assembly receipt must bind transcripts.path and transcripts.sha256')
    downloaded=json.loads(Path(PORTABLE.input('download_receipt')['path']).read_text())
    if downloaded.get('status')!='COMPLETE' or downloaded.get('pairs_done')!=28665951:raise ValueError('Full raw input validation receipt is required')
    for item in e['busco_sources'].values():
        if sha(item['path'])!=item['sha256']:raise ValueError('BUSCO installed source changed')
    for key in ['summary','main']:
        if sha(e['baseline'][key+'_path'])!=e['baseline'][key+'_sha256']:raise ValueError('Baseline matrix/summary changed')
    original=Path(e['lineage']['path']);check_lineage(original,e)
    for row in e['baseline']['markers']:
        if sha(row['raw_path'])!=row['raw_sha256']:raise ValueError('Baseline raw marker changed')
    if shutil.disk_usage(ROOT).free<8*1024**3:raise RuntimeError('Require 8 GiB free for BUSCO and copied lineage')
    environment=os.environ.copy();environment.update(PATH=os.pathsep.join(dict.fromkeys([str(Path(TOOLS['busco']).parent),str(Path(TOOLS['metaeuk']).parent)]))+os.pathsep+environment.get('PATH',''),PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for command,wanted in [([TOOLS['busco'],'--version'],'BUSCO 5.8.3'),([TOOLS['metaeuk'],'version'],'7.bba0d80')]:
        p=subprocess.run(command,capture_output=True,text=True,check=True,env=environment)
        if wanted not in p.stdout+p.stderr:raise ValueError('Tool version drift')
    with (ROOT/'MARKER_RECOVERY.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=run_directory(PORTABLE,ROOT) or ROOT/('marker_recovery_'+stamp)
        out.mkdir(parents=True,exist_ok=False);lineage=out/'lineage/embryophyta_odb10';shutil.copytree(original,lineage);check_lineage(lineage,e)
        run=out/'busco/Dioscorea_panthaica_DIPA1'
        busco=[TOOLS['busco_python'],TOOLS['busco'],'-i',str(transcript),'-m','transcriptome','--metaeuk','-l',str(lineage),'--offline','--download_path',str(out/'busco_downloads'),'-o',run.name,'--out_path',str(run.parent),'-c','4','--opt-out-run-stats','--metaeuk_parameters=--split-memory-limit=6G','--metaeuk_rerun_parameters=--split-memory-limit=6G']
        unit='dt-panthaica-busco-'+stamp.lower()
        guard=Path(__file__).with_name('exec_verified_scope.py')
        effective_path=out/'BUSCO.effective_resources.json'
        cmd=[TOOLS['systemd_run'],'--user','--scope','--quiet','--unit='+unit,'-p','MemoryMax=8G','-p','MemorySwapMax=0',sys.executable,str(guard),'--memory-bytes',str(8*1024**3),'--cpus','4','--receipt',str(effective_path),'--',TOOLS['time'],'-v']+busco
        state={'status':'RUNNING','started_utc':now(),'command':cmd,'worker_pid':os.getpid(),'process_start_ticks':Path('/proc/self/stat').read_text().split()[21],'scope_guard_sha256':sha(guard),'effective_resource_receipt':str(effective_path),'cpu_limit_method':'inherited OS affinity','transcripts_path':str(transcript),'transcripts_sha256':a.transcripts_sha256,'assembly_receipt_path':str(assembly_receipt),'assembly_receipt_sha256':a.assembly_receipt_sha256,'busco_run':str(run),'readiness_evidence_sha256':PORTABLE.input('marker_evidence',EVIDENCE_SHA)['sha256'],'driver_sha256':sha(__file__),'lineage_manifest_sha256':e['lineage']['manifest_sha256'],'threads':4,'busco_cgroup_memory_max_bytes':8*1024**3,'wall_time_limit':None,'out':str(out)}
        def stop(sig,frame):
            subprocess.run([TOOLS['systemctl'],'--user','stop',unit+'.scope'],capture_output=True)
            raise KeyboardInterrupt('Stopped only owned BUSCO scope')
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        try:
            save(out/'STATUS.json',state);save(ROOT/'MARKER_RECOVERY_CURRENT.json',state)
            with (out/'busco.launcher.log').open('w') as log:
                child=subprocess.Popen(cmd,cwd=out,env=environment,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
                while child.poll() is None:
                    state['updated_utc']=now();save(out/'STATUS.json',state);time.sleep(5)
            state['exit_code']=child.returncode
            if child.returncode:raise RuntimeError('BUSCO failed; preserve scope/log evidence and inspect before any resource change')
            effective=json.loads(effective_path.read_text())
            if effective.get('status')!='VERIFIED_BEFORE_EXEC' or effective.get('memory_max_bytes')!=8*1024**3 or effective.get('swap_max_bytes')!=0 or len(effective.get('logical_cpu_ids',[]))!=4:
                raise ValueError('BUSCO effective resources differ from declared bounds')
            state['effective_resource_receipt_sha256']=sha(effective_path)
            result=extract(run,e,out);state.update(candidate_receipt_path=str(out/'CANDIDATE_RECEIPT.json'),candidate_receipt_sha256=sha(out/'CANDIDATE_RECEIPT.json'),status='COMPLETE',finished_utc=now(),accepted_candidates=result['accepted_candidates'],denominator=173,biological_interpretation='CANDIDATE_MARKERS_ONLY_NOT_PHYLOGENY_VALIDATED')
            state['receipt']=str(out/'MARKER_RECOVERY_RECEIPT.json')
            save(out/'STATUS.json',state);save(out/'MARKER_RECOVERY_RECEIPT.json',state);save(ROOT/'MARKER_RECOVERY_CURRENT.json',state)
        except BaseException as ex:
            state.update(status='INTERRUPTED' if isinstance(ex,KeyboardInterrupt) else 'FAILED',finished_utc=now(),error=type(ex).__name__+': '+str(ex));save(out/'STATUS.json',state);save(ROOT/'MARKER_RECOVERY_CURRENT.json',state);raise
        print(json.dumps(state,indent=2))


def plan(config_path):
    c=Config(config_path);root=c.path('root');fixed=run_directory(c,root)
    return {'script':'run_marker_recovery_draft',
      'inputs':[c.input('marker_evidence',EVIDENCE_SHA),c.input('transcripts'),c.input('assembly_receipt'),c.input('download_receipt')],
      'outputs':[str(fixed/n) for n in ['MARKER_RECOVERY_RECEIPT.json','CANDIDATE_RECEIPT.json','STATUS.json']] if fixed else [], 'output_directory':str(fixed or root),'required_threads':4,
      'dynamic_outputs':{'pattern':'marker_recovery_<UTC timestamp>/','files':['STATUS.json','MARKER_RECOVERY_RECEIPT.json','CANDIDATE_RECEIPT.json','candidates/<marker>.faa','busco/Dioscorea_panthaica_DIPA1/run_embryophyta_odb10/full_table.tsv']},
      'tools':{**{k:c.tool(k) for k in ('busco','busco_python','metaeuk')},**{k:c.tool(k,d) for k,d in [('systemd_run','systemd-run'),('systemctl','systemctl'),('time','time')]}},
      'python_packages':['Bio'],'scientific_contract':'173 fixed markers; 1614 BUSCO statuses; full-read assembly only; 4 CPUs / 8 GiB / no swap',
      'transitive_inputs':'marker_evidence binds baseline main/summary/raw loci, BUSCO sources, and the complete lineage file manifest'}

def main():
    global ROOT,ENV,PORTABLE,TOOLS
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True,type=Path);p.add_argument('--execute',action='store_true');a=p.parse_args()
    PORTABLE=Config(a.config);spec=plan(a.config)
    if not a.execute:return report(spec)
    require_preflight(spec);ROOT=PORTABLE.path('root');TOOLS=spec['tools'];ENV=Path(TOOLS['busco']).parent
    for role in ('transcripts','assembly_receipt'):
        row=PORTABLE.input(role)
        if not row['sha256']:raise ValueError(role+' requires an explicit SHA-256')
        setattr(a,role,row['path']);setattr(a,role+'_sha256',row['sha256'])
    execute(a);return 0
if __name__=='__main__':raise SystemExit(main())
