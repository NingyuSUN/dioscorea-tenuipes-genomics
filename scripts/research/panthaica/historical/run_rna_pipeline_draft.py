#!/usr/bin/env python3
"""Prepare-only default. Explicit independent qc / pilot / full stages for CRR515167.
No automatic cascade, rerun, assembly restart, strand inference, or marker recovery.
"""
import argparse
from collections import Counter
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import time

ROOT=None
PORTABLE=None
TOOLS={}
FILES=[('CRR515167_f1.fq.gz',2040518277,'bd4cf0ecacc027cf691180dd5d018b33'),('CRR515167_r2.fq.gz',2084070558,'21f081f1130ee01b5fd3fd523c45bb3b')]
RAW_PAIRS=28665951
CONFIG=None
CURRENT='RNA_WORKFLOW_CURRENT.json'
CONFIG_EXPECTED_SHA='ed19abb14b85a9d10a626f2b43146a1a55bf565da038b4b07b0f0c4996a1611c'


try:
    from .portable_inputs import Config, report, require_preflight, run_directory
except ImportError:
    from portable_inputs import Config, report, require_preflight, run_directory

def now():return dt.datetime.now(dt.timezone.utc).isoformat()
def hashfile(path,kind='sha256'):
    h=hashlib.new(kind)
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def save(path,obj):
    temp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    with temp.open('w') as f:json.dump(obj,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temp,path)
def signature(p):
    s=p.stat();return {'size':s.st_size,'mtime_ns':s.st_mtime_ns,'ctime_ns':s.st_ctime_ns,'inode':s.st_ino,'device':s.st_dev}
def bound(path):
    # This helper only resolves read-only provenance inputs, never output paths.
    return PORTABLE.resolve(path) if PORTABLE else Path(path).resolve()

def raw_path(name):
    directory=PORTABLE.path('raw_dir',ROOT/'raw') if PORTABLE else ROOT/'raw'
    return directory/name

def read_record(f,mate,minlen=1,maxlen=150):
    lines=[f.readline(1000001) for _ in range(4)]
    if not lines[0]:
        if any(lines[1:]):raise ValueError('Malformed FASTQ ending')
        return None
    if any(not s or len(s)>=1000001 for s in lines):raise ValueError('Incomplete/oversized FASTQ record')
    header,seq,plus,qual=[s.rstrip(b'\r\n') for s in lines]
    if not header.startswith(b'@') or not plus.startswith(b'+'):raise ValueError('FASTQ header/plus format')
    fields=header[1:].split();ident=fields[0] if fields else b''
    if not ident:raise ValueError('Empty read ID')
    suffix=None
    if ident.endswith((b'/1',b'/2')):suffix=int(ident[-1:]);ident=ident[:-2]
    if suffix is not None and suffix!=mate:raise ValueError('Incorrect mate suffix')
    if len(fields)>1 and fields[1].startswith((b'1:',b'2:')) and int(fields[1][:1])!=mate:raise ValueError('Incorrect CASAVA mate')
    if not minlen<=len(seq)<=maxlen or len(seq)!=len(qual):raise ValueError('Sequence/quality length failure')
    if not set(seq.upper())<=set(b'ACGTN') or min(qual)<33 or max(qual)>126:raise ValueError('Invalid bases/quality')
    return ident,seq.upper(),qual,b'\n'.join([header,seq,plus,qual])+b'\n'

def pairs(paths,limit=None,minlen=1,maxlen=150):
    with gzip.open(paths[0],'rb') as a,gzip.open(paths[1],'rb') as b:
        count=0
        while limit is None or count<limit:
            one=read_record(a,1,minlen,maxlen);two=read_record(b,2,minlen,maxlen)
            if one is None or two is None:
                if one is not None or two is not None:raise ValueError('Mates have different record counts')
                break
            if one[0]!=two[0]:raise ValueError('Pair ID mismatch')
            count+=1;yield one,two
        if limit is not None and count!=limit:raise ValueError('Insufficient prefix pairs')

def diagnostic(paths,limit,config,tick=lambda **kw:None,minlen=1):
    motifs={'NEB_R1_prefix16':config['adapter_r1'][:16].encode(),'NEB_R2_prefix16':config['adapter_r2'][:16].encode()}
    stats=[{'reads':0,'bases':0,'q20_bases':0,'q30_bases':0,'sum_phred33':0,'min_ascii':127,'max_ascii':0,'lengths':Counter(),'motif_read_counts':Counter(),'motif_positions':{n:Counter() for n in motifs},'terminal_homopolymer_ge10':Counter(),'cycles':[]} for _ in range(2)]
    for count,pair in enumerate(pairs(paths,limit,minlen),1):
        for m,record in enumerate(pair):
            _,seq,qual,_=record;s=stats[m];s['reads']+=1;s['bases']+=len(seq);s['lengths'][len(seq)]+=1
            for name,motif in motifs.items():
                at=seq.find(motif)
                if at>=0:s['motif_read_counts'][name]+=1;s['motif_positions'][name][at+1]+=1
            for base in b'ACGT':
                if seq.endswith(bytes([base])*10):s['terminal_homopolymer_ge10'][chr(base)]+=1
            while len(s['cycles'])<len(seq):s['cycles'].append({'n':0,'sum_phred33':0,'q20':0,'q30':0,'bases':Counter()})
            for i,(base,qbyte) in enumerate(zip(seq,qual)):
                q=qbyte-33;c=s['cycles'][i];c['n']+=1;c['sum_phred33']+=q;c['q20']+=q>=20;c['q30']+=q>=30;c['bases'][chr(base)]+=1
                s['sum_phred33']+=q;s['q20_bases']+=q>=20;s['q30_bases']+=q>=30
            s['min_ascii']=min(s['min_ascii'],min(qual));s['max_ascii']=max(s['max_ascii'],max(qual))
        if count%1000==0:tick(stage='PREFIX_DIAGNOSTIC',pairs_done=count,pairs_total=limit)
    for s in stats:
        s['mean_phred33']=s['sum_phred33']/s['bases'];s['q20_fraction']=s['q20_bases']/s['bases'];s['q30_fraction']=s['q30_bases']/s['bases']
    return {'prefix_pairs':limit,'random_sample':False,'scope':'Ordered prefix diagnostic, not full-run QC or unbiased library estimate','phred_assumption':33,'motif_rule':'Exact first16 bases of documented NEBNext adapters; counts indicate matches, not validated adapter origin by themselves','tail_rule':'Terminal A/C/G/T runs recorded only; never automatically removed as adapters','mates':stats}

def validate_trimmed(paths,expected,tick):
    n=0;bases=[0,0];lengths=[Counter(),Counter()]
    for n,pair in enumerate(pairs(paths,None,50,150),1):
        for j,r in enumerate(pair):bases[j]+=len(r[1]);lengths[j][len(r[1])]+=1
        if n%100000==0:tick(stage='VALIDATE_TRIMMED_FULL_STREAM',pairs_done=n,pairs_total=expected)
    if n!=expected:raise ValueError('Cutadapt JSON/output pair count mismatch')
    return {'pairs':n,'bases_per_mate':bases,'lengths_per_mate':lengths,'all_ids_lengths_gzip_crc_eof_verified':True}

def check_download():
    p=(Path(PORTABLE.input('download_receipt')['path']) if PORTABLE else ROOT/'DOWNLOAD_RECEIPT.json');r=json.loads(p.read_text())
    if r.get('status')!='COMPLETE' or r.get('pairs_done')!=RAW_PAIRS:raise ValueError('Require COMPLETE full raw validation receipt')
    pv=r.get('paired_fastq_validation',{})
    if pv.get('pairs')!=RAW_PAIRS or pv.get('bases')!=RAW_PAIRS*300 or pv.get('bases_per_mate_record')!=150 or any(pv.get(k) is not True for k in ['mate_ids_all_match','record_structure_all_pass','sequence_quality_lengths_all_pass','gzip_crc_and_complete_eof_both_pass']):raise ValueError('Complete paired FASTQ/gzip validation evidence is missing')
    byname={x['name']:x for x in r['files']}
    for name,size,md5 in FILES:
        row=byname.get(name,{})
        if row.get('bytes')!=size or row.get('md5')!=md5 or row.get('official_md5_match') is not True:raise ValueError('Official file metadata mismatch in receipt')
        pth=raw_path(name)
        if pth.is_symlink() or pth.stat().st_size!=size:raise ValueError('Raw path/size changed')
    return r,hashfile(p)

def verify_raw_hashes(tick):
    rows=[]
    for name,size,md5 in FILES:
        p=raw_path(name);before=signature(p);h=hashlib.md5();s=hashlib.sha256();done=0
        with p.open('rb') as f:
            for block in iter(lambda:f.read(8*1024**2),b''):
                h.update(block);s.update(block);done+=len(block);tick(stage='VERIFY_OFFICIAL_MD5',file=name,bytes_done=done,bytes_total=size)
        if h.hexdigest()!=md5 or signature(p)!=before:raise ValueError('Raw MD5 mismatch or file changed during hash')
        rows.append({'path':str(p),'name':name,'bytes':size,'md5':md5,'sha256':s.hexdigest(),'signature':before})
    return rows

def receipt(path,phase,configsha):
    p=bound(path);r=json.loads(p.read_text())
    if r.get('status')!='COMPLETE' or r.get('phase')!=phase or r.get('config_sha256')!=configsha:raise ValueError('Required preceding stage/config receipt mismatch')
    return r,hashfile(p)

def verify_raw_signature(rows):
    if len(rows)!=2:raise ValueError('Expected two verified raw file signatures')
    for row,(name,size,md5) in zip(rows,FILES):
        if row['name']!=name or row['md5']!=md5 or row['bytes']!=size or signature(raw_path(name))!=row['signature']:raise ValueError('Raw file differs from full-MD5 QC checkpoint')

def resource_read(unit):
    p=subprocess.run([TOOLS['systemctl'],'--user','show',unit+'.scope','--property=ControlGroup','--value'],capture_output=True,text=True)
    path=p.stdout.strip()
    if not path.startswith('/user.slice/'):return {}
    root=Path('/sys/fs/cgroup')/path.lstrip('/');d={'cgroup_path':str(root)}
    for name in ['memory.current','memory.peak','memory.events','cpu.stat']:
        f=root/name
        try:d[name]=f.read_text().strip()
        except FileNotFoundError:pass
    return d

class Runner:
    def __init__(self,out,phase,state,current):self.out=out;self.phase=phase;self.state=state;self.current=current;self.last=0;self.unit=None
    def tick(self,force=False,**kw):
        self.state.update(kw)
        if force or time.monotonic()-self.last>=10:
            self.state['updated_utc']=now();save(self.out/'STATUS.json',self.state)
            self.current[self.phase]={k:self.state.get(k) for k in ['status','stage','updated_utc','worker_pid','process_start_ticks','run_directory','receipt']}
            save(ROOT/CURRENT,self.current);self.last=time.monotonic()
            print(json.dumps({k:self.state.get(k) for k in ['status','stage','updated_utc','pairs_done','pairs_total','bytes_done','bytes_total','sampled_memory_peak_bytes','error']}),flush=True)
    def stop(self,sig,frame):
        if self.unit:subprocess.run([TOOLS['systemctl'],'--user','stop',self.unit+'.scope'],capture_output=True)
        raise KeyboardInterrupt('Stopped owned scope; files preserved')
    def command(self,args,label,memory_gib,threads):
        unit='dt-panthaica-'+self.phase+'-'+label+'-'+str(os.getpid())
        guard=Path(__file__).with_name('exec_verified_scope.py')
        guard_receipt=self.out/(label+'.effective_resources.json')
        cmd=[TOOLS['systemd_run'],'--user','--scope','--quiet','--unit='+unit,'-p','MemoryMax='+str(memory_gib)+'G','-p','MemorySwapMax=0',sys.executable,str(guard),'--memory-bytes',str(memory_gib*1024**3),'--cpus',str(threads),'--receipt',str(guard_receipt),'--',TOOLS['time'],'-v']+args
        entry={'label':label,'argv':cmd,'started_utc':now(),'memory_max_bytes':memory_gib*1024**3,'threads':threads,'scope_guard_sha256':hashfile(guard),'cpu_limit_method':'inherited OS affinity, not unverified CPUQuota','resource_scope':'tool command and descendants; coordinator separately bounded','sampled_memory_peak_bytes':None,'wall_time_limit':None}
        self.state.setdefault('commands',[]).append(entry);self.unit=unit
        env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        log=self.out/(label+'.log');samples=self.out/(label+'.resources.jsonl')
        self.tick(force=True,stage=label,log=str(log),pairs_done=None,pairs_total=None,bytes_done=None,bytes_total=None)
        start=time.monotonic()
        with log.open('w') as f,samples.open('w') as sf:
            child=subprocess.Popen(cmd,cwd=self.out,env=env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
            entry['launcher_pid']=child.pid
            while child.poll() is None:
                d=resource_read(unit)
                if d:
                    d['utc']=now();sf.write(json.dumps(d)+'\n');sf.flush()
                    peak=d.get('memory.peak')
                    if peak and peak.isdigit():entry['sampled_memory_peak_bytes']=max(entry['sampled_memory_peak_bytes'] or 0,int(peak))
                if shutil.disk_usage(self.out).free<2*1024**3:
                    subprocess.run([TOOLS['systemctl'],'--user','stop',unit+'.scope'],capture_output=True)
                    raise RuntimeError('Disk reserve below2GiB; stopped owned scope')
                self.tick(sampled_memory_peak_bytes=entry['sampled_memory_peak_bytes']);time.sleep(2)
        self.unit=None;entry.update(exit_code=child.returncode,finished_utc=now(),elapsed_seconds=round(time.monotonic()-start,3),log_sha256=hashfile(log))
        self.tick(force=True)
        if child.returncode:raise RuntimeError(label+' failed; inspect preserved log and cgroup samples before any new attempt')
        effective=json.loads(guard_receipt.read_text())
        if effective.get('status')!='VERIFIED_BEFORE_EXEC' or effective['memory_max_bytes']!=memory_gib*1024**3 or effective['swap_max_bytes']!=0 or len(effective['logical_cpu_ids'])!=threads:
            raise ValueError('Effective command resources do not match the requested bounds')
        entry['effective_resource_receipt_sha256']=hashfile(guard_receipt)
        return entry

def summarize_fasta(p):
    lengths=[];length=0;n=0;ambiguous=0
    with p.open() as f:
        for line in f:
            if line.startswith('>'):
                if n:lengths.append(length)
                n+=1;length=0
            else:
                s=line.strip().upper()
                if s and (not n or not set(s)<=set('ACGTN')):raise ValueError('Malformed transcript FASTA')
                length+=len(s);ambiguous+=s.count('N')
    if n:lengths.append(length)
    if not lengths or not sum(lengths):raise ValueError('No assembled transcript sequence')
    total=sum(lengths);accum=0;n50=0
    for size in sorted(lengths,reverse=True):
        accum+=size
        if accum>=total/2:n50=size;break
    return {'records':len(lengths),'bases':total,'N_bases':ambiguous,'N50':n50,'min_length':min(lengths),'max_length':max(lengths),'meaning':'Transcript contig summaries; not genome completeness or isoform correctness'}

def execute(a,c):
    if ROOT is None or ROOT.resolve()!=ROOT or not ROOT.is_dir():raise RuntimeError('Configured task root is not an existing absolute directory')
    download,download_sha=check_download();configsha=hashfile(CONFIG)
    with (ROOT/'RNA_WORKFLOW.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        cp=ROOT/CURRENT;current=json.loads(cp.read_text()) if cp.exists() else {}
        if a.phase in current and not a.new_attempt:
            raise RuntimeError('Phase already has a current record; inspect it. New attempt requires explicit --new-attempt and preserves old outputs')
        qc=None;pilot=None
        if a.phase!='qc':
            qc,qcsha=receipt(a.qc_receipt,'qc',configsha);verify_raw_signature(qc['raw_files'])
            if qc['download_receipt_sha256']!=download_sha:raise ValueError('Download receipt changed after QC checkpoint')
        if a.phase=='full':
            pilot,pilotsha=receipt(a.pilot_receipt,'pilot',configsha)
            if pilot.get('formal_marker_eligible') is not False or pilot['qc_receipt_sha256']!=qcsha:raise ValueError('Require matching non-formal prefix pilot receipt')
        stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        out=run_directory(PORTABLE,ROOT) or ROOT/('qc' if a.phase=='qc' else 'assembly')/(a.phase+'_'+stamp);out.mkdir(parents=True,exist_ok=False)
        state={'status':'RUNNING','phase':a.phase,'stage':'PREFLIGHT','started_utc':now(),'worker_pid':os.getpid(),'process_start_ticks':Path('/proc/self/stat').read_text().split()[21],'driver_sha256':hashfile(__file__),'portable_config_sha256':hashfile(PORTABLE.filename),'config_sha256':configsha,'download_receipt_sha256':download_sha,'run_directory':str(out),'wall_time_limit':None,'auto_next_stage':False,'auto_restart':False}
        run=Runner(out,a.phase,state,current);signal.signal(signal.SIGTERM,run.stop);signal.signal(signal.SIGINT,run.stop);run.tick(force=True)
        try:
            paths=[raw_path(x[0]) for x in FILES]
            if a.phase=='qc':
                state['raw_files']=verify_raw_hashes(run.tick)
                diag=diagnostic(paths,c['diagnostic_prefix_pairs'],c,run.tick,minlen=150);save(out/'PREFIX_DIAGNOSTIC.json',diag)
                state.update(input_scope='PREFIX_100K_DIAGNOSTIC',source_raw_pairs=RAW_PAIRS,diagnostic_path=str(out/'PREFIX_DIAGNOSTIC.json'),diagnostic_sha256=hashfile(out/'PREFIX_DIAGNOSTIC.json'),formal_marker_eligible=False)
                receipt_path=out/'QC_RECEIPT.json'
            else:
                needed=(50 if a.phase=='full' else 8)*1024**3
                if shutil.disk_usage(out).free<needed:raise RuntimeError('Insufficient disk for declared '+a.phase+' preflight reserve')
                state.update(raw_files=qc['raw_files'],qc_receipt_sha256=qcsha,source_raw_pairs=RAW_PAIRS)
                if pilot:state['pilot_receipt_sha256']=pilotsha
                expected=RAW_PAIRS
                if a.phase=='pilot':
                    expected=c['pilot_prefix_pairs'];prefix=[out/'prefix_raw_R1.fq.gz',out/'prefix_raw_R2.fq.gz']
                    run.tick(force=True,stage='EXTRACT_NONRANDOM_PREFIX',pairs_done=0,pairs_total=expected)
                    with gzip.GzipFile(filename=str(prefix[0]),mode='wb',compresslevel=1,mtime=0) as f1,gzip.GzipFile(filename=str(prefix[1]),mode='wb',compresslevel=1,mtime=0) as f2:
                        for n,pair in enumerate(pairs(paths,expected,150),1):
                            f1.write(pair[0][3]);f2.write(pair[1][3])
                            if n%10000==0:run.tick(pairs_done=n)
                    paths=prefix
                version=subprocess.run([c['cutadapt_path'],'--version'],capture_output=True,text=True,check=True).stdout.strip()
                if version!=c['cutadapt_version']:raise ValueError('Cutadapt version drift')
                state['cutadapt_tool']={'path':c['cutadapt_path'],'version':version,'launcher_sha256':hashfile(c['cutadapt_path'])}
                trimmed=[out/'trimmed_R1.fq.gz',out/'trimmed_R2.fq.gz'];j=out/'cutadapt.json'
                command=[c['cutadapt_path'],'-j',str(c['trimming_threads']),'-q','20','--quality-base','33','-m','50','--pair-filter=any','-a',c['adapter_r1'],'-A',c['adapter_r2'],'-O','12','-e','0.1','--no-indels','--json',str(j),'-o',str(trimmed[0]),'-p',str(trimmed[1]),str(paths[0]),str(paths[1])]
                run.command(command,'CUTADAPT',c['trimming_cgroup_gib'],c['trimming_threads'])
                report=json.loads(j.read_text());counts=report['read_counts']
                if counts['input']!=expected or counts['output']<=0:raise ValueError('Cutadapt input count differs or no reads retained')
                state['trimming']={'input_pairs':expected,'output_pairs':counts['output'],'retained_fraction':counts['output']/expected,'filtered':counts.get('filtered'),'report_path':str(j),'report_sha256':hashfile(j)}
                state['trimmed_validation']=validate_trimmed(trimmed,counts['output'],run.tick)
                save(out/'POST_TRIM_PREFIX_DIAGNOSTIC.json',diagnostic(trimmed,min(100000,counts['output']),c,run.tick,minlen=50))
                state['trimmed_files']=[{'path':str(p),'bytes':p.stat().st_size,'sha256':hashfile(p)} for p in trimmed]
                launcher=c['rna_spades_launcher'];v=subprocess.run([TOOLS['spades_python'],launcher,'--version'],capture_output=True,text=True,check=True).stdout
                if 'v4.3.0' not in v or 'rnaSPAdes mode' not in v:raise ValueError('SPAdes version/mode mismatch')
                install=json.loads((Path(PORTABLE.input('spades_install_receipt')['path'])).read_text())
                if install.get('status')!='COMPLETE' or install.get('asset_sha256')!=c['rna_spades_source_asset_sha256']:raise ValueError('Verified RNA smoke/install receipt required')
                state['rna_spades_tool']={'launcher':launcher,'version_stdout':v.strip(),'launcher_sha256':hashfile(launcher),'install_receipt_sha256':hashfile(Path(PORTABLE.input('spades_install_receipt')['path']))}
                threads=c[a.phase+'_threads'];mem=c[a.phase+'_spades_memory_gb'];cap=c[a.phase+'_cgroup_gib'];assembly=out/'rnaspades'
                cmd=[TOOLS['spades_python'],launcher,'-1',str(trimmed[0]),'-2',str(trimmed[1]),'-t',str(threads),'-m',str(mem),'--phred-offset','33','-o',str(assembly)]
                run.command(cmd,'RNASPADES',cap,threads)
                transcript=assembly/'transcripts.fasta'
                state.update(input_scope='FULL_CRR515167' if a.phase=='full' else 'NONRANDOM_PREFIX_1M',assembly_input_pairs=counts['output'],formal_marker_eligible=a.phase=='full',transcripts={'path':str(transcript),'bytes':transcript.stat().st_size,'sha256':hashfile(transcript)},assembly_statistics=summarize_fasta(transcript))
                state['graph_files']=[{'path':str(p),'bytes':p.stat().st_size} for p in assembly.glob('*') if p.suffix in ['.gfa','.fastg','.paths']]
                state['limits']=['RNA expression and assembly isoforms do not equal genome copy number.','Pilot is a nonrandom prefix and may only diagnose resource/workflow behavior.','No linear full-run time or memory extrapolation is justified.','Formal marker eligibility means full input processing only, not validated orthology or phylogeny.']
                receipt_path=out/('ASSEMBLY_RECEIPT.json' if a.phase=='full' else 'PILOT_RECEIPT.json')
            verify_raw_signature(state['raw_files'])
            state.update(status='COMPLETE',stage='COMPLETE',finished_utc=now(),receipt=str(receipt_path),error=None)
            save(receipt_path,state)
            run.tick(force=True)
        except BaseException as ex:
            run.tick(force=True,status='INTERRUPTED' if isinstance(ex,KeyboardInterrupt) else 'FAILED',finished_utc=now(),error=type(ex).__name__+': '+str(ex));raise


def plan(config_path):
    c=Config(config_path);root=c.path('root');phase=c.data.get('phase','qc')
    if phase not in ('qc','pilot','full'):raise ValueError('phase must be qc, pilot, or full')
    inputs=[c.input('workflow_settings',CONFIG_EXPECTED_SHA),c.input('download_receipt')]
    inputs += [{'role':'raw_'+str(i+1),'path':str(c.path('raw_dir',root/'raw')/name),'sha256':None} for i,(name,_,_) in enumerate(FILES)]
    tools={}
    if phase!='qc':
        inputs += [c.input('qc_receipt'),c.input('spades_install_receipt')]
        tools={k:c.tool(k) for k in ('cutadapt','rna_spades_launcher','spades_python')}
        tools.update({k:c.tool(k,d) for k,d in [('systemd_run','systemd-run'),('systemctl','systemctl'),('time','time')]})
    if phase=='full':inputs.append(c.input('pilot_receipt'))
    folder=root/('qc' if phase=='qc' else 'assembly');fixed=run_directory(c,root)
    receipt_name={'qc':'QC_RECEIPT.json','pilot':'PILOT_RECEIPT.json','full':'ASSEMBLY_RECEIPT.json'}[phase]
    return {'script':'run_rna_pipeline_draft','phase':phase,'inputs':inputs,
      'outputs':([str(fixed/receipt_name),str(fixed/'STATUS.json')]+([] if phase=='qc' else [str(fixed/'rnaspades/transcripts.fasta')]) if fixed else []),'output_directory':str(fixed or folder),'required_threads':1 if phase=='qc' else 4,
      'dynamic_outputs':{'pattern':phase+'_<UTC timestamp>/','files':['STATUS.json',{'qc':'QC_RECEIPT.json','pilot':'PILOT_RECEIPT.json','full':'ASSEMBLY_RECEIPT.json'}[phase]],
                         'assembly_files':[] if phase=='qc' else ['trimmed_R1.fq.gz','trimmed_R2.fq.gz','cutadapt.json','rnaspades/transcripts.fasta']},
      'tools':tools,'non_executable_tools':['rna_spades_launcher'],'python_packages':[],'scientific_contract':'CRR515167: 28665951 pairs; fixed MD5/size binding, prefix pilot never formal marker input'}

def main():
    global ROOT,CONFIG,PORTABLE,TOOLS
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True,type=Path);p.add_argument('--execute',action='store_true')
    p.add_argument('--new-attempt',action='store_true')
    a=p.parse_args();PORTABLE=Config(a.config);spec=plan(a.config)
    if not a.execute:return report(spec)
    require_preflight(spec);ROOT=PORTABLE.path('root');a.phase=spec['phase'];TOOLS=spec['tools']
    CONFIG=Path(PORTABLE.input('workflow_settings',CONFIG_EXPECTED_SHA)['path'])
    c=json.loads(CONFIG.read_text())
    # Paths are deployment settings; scientific thresholds remain frozen in the hashed settings file.
    if a.phase!='qc':
        c['cutadapt_path']=TOOLS['cutadapt'];c['rna_spades_launcher']=TOOLS['rna_spades_launcher']
    a.qc_receipt=PORTABLE.input('qc_receipt')['path'] if a.phase!='qc' else None
    a.pilot_receipt=PORTABLE.input('pilot_receipt')['path'] if a.phase=='full' else None
    execute(a,c);return 0
if __name__=='__main__':raise SystemExit(main())
