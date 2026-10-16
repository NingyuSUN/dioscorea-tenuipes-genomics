#!/usr/bin/env python3
"""Dry-run default. Add qualified full-run panthaica AA to fixed 173 nuclear profiles.
Preserve original profiles and trim masks; create paired matrices only, never trees.
"""
import argparse,csv,datetime as dt,fcntl,hashlib,json,os,re,signal,subprocess,time
from pathlib import Path

ROOT=None
PORTABLE=None
NEW='Dioscorea_panthaica_DIPA1_transcriptome'
AA=set('ACDEFGHIKLMNPQRSTVWY')
EVIDENCE_SHA='ccbfb51d4c14a9853d36ce769dfc48f17b8a2c4d931409005ab1423fcbe3409d'
MARKER_EVIDENCE_SHA='de212e8f6f5502f5d75e68f8d509433c7c9e9cc6638cae304b6ce16ef456d660'
BASELINE_SHA='982addd02b6a6d0a188aa0de9b8e688b95bec926703094ab7f9238a8a91cc01e'
MAFFT=None


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
def save(p,d):
    q=p.with_name(p.name+'.tmp.'+str(os.getpid()));q.write_text(json.dumps(d,indent=2)+'\n');os.replace(q,p)
def read_fasta(p):
    result={};ident=None;seq=[]
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'):
            if ident is not None:result[ident]=''.join(seq)
            ident=line[1:].split()[0]
            if ident in result:raise ValueError('Duplicate FASTA ID')
            seq=[]
        elif line.strip():
            if ident is None:raise ValueError('Sequence before header')
            seq.append(line.strip())
    if ident is not None:result[ident]=''.join(seq)
    return result

def fasta_bytes(rows):
    return ''.join('>'+t+'\n'+'\n'.join(s[i:i+60] for i in range(0,len(s),60))+'\n' for t,s in rows.items()).encode()
def check_old_rows(reference,added,new_id=NEW):
    if set(added)!=set(reference)|{new_id}:raise ValueError('MAFFT taxon set changed')
    for taxon,old in reference.items():
        if added[taxon]!=old:raise ValueError('MAFFT changed old row or columns: '+taxon)
    lengths={len(x) for x in reference.values()}
    if len(lengths)!=1 or len(added[new_id])!=next(iter(lengths)):raise ValueError('Profile length changed')

def project(profile,taxa,columns):
    if not columns or columns!=sorted(set(columns)):raise ValueError('Frozen trim mask is empty/non-monotonic')
    n=len(next(iter(profile.values())))
    if columns[0]<0 or columns[-1]>=n:raise ValueError('Trim mask outside profile')
    if any(len(x)!=n for x in profile.values()):raise ValueError('Unequal profile rows')
    return {t:''.join(profile[t][i] if profile[t][i] in AA else '-' for i in columns) if t in profile else '-'*len(columns) for t in taxa}

def validate_map(text,sequence,aligned):
    """Check every source residue, insertion deletion and resulting profile column."""
    rows=[];headers=[]
    for line in text.splitlines():
        if line.startswith('>'):headers.append(line[1:].strip());continue
        if not line.strip() or line.startswith('#'):continue
        fields=[x.strip() for x in next(csv.reader([line]))]
        if len(fields)!=3:raise ValueError('Unexpected MAFFT map schema')
        letter,origin,target=fields;rows.append((letter,int(origin),None if target=='-' else int(target)))
    if headers!=[NEW] or [r[1] for r in rows]!=list(range(1,len(sequence)+1)):raise ValueError('Map does not cover exact candidate positions')
    reconstructed=['-']*len(aligned);targets=[]
    for letter,pos,target in rows:
        if letter!=sequence[pos-1]:raise ValueError('Map residue differs from candidate')
        if target is not None:
            if not 1<=target<=len(aligned) or reconstructed[target-1]!='-':raise ValueError('Invalid/duplicate map target')
            reconstructed[target-1]=letter;targets.append(target)
    if targets!=sorted(targets) or ''.join(reconstructed)!=aligned:raise ValueError('Map does not reconstruct added row')
    return rows

def concatenate(loci,taxa,ordered):
    return {t:''.join(loci[g][t] for g in ordered) for t in taxa}

def mafft_add(profile,sequence,folder):
    folder.mkdir(parents=True,exist_ok=False)
    ref=folder/'reference.faa';query=folder/'panthaica.faa';output=folder/'added_untrimmed.faa'
    ref.write_bytes(fasta_bytes(profile));query.write_bytes(fasta_bytes({NEW:sequence}))
    cmd=[MAFFT,'--amino','--localpair','--maxiterate','0','--thread','1','--inputorder','--add',str(query),'--keeplength','--mapout',str(ref)]
    start=time.monotonic()
    with output.open('w') as out,(folder/'mafft.log').open('w') as log:
        child=subprocess.Popen(cmd,stdout=out,stderr=log,cwd=folder,start_new_session=True)
        try:
            code=child.wait()
            if code:raise subprocess.CalledProcessError(code,cmd)
        except BaseException:
            if child.poll() is None:os.killpg(child.pid,signal.SIGTERM);child.wait()
            raise
    added=read_fasta(output);check_old_rows(profile,added)
    mapping=query.with_name(query.name+'.map')
    rows=validate_map(mapping.read_text(),sequence,added[NEW])
    return added,rows,{'command':cmd,'elapsed_seconds':round(time.monotonic()-start,3),'output_sha256':sha(output),'map_path':str(mapping),'map_sha256':sha(mapping),'reference_sha256':sha(ref),'candidate_sha256':sha(query)}

def tool_self_test(out):
    # Artificial AA strings only. Required at future execution before real profile processing.
    base='ACDEFGHIKLMNPQRSTVWYACGTHDMEKIPQRSVLNACTEFGHIKLMNPQRSTVWYACDEF'
    old1=base[:10]+'--'+base[10:40]+'-'+base[40:]
    other=base[:5]+'A'+base[6:]
    old2=other[:10]+'--'+other[10:40]+'-'+other[40:]
    candidate=base[:20]+'WWWWWW'+base[20:-10]
    added,mapping,receipt=mafft_add({'old1':old1,'old2':old2},candidate,out/'MAFFT_SYNTHETIC_SELF_TEST')
    if not any(r[2] is None for r in mapping) or '-' not in added[NEW]:raise ValueError('Synthetic insertion/deletion case not exercised')
    receipt.update(artificial_sequences_only=True,old_rows_exact=True,insertions_deleted_and_mapped=True,missing_positions_represented_by_gaps=True)
    return receipt

def read_bound(p,expected):
    if sha(p)!=expected:raise ValueError('Source SHA256 mismatch: '+str(p))
    return read_fasta(p)
def task_path(p):
    # Sources may live outside the output root; immutable receipts bind their bytes.
    return PORTABLE.resolve(p) if PORTABLE else Path(p).resolve()

def execute(a):
    def interrupted(sig,frame):raise KeyboardInterrupt('Alignment preparation interrupted; owned child group will stop')
    signal.signal(signal.SIGTERM,interrupted)
    if ROOT.resolve()!=ROOT:raise ValueError('Task root path guard')
    panel_row=PORTABLE.input('panel_evidence',EVIDENCE_SHA)
    ep=Path(panel_row['path']);panel_sha=panel_row['sha256']
    if sha(ep)!=panel_sha:raise ValueError('Frozen selected-panel evidence changed')
    e=PORTABLE.relocated(json.loads(ep.read_text()));taxa=e['taxa'];base=Path(e['baseline_root'])
    panel_id=e.get('panel_id','original_genome_main_no2661')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',panel_id):raise ValueError('Invalid panel ID')
    if len(taxa)!=12 or len(set(taxa))!=12 or NEW in taxa or e['genes']!=173 or e['columns']<=0:raise ValueError('Expected selected panel with 12 old taxa and original 173 loci')
    if not a.panel_evidence and (e['columns']!=68796 or e['main_sha256']!=BASELINE_SHA):raise ValueError('Original baseline default changed')
    marker_row=PORTABLE.input('marker_evidence',MARKER_EVIDENCE_SHA)
    mep=Path(marker_row['path'])
    if sha(mep)!=marker_row['sha256']:raise ValueError('Original fixed173 candidate-quality anchor changed')
    marker_anchor={r['marker']:r for r in json.loads(mep.read_text())['baseline']['markers']}
    if len(e['markers'])!=173 or {r['marker'] for r in e['markers']}!=set(marker_anchor):raise ValueError('Selected panel must contain exactly original fixed173 marker IDs')
    if a.panel_evidence:
        provenance=Path(e['provenance_receipt_path'])
        if sha(provenance)!=e['provenance_receipt_sha256']:raise ValueError('Selected replacement-panel provenance changed')
    baseline=read_bound(base/'main.faa',e['main_sha256'])
    if {len(s) for s in baseline.values()}!={e['columns']}:raise ValueError('Selected-panel concatenated length mismatch')
    if list(baseline)!=taxa or sha(base/'summary.json')!=e['summary_sha256']:raise ValueError('Baseline taxon order or summary changed')
    summary=json.loads((base/'summary.json').read_text())
    if len(summary['partitions'])!=173 or {p['gene'] for p in summary['partitions']}!=set(marker_anchor):raise ValueError('Selected-panel partition set differs from original173')
    marker_path=task_path(a.marker_receipt)
    if sha(marker_path)!=a.marker_receipt_sha256:raise ValueError('Marker receipt hash mismatch')
    marker=PORTABLE.relocated(json.loads(marker_path.read_text()))
    if marker.get('status')!='COMPLETE' or marker.get('readiness_evidence_sha256')!=marker_row['sha256']:raise ValueError('Require complete lineage-bound marker run')
    assembly_path=task_path(marker['assembly_receipt_path'])
    if sha(assembly_path)!=marker['assembly_receipt_sha256']:raise ValueError('Assembly receipt changed')
    assembly=PORTABLE.relocated(json.loads(assembly_path.read_text()))
    if assembly.get('status')!='COMPLETE' or assembly.get('input_scope')!='FULL_CRR515167' or assembly.get('source_raw_pairs')!=28665951 or assembly.get('formal_marker_eligible') is not True:raise ValueError('Prefix/pilot/unknown assembly rejected')
    if assembly['transcripts']['sha256']!=marker['transcripts_sha256']:raise ValueError('Assembly/marker input mismatch')
    candidate_path=task_path(marker['candidate_receipt_path'])
    if sha(candidate_path)!=marker['candidate_receipt_sha256']:raise ValueError('Candidate receipt changed')
    candidates=PORTABLE.relocated(json.loads(candidate_path.read_text()));audit={r['marker']:r for r in candidates['candidate_audit']}
    if len(candidates['candidate_audit'])!=173 or len(audit)!=173 or set(audit)!={r['marker'] for r in e['markers']}:raise ValueError('Candidate/fixed173 marker set mismatch')
    run=task_path(marker['busco_run']);table=run/'run_embryophyta_odb10/full_table.tsv'
    if sha(table)!=candidates['full_table_sha256']:raise ValueError('BUSCO full table changed')
    table_rows={}
    for line in table.read_text().splitlines():
        if line and not line.startswith('#'):
            fields=line.split('\t');table_rows.setdefault(fields[0],[]).append(fields)
    if sha(MAFFT)!=e['mafft']['launcher_sha256']:raise ValueError('MAFFT launcher changed')
    v=subprocess.run([MAFFT,'--version'],capture_output=True,text=True,check=True)
    if v.stdout+v.stderr!=e['mafft']['version_stdout']+e['mafft']['version_stderr']:raise ValueError('MAFFT version changed')
    with (ROOT/'ALIGNMENT_CONTROL.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        out=run_directory(PORTABLE,ROOT) or ROOT/('alignment_controls_'+panel_id+'_'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'));out.mkdir(parents=True,exist_ok=False)
        state={'status':'RUNNING','stage':'SELF_TEST','started_utc':now(),'worker_pid':os.getpid(),'process_start_ticks':Path('/proc/self/stat').read_text().split()[21],'driver_sha256':sha(__file__),'panel_id':panel_id,'panel_root':str(base),'selected_panel_main_sha256':e['main_sha256'],'profile_evidence_path':str(ep),'profile_evidence_sha256':panel_sha,'marker_receipt_sha256':a.marker_receipt_sha256,'assembly_receipt_sha256':marker['assembly_receipt_sha256'],'tree_run':False,'alignment_threads':1,'source_type':'GENOME_PLUS_TRANSCRIPTOME'}
        save(out/'STATUS.json',state)
        try:
            state['mafft_self_test']=tool_self_test(out)
            loci={};details=[];L=[];taxa13=taxa+[NEW];ordered=[p['gene'] for p in summary['partitions']]
            for number,ref in enumerate(e['markers'],1):
                g=ref['marker'];raw=read_bound(base/'raw'/(g+'.faa'),ref['raw_sha256']);profile=read_bound(base/'aligned'/(g+'.faa'),ref['aligned_sha256']);trim=read_bound(base/'trimmed'/(g+'.faa'),ref['trimmed_sha256']);cols=ref['retained_column_indices_0based']
                if any(profile[t].replace('-','')!=s for t,s in raw.items()) or project(profile,taxa,cols)!=trim:raise ValueError('Baseline raw/profile/mask binding failed')
                r=audit[g];row={'marker':g,'baseline_columns':len(cols),'raw_sha256':ref['raw_sha256'],'aligned_sha256':ref['aligned_sha256'],'trimmed_sha256':ref['trimmed_sha256'],'candidate_admitted':r['reason'] is None,'matrix_admitted':False,'reason':r['reason'],'coverage':0.0};newrow='-'*len(cols)
                if r['reason'] is None:
                    if r['busco_status']!='Complete' or len(table_rows[g])!=1 or table_rows[g][0][1]!='Complete' or table_rows[g][0][2]!=r['id']:raise ValueError('Candidate is not bound complete single-copy')
                    src=task_path(r['faa_path'])
                    if src!=(run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences'/(g+'.faa')).resolve():raise ValueError('FAA does not map to its exact single-copy marker path')
                    source=read_bound(src,r['faa_sha256'])
                    median_length=marker_anchor[g]['baseline_raw_median_aa']
                    if median_length!=r['baseline_raw_median_aa']:raise ValueError('Candidate-quality median differs from frozen original173 anchor')
                    if len(source)!=1 or list(source)!=[r['id']]:raise ValueError('BUSCO source ID/count mismatch')
                    seq=source[r['id']].upper();seq=seq[:-1] if seq.endswith('*') else seq;sequence=''.join(c if c in AA else 'X' for c in seq)
                    if '*' in seq or len(seq)<50 or sum(c in AA for c in seq)/len(seq)<.95 or not .70<=len(seq)/r['baseline_raw_median_aa']<=1.50:raise ValueError('Candidate quality rules do not reproduce')
                    candidate=candidate_path.parent/'candidates'/(g+'.faa')
                    if read_fasta(candidate)!={NEW:sequence}:raise ValueError('Exported FAA differs from bound BUSCO source')
                    added,mapping,receipt=mafft_add(profile,sequence,out/'loci'/g)
                    projected=project(added,taxa13,cols)
                    if {t:projected[t] for t in taxa}!=trim:raise ValueError('Old trimmed rows changed after addition')
                    (out/'loci'/g/'projected_observed.faa').write_bytes(fasta_bytes(projected))
                    coverage=sum(c in AA for c in projected[NEW])/len(cols)
                    kept={i+1 for i in cols};insertions=sum(pos is None for _,_,pos in mapping);outside=sum(pos is not None and pos not in kept for _,_,pos in mapping)
                    row.update(candidate_source_id=r['id'],candidate_faa_sha256=sha(candidate),candidate_length=len(sequence),alignment=receipt,coverage=coverage,inserted_residues_deleted=insertions,mapped_residues_outside_frozen_mask=outside,projected_sha256=sha(out/'loci'/g/'projected_observed.faa'),reason=None if coverage>=.5 else 'panthaica_under_50_percent_frozen_locus_coverage')
                    if coverage>=.5:newrow=projected[NEW];L.append(g);row['matrix_admitted']=True
                loci[g]={**trim,NEW:newrow};details.append(row)
                state.update(stage='ADD_AND_PROJECT',processed_markers=number,total_markers=173,accepted_loci=len(L),updated_utc=now());save(out/'STATUS.json',state)
            full=concatenate(loci,taxa13,ordered)
            if {t:full[t] for t in taxa}!=baseline:raise ValueError('Full173 old rows do not match original main')
            body=(base/'main.faa').read_bytes()
            if not body.endswith(b'\n'):raise ValueError('Unexpected baseline FASTA final newline')
            without=out/'fixed173_without_panthaica.faa';without.write_bytes(body)
            with_pan=out/'fixed173_with_panthaica.faa';with_pan.write_bytes(body+fasta_bytes({NEW:full[NEW]}))
            if sha(without)!=e['main_sha256'] or not with_pan.read_bytes().startswith(body):raise ValueError('Full baseline byte preservation failed')
            if read_fasta(with_pan)!=full:raise ValueError('Written full matrix mismatch')
            selected=[g for g in ordered if g in set(L)];matched=concatenate(loci,taxa13,selected)
            parts=[];column_map=[];offset=1
            oldparts={p['gene']:p for p in summary['partitions']}
            for g in selected:
                p=oldparts[g];parts.append({'gene':g,'start_1based':offset,'end_1based':offset+p['length']-1,'length':p['length'],'original_fixed173_start_1based':p['start_1based'],'original_fixed173_end_1based':p['end_1based']});offset+=p['length']
                column_map.extend(range(p['start_1based'],p['end_1based']+1))
            if selected:
                control={t:matched[t] for t in taxa};(out/'matchedL_without_panthaica.faa').write_bytes(fasta_bytes(control));(out/'matchedL_with_panthaica.faa').write_bytes(fasta_bytes(matched))
                if not (out/'matchedL_with_panthaica.faa').read_bytes().startswith((out/'matchedL_without_panthaica.faa').read_bytes()):raise ValueError('Old12 matched matrix bytes changed')
                for t in taxa:
                    projected=''.join(baseline[t][i-1] for i in column_map)
                    if matched[t]!=projected:raise ValueError('Same-locus control differs from original column projection')
                if read_fasta(out/'matchedL_without_panthaica.faa')!={t:s for t,s in read_fasta(out/'matchedL_with_panthaica.faa').items() if t!=NEW}:raise ValueError('Paired matchedL controls differ beyond one row')
            save(out/'partitions_fixed173.json',summary['partitions']);save(out/'partitions_matchedL.json',parts)
            (out/'matchedL_column_map.tsv').write_text('matched_column_1based\toriginal_fixed173_column_1based\n'+''.join(str(i)+'\t'+str(j)+'\n' for i,j in enumerate(column_map,1)))
            full_cov=sum(c in AA for c in full[NEW])/len(full[NEW]);matched_cov=sum(c in AA for c in matched[NEW])/len(matched[NEW]) if selected else None
            result={'status':'MATRICES_PREPARED_SCIENTIFIC_REVIEW_PENDING','loci_fixed_denominator':173,'matched_loci':len(selected),'matched_columns':len(column_map),'matched_marker_ids':selected,'panthaica_coverage_full173':full_cov,'panthaica_coverage_matchedL':matched_cov,'full173_panthaica_meets_50_percent_coverage':full_cov>=.5,'old12_full173_exact_sequence_match':True,'selected_panel_main_file_bytes_preserved':True,'selected_panel_id':panel_id,'selected_panel_main_sha256':e['main_sha256'],'candidate_quality_anchor':'original173_frozen_raw_median_independent_of_panel','matchedL_with_without_only_row_difference':bool(selected),'site_filter_recomputed':False,'tree_run':False,'locus_audit':details,'limits':['Fixed-reference insertions in the new taxon are deliberately discarded; raw FAA and mapout retain their provenance.','Adding to a fixed profile is not proof of orthology or absence of paralogs.','No topology-based locus or isoform selection occurs.','Different taxon-set matrices cannot share one AU likelihood comparison; competing topology tests must hold taxa, loci, columns and model fixed.']}
            save(out/'ALIGNMENT_AUDIT.json',result)
            outputs=[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in out.glob('*.faa')]
            state.update(status='COMPLETE',stage='COMPLETE',finished_utc=now(),matched_loci=len(selected),matched_columns=len(column_map),audit_sha256=sha(out/'ALIGNMENT_AUDIT.json'),outputs=outputs,scientific_status='MATRIX_PREPARATION_ONLY_NO_PHYLOGENY')
            save(out/'ALIGNMENT_CONTROL_RECEIPT.json',state);save(out/'STATUS.json',state)
        except BaseException as exc:
            state.update(status='INTERRUPTED' if isinstance(exc,KeyboardInterrupt) else 'FAILED',error=type(exc).__name__+': '+str(exc),finished_utc=now());save(out/'STATUS.json',state);raise
        print(json.dumps(state,indent=2))


def plan(config_path):
    c=Config(config_path);root=c.path('root');fixed=run_directory(c,root)
    return {'script':'prepare_alignment_controls','inputs':[c.input('panel_evidence',EVIDENCE_SHA),c.input('marker_evidence',MARKER_EVIDENCE_SHA),c.input('marker_receipt')],
      'outputs':[str(fixed/n) for n in ['ALIGNMENT_CONTROL_RECEIPT.json','ALIGNMENT_AUDIT.json','fixed173_without_panthaica.faa','fixed173_with_panthaica.faa','partitions_fixed173.json','partitions_matchedL.json','matchedL_column_map.tsv']] if fixed else [], 'output_directory':str(fixed or root),'required_threads':1,
      'dynamic_outputs':{'pattern':'alignment_controls_<panel_id>_<UTC timestamp>/',
        'files':['ALIGNMENT_CONTROL_RECEIPT.json','ALIGNMENT_AUDIT.json','STATUS.json','fixed173_without_panthaica.faa','fixed173_with_panthaica.faa','partitions_fixed173.json','partitions_matchedL.json','matchedL_column_map.tsv'],
        'conditional_files':['matchedL_without_panthaica.faa','matchedL_with_panthaica.faa']},
      'tools':{'mafft':c.tool('mafft')},'python_packages':[],
      'transitive_inputs':'panel evidence binds baseline main.faa, summary.json, raw/aligned/trimmed per-locus FASTA and frozen masks; marker receipt binds assembly, BUSCO table and candidates',
      'scientific_contract':'12 old taxa, fixed 173 loci; keeplength/mapout; 50% locus coverage; old sequence bytes preserved; no tree inference'}

def main():
    global ROOT,PORTABLE,MAFFT
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True,type=Path);p.add_argument('--execute',action='store_true');a=p.parse_args()
    PORTABLE=Config(a.config);spec=plan(a.config)
    if not a.execute:return report(spec)
    require_preflight(spec);ROOT=PORTABLE.path('root');MAFFT=spec['tools']['mafft']
    row=PORTABLE.input('marker_receipt')
    if not row['sha256']:raise ValueError('marker_receipt requires an explicit SHA-256')
    a.marker_receipt=row['path'];a.marker_receipt_sha256=row['sha256']
    a.panel_evidence=None if PORTABLE.data.get('original_panel',True) else PORTABLE.input('panel_evidence')['path']
    a.panel_evidence_sha256=PORTABLE.input('panel_evidence',EVIDENCE_SHA)['sha256']
    execute(a);return 0
if __name__=='__main__':raise SystemExit(main())
