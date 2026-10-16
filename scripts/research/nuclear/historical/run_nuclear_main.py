#!/usr/bin/env python3
"""Ready genomes + new published annotation + completed gracillima; fresh supported tree."""
from pathlib import Path
import datetime,fcntl,hashlib,json,os,subprocess,math
from Bio import SeqIO,Phylo
R=Path(__file__).resolve().parents[1]; B=Path('workspace/miniconda3/envs/busco/bin')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
def main():
 lock=(R/'main.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 rec=R/'NUCLEAR_RUN.json'
 if rec.exists():raise RuntimeError('Existing nuclear run: inspect status, do not duplicate')
 d={'status':'RUNNING','stage':'preflight','started_utc':now(),'pid':os.getpid(),'commands':[],'source_panel_sha256':sha(R/'metadata/MAIN_TAXA.json'),'script_sha256':sha(Path(__file__)),'excluded_orthogroups':['2661at3193'],'runtime_timeout':None}
 put(rec,d);logs=R/'logs';logs.mkdir(exist_ok=True)
 env=dict(os.environ);env['PATH']=str(B)+':'+env.get('PATH','')
 def run(stage,args):
  d['stage']=stage;c={'stage':stage,'argv':list(map(str,args)),'started_utc':now(),'exit_code':None};d['commands'].append(c);put(rec,d)
  with (logs/(stage+'.log')).open('w') as f:p=subprocess.run(c['argv'],stdout=f,stderr=subprocess.STDOUT,env=env)
  c.update(exit_code=p.returncode,ended_utc=now());put(rec,d)
  if p.returncode:raise RuntimeError(f'{stage} failed {p.returncode}; inspect {logs/(stage+".log")}')
 try:
  meta=json.loads((R/'metadata/WT42_DOWNLOAD.json').read_text());inputs=R/'inputs';inputs.mkdir(exist_ok=True)
  fasta=inputs/meta['name'];part=Path(str(fasta)+'.part')
  if not fasta.exists():
   run('download_WT42',['curl','--fail','--location','--retry','5','--connect-timeout','30','--continue-at','-','--output',part,meta['download_url']])
   assert part.stat().st_size==meta['size']
   assert hashlib.md5(part.read_bytes()).hexdigest()==meta['computed_md5'];part.replace(fasta)
  assert fasta.stat().st_size==meta['size'] and hashlib.md5(fasta.read_bytes()).hexdigest()==meta['computed_md5']
  seqs=list(SeqIO.parse(fasta,'fasta'));assert seqs and len({x.id for x in seqs})==len(seqs)
  put(R/'metadata/WT42_VERIFIED.json',{'status':'DOWNLOAD_MD5_AND_FASTA_PASS','source':meta,'sha256':sha(fasta),'protein_records':len(seqs),'source_compartment':'genome annotation','paper_name':'Dioscorea opposita WT42','deposited_name':'Dioscorea oppositifolia','accession':'GCA_056531125.1'})
  run('busco_WT42',[B/'python',B/'busco','-i',fasta,'-m','proteins','-l','workspace/wgs/busco_downloads/lineages/embryophyta_odb10','--offline','--download_path','workspace/wgs/busco_downloads','-o','Dioscorea_oppositifolia_WT42','--out_path',R/'busco','-c','8','--opt-out-run-stats'])
  run('build_matrix',[B/'python',R/'workflow/build_busco_matrices.py','--panel','genome_main_no2661','--taxa-json',R/'metadata/MAIN_TAXA.json','--threads','8'])
  matrix=R/'matrices/genome_main_no2661/main.faa';s=json.loads((matrix.parent/'summary.json').read_text())
  assert {'Dioscorea_tenuipes_historical','Dioscorea_gracillima'}<=set(s['taxa'])
  assert '2661at3193' not in {x['gene']for x in s['partitions']}
  assert s['main_sha256']==sha(matrix)
  put(R/'metadata/FINAL_NUCLEAR_PANEL.json',{'status':'FROZEN_BEFORE_INFERENCE','utc':now(),'taxa':s['taxa'],'excluded_input_taxa':s['excluded_taxa'],'genes':s['genes'],'columns':s['columns'],'matrix_sha256':sha(matrix),'excluded_orthogroups':['2661at3193'],'required_focal_taxa_per_locus':s['required_focal_taxa_per_locus'],'source_scope':'Existing genomic sources + gracillima genomic WGS + published WT42 genome annotation; broader WGS candidates deferred explicitly'})
  result=R/'results/genome_main_no2661';result.mkdir(parents=True,exist_ok=True);pre=result/'main'
  run('IQTREE_main',['workspace/wgs/dt_nuclear_expanded_20260915/tools/bin/iqtree2','-s',matrix,'-seed','20260928','-T','8','-mem','8000M','--prefix',pre,'-m','MFP','-mset','LG,WAG,JTT','-mrate','G4','-mfreq','F','-B','1000','-bnni','-alrt','1000','--boot-trees'])
  tree=Phylo.read(str(pre)+'.treefile','newick');names=[x.name for x in tree.get_terminals()]
  assert len(names)==len(set(names)) and set(names)==set(s['taxa'])
  assert all(x.branch_length is None or math.isfinite(x.branch_length) and x.branch_length>=0 for x in tree.find_clades())
  assert sum(bool(x.strip())for x in Path(str(pre)+'.ufboot').open())==1000
  assert 'Analysis results written to' in Path(str(pre)+'.log').read_text()
  d.update(status='INFERENCE_COMPLETE_SCIENTIFIC_REVIEW_AND_FIGURE_PENDING',stage='done',finished_utc=now(),taxa=s['taxa'],genes=s['genes'],columns=s['columns'],treefile_sha256=sha(Path(str(pre)+'.treefile')));put(rec,d)
 except Exception as e:
  d.update(status='FAILED',error=type(e).__name__+': '+str(e),finished_utc=now());put(rec,d);raise
if __name__=='__main__':main()
