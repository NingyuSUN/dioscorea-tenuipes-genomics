#!/usr/bin/env python3
"""Ready genomes + new published annotation + completed gracillima; fresh supported tree."""
from pathlib import Path
import sys
if not __debug__:
    raise RuntimeError('Optimized Python disables required scientific assertions; run without -O.')
_HELPER_DIR = str(Path(__file__).resolve().parent)
if _HELPER_DIR not in sys.path:
    sys.path.insert(0, _HELPER_DIR)
import datetime,hashlib,json,os,subprocess,math
from Bio import SeqIO,Phylo
R=CONFIG=TOOLS=TAXA=ANNOTATION=ANNOTATION_META=None
BUSCO_SOURCES={}
BUILDER=Path(__file__).resolve().parents[1]/'build_busco_matrices.py'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,d):
 t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
def execute():
 import fcntl
 R.mkdir(parents=True,exist_ok=False);(R/'metadata').mkdir()
 stage_busco_sources()
 put(R/'metadata/MAIN_TAXA.json',TAXA)
 put(R/'metadata/WT42_DOWNLOAD.json',ANNOTATION_META)
 lock=(R/'main.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 rec=R/'NUCLEAR_RUN.json'
 if rec.exists():raise RuntimeError('Existing nuclear run: inspect status, do not duplicate')
 d={'status':'RUNNING','stage':'preflight','started_utc':now(),'pid':os.getpid(),'commands':[],'source_panel_sha256':sha(R/'metadata/MAIN_TAXA.json'),'script_sha256':sha(Path(__file__)),'excluded_orthogroups':['2661at3193'],'runtime_timeout':None}
 put(rec,d);logs=R/'logs';logs.mkdir(exist_ok=True)
 env=dict(os.environ);env['PATH']=str(TOOLS['busco_python'].parent)+os.pathsep+env.get('PATH','')
 def run(stage,args):
  d['stage']=stage;c={'stage':stage,'argv':list(map(str,args)),'started_utc':now(),'exit_code':None};d['commands'].append(c);put(rec,d)
  with (logs/(stage+'.log')).open('w') as f:p=subprocess.run(c['argv'],stdout=f,stderr=subprocess.STDOUT,env=env)
  c.update(exit_code=p.returncode,ended_utc=now());put(rec,d)
  if p.returncode:raise RuntimeError(f'{stage} failed {p.returncode}; inspect {logs/(stage+".log")}')
 try:
  meta=ANNOTATION_META;fasta=ANNOTATION
  assert fasta.stat().st_size==meta['size'] and hashlib.md5(fasta.read_bytes()).hexdigest()==meta['computed_md5']
  with fasta.open() as handle:seqs=list(SeqIO.parse(handle,'fasta'))
  assert seqs and len({x.id for x in seqs})==len(seqs)
  put(R/'metadata/WT42_VERIFIED.json',{'status':'DOWNLOAD_MD5_AND_FASTA_PASS','source':meta,'sha256':sha(fasta),'protein_records':len(seqs),'source_compartment':'genome annotation','paper_name':'Dioscorea opposita WT42','deposited_name':'Dioscorea oppositifolia','accession':'GCA_056531125.1'})
  run('busco_WT42',[TOOLS['busco_python'],TOOLS['busco_script'],'-i',fasta,'-m','proteins','-l',TOOLS['lineage'],'--offline','--download_path',TOOLS['download_path'],'-o','Dioscorea_oppositifolia_WT42','--out_path',R/'busco','-c',str(CONFIG['threads']),'--opt-out-run-stats'])
  run('build_matrix',[TOOLS['busco_python'],BUILDER,'--work-root',R,'--mafft',TOOLS['mafft'],'--panel','genome_main_no2661','--taxa-json',R/'metadata/MAIN_TAXA.json','--threads',str(CONFIG['threads'])])
  matrix=R/'matrices/genome_main_no2661/main.faa';s=json.loads((matrix.parent/'summary.json').read_text())
  assert {'Dioscorea_tenuipes_historical','Dioscorea_gracillima'}<=set(s['taxa'])
  assert '2661at3193' not in {x['gene']for x in s['partitions']}
  assert s['main_sha256']==sha(matrix)
  put(R/'metadata/FINAL_NUCLEAR_PANEL.json',{'status':'FROZEN_BEFORE_INFERENCE','utc':now(),'taxa':s['taxa'],'excluded_input_taxa':s['excluded_taxa'],'genes':s['genes'],'columns':s['columns'],'matrix_sha256':sha(matrix),'excluded_orthogroups':['2661at3193'],'required_focal_taxa_per_locus':s['required_focal_taxa_per_locus'],'source_scope':'Existing genomic sources + gracillima genomic WGS + published WT42 genome annotation; broader WGS candidates deferred explicitly'})
  result=R/'results/genome_main_no2661';result.mkdir(parents=True,exist_ok=True);pre=result/'main'
  run('IQTREE_main',[TOOLS['iqtree'],'-s',matrix,'-seed','20260928','-T',str(CONFIG['threads']),'-mem','8000M','--prefix',pre,'-m','MFP','-mset','LG,WAG,JTT','-mrate','G4','-mfreq','F','-B','1000','-bnni','-alrt','1000','--boot-trees'])
  tree=Phylo.read(str(pre)+'.treefile','newick');names=[x.name for x in tree.get_terminals()]
  assert len(names)==len(set(names)) and set(names)==set(s['taxa'])
  assert all(x.branch_length is None or math.isfinite(x.branch_length) and x.branch_length>=0 for x in tree.find_clades())
  with Path(str(pre)+'.ufboot').open() as handle:assert sum(bool(x.strip())for x in handle)==1000
  assert 'Analysis results written to' in Path(str(pre)+'.log').read_text()
  d.update(status='INFERENCE_COMPLETE_SCIENTIFIC_REVIEW_AND_FIGURE_PENDING',stage='done',finished_utc=now(),taxa=s['taxa'],genes=s['genes'],columns=s['columns'],treefile_sha256=sha(Path(str(pre)+'.treefile')));put(rec,d)
 except Exception as e:
  d.update(status='FAILED',error=type(e).__name__+': '+str(e),finished_utc=now());put(rec,d);raise
 finally:
  lock.close()

def configure(config, base):
 from _nuclear_cli import path_value,require_file,require_directory,fresh_output,report_plan
 global R,CONFIG,TOOLS,TAXA,ANNOTATION,ANNOTATION_META,BUSCO_SOURCES
 CONFIG=dict(config);CONFIG.setdefault('threads',8);BUSCO_SOURCES={}
 if not isinstance(CONFIG['threads'],int) or CONFIG['threads']<1:raise ValueError('threads must be a positive integer')
 TOOLS={key:require_file(path_value(config,key,base),executable=key!='busco_script') for key in ['busco_python','busco_script','mafft','iqtree']}
 for key in ['lineage','download_path']:TOOLS[key]=require_directory(path_value(config,key,base))
 if TOOLS['lineage'].name!='embryophyta_odb10':raise ValueError('The frozen builder requires embryophyta_odb10')
 ANNOTATION=require_file(path_value(config,'annotation_fasta',base))
 meta_path=require_file(path_value(config,'annotation_metadata',base));ANNOTATION_META=json.loads(meta_path.read_text())
 assert ANNOTATION.stat().st_size==ANNOTATION_META['size'] and hashlib.md5(ANNOTATION.read_bytes()).hexdigest()==ANNOTATION_META['computed_md5']
 taxa_path=require_file(path_value(config,'taxa_json',base));TAXA=json.loads(taxa_path.read_text())
 if not isinstance(TAXA,list) or not TAXA:raise ValueError('taxa_json must contain a nonempty list')
 R=path_value(config,'work_root',base)
 wt42='Dioscorea_oppositifolia_WT42'
 if sum(row['taxon']==wt42 for row in TAXA)!=1:raise ValueError('Exactly one WT42 taxon row is required')
 if len({row['taxon'] for row in TAXA})!=len(TAXA):raise ValueError('Duplicate taxon label')
 required={'Dioscorea_tenuipes_historical','Dioscorea_gracillima','Oryza_sativa','Arabidopsis_thaliana'}
 if not required<={row['taxon'] for row in TAXA}:raise ValueError('Required focal taxa/outgroups are absent')
 inputs=[ANNOTATION,meta_path,taxa_path,*TOOLS.values()]
 for row in TAXA:
  if row['taxon']==wt42:row['busco_run']=str(R/'busco'/wt42)
  else:
   value=Path(row['busco_run']).expanduser();run=(value if value.is_absolute() else taxa_path.parent/value).resolve();require_directory(run)
   table=require_file(run/'run_embryophyta_odb10/full_table.tsv')
   summaries=list(run.glob('short_summary*.json'))
   if len(summaries)!=1:raise ValueError('Expected one BUSCO summary: '+str(run))
   seqdir=require_directory(run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences')
   files=[require_file(summaries[0]),table,*sorted(seqdir.glob('*.faa'))]
   if 'Gene Start\tGene End' in table.read_text():files.append(require_file(run/'logs/miniprot_align_embryophyta_odb10_out.log'))
   if Path(row['taxon']).name!=row['taxon'] or row['taxon'] in ('.','..'):raise ValueError('Taxon must be a simple nonempty label')
   dest=R/'input_busco'/row['taxon'];BUSCO_SOURCES[row['taxon']]={'source':run,'files':files,'destination':dest}
   row['busco_run']=str(dest);inputs.extend(files)
 fresh_output(R,[*inputs,*[item['source'] for item in BUSCO_SOURCES.values()]]);require_file(BUILDER)
 consumed={'annotation':ANNOTATION,'annotation_metadata':meta_path,'taxa_json':taxa_path,**TOOLS,'matrix_builder':BUILDER,'stop_audit':BUILDER.parent/'audit_miniprot_stops.py'}
 for taxon,source in BUSCO_SOURCES.items():
  consumed.update({taxon+':'+str(path.relative_to(source['source'])):path for path in source['files']})
 return report_plan('run_nuclear_main.py',consumed,
                    {'receipt':R/'NUCLEAR_RUN.json','staging':R/'metadata/INPUT_STAGING.json','panel':R/'metadata/FINAL_NUCLEAR_PANEL.json','matrix':R/'matrices/genome_main_no2661/main.faa','summary':R/'matrices/genome_main_no2661/summary.json','tree':R/'results/genome_main_no2661/main.treefile','bootstrap':R/'results/genome_main_no2661/main.ufboot'},
                    ['biopython','BUSCO','MAFFT','IQ-TREE','Linux fcntl'],threads=CONFIG['threads'],
                    note='Annotation is supplied locally and MD5/size-verified; this driver performs no network download.')

def stage_busco_sources():
 """Keep the unchanged builder's audit sidecars inside this fresh output root."""
 import shutil
 from _nuclear_cli import sha256
 records=[]
 for taxon,item in BUSCO_SOURCES.items():
  for source in item['files']:
   relative=source.relative_to(item['source']);destination=item['destination']/relative
   destination.parent.mkdir(parents=True,exist_ok=True)
   before=sha256(source);shutil.copyfile(source,destination)
   if sha256(source)!=before or sha256(destination)!=before:raise RuntimeError('BUSCO staging identity changed: '+str(source))
   records.append({'taxon':taxon,'source':str(source),'staged_path':str(destination.relative_to(R)),'sha256':before})
  # An empty single-copy directory is still meaningful input evidence.
  (item['destination']/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences').mkdir(parents=True,exist_ok=True)
 put(R/'metadata/INPUT_STAGING.json',{'status':'EXACT_INPUT_COPIES_VERIFIED','files':records,'original_inputs_modified':False})


def plan(config_path):
    """Return a read-only, JSON-compatible artifact and dependency plan."""
    config_path=Path(config_path).resolve();return configure(json.loads(config_path.read_text()),config_path.parent)


def main(argv=None):
 from _nuclear_cli import read_config
 config,base,run=read_config(__doc__,argv);print(json.dumps(configure(config,base),indent=2))
 if run:execute()
 return 0

if __name__=='__main__':raise SystemExit(main())
