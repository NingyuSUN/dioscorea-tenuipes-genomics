#!/usr/bin/env python3
"""Strict receipt-gated fixed-column CP matrix preparation. Never runs a tree.
Default is read-only validation. --execute writes only after all selected inputs pass.
No missing sample is silently removed; --subset is explicit and changes output names.
"""
import argparse,csv,hashlib,json,os,sys
from pathlib import Path

class AdmissionError(ValueError):pass

def require(condition,message):
 if not condition:raise AdmissionError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
 return h.hexdigest()
def digest(data):return hashlib.sha256(data).hexdigest()
def read_json(path):return json.loads(Path(path).read_text())
def bound_file(path,expected):
 q=Path(path);require(q.is_file(),'MISSING_INPUT: '+str(q));require(sha(q)==expected,'HASH_MISMATCH: '+str(q));return q

def fasta_bytes(data,columns=None,alphabet='ACGTRYSWKMBDHVN-?'):
 require(data.endswith(b'\n'),'FASTA must end with newline')
 blocks={};seqs={};order=[];name=None;buf=[];seq=[]
 def finish():
  if name is None:return
  s=''.join(seq);require(s and set(s)<=set(alphabet),'Invalid sequence alphabet/empty sequence: '+name)
  if columns is not None:require(len(s)==columns,'UNEQUAL_COLUMNS: '+name)
  blocks[name]=b''.join(buf);seqs[name]=s
 for raw in data.splitlines(keepends=True):
  if raw.startswith(b'>'):
   finish();fields=raw[1:].strip().decode('ascii').split();require(len(fields)==1,'Expected exact FASTA tip ID without description')
   name=fields[0];require(name not in order,'DUPLICATE_TIP: '+name);order.append(name);buf=[raw];seq=[]
  else:
   require(name is not None,'Sequence before FASTA header');buf.append(raw);seq.append(raw.strip().decode('ascii'))
 finish();require(order,'Empty FASTA');return order,seqs,blocks

def choose_samples(config,subset):
 names=[x['name'] for x in config['samples']]
 if subset is None:return config['samples'],False
 requested=subset.split(',');require(requested and len(requested)==len(set(requested)),'Subset must be nonempty and unique')
 require(set(requested)<=set(names),'Unknown subset sample')
 # Use frozen biological identity order, not command-line order.
 chosen=[x for x in config['samples'] if x['name'] in requested]
 return chosen,len(chosen)!=len(names)

def sample_identity(record,spec,label):
 require(record.get('run')==spec['run'] and record.get('sample_accession')==spec['sample_accession'],label+' sample identity mismatch')

def validate_ddbj_yuga_concordance(audit,prefix,spec):
 if spec['run']!='DRR351755':return
 ddbj='ddbj' in Path(prefix).name.lower() or any('ddbj.nig.ac.jp' in r.get('url','') for r in audit.get('inputs',[]))
 if not ddbj:return
 c=audit.get('cross_source_R1_full_comparison',{})
 require(c.get('status')=='PASS_FULL_R1_ID_SEQUENCE_QUALITY_EQUAL' and c.get('records_compared')==1000000 and all(c.get(k)==0 for k in ['ID_differences','sequence_differences','quality_differences']),'Yuga DDBJ requires full1M R1 cross-source concordance; preview-only is insufficient')
 r1=[r for r in audit.get('inputs',[]) if r.get('mate')==1]
 require(len(r1)==1 and isinstance(c.get('DDBJ_R1_sha256'),str) and len(c['DDBJ_R1_sha256'])==64 and c['DDBJ_R1_sha256']==r1[0].get('file_sha256'),'Yuga cross-source receipt is not bound to this DDBJ R1 hash')
 cross_path=Path(prefix)/'FULL_R1_CROSS_SOURCE_COMPARISON.json'
 require(cross_path.is_file() and read_json(cross_path)==c,'Yuga standalone and embedded cross-source receipts differ or are missing')

def complete_stage(path,wanted):
 r=read_json(path);require(r.get('status')==wanted,'INCOMPLETE_STAGE: '+str(path)+' = '+str(r.get('status')));return r

def validate_known(config,profile):
 paths=config['source_profiles'][profile]
 for key,h in config['source_hashes'].items():bound_file(paths[key],h)
 lr=Path(paths['parent_ledger'])/'LEDGER_RECEIPT.json';bound_file(lr,config['parent_ledger_receipt_sha256']);l=read_json(lr)
 require(l['status']=='PASS_EXACT_LEDGER_AND_HIME2_SELF_READDITION','Frozen parent ledger not complete')
 for name,h in config['geometry_hashes'].items():bound_file(Path(paths['parent_ledger'])/name,h)
 k=read_json(paths['kita_receipt']);h=read_json(paths['hachi_receipt'])
 require(k['status']=='HISTORICAL_2022_PREPARED_SOURCE_SENSITIVITY_NO_TREE','Kita source receipt state mismatch')
 require(h['status']=='PREPARED_94_TIP_PARTIAL_COORDINATE_PANEL_NO_TREE' and h['new_tip']==config['hachi_tip'],'Hachi receipt state/identity mismatch')
 for kind,receipt in [('kita',k),('hachi',h)]:
  basename=Path(paths[kind]).name
  require(receipt['files'][basename]['sha256']==sha(paths[kind]),kind+' alias not bound to frozen receipt')
 base=Path(paths['original93']).read_bytes();order,rows,blocks=fasta_bytes(base,config['columns'])
 require(order==config['original_order'] and len(order)==93,'Frozen original93 order changed')
 extra=[]
 for kind in ['kita','hachi']:
  data=Path(paths[kind]).read_bytes();tips,seqs,_=fasta_bytes(data,config['columns'],'ACGTN')
  require(tips==[config[kind+'_tip']],kind+' tip mismatch');extra.append(data)
 full95=base+b''.join(extra);require(digest(full95)==config['frozen95_sha256'],'Reconstructed95 differs from verified existing95')
 return paths,base,full95

def relative_bound(record_path,sample_root,original_root):
 p=Path(record_path)
 try:rel=p.relative_to(original_root)
 except ValueError:raise AdmissionError('Sample provenance path escapes expected run root: '+str(p))
 require('..' not in rel.parts,'Parent traversal in sample provenance')
 return sample_root/rel

def validate_evidence(fixed,normalized,ledger_path,evidence_path,normalized_path):
 # This complements hash checks: N cannot become a match/called base at a gap or failed position.
 with Path(normalized_path).open() as handle:nr=list(csv.DictReader(handle,delimiter='\t'))
 require(len(nr)==len(normalized),'Normalized evidence row count mismatch')
 for i,e in enumerate(nr,1):
  require(int(e['normalized_singleIR_position_1based'])==i and e['projected_base']==normalized[i-1],'Normalized evidence/base mismatch')
  if normalized[i-1] in 'ACGT':
   tokens=e['query_evidence'].split(';') if e['query_evidence'] else []
   require(e['status']=='PASS' and tokens and int(e['unique_query_bases'])==len(tokens),'Callable base lacks PASS query evidence')
   require(all(x.rsplit(':',1)[1]==normalized[i-1] for x in tokens),'Query evidence base mismatch')
  else:require(e['status']!='PASS','N cannot be PASS')
 with Path(ledger_path).open() as lf,Path(evidence_path).open() as ef:
  lrows=csv.DictReader(lf,delimiter='\t');erows=csv.DictReader(ef,delimiter='\t');count=0
  from itertools import zip_longest
  for count,(l,e) in enumerate(zip_longest(lrows,erows),1):
   require(l is not None and e is not None,'Fixed evidence/ledger row count mismatch')
   require(int(l['frozen_column_1based'])==count==int(e['frozen_column_1based']),'Column order mismatch')
   pos=l['normalized_singleIR_position_1based'];base=normalized[int(pos)-1] if pos else 'N'
   require(e['normalized_singleIR_position_1based']==pos and e['region']==l['region'],'Frozen coordinate geometry mismatch')
   require(base==fixed[count-1]==e['edo_base'],'REFERENCE_FILL_OR_FIXED_PROJECTION_MISMATCH')
   if not pos:require(base=='N' and e['status']=='hime2_gap_no_coordinate' and not e['query_evidence'],'Reference-gap position must remain unknown')
   else:require(e['query_evidence']==nr[int(pos)-1]['query_evidence'] and e['status']==nr[int(pos)-1]['status'],'Fixed evidence differs from normalized evidence')
 require(count==len(fixed),'Fixed evidence length mismatch')

def validate_sample(config,paths,sample_parent,spec,original93):
 root=Path(sample_parent)/spec['run'];proj=root/'fixed_projection';ledger=root/'fixed_ledger_v2';orig=Path(config['default_sample_parent'])/spec['run']
 for p in [proj/'SAMPLE_BINDING.json',proj/'PROJECTION_RECEIPT.json',ledger/'LEDGER_RECEIPT.json',root/'denovo_readback/STATUS.json',root/'pilot_1m_qualityonly/STATUS.json',root/'denovo_1m/STATUS.json']:
  require(p.is_file(),'MISSING_INPUT: '+str(p))
 b=read_json(proj/'SAMPLE_BINDING.json');p=read_json(proj/'PROJECTION_RECEIPT.json');lr=ledger/'LEDGER_RECEIPT.json';l=read_json(lr)
 sample_identity(b,spec,'Projection binding');require(b['actual_tip']==spec['tip'],'Projection actual tip mismatch')
 require(b.get('all_original93_unchanged') is True and b.get('no_tree_run') is True,'Sample binding changed original93/tree status')
 require(b['parent_ledger_receipt_sha256']==config['parent_ledger_receipt_sha256'],'Parent ledger changed')
 require(b['fresh_ledger_receipt_sha256']==sha(lr)==p['ledger_receipt_sha256'],'Fresh ledger not bound')
 require(l['status']=='PASS_EXACT_LEDGER_AND_HIME2_SELF_READDITION' and l['frozen_columns']==config['columns'] and l['frozen_taxa']==93,'Fresh ledger state changed')
 require(p['status']=='PREPARED_94_TIP_PARTIAL_COORDINATE_PANEL_NO_TREE' and p['new_tip']==spec['tip'],'Completed projection schema/identity mismatch')
 require(p['source93_rows_unchanged'] is True and p['control93_byte_identical_to_frozen_input'] is True and p['occupancy_recomputed'] is False,'Original rows or occupancy changed')
 require(p['columns_both']==config['columns'] and b['fixed_denominator']==config['columns'] and b['normalized_denominator']==128097,'Projection denominator mismatch')
 for name,r in l['files'].items():
  require(Path(name).name==name,'Unsafe ledger filename');bound_file(ledger/name,r['sha256'])
 for name,h in config['geometry_hashes'].items():bound_file(ledger/name,h)
 for name,r in p['files'].items():
  require(Path(name).name==name,'Unsafe projection filename');bound_file(proj/name,r['sha256'])
 require((proj/'control93.fixed_columns.fa').read_bytes()==original93,'Old93 bytes changed in sample projection')
 qc=complete_stage(root/'pilot_1m_qualityonly/STATUS.json','COMPLETE_REVIEW_PENDING')
 assembly=complete_stage(root/'denovo_1m/STATUS.json','COMPLETE_REVIEW_PENDING')
 readback=complete_stage(root/'denovo_readback/STATUS.json','COMPLETE_READBACK_PARTIAL_CONTIGS')
 binding_path=root/'pilot_1m_qualityonly/INPUT_OUTPUT_BINDING.json';binding=read_json(binding_path);sample_identity(binding,spec,'Quality binding')
 require(b['quality_binding_sha256']==sha(binding_path)==assembly['quality_binding_sha256'],'Quality-to-assembly/projection binding mismatch')
 raw_paths=[Path(n) for n in binding['raw_inputs']]
 require(len(raw_paths)==2 and {q.name for q in raw_paths}=={'prefix_R1.fastq.gz','prefix_R2.fastq.gz'} and len({q.parent for q in raw_paths})==1,'Raw paired-prefix path set mismatch')
 prefix_original=raw_paths[0].parent;prefix=relative_bound(prefix_original,root,orig)
 prefix_status=complete_stage(prefix/'STATUS.json','COMPLETE_REVIEW_PENDING');sample_identity(prefix_status,spec,'Prefix status')
 pa=read_json(prefix/'PAIRED_PREFIX_AUDIT.json');sample_identity(pa,spec,'Prefix audit')
 validate_ddbj_yuga_concordance(pa,prefix,spec)
 require(pa['matched_ids'] is True and pa['records_per_mate']==1000000 and prefix_status['records_each']==1000000,'Prefix pair identity/count not complete')
 require(binding['prefix_audit_sha256']==sha(prefix/'PAIRED_PREFIX_AUDIT.json'),'Prefix audit changed')
 require(len(pa['inputs'])==2 and {r['mate'] for r in pa['inputs']}=={1,2} and len({r['normalized_id_sha256'] for r in pa['inputs']})==1,'Prefix mate ID hashes differ')
 for r in pa['inputs']:
  require(r['records']==1000000 and r['local_recompressed_gzip_EOF_CRC_verified'] is True,'Prefix FASTQ validation incomplete')
  bound_file(prefix/('prefix_R'+str(r['mate'])+'.fastq.gz'),r['file_sha256'])
 frozen=read_json(root/'denovo_readback/FROZEN_PARAMETERS.json')
 for name,r in {**binding['raw_inputs'],**binding['clean_outputs']}.items():bound_file(relative_bound(name,root,orig),r['sha256'])
 for name,r in binding['clean_outputs'].items():require(assembly['source_input_hashes'][name]==r==frozen['source_input_hashes'][name],'Clean input differs across stages')
 contig_path=root/'denovo_1m/assembly/final.contigs.fa';bound_file(contig_path,assembly['contigs_sha256'])
 require(frozen['source_input_hashes'][str(orig/'denovo_1m/assembly/final.contigs.fa')]['sha256']==assembly['contigs_sha256'],'Readback assembly binding mismatch')
 for record_path,h in p['readback_hashes'].items():bound_file(relative_bound(record_path,root,orig),h)
 cfg=read_json(ledger/'CONFIG.json');sample_identity(cfg['new_sample_identity'],spec,'Ledger identity');require(cfg['projection_parameters']['edo_tip']==spec['tip'] and cfg['projection_parameters']['indel_policy_version']==2,'Projection tip/policy changed')
 # Local copied runs may relocate only run-owned source files; frozen external source paths retain their provenance.
 for key,record in cfg['sources'].items():
  q=Path(record['path'])
  if q.is_relative_to(orig):q=relative_bound(q,root,orig)
  bound_file(q,record['sha256'])
 alias=relative_bound(b['actual_output'],root,orig);bound_file(alias,b['actual_output_sha256'])
 raw=alias.read_bytes();names,seqs,_=fasta_bytes(raw,config['columns'],'ACGTN');require(names==[spec['tip']],'Alias identity mismatch')
 fixed=seqs[spec['tip']];norm_names,norm_rows,_=fasta_bytes((proj/'Edo.normalized_singleIR.partial.fa').read_bytes(),128097,'ACGTN')
 require(norm_names==names,'Normalized/fixed sample mismatch');normalized=norm_rows[spec['tip']]
 require(fasta_bytes((proj/'Edo.projected_fixed124448.fa').read_bytes(),config['columns'],'ACGTN')[1]==seqs,'Legacy-named output differs from alias')
 require(p['projection_stats']==b['projection_stats'],'Projection statistics differ from binding')
 require(sum(c in 'ACGT' for c in normalized)==p['projection_stats']['normalized_callable_ACGT'] and normalized.count('N')==p['projection_stats']['normalized_unknown_N'],'Normalized callable/unknown counts mismatch')
 require(sum(c in 'ACGT' for c in fixed)==p['projection_stats']['frozen_callable_ACGT'] and fixed.count('N')==p['projection_stats']['frozen_unknown_N'],'Callable/unknown counts mismatch')
 require(any(c in 'ACGT' for c in fixed),'All-missing taxon rejected by existing projector rule')
 validate_evidence(fixed,normalized,ledger/'FROZEN_COLUMN_LEDGER.tsv',proj/'FROZEN_QUERY_EVIDENCE.tsv',proj/'NORMALIZED_QUERY_EVIDENCE.tsv')
 # The stored 94-row sample matrix must be exactly old93 plus this same query row.
 require((proj/'panel94.fixed_columns.fa').read_bytes()==original93+raw,'Sample94 is not byte-exact old93 plus alias')
 return raw,{'name':spec['name'],'run':spec['run'],'sample_accession':spec['sample_accession'],'tip':spec['tip'],'source_scope':'nonrandom first1M pairs, readback-supported partial coordinate projection','sample_root':str(root),'projection_path':str(alias),'projection_sha256':sha(alias),'projection_receipt_path':str(proj/'PROJECTION_RECEIPT.json'),'projection_receipt_sha256':sha(proj/'PROJECTION_RECEIPT.json'),'sample_binding_sha256':sha(proj/'SAMPLE_BINDING.json'),'normalized_path':str(proj/'Edo.normalized_singleIR.partial.fa'),'normalized_sha256':sha(proj/'Edo.normalized_singleIR.partial.fa'),'ledger_path':str(ledger),'ledger_receipt_sha256':sha(lr),'callable_ACGT':sum(c in 'ACGT' for c in fixed),'unknown_N':fixed.count('N'),'author_alias':spec['author_alias'],'author829_comparison_performed_here':False}

def matrix_names(samples,partial):
 n=len(samples);suffix='_partial_'+'_'.join(x['name'] for x in samples) if partial else ''
 return {'main':'main'+str(95+n)+suffix,'without_new':'control95_without_new_tenuipes','without_hachi':'control'+str(94+n)+'_without_Hachi'+suffix}

def build_matrices(base95,extras,hachi_tip,columns):
 main=base95+b''.join(extras);order,seqs,blocks=fasta_bytes(main,columns)
 require(hachi_tip in blocks,'Hachi missing from primary');without_hachi=b''.join(blocks[n] for n in order if n!=hachi_tip)
 return {'main':main,'without_new':base95,'without_hachi':without_hachi}

def write_outputs(out,matrices,names,regions,focal,result,partial):
 require(not out.exists(),'Output already exists; never overwrite')
 out.mkdir(parents=True);result.update(output_directory_created=True,output_directory=str(out),output_files_created=[],output_files_completed=[],matrix_files_created=[],matrix_files_completed=[])
 def write_file(name,data,is_matrix=False):
  with (out/name).open('xb') as stream:
   result['output_files_created'].append(name)
   if is_matrix:result['matrix_files_created'].append(name)
   stream.write(data);stream.flush();os.fsync(stream.fileno())
  result['output_files_completed'].append(name)
  if is_matrix:result['matrix_files_completed'].append(name)
 for kind,data in matrices.items():write_file(names[kind]+'.fixed_columns.fa',data,True)
 write_file('fixed_regions.nex',regions)
 write_file('FOCAL_TIP_PROVENANCE.json',(json.dumps(focal,indent=2)+'\n').encode())
 result.update(status='PARTIAL_MATRICES_PREPARED_REVIEW_PENDING' if partial else 'MATRICES98_95_97_PREPARED_REVIEW_PENDING',matrices_written=True)
 # A successful ENTRY_RECEIPT becomes visible only after every other artifact is complete.
 pending=out/'ENTRY_RECEIPT.json.pending';final=out/'ENTRY_RECEIPT.json'
 with pending.open('x') as stream:
  json.dump(result,stream,indent=2);stream.write('\n');stream.flush();os.fsync(stream.fileno())
 require(not final.exists(),'ENTRY_RECEIPT already exists; never overwrite')
 os.replace(pending,final)

def record_failure(result,exc):
 created=result.get('output_directory_created',False)
 result.update(status='PARTIAL_OUTPUT_PRESERVED' if created else 'REJECTED_MISSING_OR_INVALID_INPUTS',error=type(exc).__name__+': '+str(exc),matrices_written=bool(result.get('matrix_files_created')),successful_entry_receipt_written=False)
 return result

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,default=Path(__file__).with_name('ENTRY_CONFIG.json'));ap.add_argument('--source-profile',choices=['local','remote'],default='local');ap.add_argument('--sample-parent',type=Path);ap.add_argument('--subset');ap.add_argument('--execute',action='store_true');ap.add_argument('--out',type=Path);a=ap.parse_args()
 c=read_json(a.config);samples,partial=choose_samples(c,a.subset);parent=a.sample_parent or Path(c['default_sample_parent']);names=matrix_names(samples,partial)
 result={'status':'INPUTS_NOT_YET_VERIFIED','tree_run':False,'scientific_review_passed':False,'config_sha256':sha(a.config),'selected_samples':[x['name'] for x in samples],'omitted_samples':[x['name'] for x in c['samples'] if x not in samples],'explicit_partial_mode':partial,'panels':names,'columns':c['columns'],'unknown_N_is_never_reference_filled':True,'all3_required_by_default':True,'source_profile':a.source_profile,'sample_parent':str(parent)}
 try:
  paths,base,base95=validate_known(c,a.source_profile);extras=[];provenance=[]
  for spec in samples:
   raw,proof=validate_sample(c,paths,parent,spec,base);extras.append(raw);provenance.append(proof)
  matrices=build_matrices(base95,extras,c['hachi_tip'],c['columns'])
  require(all(data.startswith(base) for data in matrices.values()),'Original93 bytes changed')
  result.update(status='ALL_SELECTED_INPUTS_VALIDATED_MATRICES_NOT_WRITTEN',old93_bytes_unchanged=True,control95_sha256=digest(matrices['without_new']),sample_provenance=provenance,matrix_hashes={names[k]:digest(v) for k,v in matrices.items()})
  focal=[{'tip':'DT_local_path2','name':'hime2','source':'frozen original93','row_from':paths['original93']},{'tip':c['hachi_tip'],'name':'Hachi-1','source':'DRR351751 / SAMD00444457','row_from':paths['hachi']},{'tip':c['kita_tip'],'name':'Kita1','source':'Natsume2022 author-declared public CP','row_from':paths['kita']},{'tip':'OQ525999.1','name':'tokoro OQ525999.1','source':'Hu2023 public CP','row_from':paths['original93']}]
  focal += [{'tip':x['tip'],'name':x['name'],'source':x['run']+' / '+x['sample_accession'],'included':x in samples,'row_from':next((p['projection_path'] for p in provenance if p['name']==x['name']),None),'author_alias':x['author_alias']} for x in c['samples']]
  result['focal_tip_provenance']=focal
  if a.execute:
   require(a.out is not None,'--execute requires a fresh --out directory');require(not a.out.exists(),'Output already exists; never overwrite')
   regions=(Path(paths['parent_ledger'])/'regions.nex').read_bytes()
   write_outputs(a.out,matrices,names,regions,focal,result,partial)
 except (OSError,ValueError,KeyError) as exc:
  record_failure(result,exc)
  print(json.dumps(result,indent=2));raise SystemExit(2)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
