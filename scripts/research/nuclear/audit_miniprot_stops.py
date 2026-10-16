"""Recover stop-codon evidence before BUSCO removes '*' from predicted proteins."""
from pathlib import Path
import json,hashlib
from Bio import SeqIO

def audit(run):
 run=Path(run);seqdir=run/'run_embryophyta_odb10/busco_sequences/single_copy_busco_sequences'
 identifiers={}
 for p in seqdir.glob('*.faa'):
  records=list(SeqIO.parse(p,'fasta'));assert len(records)==1
  assert records[0].id not in identifiers
  identifiers[records[0].id]=p.stem
 evidence={};raw=run/'logs/miniprot_align_embryophyta_odb10_out.log'
 with raw.open() as f:
  for line in f:
   if not line.startswith('##PAF'):continue
   fields=line.rstrip().split('\t')[1:]
   if fields[5]=='*':continue
   gene='{}|{}:{}-{}|{}'.format(fields[0],fields[5],fields[7],fields[8],fields[4])
   if gene not in identifiers:continue
   sta=next(f);assert sta.startswith('##STA\t'),sta[:80]
   seq=sta.strip().split('\t')[1].upper()
   row={'marker':identifiers[gene],'gene_id':gene,'internal_stop_codons':seq.strip('*').count('*'),'pre_strip_length':len(seq)}
   if gene in evidence:assert evidence[gene]==row
   evidence[gene]=row
 unbound=set(identifiers)-set(evidence)
 with raw.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
 result={'raw_alignment_sha256':sha,'single_copy_genes_checked':len(identifiers),
         'excluded_markers':sorted({r['marker'] for r in evidence.values() if r['internal_stop_codons']}),
         'excluded_unbound_markers':sorted({identifiers[i] for i in unbound}),
         'unbound_postprocessed_gene_ids':sorted(unbound),
         'records':list(evidence.values()),'status':'RAW_TRANSLATION_STOP_AUDIT_COMPLETE'}
 (run/'miniprot_stop_audit.json').write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':
 import sys
 d=audit(sys.argv[1]);print(json.dumps({k:v for k,v in d.items() if k!='records'}))
