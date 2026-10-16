#!/usr/bin/env python3
"""Public accession discovery through a user-supplied Entrez API adapter and ENA.
No raw reads downloaded. Stores compact API responses and an accession inventory.
"""
import argparse,datetime,json,subprocess,urllib.parse,urllib.request,csv
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--entrez-script',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True);queries=[]
def call(name,q):
 r=subprocess.run(['python3',str(a.entrez_script)],input=json.dumps(q),capture_output=True,text=True,check=True);data=json.loads(r.stdout);(a.out/(name+'.json')).write_text(r.stdout)
 queries.append({'name':name,'request':q,'ok':data.get('ok'),'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});assert data.get('ok'),data
 return data
for db in ['taxonomy','sra','nuccore','assembly']:
 term='Dioscorea tenuipes[Scientific Name]' if db=='taxonomy' else 'txid64476[Organism:exp]'
 r=call(db+'_search',{'endpoint':'esearch','params':{'db':db,'term':term,'retmode':'json','retmax':10},'record_path':'esearchresult','max_items':20})
 if db=='nuccore':
  count=int(r['summary']['count']);ids=r['summary']['idlist']
  for start in range(10,count,10):
   z=call('nuccore_page_'+str(start),{'endpoint':'esearch','params':{'db':db,'term':term,'retmode':'json','retstart':start,'retmax':10},'record_path':'esearchresult','max_items':20});ids+=z['summary']['idlist']
  d=call('nucleotide_summary',{'endpoint':'esummary','params':{'db':'nuccore','id':','.join(ids),'retmode':'json'},'record_path':'result','max_items':max(60,count+1),'max_depth':8})['summary']
  rows=[{'accession':d[i]['accessionversion'],'length_bp':d[i]['slen'],'title':d[i]['title'],'sample_annotation':d[i].get('subname','')} for i in d['uids']]
  with (a.out/'nucleotide_inventory.tsv').open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
call('new_wgs_biosample',{'endpoint':'efetch','params':{'db':'biosample','id':'SAMN48541647','retmode':'xml'},'response_format':'xml','max_items':10,'max_depth':12})
call('PZ211775_source',{'endpoint':'efetch','params':{'db':'nuccore','id':'PZ211775.1','rettype':'gb','retmode':'xml'},'response_format':'xml','record_path':'GBSet.GBSeq','max_items':30,'max_depth':10})
for name,term in [('sra_free_text','"Dioscorea tenuipes"'),('PZ211775_sra','"241029-ANU-P002" OR "PZ211775"')]:
 call(name,{'endpoint':'esearch','params':{'db':'sra','term':term,'retmode':'json','retmax':10},'record_path':'esearchresult','max_items':20})
fields='run_accession,study_accession,sample_accession,secondary_sample_accession,experiment_accession,scientific_name,library_strategy,library_source,library_selection,instrument_model,read_count,base_count,fastq_ftp'
u='https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode({'result':'read_run','query':'tax_eq(64476)','fields':fields,'format':'tsv','limit':0})
(a.out/'ena_runs.tsv').write_bytes(urllib.request.urlopen(u,timeout=60).read());queries.append({'name':'ENA_read_run','query_url':u,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
(a.out/'query_manifest.json').write_text(json.dumps(queries,indent=2)+'\n')
