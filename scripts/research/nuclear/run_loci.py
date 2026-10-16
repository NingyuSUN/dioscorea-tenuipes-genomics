"""Bounded concurrent locus inference. Finished receipt required for reuse; no automatic retries."""
import argparse,concurrent.futures,fcntl,hashlib,json,os,re,subprocess,time,traceback
from pathlib import Path
from datetime import datetime,timezone
from Bio import Phylo
R=Path(__file__).resolve().parents[1];IQ='iqtree2';stamp=lambda:datetime.now(timezone.utc).isoformat();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def run_one(row):
 g=row['gene'];dest=R/'gene_trees'/g;receipt=dest/'RECEIPT.json'
 if receipt.exists():
  old=json.loads(receipt.read_text());assert old['status']=='COMPLETE' and old['input_sha256']==row['input_sha256'] and sha(dest/'gene.treefile')==old['tree_sha256'];return old
 if dest.exists():raise RuntimeError('Incomplete or occupied locus output requires inspection: '+g)
 dest.mkdir(parents=True);inp=R/'loci'/(g+'.faa');assert sha(inp)==row['input_sha256'];prefix=dest/'gene';cmd=[IQ,'-s',str(inp),'-m','MFP','-mset','JTT,WAG,LG','-mfreq','F','-mrate','G4','-B','1000','-bnni','-alrt','1000','-T','1','-seed',str(row['seed']),'-safe','--prefix',str(prefix),'--boot-trees'];started=time.monotonic();s={'gene':g,'status':'RUNNING','input_sha256':sha(inp),'command':cmd,'started_utc':stamp()};receipt.write_text(json.dumps(s,indent=2)+'\n')
 with (dest/'driver.log').open('w') as log:
  child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT);s['pid']=child.pid;receipt.write_text(json.dumps(s,indent=2)+'\n');code=child.wait()
 s.update(exit_code=code,finished_utc=stamp(),elapsed_seconds=time.monotonic()-started,status='FAILED' if code else 'VALIDATION_PENDING');receipt.write_text(json.dumps(s,indent=2)+'\n')
 if code:raise RuntimeError(g+' inference failed')
 tree=Phylo.read(dest/'gene.treefile','newick');tips=[t.name for t in tree.get_terminals()];assert len(tips)==len(set(tips)) and set(tips)==set(row['retained_taxa']);text=(dest/'gene.log').read_text();assert 'Total wall-clock time' in text;report=(dest/'gene.iqtree').read_text();model=re.search(r'Best-fit model according to BIC:\s*(\S+)',report);assert model
 s.update(status='COMPLETE',tree_sha256=sha(dest/'gene.treefile'),model=model.group(1),tip_count=len(tips),log_warnings=[l for l in text.splitlines() if re.search(r'WARNING|ERROR|underflow|overflow',l,re.I)]);receipt.write_text(json.dumps(s,indent=2)+'\n');return s
def main():
 global R,IQ
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--iqtree',default='iqtree2');p.add_argument('--workers',type=int,default=3);p.add_argument('--mode',choices=['pilot','all'],required=True);a=p.parse_args();R=a.root.resolve();IQ=a.iqtree
 if a.workers < 1:p.error('--workers must be at least 1')
 lock=(R/'loci.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);manifest=json.loads((R/'metadata/LOCUS_MANIFEST.json').read_text());rows=manifest['loci'];rows=[r for r in rows if a.mode=='all' or r['gene'] in manifest['pilot_genes']];statuspath=R/('LOCUS_'+a.mode.upper()+'_STATUS.json');assert not statuspath.exists();s={'status':'RUNNING','mode':a.mode,'pid':os.getpid(),'started_utc':stamp(),'requested_loci':len(rows),'completed':[],'failed':[],'workers':a.workers,'wall_time_limit':None}
 def save():
  tmp=statuspath.with_suffix('.tmp');tmp.write_text(json.dumps(s,indent=2)+'\n');tmp.replace(statuspath)
 save()
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  future={pool.submit(run_one,row):row['gene'] for row in rows}
  for f in concurrent.futures.as_completed(future):
   try:r=f.result();s['completed'].append({'gene':r['gene'],'seconds':r['elapsed_seconds'],'model':r['model']})
   except Exception as e:s['failed'].append({'gene':future[f],'error':str(e),'traceback':traceback.format_exc()})
   save()
 s.update(status='COMPLETE_REVIEW_PENDING' if not s['failed'] else 'PARTIAL_FAILURE_REVIEW_REQUIRED',finished_utc=stamp());save();print(json.dumps(s,indent=2))
if __name__=='__main__':main()
