"""Dry-run only: create marked non-data placeholders outside the repository."""
import importlib.util, json, tempfile, subprocess, os
from pathlib import Path
R=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('stage_engine',R/'workflow/portable/stages.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
cat=json.loads((R/'workflow/research/PORTABLE_STAGES.json').read_text())
assert len(cat['stages'])==22 and len({x['script'] for x in cat['stages']})==22
results=[]
for name in cat['workflow_templates']:
 original=R/name;d=json.loads(original.read_text());base=original.parent
 with tempfile.TemporaryDirectory(prefix='portable manifest plan ') as temp:
  temp=Path(temp)
  def relocate(v):
   p=(base/v).resolve()
   if p.is_relative_to(R/'config'):return str(p)
   if p.is_relative_to(R/'scripts'):return str(p)
   return str(temp/p.relative_to(R))
  d['receipts']=relocate(d['receipts']);d['frozen_dirs']=[relocate(p) for p in d['frozen_dirs']]
  for s in d['stages'].values():
   s['config']=str((base/s['config']).resolve())
   for k in ['inputs','outputs','tools','input_directories','generated_input_directories','output_directories']:s[k]=[relocate(p) for p in s.get(k,[])]
   for k in ['bind_inputs','bind_fields']:s[k]={role:relocate(p) for role,p in s.get(k,{}).items()}
  generated={p for s in d['stages'].values() for p in s['outputs']}
  for s in d['stages'].values():
   for folder in s['input_directories']:Path(folder).mkdir(parents=True,exist_ok=True)
   for p in s['inputs']+s['tools']:
    if p in generated or Path(p).is_file():continue
    target=Path(p);target.parent.mkdir(parents=True,exist_ok=True);target.write_text('WORKFLOW PLANNING PLACEHOLDER ONLY. NOT SCIENTIFIC DATA.\n')
  f=temp/'manifest.json';f.write_text(json.dumps(d));loaded=m.load_manifest(f)
  assert len(loaded['stages'])==len(d['stages'])
  check=subprocess.run(['snakemake','--snakefile',str(R/'workflow/portable/Snakefile'),'--config','manifest='+str(f),'--cores','8','--dry-run','--quiet'],cwd=temp,text=True,capture_output=True,env={**os.environ,'XDG_CACHE_HOME':str(temp/'cache')})
  if check.returncode:raise RuntimeError(check.stdout+'\n'+check.stderr)
  assert not Path(d['receipts']).exists()
  assert not any(Path(p).exists() for p in generated)
  results.append({'template':name,'stage_instances':len(loaded['stages']),'status':'DAG_STRUCTURE_PASSED_WITH_EXTERNAL_PLACEHOLDERS_ONLY','snakemake_dry_run_exit':check.returncode,'scientific_preflight':'NOT_RUN','scientific_execution':'NOT_RUN'})
report={'catalog_unique_scripts':22,'configuration_variants':24,'templates':results,'all_passed':True,'scope':'Manifest loading and DAG validation only. No stage execution; no research inputs created or altered. Placeholder data exist only in temporary directories and are removed.'}
print(json.dumps(report,indent=2))
