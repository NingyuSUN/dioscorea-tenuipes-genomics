from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import copy,hashlib,importlib.util,json,sys

P=Path('workspace/historical_research/ils_execution_01/msc_model_03/inference')
spec=importlib.util.spec_from_file_location('coordinator',P/'reestimate_panel.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
projection=m.load_vendor('original_quartet_projector');parser=m.load_vendor('engineering_parser')
checks=[]
def record(name):checks.append({'check':name,'pass':True})
def reject(call):
    try:call()
    except (ValueError,AssertionError):return
    raise AssertionError('Expected rejection')
with TemporaryDirectory(prefix='dt_msc03_artificial_') as td:
    root=Path(td);fa=root/'toy.faa';true=root/'truth.nwk';seal=root/'toy_provenance.json';seal.write_text('{"artificial":true}\n')
    taxa=projection.TAXA+['Toy_'+str(i) for i in range(8)];h,t,g,z=taxa[:4];tail=','.join(taxa[4:])
    nwk=f'(({h},{t}),({g},{z}),({tail}));';other=f'(({h},{g}),({t},{z}),({tail}));';true.write_text(nwk)
    fa.write_text(''.join('>'+n+'\nACDE\n' for n in taxa))
    names=list(m.FIRST)+[str(1000+i)+'at3193' for i in range(170)]
    base=['/not-executable/pinned-iqtree-placeholder','-s',str(fa),'-m','MFP','-mset','JTT,WAG,LG','-mfreq','F','-mrate','G4','-B','1000','-bnni','-alrt','1000','-T','1','-seed','123','-safe','--prefix',str(root/'old-prefix'),'--boot-trees']
    jobs=[dict(gene=gene,input_fasta=str(fa),input_sha256=m.sha(fa),expected_inference_taxa=taxa,original_full_length=4,argv=base,original_argv=base,seed=123,output_prefix=str(root/'old-prefix'),true_tree_path=str(true),true_tree_sha256=m.sha(true),source_bindings=[dict(role='artificial_fixture',path=str(seal),sha256=m.sha(seal))]) for gene in names]
    data=dict(gene_order=names,taxa12=taxa,jobs=jobs);jobfile=root/'jobs.json'
    def write(d):jobfile.write_text(json.dumps(d))
    write(data)
    with patch.object(m.subprocess,'Popen',side_effect=AssertionError('No subprocess allowed in local tests')) as popen:
        _,planned,_=m.prepare(jobfile,root/'outputs');assert len(planned)==173 and not popen.called
    record('synthetic_plan_validates_without_any_subprocess')
    fa.write_text(fa.read_text().replace('ACDE','ACDF',1));reject(lambda:m.prepare(jobfile,root/'outputs'));fa.write_text(''.join('>'+n+'\nACDE\n' for n in taxa));record('single_residue_source_tampering_rejected')
    bad=copy.deepcopy(data);bad['jobs'][2]['expected_inference_taxa']=taxa[:-1];write(bad);reject(lambda:m.prepare(jobfile,root/'outputs'));record('wrong_leaf_set_rejected')
    bad=copy.deepcopy(data);bad['jobs'][0]['argv'][bad['jobs'][0]['argv'].index('-B')+1]='100';write(bad);reject(lambda:m.prepare(jobfile,root/'outputs'));record('shortened_bootstrap_and_original_command_drift_rejected')
    bad=copy.deepcopy(data);bad['jobs'][0]['argv']=bad['jobs'][0]['argv']+['-keep-ident'];write(bad);reject(lambda:m.prepare(jobfile,root/'outputs'));record('undeclared_command_change_rejected_even_if_original_settings_retained')
    bad=copy.deepcopy(data);bad['gene_order'][1]=bad['gene_order'][0];write(bad);reject(lambda:m.prepare(jobfile,root/'outputs'));record('duplicated_gene_identity_rejected')
    write(data);seal.write_text('{"artificial":false}\n');reject(lambda:m.prepare(jobfile,root/'outputs'));seal.write_text('{"artificial":true}\n');record('generation_seal_drift_rejected')
    _,planned,_=m.prepare(jobfile,root/'outputs');job=planned[0];prefix=Path(job['output_prefix']);prefix.parent.mkdir(parents=True)
    prefix.with_suffix('.treefile').write_text(other);prefix.with_suffix('.ufboot').write_text((other+'\n')*950+(nwk+'\n')*50)
    with patch.object(m.subprocess,'Popen',side_effect=AssertionError('No IQ-TREE used for artificial output check')):
        r=m.summarize_gene(job,projection,parser)
    assert r['true_focal_split']=='DT_TOK' and r['point_split']=='DT_GRA' and r['bootstrap_counts']['DT_GRA']==950 and r['bootstrap_supported_split_at_950']=='DT_GRA' and r['full_unrooted_RF_on_retained_taxa']==4
    record('known_truth_estimate_disagreement_recorded_not_treated_as_run_failure')
    prefix.with_suffix('.ufboot').write_text((other+'\n')*999);reject(lambda:m.summarize_gene(job,projection,parser));record('incomplete_bootstrap_output_rejected')
    prefix.with_suffix('.ufboot').write_text((other+'\n')*1000);prefix.with_suffix('.treefile').write_text(other.replace('Toy_0','Unexpected_leaf'));reject(lambda:m.summarize_gene(job,projection,parser));record('estimated_tree_wrong_label_rejected')
    missing=copy.deepcopy(job);missing['expected_inference_taxa']=[v for v in taxa if v!=t];missing['output_prefix']=str(root/'missing-focal'/'inferred');mp=Path(missing['output_prefix']);mp.parent.mkdir()
    missing_tree=f'({h},({g},{z}),({tail}));';mp.with_suffix('.treefile').write_text(missing_tree);mp.with_suffix('.ufboot').write_text((missing_tree+'\n')*1000)
    r=m.summarize_gene(missing,projection,parser)
    assert not r['focal_eligible'] and r['bootstrap_counts'] is None and r['bootstrap_denominator'] is None and r['saved_bootstrap_trees']==1000
    record('missing_focal_taxon_has_no_quartet_denominator_despite_saved_full_trees')
assert projection.fixtures()==6;record('six_original_quartet_projection_fixtures')
result=dict(status='PASS_ARTIFICIAL_INPUT_OUTPUT_CHECKS_ONLY',checks=checks,check_count=len(checks),actual_IQTREE_runs=0,actual_MSC_runs=0,actual_server_access=False,not_validated=['systemd runtime scopes','real msprime or AliSim API','full real173 panel integration','scientific model adequacy or ILS calibration'])
(P/'LOCAL_ARTIFICIAL_TESTS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
