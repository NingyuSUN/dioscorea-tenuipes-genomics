#!/usr/bin/env python3
"""Meaningful integrity/union/order/failure controls; never run IQ-TREE."""
from pathlib import Path
import json,sys,tempfile
from source_au_common import *
import run_source_au as run

def reject(fn):
    try:fn()
    except (AssertionError,ValueError,RuntimeError):return
    raise AssertionError('Expected rejection did not occur')

def main():
    w=Path(__file__).resolve().parent;inputs=w.parent/'inputs';m,c=run.verify_inputs(inputs);checks=[]
    a=fasta(inputs/'OQ_joint.fa');b=fasta(inputs/'Kita_joint.fa')
    assert list(a)==list(b)==c['tips'] and all(len(s)==122234 for s in [*a.values(),*b.values()])
    assert all(a[k]==b[k] for k in a if k!='tokoro') and sum(x!=y for x,y in zip(a['tokoro'],b['tokoro']))==683
    checks.append('actual_joint_matrices_same9_rows_and683_source_changes_PASS')
    old=json.loads((inputs/'historical.CANDIDATE_MANIFEST.json').read_text())['candidates']
    rows=[{'newick':x['newick'],'origin':str(i)} for i,x in enumerate(old)];rows.append(dict(rows[0],origin='duplicate_other_source'))
    union=freeze_union(rows,c['tips'],c['constraints']);rev=freeze_union(list(reversed(rows)),c['tips'],c['constraints'])
    assert len(union)==2 and [e['topology_sha256'] for e in union]==[e['topology_sha256'] for e in rev]
    assert [e['newick'] for e in union]==[e['newick'] for e in rev]
    checks.append('deduplicated_union_order_independent_of_source_order_PASS')
    reject(lambda:freeze_union(rows*7,c['tips'],c['constraints']))
    reject(lambda:freeze_union([dict(rows[0],newick=rows[0]['newick'].replace('hime2','wrong_tip')),rows[1]],c['tips'],c['constraints']))
    checks.extend(['candidate_budget_refused_PASS','candidate_label_mismatch_refused_PASS'])
    assert sig_radius('-2',6)==Decimal('.000005') and sig_radius('-0.002',6)==Decimal('.000000005')
    checks.append('six_significant_digit_omitted_zero_precision_PASS')
    fake=[{'hypothesis':'unconstrained','logL_search':-100.0},{'hypothesis':'HT_cherry','logL_search':-99.98}]
    assert run.rescue_needed(fake,c)[0]
    fake[1]['logL_search']=-99.995;assert not run.rescue_needed(fake,c)[0]
    checks.append('predefined_rescue_margin_trigger_PASS')
    with tempfile.TemporaryDirectory(prefix='dt_sourceau_tests_') as td:
        td=Path(td);malformed=td/'bad.sitelh';malformed.write_text('2 2\nTree1 -2 -2\nTree1 -2 -2\n')
        reject(lambda:sitelh(malformed));checks.append('candidate_sitelh_row_order_mismatch_refused_PASS')
        prefix=td/'failed_job';reject(lambda:run.run_command([sys.executable,'-c','raise SystemExit(7)'],prefix))
        receipt=json.loads(prefix.with_suffix('.runtime.json').read_text());assert receipt['status']=='FAILED_NO_AUTO_RESTART' and receipt['exit_code']==7
        reject(lambda:run.run_command([sys.executable,'-c','pass'],prefix));checks.append('child_failure_receipt_and_existing_prefix_rejected_PASS')
    # Replay the exact completed old output binding; no model calculations.
    entries=[]
    for e in old:entries.append({**e,'topology_sha256':topology_id(e['newick'])})
    hist=Path('workspace/wgs/dt_research_p57_20260930/topology/full_01/cp/au')
    audit=run.audit_au(hist,entries,124448);assert len(audit['rows'])==2
    checks.append('real_IQTREE24_trees_sitelh_order_rounding_replay_PASS')
    reject(lambda:cgroup_check(c));checks.append('unbounded_session_rejected_PASS')
    report={'status':'PASS_NO_IQTREE_EXECUTED','checks':checks,'tests':len(checks),'input_receipt_sha256':sha(inputs/'PREPARED.json'),'real_old_CP_numeric_checks':audit['rows']};save(w/'TEST_RECEIPT.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':main()
