#!/usr/bin/env python3
"""Run real preparation checks and negative input controls without running IQ-TREE."""
from pathlib import Path
import importlib.util, json, tempfile

w=Path(__file__).resolve().parent;inputs=w.parent/'inputs'
s=importlib.util.spec_from_file_location('cpdriver',w/'run_cp_trees.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
m,c=d.verify_inputs(inputs);checks=[]
for mode in ['pilot','full']:
    commands=[]
    for p in c['sequential_panels']:
        cmd=d.argv(p,inputs,w.parent/(mode+'_01'),mode,c);commands.append(cmd)
        assert '-keep-ident' in cmd and cmd[cmd.index('-T')+1]=='4'
        if mode=='pilot':assert '-B' not in cmd and '-alrt' not in cmd and '-bnni' not in cmd and cmd[cmd.index('-n')+1]=='10'
        else:assert cmd[cmd.index('-B')+1]=='1000' and cmd[cmd.index('-alrt')+1]=='1000' and '-bnni' in cmd
    (w/('DRY_RUN_'+mode+'.json')).write_text(json.dumps({'status':'NO_INFERENCE_COMMANDS_ONLY','commands':commands},indent=2)+'\n')
    checks.append(mode+'_commands_and_support_policy_PASS')
try:d.cgroup_check(c)
except AssertionError:checks.append('unbounded_session_refused_PASS')
else:raise AssertionError('Run this negative control from the ordinary unbounded SSH session')
with tempfile.TemporaryDirectory(prefix='dt_cp9594_input_tamper_') as td:
    q=Path(td)
    for f in inputs.iterdir():
        if f.name=='CONFIG.json':q.joinpath(f.name).write_text(f.read_text()+' ')
        else:q.joinpath(f.name).symlink_to(f)
    try:d.verify_inputs(q)
    except AssertionError:checks.append('input_CONFIG_tamper_rejected_PASS')
    else:raise AssertionError('Tampered input accepted')
from Bio import Phylo
checks.append('Bio_Phylo_runtime_available_PASS')
receipt={'status':'PREPARATION_TESTS_PASS_NO_IQTREE_EXECUTION','checks':checks,
         'identical_groups':{p:v['identical_sequence_groups'] for p,v in m['panels'].items()},
         'new_tip_callable':{p:v['new_tip_callable'] for p,v in m['panels'].items()}}
(w/'PREPARATION_TESTS.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
