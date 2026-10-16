#!/usr/bin/env python3
"""Frozen full12 conditional S0 fitter. Default validates sources; never simulates.
Only --execute invokes scipy optimization, after an actual Linux cgroup guard.
No population-history, biological rate, significance, or ILS-exclusion inference.
"""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import sys
import time

SCHEMA = 'dioscorea.conditional_s0.v1'
AA = 'ARNDCQEGHILKMFPSTWYV'
L = 'workspace/historical_research'
REPO = 'workspace/projects/dioscorea-tenuipes-reference-genome'
CONTRACT_SHA = 'b1e308c458636878d268a8d9b1ab83fee0236197991a58f530ba77eed4838fbf'
SOURCE_PINS = {
 'snapshot_manifest': (L+'/ils_forward_preparation/SOURCE_SNAPSHOT_MANIFEST.json', '24e07a9d548a4d063fefdffa8311acfe6e052131af445648c9a5f081c7e305c0'),
 'locus_manifest': (REPO+'/docs/evidence/2026-09-29/research_p50/inputs/LOCUS_MANIFEST.json', '84c3fcd80ce9d021e03e2ee925f1dc9fc77f852757cd971683d0fb4fd3f0525d'),
 'model_table': (L+'/ils_forward_preparation/PER_LOCUS_INPUT_CAPABILITY.tsv', '5b3b881ab39a3df05c4434c37502fff746ea830d978d0785f4d24aa194f06492')}
PROFILE_POINTS = 41
NUMERIC = {'CU_starts': [0.1, 1., 3.], 'CU_bounds': [0.,10.],
 'CU_optimizer': 'L-BFGS-B', 'CU_maxiter': 2000, 'CU_ftol': 1e-12, 'CU_gtol': 1e-8,
 'time_starts_s': [.25,1.,4.], 'time_bounds_s': [1/16,16.],
 'time_optimizer': 'L-BFGS-B', 'time_maxiter': 2000, 'time_ftol': 1e-12,
 'profile_points_per_axis': PROFILE_POINTS, 'boundary_flag_tolerance': 1e-7,
 'equivalent_objective_tolerance': 1e-8,
 'profile_near_optimum_delta_objective': 1e-6,
 'near_optimum_log_axis_span_nonidentified': math.log(2)}


def sha(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1048576), b''): h.update(b)
    return h.hexdigest()


def save(path, obj):
    with open(path,'x') as f: json.dump(obj,f,indent=2,allow_nan=False); f.write('\n')


def fasta(path):
    d = {}; name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in d: raise ValueError('Duplicate FASTA label '+name)
            d[name] = ''
        elif line.strip():
            if name is None: raise ValueError('Sequence before label')
            d[name] += line.strip()
    if not d or len({len(x) for x in d.values()}) != 1: raise ValueError('Nonrectangular/empty alignment')
    return d


class Tree:
    """Small explicit unrooted Newick graph; branch lengths remain AA distances.
    A degree-two Newick display root is suppressed, summing its two lengths.
    Support labels (including slash labels) are ignored, never used to filter.
    """
    def __init__(self, text):
        tokens = re.findall(r"[^\s(),:;]+|[(),:;]",text)
        self.adj = {}; self.labels = {}; self.rounding = {}; p = 0
        def node():
            nonlocal p
            n = len(self.adj); self.adj[n] = {}; children = []
            if tokens[p] == '(':
                p += 1
                while True:
                    c,b,err = node(); children.append((c,b,err))
                    if tokens[p] == ',': p += 1; continue
                    if tokens[p] != ')': raise ValueError('Malformed Newick')
                    p += 1; break
                if p < len(tokens) and tokens[p] not in [',',')',':',';']: p += 1
            else:
                label = tokens[p]; p += 1
                if label in self.labels.values(): raise ValueError('Duplicate tree tip')
                self.labels[n] = label
            b = 0.; rounding_error = 0.
            if p < len(tokens) and tokens[p] == ':':
                p += 1; printed = tokens[p]; b = float(printed); p += 1
                if not math.isfinite(b) or b < 0: raise ValueError('Invalid tree branch')
                rounding_error = .5 * 10. ** Decimal(printed).as_tuple().exponent
            for c,v,err in children:
                self.adj[n][c] = v; self.adj[c][n] = v
                self.rounding[frozenset((n,c))] = err
            return n,b,rounding_error
        root,_,_ = node()
        if p != len(tokens)-1 or tokens[p] != ';': raise ValueError('Trailing/missing Newick tokens')
        if len(self.adj[root]) == 2 and root not in self.labels:
            (a,x),(b,y) = self.adj[root].items()
            del self.adj[a][root]; del self.adj[b][root]; del self.adj[root]
            self.adj[a][b] = x+y; self.adj[b][a] = x+y
            self.rounding[frozenset((a,b))] = self.rounding.pop(frozenset((a,root))) + self.rounding.pop(frozenset((b,root)))
        self.tipnodes = {v:k for k,v in self.labels.items()}
        self.taxa = set(self.tipnodes)

    def side(self, a, b):
        out = set(); stack = [(a,b)]
        while stack:
            n,p = stack.pop()
            if n in self.labels: out.add(self.labels[n])
            stack.extend((c,n) for c in self.adj[n] if c != p)
        return frozenset(out)

    def edges(self):
        for a in sorted(self.adj):
            for b in sorted(self.adj[a]):
                if a < b: yield a,b,self.adj[a][b],self.side(a,b)

    def distance(self, ta, tb):
        target = self.tipnodes[tb]; stack = [(self.tipnodes[ta],None,0.)]
        while stack:
            n,p,d = stack.pop()
            if n == target: return d
            stack.extend((c,n,d+b) for c,b in self.adj[n].items() if c != p)
        raise ValueError('Disconnected tree')


    def distance_rounding_bound(self, ta, tb):
        target=self.tipnodes[tb]; stack=[(self.tipnodes[ta],None,0.)]
        while stack:
            n,p,err=stack.pop()
            if n==target:return err
            stack.extend((c,n,err+self.rounding[frozenset((n,c))]) for c in self.adj[n] if c!=p)
        raise ValueError('Disconnected tree')


def canonical_split(side, taxa):
    a,b = tuple(sorted(side)),tuple(sorted(set(taxa)-set(side)))
    return min(a,b)


def quartet_key(side, q):
    x = set(side)&set(q)
    return canonical_split(x,q) if len(x)==2 else None


def project_quartet(tree, q):
    if not set(q)<=tree.taxa: return None,'missing'
    observed = {quartet_key(s,q) for a,b,l,s in tree.edges() if l>0 and quartet_key(s,q) is not None}
    if len(observed)>1: raise ValueError('Incompatible induced quartet')
    return (next(iter(observed)),'resolved') if observed else (None,'unresolved')


def topology_edges(tree):
    edges = [(canonical_split(s,tree.taxa),a,b) for a,b,l,s in tree.edges() if 1<len(s)<len(tree.taxa)-1]
    return sorted(edges)


def bind_file(path, expected, path_map):
    actual = path_map.get(path, path)
    if not isinstance(actual, str) or not Path(actual).is_absolute():
        raise ValueError('Nonabsolute path mapping')
    got = sha(actual)
    if got != expected: raise ValueError('SHA mismatch '+path)
    return {'path':path, 'resolved_path':actual, 'sha256':got}


def report_tree_binding(full_newick, report_newick, source_sequences):
    """Only exact report trees or proved identical-sequence reductions qualify.
    Full receipt-bound tree remains authoritative; no source tree is rewritten.
    """
    full=Tree(full_newick);report=Tree(report_newick)
    if full.taxa != set(source_sequences): raise ValueError('Full tree/source AA labels differ')
    if not report.taxa <= full.taxa: raise ValueError('Report tree contains extra labels')
    if len(report.taxa)<4: raise ValueError('Report reduction has fewer than four tips')
    same=re.sub(r'\s+','',full_newick)==re.sub(r'\s+','',report_newick)
    missing=sorted(full.taxa-report.taxa)
    if not same and not missing: raise ValueError('Changed full-label report tree is not an identical-sequence reduction')
    maps=[]
    for absent in missing:
        identical=[t for t in sorted(report.taxa) if source_sequences[t]==source_sequences[absent]]
        if not identical: raise ValueError('Report missing leaf is not AA-identical to any retained leaf: '+absent)
        # IQ-TREE insertion must itself be at zero distance, within printed precision.
        zero=[t for t in identical if full.distance(t,absent)<=full.distance_rounding_bound(t,absent)+1e-12]
        if not zero: raise ValueError('AA-identical missing leaf has nonzero full-tree insertion distance: '+absent)
        maps.append({'missing_taxon':absent,'identical_representative':zero[0],
          'all_exact_AA_identical_report_taxa':identical,'representatives_at_zero_distance':zero,
          'AA_sequence_sha256':hashlib.sha256(source_sequences[absent].encode()).hexdigest(),
          'AA_length':len(source_sequences[absent]),'full_tree_insertion_distance':full.distance(absent,zero[0])})
    def induced_splits(tree):
        observed=set()
        for a,b,length,side in tree.edges():
            restricted=set(side)&report.taxa
            if 1<len(restricted)<len(report.taxa)-1:
                observed.add(canonical_split(restricted,report.taxa))
        return observed
    if induced_splits(full)!=induced_splits(report): raise ValueError('Report reduction unrooted splits differ')
    comparisons=[]
    for a,b in itertools.combinations(sorted(report.taxa),2):
        fd,rd=full.distance(a,b),report.distance(a,b)
        bound=full.distance_rounding_bound(a,b)+report.distance_rounding_bound(a,b)+1e-12
        difference=abs(fd-rd)
        if difference>bound: raise ValueError('Report reduction patristic distance differs: '+a+' / '+b)
        comparisons.append({'taxon_a':a,'taxon_b':b,'full_restricted_distance':fd,'report_distance':rd,
          'absolute_difference':difference,'combined_printed_rounding_bound':bound})
    return {'status':'EXACT_WHITESPACE_NORMALIZED_TREE' if same else 'PROVEN_IDENTICAL_SEQUENCE_REDUCTION',
      'full_tree_sha256_LF':hashlib.sha256((full_newick.strip()+'\n').encode()).hexdigest(),
      'report_tree_sha256_LF':hashlib.sha256((report_newick.strip()+'\n').encode()).hexdigest(),
      'missing_to_identical_representative':maps,'unrooted_restricted_splits_equal':True,
      'patristic_pair_count':len(comparisons),'max_absolute_patristic_difference':max(x['absolute_difference'] for x in comparisons),
      'rounding_rule':'sum of half-unit-in-last-printed-place along BOTH paths, plus1e-12 floating tolerance',
      'patristic_comparisons':comparisons if missing else [],
      'mechanism_interpretation':'compatible with identical-sequence deletion/reinsertion; exact backend event not asserted',
      'full_receipt_bound_tree_retained':True}


def validate_model_distribution(freq, pair_count, alpha):
    """Preserve printed zeros; do not impute unprinted positive frequencies.
    Four-decimal reported values have a total rounding bound of20*0.00005.
    IQ-TREE2.4.0 readStateFreq accepts zero; its actual AliSim run is a separate gate.
    """
    if pair_count != 20 or set(freq) != set(AA):
        raise ValueError('Invalid frequency state labels/count')
    if any(not math.isfinite(v) or v < 0 or v > 1 for v in freq.values()):
        raise ValueError('Frequency must be finite in [0,1]')
    if sum(v > 0 for v in freq.values()) < 2:
        raise ValueError('At least two positive-frequency states required')
    if not math.isfinite(alpha) or alpha <= 0:
        raise ValueError('Alpha must be finite and positive')
    total = sum(freq.values())
    if total <= 0 or abs(total - 1) > 0.001 + 1e-12:
        raise ValueError('Frequency rounding exceeds frozen0.001 plus floating-point tolerance')
    return total


def bind_sources(contract, path_map):
    """Only exact mapping keys; no basename, prefix, cwd or remote guesses."""
    ledger = []
    def bind(path, expected):
        item = bind_file(path, expected, path_map)
        ledger.append(item); return item
    pins = {k:bind(*v) for k,v in SOURCE_PINS.items()}
    snapshot = json.loads(Path(pins['snapshot_manifest']['resolved_path']).read_text())
    roles = {}
    for row in snapshot['files']:
        key = (row['role'],row.get('gene'))
        if key in roles: raise ValueError('Duplicate manifest role/gene')
        roles[key] = row
    def role(r,g=None):
        x = roles[(r,g)]; return bind(x['local_path'],x['sha256'])
    main = role('matrix'); summaryb = role('summary')
    if main['path'] != contract['panel']['matrix_path'] or main['sha256'] != contract['panel']['matrix_sha256']: raise ValueError('Matrix contract mismatch')
    if summaryb['path'] != contract['sources']['mask_summary']: raise ValueError('Summary contract mismatch')
    matrix = fasta(main['resolved_path']); taxa = list(matrix)
    if len(taxa)!=12 or len(next(iter(matrix.values())))!=68796 or any(set(v)-set(AA+'-') for v in matrix.values()): raise ValueError('Original12 matrix shape/alphabet')
    summary = json.loads(Path(summaryb['resolved_path']).read_text())
    original = json.loads(Path(pins['locus_manifest']['resolved_path']).read_text())['loci']
    parts = summary['partitions']; genes = [x['gene'] for x in original]
    if len(genes)!=173 or len(set(genes))!=173 or genes!=[x['gene'] for x in parts]: raise ValueError('Original173 order mismatch')
    if '2661at3193' in genes: raise ValueError('Globally excluded marker present')
    sources = {k:bind(v['path'],v['sha256']) for k,v in contract['sources'].items() if isinstance(v,dict) and 'sha256' in v}
    tree_lines = Path(sources['gene_tree_corpus']['resolved_path']).read_text().splitlines()
    if len(tree_lines)!=173: raise ValueError('173 tree lines required')
    modelrows = list(csv.DictReader(open(pins['model_table']['resolved_path']),delimiter='\t'))
    if [x['gene'] for x in modelrows]!=genes: raise ValueError('Model table gene order mismatch')
    entries = []; trees = []; previous_end = 0
    for i,(g,orig,part,mt,line) in enumerate(zip(genes,original,parts,modelrows,tree_lines)):
        loc,report,receipt = role('locus_AA',g),role('model_report',g),role('locus_receipt',g)
        rc = json.loads(Path(receipt['resolved_path']).read_text())
        if rc['gene']!=g or rc['status']!='COMPLETE' or rc['exit_code']!=0: raise ValueError('Original receipt not complete')
        if loc['sha256']!=orig['input_sha256'] or rc['input_sha256']!=loc['sha256']: raise ValueError('Original input binding')
        th = hashlib.sha256((line+'\n').encode()).hexdigest()
        if th != rc['tree_sha256']: raise ValueError('Gene tree line receipt mismatch '+g)
        a,b = part['start_1based'],part['end_1based']; length = b-a+1
        if a!=previous_end+1 or length!=orig['length'] or length!=part['length']: raise ValueError('Noncontiguous partition')
        if part['sha256']!=orig['source_sha256'] or int(mt['length_AA'])!=length or int(mt['matrix_start_1based'])!=a or int(mt['matrix_end_1based'])!=b: raise ValueError('Partition/model table source coordinates differ')
        previous_end=b
        seq = fasta(loc['resolved_path']); full = {t:s[a-1:b] for t,s in matrix.items()}
        present = {t:s for t,s in full.items() if set(s)!={'-'}}
        if seq!=present or set(seq)!=set(orig['retained_taxa']) or set(taxa)-set(seq)!=set(orig['all_missing_taxa_omitted']): raise ValueError('Exact original mask/source mismatch '+g)
        tr = Tree(line)
        if tr.taxa != set(seq): raise ValueError('Tree tip/AA mismatch')
        rep = Path(report['resolved_path']).read_text()
        report_trees = re.findall(r'Tree in newick format:\s*([^;]+;)',rep)
        if len(report_trees)!=1: raise ValueError('Ambiguous model report tree '+g)
        report_binding=report_tree_binding(line, report_trees[0], seq)
        model = re.search(r'Model of substitution:\s*(\S+)',rep).group(1)
        if model!=rc['model'] or model!=mt['model'] or not re.fullmatch(r'(JTT|WAG|LG)\+F\+G4',model): raise ValueError('Original model mismatch')
        alpha = float(re.search(r'Gamma shape alpha:\s*([\d.eE+-]+)',rep).group(1))
        pairs = re.findall(r'pi\(([A-Z])\)\s*=\s*([\d.eE+-]+)',rep)
        freq = {x:float(y) for x,y in pairs}
        total = validate_model_distribution(freq, len(pairs), alpha)
        if float(mt['gamma_alpha_reported'])!=alpha or json.loads(mt['frequencies_reported'])!=freq: raise ValueError('Model table parameter mismatch')
        if mt['model_report_sha256']!=report['sha256'] or mt['receipt_sha256']!=receipt['sha256'] or mt['locus_AA_sha256']!=loc['sha256']: raise ValueError('Model table SHA mismatch')
        mask = dict(main,matrix_start_1based=a,matrix_end_1based=b)
        entries.append({'gene':g,'length':length,'model_matrix':model.split('+')[0], 'alpha':alpha,
          'frequency_order':AA,'ordered_normalized_frequencies':[freq[x]/total for x in AA],
          'frequency_reported_sum':total,'frequency_normalization':'printed AA-labelled frequency / printed sum; zeros preserved; no pseudocount',
          'printed_zero_states':[x for x in AA if freq[x]==0],
          'frequency_report_precision_decimal_places':4,
          'frequency_precision_limitation':'Printed0.0000 may conceal a small positive fitted frequency; retained exactly as reported; not biological absence. Actual AliSim zero-frequency API remains to be verified on Linux.',
          'source_locus':loc,'source_model':report,'source_report':report,'source_mask':mask,
          'source_mask_summary':summaryb,'source_receipt':receipt,
          'source_tree':dict(sources['gene_tree_corpus'],line_1based=i+1,line_sha256=th),
          'report_tree_binding':report_binding,
          'retained_taxa':list(seq),'all_missing_taxa':[t for t in taxa if t not in seq],
          'mask_missing_cells':sum(v.count('-') for v in full.values()),
          'original_inference_argv':rc['command'], 'source_warning_records':rc.get('log_warnings',[])})
        trees.append(tr)
    if previous_end!=68796: raise ValueError('Partition total mismatch')
    fixed_text = Path(sources['fixed_species_topology_source']['resolved_path']).read_text()
    found = re.findall(r'Tree in newick format:\s*([^;]+;)',fixed_text)
    if len(found)!=1: raise ValueError('Ambiguous fixed species tree extraction')
    species = Tree(found[0])
    if species.taxa!=set(taxa) or len(topology_edges(species))!=9 or any(len(species.adj[n])!=3 for n in species.adj if n not in species.labels): raise ValueError('Fixed full12 topology shape')
    return {'taxa':taxa,'gene_entries':entries,'gene_trees':trees,'species':species,'species_newick':found[0],
     'source_manifest_hashes':{k:v['sha256'] for k,v in pins.items()},'source_bindings':ledger}


def quartet_observations(species, gene_trees, gene_names):
    import numpy as np
    qs = list(itertools.combinations(sorted(species.taxa),4)); edges = topology_edges(species)
    B = np.zeros((len(qs),len(edges))); concord = []
    choices = []
    for qi,q in enumerate(qs):
        options = [tuple((q[0],x)) for x in q[1:]]; choices.append(options)
        induced = set()
        for ei,(s,a,b) in enumerate(edges):
            split = quartet_key(s,q)
            if split is not None: B[qi,ei]=1; induced.add(split)
        if len(induced)!=1: raise ValueError('Species quartet not uniquely binary')
        concord.append(options.index(next(iter(induced))))
    weighted = np.zeros((len(qs),3)); raw = np.zeros((len(qs),3),dtype=int)
    rows=[]; gene_audit=[]
    for g,t in zip(gene_names,gene_trees):
        result=[project_quartet(t,q) for q in qs]; m=sum(s=='resolved' for _,s in result)
        audit={'gene':g,'resolved_quartets':m,'unresolved_quartets':sum(s=='unresolved' for _,s in result),
          'missing_quartets':sum(s=='missing' for _,s in result),'per_resolved_weight':1/m if m else None,'total_weight':1 if m else 0}
        gene_audit.append(audit)
        for qi,(split,status) in enumerate(result):
            oi=choices[qi].index(split) if status=='resolved' else None
            if oi is not None: raw[qi,oi]+=1; weighted[qi,oi]+=1/m
            rows.append((g,qi,status,oi,1/m if oi is not None else 0))
    c=np.array([weighted[i,k] for i,k in enumerate(concord)])
    d=weighted.sum(axis=1)-c
    return {'quartets':qs,'edges':edges,'B':B,'concord_indices':concord,'raw':raw,'weighted':weighted,
      'c':c,'d':d,'rows':rows,'gene_audit':gene_audit}


def cu_objective(x,B,c,d):
    import numpy as np
    t=B@x; e=np.exp(-t); p=1-2*e/3
    value=-np.sum(c*np.log(p)+d*(-t-math.log(3)))
    grad=B.T@(d-c*(2*e/3)/p)
    return float(value),grad


def optimize_cu(obs):
    import numpy as np
    from scipy.optimize import minimize
    records=[]
    for s in NUMERIC['CU_starts']:
        r=minimize(lambda x:cu_objective(x,obs['B'],obs['c'],obs['d']),np.full(len(obs['edges']),s),jac=True,
          method='L-BFGS-B',bounds=[(0,10)]*len(obs['edges']),
          options={'maxiter':2000,'ftol':1e-12,'gtol':1e-8})
        records.append({'start':s,'success':bool(r.success),'message':str(r.message),'objective':float(r.fun),
          'x':r.x.tolist(),'nit':int(r.nit),'gradient':r.jac.tolist()})
    good=[r for r in records if r['success'] and math.isfinite(r['objective'])]
    if not good: raise ValueError('No CU optimizer start converged')
    best=min(good,key=lambda r:(r['objective'],r['x']))
    return best['x'],records


def rooted_structure(species, edges, values, delta, rho):
    """Root on Arabidopsis terminal edge. Collapse exactly-zero internal edges.
    A tolerance marks boundaries only; it NEVER creates nonzero population time.
    """
    out='Arabidopsis_thaliana'; tip=species.tipnodes[out]; V=next(iter(species.adj[tip]))
    lengths={frozenset((a,b)):float(x) for (_,a,b),x in zip(edges,values)}
    internal=set(species.adj)-set(species.labels); parent={V:None}; dist={V:0.}; order=[V]
    for n in order:
        for c in sorted(species.adj[n]):
            if c not in internal or c==parent[n]: continue
            parent[c]=n;dist[c]=dist[n]+lengths[frozenset((n,c))];order.append(c)
    if set(order)!=internal: raise ValueError('Internal orientation failed')
    dmax=max(dist.values()); H=dmax+delta
    representative={}
    for n in order:
        p=parent[n]
        representative[n]=representative[p] if p is not None and lengths[frozenset((p,n))]==0 else n
    groups={r:[] for r in representative.values()}
    for n,r in representative.items():groups[r].append(n)
    ids={r:'pop_%03d'%i for i,r in enumerate(sorted(groups))}
    nodes={}
    for r,old in groups.items():
        p=parent[r]; pid='pop_ROOT' if p is None else ids[representative[p]]
        nodes[ids[r]]={'id':ids[r],'parent_id':pid,'children':[],'descendant_taxa':[],'time_CU':H-dist[r]}
    nodes['pop_ROOT']={'id':'pop_ROOT','parent_id':None,'children':[],'descendant_taxa':[],'time_CU':H+rho}
    for n,t in species.labels.items():
        pid='pop_ROOT' if t==out else ids[representative[next(iter(species.adj[n]))]]
        nodes[t]={'id':t,'parent_id':pid,'children':[],'descendant_taxa':[t],'time_CU':0.}
    for n in nodes.values():
        if n['parent_id'] is not None: nodes[n['parent_id']]['children'].append(n['id'])
    def descend(n):
        x=nodes[n]
        if x['children']:x['descendant_taxa']=sorted(itertools.chain.from_iterable(descend(c) for c in sorted(x['children'])))
        return x['descendant_taxa']
    descend('pop_ROOT')
    for n in nodes.values():
        n['children'].sort()
        if n['parent_id'] is not None and nodes[n['parent_id']]['time_CU']<=n['time_CU']:raise ValueError('Nonpositive edge after zero-edge collapse')
    return list(nodes.values()),{'V_original_node':V,'H_CU':H,'dmax_CU':dmax,'collapsed_edges':sum(v==0 for v in values),
      'original_internal_node_to_population':{str(n):ids[r] for n,r in representative.items()}}


def mrca_age(nodes,a,b):
    d={x['id']:x for x in nodes}; ancestors=set(); n=a
    while n is not None:ancestors.add(n);n=d[n]['parent_id']
    n=b
    while n not in ancestors:n=d[n]['parent_id']
    return d[n]['time_CU']


def pair_data(species,edges,values,trees,genes):
    import numpy as np
    # Unit offsets here are only algebraic templates; no fitted value is claimed.
    nodes,_=rooted_structure(species,edges,values,1.,1.)
    pairs=[];invalid=[]
    for g,t in zip(genes,trees):
        rows=[]
        for a,b in itertools.combinations(sorted(t.taxa),2):
            D=t.distance(a,b)
            if not math.isfinite(D) or D<=0:invalid.append({'gene':g,'a':a,'b':b,'D':D});continue
            root_pair='Arabidopsis_thaliana' in (a,b)
            base=mrca_age(nodes,a,b)-1.-int(root_pair)
            rows.append((a,b,math.log(D),base,int(root_pair)))
        if not rows:raise ValueError('No positive pair distance: cannot assign source-calibrated rate '+g)
        pairs.append({'gene':g,'labels':[(r[0],r[1]) for r in rows],
          'log_D':np.array([r[2] for r in rows]),'base':np.array([r[3] for r in rows]),
          'root_pair':np.array([r[4] for r in rows])})
    return pairs,invalid


def clock_profile(delta,rho,pairs,details=False):
    import numpy as np
    f=[];resids=[];loss=[]
    for p in pairs:
        tau=p['base']+delta+p['root_pair']*rho
        y=p['log_D']-np.log(2*(tau+1))
        fj=float(y.mean());r=y-fj;f.append(fj);resids.append(r);loss.append(float(np.mean(r*r)))
    lk=float(np.mean(f));value=float(np.mean(loss))
    if not details:return value
    return {'objective':value,'kappa':math.exp(lk),'log_kappa':lk,
      'gene_relative_rates':{p['gene']:math.exp(v-lk) for p,v in zip(pairs,f)},
      'gene_log_scale':{p['gene']:v for p,v in zip(pairs,f)},'residuals':resids,'gene_MSE':loss}


def optimize_clock(pairs,values):
    import numpy as np
    from scipy.optimize import minimize
    positive=[x for x in values if x>0 and x<10-NUMERIC['boundary_flag_tolerance']]
    s=float(np.median(positive)) if positive else 1.
    lo,hi=s/16,s*16; bounds=[(math.log(lo),math.log(hi))]*2; records=[]
    for a,b in itertools.product([.25,1.,4.],repeat=2):
        r=minimize(lambda x:clock_profile(*np.exp(x),pairs),np.log([s*a,s*b]),method='L-BFGS-B',bounds=bounds,
          options={'maxiter':2000,'ftol':1e-12,'gtol':1e-8})
        records.append({'start_multipliers':[a,b],'success':bool(r.success),'message':str(r.message),
          'objective':float(r.fun),'delta':float(math.exp(r.x[0])),'rho':float(math.exp(r.x[1])),'nit':int(r.nit)})
    good=[r for r in records if r['success'] and math.isfinite(r['objective'])]
    if not good:raise ValueError('No clock optimizer start converged')
    minobj=min(r['objective'] for r in good); equal=[r for r in good if r['objective']<=minobj+1e-8]
    best=min(equal,key=lambda r:(r['delta'],r['rho']))
    grid=np.geomspace(lo,hi,PROFILE_POINTS); profile=[]
    for d,r in itertools.product(grid,repeat=2):profile.append({'delta':float(d),'rho':float(r),'objective':clock_profile(d,r,pairs)})
    near=[r for r in profile+good if r['objective']<=minobj+1e-6]
    spans={k:max(math.log(r[k]) for r in near)-min(math.log(r[k]) for r in near) for k in ['delta','rho']}
    flags={'boundary':{k:abs(math.log(best[k]/lo))<1e-7 or abs(math.log(best[k]/hi))<1e-7 for k in ['delta','rho']},
      'near_optimum_log_spans':spans,'nonidentified_profile_diagnostic':any(x>=math.log(2) for x in spans.values()),
      'nonunique_start_optima':any(abs(math.log(r[k]/best[k]))>1e-4 for r in equal for k in ['delta','rho']),
      'profile_flags_are_not_confidence_intervals':True}
    return best,{'anchor_s':s,'anchor_is_fallback_unit_scenario':not bool(positive),'bounds':[lo,hi],
      'starts':records,'selection':'minimum objective within 1e-8, then delta/rho ascending','flags':flags},profile


def guard(receipt):
    if not receipt or sys.platform!='linux':raise ValueError('Actual fitting requires Linux --scope-receipt')
    r=json.loads(Path(receipt).read_text()); rel=next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    cg=Path('/sys/fs/cgroup')/rel.lstrip('/'); mem=(cg/'memory.max').read_text().strip();swap=(cg/'memory.swap.max').read_text().strip()
    cpus=sorted(os.sched_getaffinity(0))
    if mem!='4294967296' or swap!='0' or len(cpus)!=1:raise ValueError('Actual4GiB/no-swap/oneCPU guard failed')
    if r.get('status')!='VERIFIED_BEFORE_EXEC' or r.get('cgroup')!=str(cg) or r.get('memory_max_bytes')!=4294967296 or r.get('swap_max_bytes')!=0 or r.get('logical_cpu_ids')!=cpus:raise ValueError('Guard receipt differs from current process')
    return {'receipt_path':str(receipt),'receipt_sha256':sha(receipt),'actual_cgroup':str(cg),'memory_bytes':int(mem),'swap_bytes':0,'cpus':cpus,
      'does_not_establish_human_authorization':True}


def write_tsv(path,head,rows):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'xt') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(head);w.writerows(rows)


def execute(data,out):
    import numpy as np
    names=[g['gene'] for g in data['gene_entries']]
    obs=quartet_observations(data['species'],data['gene_trees'],names)
    save(out/'QUARTET_GENE_AUDIT.json',obs['gene_audit'])
    save(out/'QUARTET_COUNTS.json',{'quartets':obs['quartets'],'split_choices':'pair containing lexicographic-first taxon; pair with next3 taxa in order',
      'concord_indices':obs['concord_indices'],'B_internal_path':obs['B'].astype(int).tolist(),
      'raw_integer_counts':obs['raw'].tolist(),'equal_gene_weighted_counts':obs['weighted'].tolist()})
    write_tsv(out/'PER_GENE_QUARTETS.tsv.gz',['gene','quartet_index_0based','status','observed_split_index','weight'],obs['rows'])
    values,starts=optimize_cu(obs)
    pairs,invalid=pair_data(data['species'],obs['edges'],values,data['gene_trees'],names)
    best,clockmeta,profile=optimize_clock(pairs,values)
    detail=clock_profile(best['delta'],best['rho'],pairs,True)
    nodes,orient=rooted_structure(data['species'],obs['edges'],values,best['delta'],best['rho'])
    save(out/'CLOCK_PROFILE.json',{'grid':'41x41 logarithmically equally spaced including bounds','points':profile})
    rows=[];taxon_residuals={t:[] for t in data['taxa']}
    for p,resid,mse in zip(pairs,detail['residuals'],detail['gene_MSE']):
        for (a,b),r in zip(p['labels'],resid):rows.append((p['gene'],a,b,float(r),1/len(resid),mse))
        for t in data['taxa']:
            rr=[float(r) for (a,b),r in zip(p['labels'],resid) if t in (a,b)]
            if rr:taxon_residuals[t].append(float(np.mean(rr)))
    write_tsv(out/'CLOCK_PAIR_RESIDUALS.tsv',['gene','taxon_a','taxon_b','log_distance_residual','within_gene_pair_weight','gene_MSE'],rows)
    save(out/'CLOCK_MISFIT.json',{'equal_gene_mean_MSE':detail['objective'],'invalid_pairs':invalid,
      'by_taxon_equal_gene_mean_log_residual':{t:float(np.mean(r)) if r else None for t,r in taxon_residuals.items()},
      'interpretation':'conditional strict-clock mismatch; not rejection of ILS'})
    eout=[{'edge_id':'CU_%02d'%i,'split_side':list(e[0]),'CU':v,'exact_zero':v==0,
      'near_lower_boundary':v<=1e-7,'right_censored':v>=10-1e-7} for i,(e,v) in enumerate(zip(obs['edges'],values))]
    result={'schema_version':SCHEMA,'scenario':'S0','status':'FITTED_CONDITIONAL_SCENARIO',
      'taxa':data['taxa'],'rooted_population_nodes':nodes,'kappa':detail['kappa'],
      'kappa_units':'AA substitutions per site per coalescent-time unit',
      'gene_relative_rates':detail['gene_relative_rates'],'rate_normalization':'geometric mean across original173 genes equals1',
      'gene_entries':data['gene_entries'],'CU_edges':eout,'CU_fit_starts':starts,
      'delta':best['delta'],'rho':best['rho'],'root_time_fit':clockmeta,'root_orientation':orient,
      'population_unit':{'mean_two_lineage_coalescence_time':1,'global_ploidy':2,'initial_size':.5,'samples_per_taxon':1,'sample_ploidy':1},
      'clock_profile':{'path':str(out/'CLOCK_PROFILE.json'),'sha256':sha(out/'CLOCK_PROFILE.json')},
      'clock_misfit':{'path':str(out/'CLOCK_MISFIT.json'),'sha256':sha(out/'CLOCK_MISFIT.json')},
      'simulation_admissibility':{'structural_parameter_checks_passed':True,'conditional_S0_only':True,
        'actual_execution_authorization_established':False,'engineering_gate_acceptance_established':False,
        'required_cap20_sensitivity':any(e['right_censored'] for e in eout),
        'root_time_boundary_or_nonidentification':clockmeta['flags']},
      'uncalibrated_limitations':['estimated gene-tree errors uncorrected','conditional fixed topology/root/constantNe/strict within-gene clock',
        'overlapping quartets and unknown linkage; no standard errors or P values','no physical Ne, generations or mutation rate estimated',
        'no indirect-inference correction or parameter-bootstrap calibration','S0 output cannot exclude general ILS or establish gene flow']}
    if abs(sum(math.log(v) for v in result['gene_relative_rates'].values())/len(names))>1e-10:raise ValueError('Rate normalization failed')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--contract',required=True);p.add_argument('--out',required=True)
    p.add_argument('--path-map',help='JSON object of exact original absolute path to absolute actual path; all SHA remain mandatory')
    p.add_argument('--scope-receipt');p.add_argument('--execute',action='store_true');a=p.parse_args()
    # Guard before data loading or optimizer dependency import. Setting thread env cannot replace actual affinity.
    resource_guard=guard(a.scope_receipt) if a.execute else None
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
    out=Path(a.out).absolute();out.mkdir(parents=True,exist_ok=False);started=time.time()
    try:
        if sha(a.contract)!=CONTRACT_SHA:raise ValueError('Frozen contract SHA mismatch')
        contract=json.loads(Path(a.contract).read_text());path_map=json.loads(Path(a.path_map).read_text()) if a.path_map else {}
        if not isinstance(path_map,dict) or any(not isinstance(k,str) or not isinstance(v,str) or not Path(k).is_absolute() or not Path(v).is_absolute() for k,v in path_map.items()):raise ValueError('Invalid exact path map')
        data=bind_sources(contract,path_map)
        save(out/'SOURCE_BINDINGS.json',data['source_bindings']);save(out/'NUMERICAL_RULES.json',NUMERIC)
        with open(out/'FIXED_SPECIES_TOPOLOGY.nwk','x') as f:f.write(data['species_newick']+'\n')
        base={'schema_version':SCHEMA,'scenario':'S0','status':'CHECKED_NOT_FITTED','taxa':data['taxa'],
          'source_manifest_hashes':dict(data['source_manifest_hashes'],contract=CONTRACT_SHA,source_bindings=sha(out/'SOURCE_BINDINGS.json')),
          'contract_path':str(Path(a.contract).absolute()),'path_map_sha256':sha(a.path_map) if a.path_map else None,
          'implementation_sha256':sha(__file__),'numerical_rules_sha256':sha(out/'NUMERICAL_RULES.json'),
          'gene_entries':data['gene_entries'],'rooted_population_nodes':None,'kappa':None,'gene_relative_rates':None,
          'simulation_admissibility':{'structural_parameter_checks_passed':False,'reason':'not fitted; no simulation authority'}}
        if a.execute:base.update(execute(data,out))
        save(out/('S0_FITTED.json' if a.execute else 'CHECKED_PLAN.json'),base)
        import resource
        save(out/'RECEIPT.json',{'status':'COMPLETE','mode':'fit' if a.execute else 'source_checks_only',
          'elapsed_seconds':time.time()-started,'maxrss_KiB_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          'guard':resource_guard,'argv':sys.argv,'actual_simulation_performed':False})
    except Exception as e:
        save(out/'FAILURE.json',{'status':'STOPPED','error':repr(e),'elapsed_seconds':time.time()-started,'no_automatic_retry':True})
        raise

if __name__=='__main__':main()
