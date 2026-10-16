"""Descriptive comparison of the frozen full-panel fits; never refits inputs."""
from pathlib import Path
import sys
if not __debug__:
    raise RuntimeError('Optimized Python disables required scientific assertions; run without -O.')
_HELPER_DIR = str(Path(__file__).resolve().parent)
if _HELPER_DIR not in sys.path:
    sys.path.insert(0, _HELPER_DIR)
import csv, hashlib, json
import os
# Bound native numerical libraries before importing them in the CLI process.
for _thread_env in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_thread_env] = '1'
import numpy as np
from scipy.spatial.distance import pdist
from scipy.linalg import subspace_angles
PRIVATE=None
OUT=None
PLINK_ROOT=None
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_table(name, rows):
    with (OUT/name).open('x') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)
def execute():
    OUT.mkdir(parents=True,exist_ok=False)
    r=json.loads((PRIVATE/'RUN_RECEIPT.json').read_text())
    assert r['status']=='COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT'
    assert r['selected_sites']==37566 and len(r['arms'])==12
    for x in r['output_files']: assert sha(PRIVATE/x['path'])==x['sha256']
    a=PLINK_ROOT; ar=json.loads((a/'RECEIPT.json').read_text())
    assert ar['status']=='COMPLETED_PENDING_INDEPENDENT_AUDIT'
    for x in ar['outputs']: assert sha(a/x['path'])==x['sha256']
    names=['A_PLINK',*r['arms']]; matrices={};vectors={};values={};records=[]
    for name in names:
        if name=='A_PLINK':
            c=np.loadtxt(a/'pca.rel'); lines=(a/'pca.eigenvec').read_text().splitlines()
            assert [x.split()[1] for x in lines[1:]]==r['samples']
            v=np.array([[float(x) for x in row.split()[2:]] for row in lines[1:]])
            eig=np.loadtxt(a/'pca.eigenval')
            residual=float(np.max(np.abs(c@v-v*eig)))
            assert residual<1e-4 and np.max(np.abs(v.T@v-np.eye(14)))<1e-4
            # Full-precision reconstruction of PLINK matrix is a separate diagnosis.
        else:
            raw=np.load(PRIVATE/(name+'.npz'))
            c=np.asarray(raw['covariance'],dtype=np.float64)
            eig,v=np.linalg.eigh(c);eig=eig[::-1];v=v[:,::-1]
            assert np.max(np.abs(eig-r['arms'][name]['all_14_eigenvalues']))<1e-10
            residual=float(np.max(np.abs(c@v-v*eig)))
        assert c.shape==(14,14) and np.isfinite(c).all() and eig[1]>0
        matrices[name]=c; vectors[name]=v; values[name]=eig
        records.append({'arm':name,'eigen_equation_max_absolute_residual':residual,
            'trace':float(np.trace(c)),'eigenvalues':eig.tolist(),
            'top2_trace_ratio':float(eig[:2].sum()/eig.sum()),
            'minimum_eigenvalue':float(eig[-1]),
            'trace_ratio_is_not_labeled_explained_variance':True,
            'diagonal_min':float(np.diag(c).min()),'diagonal_max':float(np.diag(c).max()),
            'covariance_asymmetry':float(np.max(np.abs(c-c.T)))})
    pairs=[('A_PLINK','B_HC_e1e6_r2','estimator_bridge')]
    for rank in (2,1,3):
        for rep,kind in [('C_PL','source_PL_vs_HC'),('B_HC_e1e4','epsilon_sensitivity'),('B_official_HC_e1e6','official_HC_bridge')]:
            pairs.append((f'{rep}_r{rank}',f'B_HC_e1e6_r{rank}',kind))
    for rep in ('B_HC_e1e6','C_PL','B_HC_e1e4','B_official_HC_e1e6'):
        for rank in (1,3): pairs.append((f'{rep}_r{rank}',f'{rep}_r2','rank_sensitivity'))
    comparisons=[]; off=~np.eye(14,dtype=bool)
    for name,ref,kind in pairs:
        c=matrices[name];rc=matrices[ref];delta=c-rc
        v=vectors[name][:,:2];rv=vectors[ref][:,:2]
        x=v-v.mean(0);base=rv-rv.mean(0);x/=np.linalg.norm(x);base/=np.linalg.norm(base)
        u,s,vt=np.linalg.svd(x.T@base);aligned=(x@u@vt)*s.sum()
        scores=v*np.sqrt(values[name][:2]);target=rv*np.sqrt(values[ref][:2])
        distances=pdist(scores);refdist=pdist(target)
        comparisons.append({'arm':name,'reference':ref,'comparison':kind,
            'covariance_max_absolute_difference':float(np.max(np.abs(delta))),
            'covariance_relative_frobenius_difference':float(np.linalg.norm(delta)/np.linalg.norm(rc)),
            'diagonal_max_absolute_difference':float(np.max(np.abs(np.diag(delta)))),
            'off_diagonal_max_absolute_difference':float(np.max(np.abs(delta[off]))),
            'top2_centered_subspace_max_angle_degrees':float(np.degrees(subspace_angles(x,base)).max()),
            'top2_centered_normalized_eigenvector_Procrustes_squared_residual':float(np.square(aligned-base).sum()),
            'top2_eigenvector_pairwise_distance_Pearson_r':float(np.corrcoef(pdist(x),pdist(base))[0,1]),
            'top2_sqrt_eigenvalue_score_pairwise_distance_Pearson_r':float(np.corrcoef(distances,refdist)[0,1]),
            'top2_score_distance_norm_ratio':float(np.linalg.norm(distances)/np.linalg.norm(refdist)),
            'top2_score_distance_relative_L2_difference':float(np.linalg.norm(distances-refdist)/np.linalg.norm(refdist)),
            'dependent_pairwise_distances':91,'inferential_test':False})
    coordrows=[];base=vectors['B_HC_e1e6_r2'][:,:2];base=base-base.mean(0);base/=np.linalg.norm(base)
    for name in names:
        v=vectors[name][:,:2];x=v-v.mean(0);x/=np.linalg.norm(x)
        u,s,vt=np.linalg.svd(x.T@base);aligned=(x@u@vt)*s.sum()
        for j,sample in enumerate(r['samples']):
            coordrows.append({'arm':name,'sample':sample,'group':sample.split('-')[0],
                'PC1_eigenvector':float(v[j,0]),'PC2_eigenvector':float(v[j,1]),
                'PC1_Procrustes_aligned':float(aligned[j,0]),'PC2_Procrustes_aligned':float(aligned[j,1])})
    result={'status':'DESCRIPTIVE_FULL_PANEL_COMPARISON','sites':37566,'samples':14,
        'called_cells':r['called_cells'],'missing_cells':r['missing_cells'],
        'reference_arm':'B_HC_e1e6_r2','arms':records,'comparisons':comparisons,
        'source_receipt_sha256':sha(PRIVATE/'RUN_RECEIPT.json'),
        'PLINK_receipt_sha256':sha(a/'RECEIPT.json'),'script_sha256':sha(Path(__file__)),
        'limits':['Same original discovered/filtered/pruned 37566 sites and 14 individuals; not unbiased all-genome sampling.',
        'Ninety-one distances share fourteen samples; no independent significance test or accuracy claim.',
        'Correlation assesses linear agreement, not absolute scale; separate score-distance norms are reported.',
        'Procrustes handles translation, scale and rotation; the plot does not preserve absolute covariance scale.',
        'A versus B changes missing-data model and covariance estimator; B versus C includes source-likelihood-specific frequency estimation.',
        'Rank is a prespecified latent dimension, not ancestral population count or an independently optimal choice.',
        'Agreement cannot establish genotype truth, diploidy, absence of reference bias, or a demographic mechanism.']}
    (OUT/'COMPARISON.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_table('COMPARISONS.tsv',comparisons);write_table('COORDINATES.tsv',coordrows)
    write_table('EIGENVALUES.tsv',[{'arm':name,'component':i+1,'eigenvalue':float(v)} for name in names for i,v in enumerate(values[name])])
    print(json.dumps({'status':result['status'],'comparisons':comparisons},indent=2))

def configure(config, base):
    from _population_cli import path_value, require_directory, fresh_output
    global PRIVATE, OUT, PLINK_ROOT
    PRIVATE = require_directory(path_value(config, 'run_root', base))
    PLINK_ROOT = require_directory(path_value(config, 'plink_root', base))
    OUT = fresh_output(path_value(config, 'outdir', base), [PRIVATE, PLINK_ROOT])


def preflight():
    from _population_cli import report_plan
    receipt = json.loads((PRIVATE / 'RUN_RECEIPT.json').read_text())
    assert receipt['status'] == 'COMPLETED_PENDING_ROOT_SCIENTIFIC_AUDIT'
    assert receipt['selected_sites'] == 37566 and len(receipt['arms']) == 12
    required_arms = {rep + '_r' + str(rank) for rep in ('B_HC_e1e6', 'C_PL', 'B_HC_e1e4', 'B_official_HC_e1e6') for rank in (2, 1, 3)}
    if set(receipt['arms']) != required_arms:
        raise ValueError('The twelve fixed representation/rank arms are required')
    bound_outputs = {item['path'] for item in receipt['output_files']}
    if not {name + '.npz' for name in required_arms} <= bound_outputs:
        raise ValueError('Covariance arrays are missing from the fit receipt bindings')
    for item in receipt['output_files']:
        assert sha(PRIVATE / item['path']) == item['sha256']
    plink = json.loads((PLINK_ROOT / 'RECEIPT.json').read_text())
    assert plink['status'] == 'COMPLETED_PENDING_INDEPENDENT_AUDIT'
    if not {'pca.rel', 'pca.eigenvec', 'pca.eigenval'} <= {item['path'] for item in plink['outputs']}:
        raise ValueError('PLINK numerical outputs lack receipt bindings')
    for item in plink['outputs']:
        assert sha(PLINK_ROOT / item['path']) == item['sha256']
    inputs = {'fit_receipt': PRIVATE / 'RUN_RECEIPT.json', 'plink_receipt': PLINK_ROOT / 'RECEIPT.json', **{'fit:' + item['path']: PRIVATE / item['path'] for item in receipt['output_files']}, **{'plink:' + item['path']: PLINK_ROOT / item['path'] for item in plink['outputs']}}
    return report_plan('compare_full_panel.py', inputs,
                       {name: OUT / name for name in ['COMPARISON.json', 'COMPARISONS.tsv', 'COORDINATES.tsv', 'EIGENVALUES.tsv']}, ['numpy', 'scipy'])


def plan(config_path):
    """Return a read-only, JSON-compatible artifact and dependency plan."""
    config_path=Path(config_path).resolve();config=json.loads(config_path.read_text());configure(config,config_path.parent);return preflight()


def main(argv=None):
    from _population_cli import read_config
    config, base, run = read_config(__doc__, argv)
    configure(config, base)
    print(json.dumps(preflight(), indent=2))
    if run:
        execute()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
