"""Synthetic-only formula/invariance/kernel checks; no private input or scientific validation."""
from pathlib import Path
import hashlib,json,time,sys,platform
import numpy as np
from pcangsd import shared_cy,covariance_cy
import explicit3 as r

HERE=Path(__file__).resolve().parent
RESULT={'scope':'Synthetic mathematical reference checks only; no DT data loaded','checks':{},'differences':[]}
def close(name,a,b,atol=1e-11,rtol=1e-11):
    a,b=np.asarray(a),np.asarray(b)
    delta=float(np.max(np.abs(a-b)))
    np.testing.assert_allclose(a,b,atol=atol,rtol=rtol)
    RESULT['checks'][name]={'status':'PASS','max_abs_difference':delta,'atol':atol,'rtol':rtol}
def reject(name,call):
    try:call()
    except ValueError as e:RESULT['checks'][name]={'status':'PASS_REJECTED','message':str(e)}
    else:raise AssertionError('Expected rejection: '+name)

def main():
    start=time.monotonic()
    rng=np.random.default_rng(20260924)
    # Flat means no genotype evidence: posterior equals prior, not confident genotype.
    f=np.array([.15,.25,.4,.65,.85]);flat=np.ones((5,4,3))
    q=r.posterior(flat,f[:,None]);expected=np.stack(((1-f)**2,2*f*(1-f),f**2),axis=1)[:,None,:]
    close('flat_posterior_is_HWE_prior',q,np.broadcast_to(expected,q.shape))
    close('flat_dosage_is_2f',r.dosage(q),np.broadcast_to(2*f[:,None],(5,4)))
    close('flat_normal_centered_zero',r.update_normal(flat,f),np.zeros((5,4)))
    close('flat_normal_covariance_identity',r.cov_normal(flat,f)['covariance'],np.eye(4))
    maf=r.em_maf(flat);close('flat_unidentifiable_f_stays_initial',maf['f'],np.full(5,.25))
    assert maf['all_cells_flat_sites'].all()
    flat_fit=r.fit(flat,rank=1,f=np.full(5,.25))
    assert flat_fit['all_sites_flat'] and flat_fit['all_cells_flat_site_mask'].all()
    close('flat_external_f_mask_and_covariance',flat_fit['covariance'],np.eye(4))
    assert len(flat_fit['all_singular_values'])==4
    RESULT['checks']['flat_fit_diagnostics_not_identifiability']={'status':'PASS','all_sites_flat':True,'pca_converged':flat_fit['pca_converged'],'all_singular_values':flat_fit['all_singular_values'].tolist(),'rank_cut_gap':flat_fit['rank_cut_singular_value_gap'],'interpretation':'Algorithmic convergence of all-flat data cannot identify population structure'}
    _,_,tie_s=r.truncated_svd(np.diag([3.,1.,1.,0.]),rank=2)
    close('rank_cut_exact_tie_full_spectrum',tie_s,np.array([3.,1.,1.,0.]))
    pi=np.broadcast_to(np.array([.01,.2,.6,.99])[None,:],(5,4));cc=r.covariance_from_pi(flat,pi,f)
    diagonal=np.mean((2*pi*(1-pi)+(2*pi-2*f[:,None])**2)/(2*f*(1-f))[:,None],axis=0)
    close('flat_structured_diagonal_keeps_posterior_variance',np.diag(cc['covariance']),diagonal)
    # One-hot GT gives deterministic dosage independent of an interior prior.
    gt=np.array([[0,0,1,2],[0,1,2,2],[2,0,1,1],[1,1,0,2],[2,2,0,0]])
    hard=np.eye(3)[gt];hf=gt.mean(axis=1)/2
    close('certain_GT_dosage',r.dosage(r.posterior(hard,pi)),gt)
    close('certain_GT_frequency_EM',r.em_maf(hard)['f'],hf)
    z=(gt-2*hf[:,None])/np.sqrt(2*hf*(1-hf))[:,None]
    close('certain_GT_covariance_equals_direct_genotype_GRM',r.covariance_from_pi(hard,pi,hf)['covariance'],z.T@z/len(gt))
    hfit=r.fit(hard,rank=2);close('certain_GT_full_fit_equals_direct_GRM',hfit['covariance'],z.T@z/len(gt))
    assert hfit['pca_converged']
    # Analytic extreme counterexample: third value is never recovered by subtraction.
    pl=np.array([0.,80.,300.]);gl=r.likelihoods((10**(-pl/10))[None,None,:]);f0=np.array([.5]);recon=np.array([[.9998]])
    exact=r.cov_pcangsd(gl,recon,f0);q0=exact['posterior'][0,0]
    assert (q0>=0).all() and abs(q0.sum()-1)<1e-15 and 0<=r.dosage(exact['posterior'])[0,0]<=2
    l2=np.ascontiguousarray(gl[...,:2].reshape(1,2).astype(np.float32));reconstructed=1-float(l2[0,0])-float(l2[0,1])
    p=np.ascontiguousarray(recon,dtype=np.float32);e=np.zeros((1,1),np.float32);diag=np.zeros(1)
    covariance_cy.covPCAngsd(l2,p,e,f0,diag)
    RESULT['checks']['PL_0_80_300_high_pi_counterexample']={'status':'PASS_SYNTHETIC_COUNTEREXAMPLE','pi_reference':float(exact['pi'][0,0]),'explicit_GL':gl[0,0].tolist(),'official_stored_two_float32':l2[0].tolist(),'official_residual_L2':reconstructed,'reference_posterior':q0.tolist(),'reference_dosage':float(r.dosage(exact['posterior'])[0,0]),'reference_diagonal_sum':exact['diagonal_sum'].tolist(),'official_standardized_dosage':e.tolist(),'official_diagonal_sum':diag.tolist(),'interpretation':'Synthetic arithmetic failure, not measured frequency or impact in DT data'}
    assert reconstructed<0 and ((diag<0).any() or (np.abs(e)>3).any() or not np.isfinite(e).all())
    # Moderate exact dyadic GL avoids conflating input quantization with formula mismatch.
    patterns=np.array([[1,2,5],[2,4,2],[5,2,1],[3,1,4],[1,6,1]],dtype=np.float64)/8
    moderate=patterns[rng.integers(0,len(patterns),(64,8))]
    f=rng.uniform(.1,.9,64);recon=np.ascontiguousarray(rng.uniform(-.2,.2,(64,8)),dtype=np.float32)
    legacy=np.ascontiguousarray(moderate[...,:2].reshape(64,16),dtype=np.float32)
    legacy_f=f.copy();shared_cy.emMAF_update(legacy,legacy_f)
    close('compiled_emMAF_update',r.em_update(moderate,f),legacy_f)
    e=np.zeros((64,8),np.float32);covariance_cy.updateNormal(legacy,e,f)
    close('compiled_updateNormal_each_dosage',r.update_normal(moderate,f),e,atol=2e-7,rtol=2e-7)
    p=recon.copy();e.fill(0);covariance_cy.updatePCAngsd(legacy,p,e,f)
    ee,pp=r.update_pcangsd(moderate,recon,f)
    close('compiled_updatePCAngsd_each_dosage',ee,e,atol=2e-7,rtol=2e-7)
    close('compiled_updatePCAngsd_each_pi',pp,p,atol=5e-8,rtol=5e-8)
    post=np.zeros((64,24),np.float32);shared_cy.post(legacy,p,post)
    close('compiled_post_each_posterior_moderate',r.posterior(moderate,p),post.reshape(64,8,3),atol=1e-7,rtol=1e-7)
    for name in ['normal','pcangsd']:
        e.fill(0);diag=np.zeros(8)
        if name=='normal':covariance_cy.covNormal(legacy,e,f,diag);expected=r.cov_normal(moderate,f)
        else:p=recon.copy();covariance_cy.covPCAngsd(legacy,p,e,f,diag);expected=r.cov_pcangsd(moderate,recon,f)
        close('compiled_cov_'+name+'_each_standardized_dosage',expected['standardized_dosage'],e,atol=3e-7,rtol=3e-7)
        close('compiled_cov_'+name+'_each_diagonal_sum',expected['diagonal_sum'],diag,atol=1e-10,rtol=1e-11)
        c=e.T@e;np.fill_diagonal(c,diag);c/=64
        close('compiled_cov_'+name+'_full_covariance',expected['covariance'],c,atol=2e-7,rtol=2e-7)
    # Asymmetric synthetic data, fixed iterations: invariance of C, not eigenvector signs.
    genotypes=rng.integers(0,3,(64,8));genotypes[:,0]=0;genotypes[:,1]=2;genotypes[:,2]=1
    moderate=r.likelihoods(.95*np.eye(3)[genotypes]+.05*rng.uniform(.05,1,(64,8,3)))
    moderate[::3,7,:]=1/3
    f=r.em_maf(moderate,max_iter=10000)['f'];ff=r.em_maf(moderate[...,::-1],initial=.75,max_iter=10000)['f']
    assert np.all((f>.05)&(f<.95)),f
    close('global_allele_flip_frequency',ff,1-f,atol=1e-10,rtol=1e-10)
    baseline=r.fit(moderate,rank=2,f=f,pca_max_iter=15,pca_tolerance=1e-30)
    flipped=r.fit(moderate[...,::-1],rank=2,f=1-f,pca_max_iter=15,pca_tolerance=1e-30)
    close('global_allele_flip_full_iterated_covariance',flipped['covariance'],baseline['covariance'],atol=2e-10,rtol=2e-10)
    close('global_allele_flip_individual_pi',flipped['pi'],1-baseline['pi'],atol=2e-10,rtol=2e-10)
    order=rng.permutation(8);permuted=r.fit(moderate[:,order,:],rank=2,f=f,pca_max_iter=15,pca_tolerance=1e-30)
    close('individual_permutation_full_iterated_covariance',permuted['covariance'],baseline['covariance'][np.ix_(order,order)],atol=2e-10,rtol=2e-10)
    close('individual_permutation_pi',permuted['pi'],baseline['pi'][:,order],atol=2e-10,rtol=2e-10)
    close('likelihood_cell_positive_scaling_invariance',r.posterior(moderate*rng.uniform(.1,10,(64,8,1)),baseline['pi']),r.posterior(moderate,baseline['pi']))
    assert np.linalg.eigvalsh(baseline['covariance']).min()>-1e-10
    RESULT['checks']['fixed_iteration_status_preserved']={'status':'PASS','pca_iterations':baseline['pca_iterations'],'pca_converged':baseline['pca_converged'],'purpose':'15 finite iterations exercise equivariance, not a convergence claim'}
    reject('reject_negative_GL',lambda:r.likelihoods(np.array([[[-1,1,1]]])))
    reject('reject_allzero_GL',lambda:r.likelihoods(np.zeros((1,1,3))))
    reject('reject_nan_GL',lambda:r.likelihoods(np.array([[[np.nan,1,1]]])))
    reject('reject_boundary_frequency_covariance',lambda:r.cov_normal(np.ones((1,2,3)),np.array([0.])))
    reject('reject_invalid_rank',lambda:r.truncated_svd(np.ones((4,3)),3))
    RESULT.update(status='PASS_SYNTHETIC_MATH_CHECKS_ONLY',elapsed_seconds=time.monotonic()-start,python=sys.version,executable=sys.executable,numpy=np.__version__,platform=platform.platform(),fixture_seed=20260924,thread_environment={'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},private_data_used=False,scientific_validation=False)
    RESULT['compiled_modules']=[{'path':m.__file__,'sha256':hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()} for m in [shared_cy,covariance_cy]]
    RESULT['reference_sha256']=hashlib.sha256((HERE/'explicit3.py').read_bytes()).hexdigest();RESULT['test_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with (HERE/'TEST_RESULTS_v2.json').open('x') as h:json.dump(RESULT,h,indent=2,allow_nan=False);h.write('\n')
    print(json.dumps(RESULT,indent=2,allow_nan=False))
if __name__=='__main__':main()
