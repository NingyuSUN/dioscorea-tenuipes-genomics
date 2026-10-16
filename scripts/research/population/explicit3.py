"""Auditable explicit-three-genotype float64 mathematical reference, not production PCA.

Equations follow pinned PCAngsd kernels and Meisner & Albrechtsen (2018).
No source files, private data, CLI, automatic site filters or MAP rank selection.
"""
import numpy as np

PI_MIN = 1e-4


def likelihoods(values):
    """Validate M x N x 3 nonnegative finite GL; normalize all three explicitly."""
    x = np.array(values, dtype=np.float64, copy=True)
    if x.ndim != 3 or x.shape[-1] != 3 or min(x.shape[:2]) < 1:
        raise ValueError('GL must have nonempty shape M x N x 3')
    if not np.isfinite(x).all() or (x < 0).any():
        raise ValueError('GL must be finite and nonnegative')
    scale = x.max(axis=2, keepdims=True)
    if (scale == 0).any():
        raise ValueError('All-zero GL is invalid; represent missing with a flat triple')
    x /= scale
    x /= x.sum(axis=2, keepdims=True)
    return x


def _frequencies(f, m, interior=False):
    f = np.asarray(f, dtype=np.float64)
    if f.shape != (m,) or not np.isfinite(f).all() or (f < 0).any() or (f > 1).any():
        raise ValueError('f must be a finite length-M frequency vector')
    if interior and ((f <= 0).any() or (f >= 1).any()):
        raise ValueError('Covariance requires 0 < f < 1; no automatic site deletion')
    return f


def posterior(gl, pi):
    """q_g proportional to L_g h_g(pi); no residual L2 or posterior flooring."""
    gl = likelihoods(gl)
    pi = np.asarray(pi, dtype=np.float64)
    try:
        pi = np.broadcast_to(pi, gl.shape[:2])
    except ValueError as e:
        raise ValueError('pi must broadcast to M x N') from e
    if not np.isfinite(pi).all() or (pi < 0).any() or (pi > 1).any():
        raise ValueError('pi must lie in [0,1]')
    prior = np.stack(((1-pi)**2, 2*pi*(1-pi), pi**2), axis=2)
    weighted = gl * prior
    denominator = weighted.sum(axis=2, keepdims=True)
    if (denominator <= 0).any() or not np.isfinite(denominator).all():
        raise ValueError('Zero/nonfinite posterior denominator; no fabricated posterior')
    return weighted / denominator


def dosage(q):
    return q[..., 1] + 2*q[..., 2]


def em_update(gl, f):
    gl = likelihoods(gl)
    f = _frequencies(f, gl.shape[0])
    return dosage(posterior(gl, f[:, None])).mean(axis=1) * 0.5


def em_maf(gl, initial=0.25, max_iter=10000, tolerance=1e-10):
    """Ordinary unaccelerated EM; returned f is encoded-allele frequency, not folded MAF."""
    gl = likelihoods(gl)
    if max_iter < 1 or tolerance <= 0:
        raise ValueError('Positive max_iter and tolerance required')
    f = np.full(gl.shape[0], initial, dtype=np.float64) if np.ndim(initial) == 0 else np.array(initial, dtype=np.float64, copy=True)
    f = _frequencies(f, gl.shape[0])
    trace = []
    for iteration in range(1, max_iter+1):
        updated = em_update(gl, f)
        delta = float(np.max(np.abs(updated-f)))
        trace.append(delta)
        f = updated
        if delta < tolerance:
            break
    return {'f': f, 'converged': delta < tolerance, 'iterations': iteration,
            'max_abs_update_trace': trace, 'initial': np.asarray(initial).tolist(),
            'all_cells_flat_sites': np.all(gl == gl[..., :1], axis=(1, 2))}


def update_normal(gl, f):
    gl = likelihoods(gl)
    f = _frequencies(f, gl.shape[0])
    return dosage(posterior(gl, f[:, None])) - 2*f[:, None]


def reconstruction_to_pi(reconstruction, f):
    r = np.asarray(reconstruction, dtype=np.float64)
    if r.ndim != 2 or not np.isfinite(r).all():
        raise ValueError('Reconstruction must be a finite M x N matrix')
    f = _frequencies(f, r.shape[0])
    return np.clip((r + 2*f[:, None])*0.5, PI_MIN, 1-PI_MIN)


def update_pcangsd(gl, reconstruction, f):
    gl = likelihoods(gl)
    if np.shape(reconstruction) != gl.shape[:2]:
        raise ValueError('Reconstruction shape differs from GL')
    f = _frequencies(f, gl.shape[0])
    pi = reconstruction_to_pi(reconstruction, f)
    q = posterior(gl, pi)
    return dosage(q) - 2*f[:, None], pi


def covariance_from_pi(gl, pi, f):
    """Replace dosage-square diagonal with posterior second centered moment."""
    gl = likelihoods(gl)
    f = _frequencies(f, gl.shape[0], interior=True)
    q = posterior(gl, pi)
    center = 2*f[:, None]
    denominator = 2*f*(1-f)
    z = (dosage(q)-center) / np.sqrt(denominator)[:, None]
    centered_genotypes = np.arange(3, dtype=np.float64)[None, None, :] - center[..., None]
    diagonal_sum = ((q * centered_genotypes**2).sum(axis=2) / denominator[:, None]).sum(axis=0)
    c = z.T @ z
    np.fill_diagonal(c, diagonal_sum)
    c /= gl.shape[0]
    return {'covariance': c, 'standardized_dosage': z, 'diagonal_sum': diagonal_sum, 'posterior': q}


def cov_normal(gl, f):
    gl = likelihoods(gl)
    f = _frequencies(f, gl.shape[0], interior=True)
    return covariance_from_pi(gl, f[:, None], f)


def cov_pcangsd(gl, reconstruction, f):
    gl = likelihoods(gl)
    if np.shape(reconstruction) != gl.shape[:2]:
        raise ValueError('Reconstruction shape differs from GL')
    pi = reconstruction_to_pi(reconstruction, f)
    result = covariance_from_pi(gl, pi, f)
    result['pi'] = pi
    return result


def truncated_svd(centered_dosage, rank):
    e = np.asarray(centered_dosage, dtype=np.float64)
    if e.ndim != 2 or not np.isfinite(e).all() or not isinstance(rank, (int, np.integer)) or not 1 <= rank < min(e.shape):
        raise ValueError('Require finite M x N and integer 1 <= rank < min(M,N)')
    u, all_s, vt = np.linalg.svd(e, full_matrices=False)
    u, s, vt = u[:, :rank], all_s[:rank], vt[:rank, :]
    # Same sign rule as covariance.py:8-13; only ordering differs (descending).
    signs = np.sign(u[np.argmax(np.abs(u), axis=0), np.arange(rank)])
    u = u * signs
    vt = vt * signs[:, None]
    return (u*s) @ vt, vt.T.copy(), all_s


def fit(gl, rank, f=None, em_max_iter=10000, em_tolerance=1e-10,
        pca_max_iter=100, pca_tolerance=1e-5):
    """Fixed-rank iteration order of covariance.py:73-104, in float64.

Returns explicit convergence status. A nonconverged PCA covariance is diagnostic,
never silently accepted. Population f stays fixed during all PCA iterations.
"""
    gl = likelihoods(gl)
    if pca_max_iter < 1 or pca_tolerance <= 0:
        raise ValueError('Positive PCA max_iter and tolerance required')
    frequency_em = None
    if f is None:
        frequency_em = em_maf(gl, max_iter=em_max_iter, tolerance=em_tolerance)
        if not frequency_em['converged']:
            raise ValueError('Population-frequency EM did not converge')
        f = frequency_em['f']
    f = _frequencies(f, gl.shape[0], interior=True)
    reconstruction, v_previous, s = truncated_svd(update_normal(gl, f), rank)
    trace = []
    for iteration in range(1, pca_max_iter+1):
        e, _ = update_pcangsd(gl, reconstruction, f)
        reconstruction, v, s = truncated_svd(e, rank)
        difference = float(np.sqrt(np.mean((v-v_previous)**2)))
        trace.append(difference)
        if difference < pca_tolerance:
            break
        v_previous = v
    result = cov_pcangsd(gl, reconstruction, f)
    result.update(f=f.copy(), reconstruction=reconstruction, rank=int(rank),
                  pca_converged=difference < pca_tolerance, pca_iterations=iteration,
                  v_rmse_trace=trace, singular_values=s[:rank].copy(),
                  all_singular_values=s.copy(),
                  rank_cut_singular_value_gap=float(s[rank-1]-s[rank]),
                  rank_cut_gap_relative_to_last_retained=(float((s[rank-1]-s[rank])/s[rank-1]) if s[rank-1] > 0 else None),
                  all_cells_flat_site_mask=np.all(gl == gl[..., :1], axis=(1, 2)),
                  all_sites_flat=bool(np.all(gl == gl[..., :1])),
                  frequency_em=frequency_em)
    return result
