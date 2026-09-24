"""Diploid IBS and unrooted topology helpers; missing genotypes stay missing."""
from pathlib import Path
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def read_bed(path,n_samples,n_snps):
    raw=Path(path).read_bytes(); width=(n_samples+3)//4
    if raw[:3]!=b'\x6c\x1b\x01' or len(raw)!=3+width*n_snps:
        raise ValueError('Invalid SNP-major BED header or size')
    a=np.frombuffer(raw[3:],dtype=np.uint8).reshape(n_snps,width)
    codes=((a[:,:,None] >> (2*np.arange(4,dtype=np.uint8)))&3).reshape(n_snps,-1)[:,:n_samples]
    return np.array([2,-1,1,0],dtype=np.int8)[codes]

def ibs_distance(g):
    g=np.asarray(g)
    if g.ndim!=2 or not np.isin(g,[-1,0,1,2]).all():raise ValueError('Invalid diploid dosages')
    n=g.shape[1];d=np.zeros((n,n));counts=np.zeros((n,n),dtype=np.int64)
    for i in range(n):
        for j in range(i,n):
            ok=(g[:,i]>=0)&(g[:,j]>=0);k=int(ok.sum())
            if k==0:raise ValueError('No jointly called sites')
            d[i,j]=d[j,i]=np.abs(g[ok,i]-g[ok,j]).sum()/(2*k)
            counts[i,j]=counts[j,i]=k
    return d,counts

def nj_tree(names,d):
    d=np.asarray(d)
    if len(set(names))!=len(names) or d.shape!=(len(names),len(names)) or not np.isfinite(d).all() or not np.allclose(d,d.T) or not np.allclose(np.diag(d),0):raise ValueError('Invalid distance matrix')
    t=DistanceTreeConstructor().nj(DistanceMatrix(list(names),[list(d[i,:i+1]) for i in range(len(names))]));t.rooted=False
    return t

def canonical_split(part,all_names):
    a=tuple(sorted(part));b=tuple(sorted(set(all_names)-set(part)))
    if min(len(a),len(b))<2:return None
    return min((a,b),key=lambda x:(len(x),x))

def tree_splits(tree):
    names={t.name for t in tree.get_terminals()}
    return {s for c in tree.find_clades() if (s:=canonical_split({t.name for t in c.get_terminals()},names)) is not None}

def residual_fit(d,estimated):
    d=np.asarray(d);estimated=np.asarray(estimated)
    if d.shape!=estimated.shape:raise ValueError('Different matrix sizes')
    sel=np.triu_indices(len(d),1);den=np.sum(d[sel]**2)
    if den==0:raise ValueError('Zero distance norm')
    return float(100*(1-np.sum((d[sel]-estimated[sel])**2)/den))

def checked_neighbor_net(names,d,cutoff=1e-10):
    """Original Neighbor-Net circular ordering + explicit SciPy NNLS split weights.

    Avoids the upstream 0.0.10 conjugate solver, which failed repeatability checks.
    Enumerates the n(n-1)/2 circular splits and checks NNLS KKT optimality.
    """
    from scipy.optimize import nnls
    from splitspy.nnet import nnet_cycle
    from splitspy.splits.basic_split import cyc_split
    d=np.asarray(d,dtype=float);n=len(names)
    cycle=nnet_cycle.compute(names,d.tolist())
    candidates=[cyc_split(cycle,i+1,j,0.) for i in range(1,n+1) for j in range(i+1,n+1)]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    design=np.array([[s.separates(i+1,j+1) for s in candidates] for i,j in pairs],dtype=float)
    target=np.array([d[i,j] for i,j in pairs]);weights,rnorm=nnls(design,target,maxiter=10000)
    grad=design.T@(design@weights-target);active=weights>cutoff
    violation=float(max(0.,-grad.min(),np.max(np.abs(grad[active])) if active.any() else 0.))
    if violation>1e-9:raise ValueError('NNLS optimality check failed')
    splits=[]
    for s,w in zip(candidates,weights):
        if w>cutoff:s.weight=float(w);splits.append(s)
    return cycle,splits,{'solver':'scipy.optimize.nnls','candidates':len(candidates),'positive_splits':len(splits),'residual_norm':float(rnorm),'kkt_violation':violation}
