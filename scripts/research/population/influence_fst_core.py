"""Vectorized diploid biallelic WC84 FST matching VCFtools v0.1.16.

Reference implementation: vcftools/vcftools tag v0.1.16,
 src/cpp/variant_file_output.cpp, output_weir_and_cockerham_fst and
 output_windowed_weir_and_cockerham_fst. This implements their global ratio
 of summed variance components, not the arithmetic mean of site FST.

Input genotype dosage is an integer matrix (sites, samples), with -1 missing
and 0/1/2 copies of either consistently encoded allele. Partial calls must be
masked before import. DP/GQ masking is upstream; this module never imputes GT.
Both alleles are evaluated, negative FST is retained, and NaN site FST is
excluded from global numerator/denominator exactly as in VCFtools.
"""
from dataclasses import dataclass
import numpy as np


def validate_genotypes(g):
    if g.ndim != 2 or not np.issubdtype(g.dtype, np.integer):
        raise ValueError('Genotypes must be a 2-D integer matrix')
    if g.size and (g.min() < -1 or g.max() > 2):
        raise ValueError('Only genotype codes -1, 0, 1, 2 are allowed')


@dataclass(frozen=True)
class PopulationCounts:
    called: np.ndarray
    allele_count: np.ndarray
    heterozygotes: np.ndarray


def population_counts(g, columns):
    """Count called diploid individuals, allele copies, and heterozygotes."""
    columns=list(columns)
    if not columns or len(columns)!=len(set(columns)):
        raise ValueError('A population requires unique nonempty sample columns')
    x=g[:,columns]
    called=(x>=0).sum(axis=1, dtype=np.int32)
    allele_count=np.maximum(x,0).sum(axis=1,dtype=np.int32)
    heterozygotes=(x==1).sum(axis=1,dtype=np.int32)
    return PopulationCounts(called,allele_count,heterozygotes)


def subtract_counts(full, g, removed_columns):
    """Reuse baseline sufficient statistics when a few individuals are removed."""
    if not removed_columns:
        return full
    removed=population_counts(g,removed_columns)
    return PopulationCounts(full.called-removed.called,
                            full.allele_count-removed.allele_count,
                            full.heterozygotes-removed.heterozygotes)


def mask_from_counts(called, allele_count, n_samples, max_missing=.1, min_mac=1):
    if n_samples<1 or not 0<=max_missing<=1 or min_mac<0:
        raise ValueError('Invalid site filter parameters')
    required=int(np.ceil(n_samples*(1-max_missing)-1e-12))
    minor=np.minimum(allele_count,2*called-allele_count)
    return (called>=required)&(minor>=min_mac)


def panel_site_mask(g, columns, max_missing=.1, min_mac=1):
    columns=list(columns)
    c=population_counts(g,columns)
    return mask_from_counts(c.called,c.allele_count,len(columns),max_missing,min_mac)


@dataclass(frozen=True)
class SiteComponents:
    numerator: np.ndarray
    denominator: np.ndarray
    fst: np.ndarray


def wc84_pair(pop1, pop2):
    """Return both-allele per-site components with VCFtools NaN semantics.

    n_i is the number of called diploid individuals in population i, not
    allele count or the nominal population size. Empty populations and
    insufficient sample size consequently yield undefined sites.
    """
    n1=pop1.called.astype(np.float64)
    n2=pop2.called.astype(np.float64)
    k1=pop1.allele_count.astype(np.float64)
    k2=pop2.allele_count.astype(np.float64)
    h1=pop1.heterozygotes.astype(np.float64)
    h2=pop2.heterozygotes.astype(np.float64)
    numerator=np.zeros_like(n1)
    denominator=np.zeros_like(n1)
    with np.errstate(invalid='ignore',divide='ignore',over='ignore'):
        nsum=n1+n2
        nbar=nsum/2.0
        nc=nsum-(n1*n1+n2*n2)/nsum
        hbar=(h1+h2)/nsum
        # Source sums reference then alternate; arithmetic is symmetric.
        for copies1,copies2 in ((2*n1-k1,2*n2-k2),(k1,k2)):
            p1=copies1/(2.0*n1)
            p2=copies2/(2.0*n2)
            pbar=(copies1+copies2)/(2.0*nsum)
            s2=(n1*(p1-pbar)*(p1-pbar)+n2*(p2-pbar)*(p2-pbar))/nbar
            a=(s2-(pbar*(1-pbar)-s2/2.0-hbar/4.0)/(nbar-1.0))*nbar/nc
            b=(pbar*(1-pbar)-s2/2.0-hbar*((2*nbar-1.0)/(4*nbar)))*nbar/(nbar-1.0)
            c=hbar/2.0
            include=~(np.isnan(a)|np.isnan(b)|np.isnan(c))
            numerator+=np.where(include,a,0.0)
            denominator+=np.where(include,a+b+c,0.0)
        fst=numerator/denominator
    return SiteComponents(numerator,denominator,fst)


@dataclass
class WeightedFstAccumulator:
    numerator: float=0.0
    denominator: float=0.0
    site_fst_sum: float=0.0
    selected_sites: int=0
    defined_sites: int=0
    negative_sites: int=0
    infinite_sites: int=0

    def update(self,components,selection=None):
        if selection is None:
            selection=np.ones(components.fst.shape,dtype=bool)
        if selection.dtype!=bool or selection.shape!=components.fst.shape:
            raise ValueError('Selection must be a Boolean vector matching sites')
        include=selection&~np.isnan(components.fst)
        self.selected_sites+=int(selection.sum())
        self.defined_sites+=int(include.sum())
        self.negative_sites+=int((include&(components.fst<0)).sum())
        self.infinite_sites+=int((include&np.isinf(components.fst)).sum())
        self.numerator+=float(components.numerator[include].sum(dtype=np.float64))
        self.denominator+=float(components.denominator[include].sum(dtype=np.float64))
        self.site_fst_sum+=float(components.fst[include].sum(dtype=np.float64))

    @property
    def weighted_fst(self):
        with np.errstate(invalid='ignore',divide='ignore'):
            return float(np.divide(self.numerator,self.denominator))

    @property
    def mean_fst(self):
        with np.errstate(invalid='ignore',divide='ignore'):
            return float(np.divide(self.site_fst_sum,self.defined_sites))

    def as_dict(self):
        return dict(weighted_fst=self.weighted_fst,numerator=self.numerator,
                    denominator=self.denominator,selected_sites=self.selected_sites,
                    defined_sites=self.defined_sites,
                    undefined_sites=self.selected_sites-self.defined_sites,
                    mean_fst=self.mean_fst,negative_sites=self.negative_sites,
                    infinite_sites=self.infinite_sites)
