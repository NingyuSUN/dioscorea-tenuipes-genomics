#!/usr/bin/env python3
"""Summarize descriptive PCA/FST influence: 14 deletions and 4 pair controls.

PCA axes can change sign, rotate or change scale after a deletion. Comparisons
therefore use shared sample IDs, distance correlation, and least-squares
translation/rotation/reflection/uniform-scale alignment. Group-centroid checks
exclude the focal sample from its own centroid; they are descriptive checks of
these same fitted samples, not held-out classification or population inference.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = {'AIC': '#d97706', 'HCO': '#2563eb', 'IZU': '#059669'}
MODE_COLORS = {'fixed': '#5b51a5', 'refilter': '#15988f'}
FST_PAIRS = ('AIC_HCO', 'AIC_IZU', 'HCO_IZU')


def read_pca(directory):
    directory = Path(directory)
    frame = pd.read_csv(directory/'pca.eigenvec', sep=r'\s+', dtype={'FID': str, 'IID': str})
    if not {'IID', 'PC1', 'PC2'}.issubset(frame.columns):
        raise ValueError(f'Expected a PLINK eigenvec header in {directory}')
    if frame['IID'].duplicated().any():
        raise ValueError(f'Duplicate sample IDs in {directory}')
    frame = frame.set_index('IID')
    if not np.isfinite(frame[['PC1', 'PC2']].to_numpy(dtype=float)).all():
        raise ValueError(f'Nonfinite PCA coordinates in {directory}')
    eigenvalues = np.atleast_1d(np.loadtxt(directory/'pca.eigenval'))
    if not np.isfinite(eigenvalues).all():
        raise ValueError(f'Nonfinite eigenvalues in {directory}')
    return frame, eigenvalues


def match_coordinates(baseline, other, expected_ids):
    ids = list(expected_ids)
    if len(ids) != len(set(ids)) or set(other.index) != set(ids):
        raise ValueError('PCA sample set does not match the declared retained samples')
    if not set(ids).issubset(baseline.index):
        raise ValueError('Retained samples absent from baseline PCA')
    return ids, baseline.loc[ids, ['PC1', 'PC2']].to_numpy(dtype=float), other.loc[ids, ['PC1', 'PC2']].to_numpy(dtype=float)


def compare_coordinates(reference, altered):
    reference, altered = np.asarray(reference, dtype=float), np.asarray(altered, dtype=float)
    if reference.shape != altered.shape or reference.ndim != 2 or reference.shape[1] != 2:
        raise ValueError('Expected equal n by 2 coordinate arrays')
    if not np.isfinite(reference).all() or not np.isfinite(altered).all():
        raise ValueError('PCA coordinates must be finite')
    result = {'shared_sample_count': len(reference), 'PC12_pairwise_distance_r': None,
              'procrustes_relative_residual': None, 'procrustes_scale': None,
              'comparison_status': 'ok'}
    if len(reference) < 3:
        result['comparison_status'] = 'fewer_than_three_samples'
        return result, None
    y = reference-reference.mean(axis=0)
    x = altered-altered.mean(axis=0)
    ny, nx = float(np.linalg.norm(y)), float(np.linalg.norm(x))
    result.update(reference_rank=int(np.linalg.matrix_rank(y)), panel_rank=int(np.linalg.matrix_rank(x)))
    if ny == 0 or nx == 0:
        result['comparison_status'] = 'zero_spread'
        return result, None
    # Normalize before SVD to avoid overflow; reflections are allowed since the
    # signs of PCA eigenvectors carry no meaning.
    u, singular, vt = np.linalg.svd((x/nx).T @ (y/ny))
    rotation = u @ vt
    scale = float(singular.sum() * ny / nx)
    aligned = scale*x@rotation + reference.mean(axis=0)
    result['procrustes_scale'] = scale
    result['procrustes_relative_residual'] = float(np.linalg.norm(aligned-reference)/ny)
    d0, d1 = pdist(y/ny), pdist(x/nx)
    if np.std(d0) > 0 and np.std(d1) > 0:
        result['PC12_pairwise_distance_r'] = float(np.corrcoef(d0, d1)[0, 1])
    else:
        result['comparison_status'] = 'constant_pairwise_distance'
    if min(result['reference_rank'], result['panel_rank']) < 2:
        result['comparison_status'] = 'rank_deficient'
    return result, aligned


def centroid_separation(ids, coords, groups):
    coords, groups = np.asarray(coords, dtype=float), np.asarray(groups)
    if coords.shape != (len(ids), 2) or len(groups) != len(ids) or not np.isfinite(coords).all():
        raise ValueError('Invalid centroid inputs')
    radius = float(np.sqrt(np.mean(np.sum((coords-coords.mean(axis=0))**2, axis=1))))
    rows = []
    for i, sample in enumerate(ids):
        own = np.flatnonzero((groups == groups[i]) & (np.arange(len(ids)) != i))
        others = sorted(set(groups)-{groups[i]})
        row = {'sample': sample, 'group': str(groups[i]), 'own_peer_count': len(own),
               'own_loo_centroid_distance': None, 'nearest_other_group': None,
               'nearest_other_centroid_distance': None, 'centroid_margin': None,
               'normalized_centroid_margin': None, 'own_centroid_nearest': None,
               'status': 'ok'}
        if not len(own) or not others or radius == 0:
            row['status'] = 'singleton_group' if not len(own) else ('no_other_group' if not others else 'zero_spread')
        else:
            own_distance = float(np.linalg.norm(coords[i]-coords[own].mean(axis=0)))
            other_distances = [(float(np.linalg.norm(coords[i]-coords[groups == g].mean(axis=0))), str(g)) for g in others]
            other_distance, other_group = min(other_distances)
            margin = other_distance-own_distance
            row.update(own_loo_centroid_distance=own_distance, nearest_other_group=other_group,
                       nearest_other_centroid_distance=other_distance, centroid_margin=margin,
                       normalized_centroid_margin=margin/radius, own_centroid_nearest=bool(margin > 0))
        rows.append(row)
    valid = [r for r in rows if r['normalized_centroid_margin'] is not None]
    summary = {'centroid_evaluable_count': len(valid),
               'own_centroid_nearest_count': sum(r['own_centroid_nearest'] for r in valid),
               'centroid_tie_count': sum(r['centroid_margin'] == 0 for r in valid),
               'min_normalized_centroid_margin': min((r['normalized_centroid_margin'] for r in valid), default=None),
               'PC12_rms_radius': radius}
    return summary, rows


def sample_groups(ids, metadata):
    if metadata is not None:
        groups = [str(metadata.loc[s, 'group']) for s in ids]
    else:
        groups = [s.split('-', 1)[0] for s in ids]
    if not set(groups).issubset(COLORS):
        raise ValueError('Only the declared AIC/HCO/IZU study groups are allowed')
    return groups


def validate_manifest(panels):
    names = [p['name'] for p in panels]
    if len(names) != 38 or len(set(names)) != 38:
        raise ValueError('Expected exactly 38 unique panels: 2 baselines + 28 LOO + 8 pair controls')
    expected = {'baseline_main', 'baseline_complete'}
    expected.update(f'loo_{i:02d}_{mode}' for i in range(14) for mode in ('fixed', 'refilter'))
    expected.update(f'pair_{i:02d}_{mode}' for i in range(4) for mode in ('fixed', 'refilter'))
    if set(names) != expected:
        raise ValueError('Missing or unexpected named perturbations')
    lookup = {p['name']: p for p in panels}
    baseline = set(lookup['baseline_main']['keep'])
    if len(baseline) != 14 or set(lookup['baseline_complete']['keep']) != baseline:
        raise ValueError('Baselines must contain the same 14 study samples')
    for p in panels:
        keep, omitted = p['keep'], p.get('omitted', [])
        expected_kind = 'baseline' if p['name'].startswith('baseline_') else p['name'].split('_', 1)[0]
        expected_mode = ('complete' if p['name'] == 'baseline_complete' else 'fixed') if expected_kind == 'baseline' else p['name'].rsplit('_', 1)[1]
        if p['kind'] != expected_kind or p['mode'] != expected_mode:
            raise ValueError(f'Panel kind/mode conflicts with its name: {p["name"]}')
        expected_n = 14 if p['kind'] == 'baseline' else (13 if p['kind'] == 'loo' else 12)
        if len(keep) != expected_n or len(set(keep)) != expected_n:
            raise ValueError(f'Incorrect retained sample count: {p["name"]}')
        if set(keep) & set(omitted) or set(keep) | set(omitted) != baseline:
            raise ValueError(f'Inconsistent keep/omitted sample partition: {p["name"]}')
    for kind, count in (('loo', 14), ('pair', 4)):
        retained_sets = []
        for i in range(count):
            if set(lookup[f'{kind}_{i:02d}_fixed']['keep']) != set(lookup[f'{kind}_{i:02d}_refilter']['keep']):
                raise ValueError('Fixed and refilter scenarios must retain the same samples')
            retained_sets.append(frozenset(lookup[f'{kind}_{i:02d}_fixed']['keep']))
        if len(set(retained_sets)) != count:
            raise ValueError('Repeated retained-sample scenarios')


def summarize_fst(root, panel_table, out):
    path = root/'fst_summary.tsv'
    if not path.exists():
        return panel_table, {'status': 'not_supplied'}
    fst = pd.read_csv(path, sep='\t')
    required = {'panel', 'pair', 'weighted_fst', 'numerator', 'denominator', 'selected_sites', 'defined_sites'}
    if not required.issubset(fst.columns) or len(fst) != 114 or fst.duplicated(['panel', 'pair']).any():
        raise ValueError('Expected 114 unique panel/pair FST rows with estimator totals')
    if set(fst['panel']) != set(panel_table['panel']) or set(fst['pair']) != set(FST_PAIRS):
        raise ValueError('FST panel or population pair mismatch')
    if 'mode' in fst.columns:
        declared_modes = fst['panel'].map(panel_table.set_index('panel')['mode'])
        if not fst['mode'].eq(declared_modes).all():
            raise ValueError('FST site-selection mode disagrees with the panel manifest')
        # The engine already records mode. Keep the manifest's validated copy
        # when joining metadata, rather than introducing mode_x/mode_y columns.
        fst = fst.drop(columns='mode')
    if not np.isfinite(fst[['weighted_fst', 'numerator', 'denominator']]).all().all() or (fst.denominator <= 0).any():
        raise ValueError('Nonfinite or undefined FST totals require an explicit handling decision')
    if not np.allclose(fst.weighted_fst, fst.numerator/fst.denominator, atol=1e-12, rtol=1e-10):
        raise ValueError('Weighted FST does not match the stored ratio of component sums')
    baseline = fst[fst.panel == 'baseline_main'].set_index('pair').weighted_fst
    fst['baseline_main_weighted_fst'] = fst.pair.map(baseline)
    fst['delta_vs_baseline_main'] = fst.weighted_fst-fst.baseline_main_weighted_fst
    fst = fst.merge(panel_table[['panel', 'kind', 'mode', 'omitted']], on='panel', validate='many_to_one')
    fst.to_csv(out/'fst_influence.tsv', sep='\t', index=False)
    values = fst.pivot(index='panel', columns='pair', values='weighted_fst')
    values.columns = [f'weighted_fst_{p}' for p in values.columns]
    values['smaller_AIC_pair_minus_HCO_IZU'] = values[['weighted_fst_AIC_HCO', 'weighted_fst_AIC_IZU']].min(axis=1)-values['weighted_fst_HCO_IZU']
    values['both_AIC_pairs_exceed_HCO_IZU'] = values['smaller_AIC_pair_minus_HCO_IZU'] > 0
    panel_table = panel_table.merge(values, left_on='panel', right_index=True, validate='one_to_one')
    ranges = []
    for (kind, mode, pair), group in fst.groupby(['kind', 'mode', 'pair'], sort=True):
        ranges.append({'kind': kind, 'mode': mode, 'pair': pair, 'panels': len(group),
                       'min_weighted_fst': float(group.weighted_fst.min()), 'max_weighted_fst': float(group.weighted_fst.max()),
                       'min_delta_vs_baseline_main': float(group.delta_vs_baseline_main.min()),
                       'max_delta_vs_baseline_main': float(group.delta_vs_baseline_main.max())})
    summary = {'status': 'complete', 'panel_pair_rows': len(fst), 'ranges': ranges,
               'both_AIC_pairs_exceed_HCO_IZU_panels': int(values.both_AIC_pairs_exceed_HCO_IZU.sum()),
               'total_panels': len(values),
               'minimum_AIC_vs_HCO_IZU_gap': float(values.smaller_AIC_pair_minus_HCO_IZU.min()),
               'interpretation': 'Weighted WC84 FST on filtered SNPs; perturbation ranges are not confidence intervals. The ordering between the two AIC contrasts is not a primary stability claim.'}
    return panel_table, summary


def save_figure(fig, out, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out/f'{name}.{ext}', dpi=210)
    plt.close(fig)


def plot_overlays(panels, baseline, aligned, groups, out):
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 13), layout='constrained')
    cases = [('baseline_main', 'Baseline: 14 study samples'),
             ('baseline_complete', 'Complete-genotype panel; aligned'),
             (('loo', 'fixed'), '14 individual deletions | fixed sites'),
             (('loo', 'refilter'), '14 individual deletions | refiltered'),
             (('pair', 'fixed'), '4 candidate-pair controls | fixed sites'),
             (('pair', 'refilter'), '4 candidate-pair controls | refiltered')]
    lookup = {p['name']: p for p in panels}
    all_points = np.concatenate([baseline[['PC1', 'PC2']].to_numpy()]+[x for x in aligned.values() if x is not None])
    lower, upper = all_points.min(axis=0), all_points.max(axis=0)
    padding = np.maximum((upper-lower)*.20, .01)
    for ax, (selector, title) in zip(axes.flat, cases):
        for group, color in COLORS.items():
            ids = [s for s in baseline.index if groups[s] == group]
            pts = baseline.loc[ids, ['PC1', 'PC2']].to_numpy()
            ax.scatter(pts[:, 0], pts[:, 1], s=64, facecolors='none', edgecolors=color, linewidths=1.3, zorder=3)
        names = [selector] if isinstance(selector, str) else [p['name'] for p in panels if (p['kind'], p['mode']) == selector]
        for name in names:
            coords = aligned.get(name)
            if coords is None:
                continue
            ids = lookup[name]['keep']
            ref = baseline.loc[ids, ['PC1', 'PC2']].to_numpy()
            for i, sample in enumerate(ids):
                color = COLORS[groups[sample]]
                ax.plot([ref[i, 0], coords[i, 0]], [ref[i, 1], coords[i, 1]], color=color, alpha=.23, linewidth=.65)
                ax.scatter(coords[i, 0], coords[i, 1], c=color, s=21 if len(names) > 1 else 35, alpha=.42 if len(names) > 1 else .9, linewidths=0, zorder=4)
            if name == 'baseline_main':
                for group, color in COLORS.items():
                    group_points = coords[[groups[s] == group for s in ids]]
                    ax.annotate(f'{group} (n={len(group_points)})', group_points.mean(axis=0),
                                xytext=(10, 12), textcoords='offset points', fontsize=9, color=color)
        ax.set(title=title, xlabel='Aligned PC1', ylabel='Aligned PC2')
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(lower[0]-padding[0], upper[0]+padding[0])
        ax.set_ylim(lower[1]-padding[1], upper[1]+padding[1])
        ax.spines[['top', 'right']].set_visible(False)
    handles = [Line2D([], [], marker='o', linestyle='', color=c, label=g) for g, c in COLORS.items()]
    handles += [Line2D([], [], marker='o', linestyle='', markerfacecolor='none', color='#666666', label='Main-panel position')]
    fig.legend(handles=handles, loc='outside lower center', ncol=4, frameon=False)
    fig.suptitle('D. tenuipes | PCA influence on retained individuals\nTranslation, rotation/reflection and uniform scale fitted to shared baseline samples', fontsize=13)
    save_figure(fig, out, 'pca_influence_overlays')


def scenario_order(panels):
    scenarios = [f'loo_{i:02d}' for i in range(14)] + [f'pair_{i:02d}' for i in range(4)]
    lookup = {p['name']: p for p in panels}
    labels = []
    for stem in scenarios:
        omitted = lookup[stem+'_fixed']['omitted']
        if stem.startswith('loo'):
            labels.append('omit '+omitted[0])
        else:
            labels.append('pair '+stem[-2:]+'\nomit '+ '\n'.join(omitted))
    return scenarios, labels


def finite_stat(series, operation):
    """Keep undefined/degenerate comparisons null in machine-readable summaries."""
    values = np.asarray(series, dtype=float)
    values = values[np.isfinite(values)]
    return float(getattr(np, operation)(values)) if len(values) else None


def plot_metrics(panels, table, out):
    scenarios, labels = scenario_order(panels)
    indexed = table.set_index('panel')
    specs = [('PC12_pairwise_distance_r', 'PC1/2 distance correlation'),
             ('procrustes_relative_residual', 'Relative Procrustes residual'),
             ('min_normalized_centroid_margin', 'Minimum centroid margin / RMS radius')]
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True, layout='constrained')
    for ax, (metric, ylabel) in zip(axes, specs):
        for mode, color in MODE_COLORS.items():
            values = [indexed.loc[s+'_'+mode, metric] for s in scenarios]
            ax.plot(range(18), values, 'o-', color=color, label=mode, markersize=4, linewidth=1.2)
        for panel, color, style, label in [('baseline_main', '#4b5563', '--', 'main baseline'), ('baseline_complete', '#c25e11', ':', 'complete genotypes')]:
            ax.axhline(indexed.loc[panel, metric], color=color, linestyle=style, linewidth=1., label=label)
        ax.axvline(13.5, color='#bbbbbb', linewidth=.8)
        ax.set(ylabel=ylabel)
        if metric == 'PC12_pairwise_distance_r':
            low = finite_stat(indexed[metric], 'min')
            ax.set_ylim(max(-1.005, (low-.02) if low is not None else -.02), 1.005)
        ax.grid(axis='y', color='#eeeeee')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(ncol=4, fontsize=9, frameon=False)
    axes[-1].set_xticks(range(18), labels, rotation=60, ha='right', fontsize=7.5)
    fig.suptitle('Descriptive PCA influence | 14 individual deletions and 4 pair controls\nPerturbations of the observed samples; no bootstrap support or confidence intervals', fontsize=13)
    save_figure(fig, out, 'pca_influence_metrics')
    if not all(f'weighted_fst_{p}' in indexed.columns for p in FST_PAIRS):
        return
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True, layout='constrained')
    for ax, pair in zip(axes, FST_PAIRS):
        metric = f'weighted_fst_{pair}'
        for mode, color in MODE_COLORS.items():
            ax.plot(range(18), [indexed.loc[s+'_'+mode, metric] for s in scenarios], 'o-', color=color, markersize=4, linewidth=1.2, label=mode)
        for panel, color, style, label in [('baseline_main', '#4b5563', '--', 'main baseline'), ('baseline_complete', '#c25e11', ':', 'complete genotypes')]:
            ax.axhline(indexed.loc[panel, metric], color=color, linestyle=style, linewidth=1, label=label)
        ax.axvline(13.5, color='#bbbbbb', linewidth=.8)
        ax.set(ylabel=pair.replace('_', ' vs ')+'\nweighted WC84 FST')
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#eeeeee')
    axes[0].legend(ncol=4, fontsize=9, frameon=False)
    axes[-1].set_xticks(range(18), labels, rotation=60, ha='right', fontsize=7.5)
    fig.suptitle('Filtered-SNP FST influence | 14 individual deletions and 4 pair controls\nFixed sites and sample-dependent refiltering shown separately; no confidence intervals', fontsize=13)
    save_figure(fig, out, 'fst_influence')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, help='Optional sample/group TSV; otherwise use validated study sample prefixes')
    args = parser.parse_args()
    panels = json.loads((args.root/'panels.json').read_text())['panels']
    validate_manifest(panels)
    args.out.mkdir(parents=True, exist_ok=False)
    metadata = pd.read_csv(args.metadata, sep='\t').set_index('sample') if args.metadata else None
    baseline, baseline_eigenvalues = read_pca(args.root/'panels'/'baseline_main')
    lookup = {p['name']: p for p in panels}
    if set(baseline.index) != set(lookup['baseline_main']['keep']):
        raise ValueError('Baseline PCA and manifest disagree')
    all_groups = dict(zip(baseline.index, sample_groups(baseline.index, metadata)))
    rows, centroid_rows, aligned_rows, aligned_by_panel, raw_eigenvalues = [], [], [], {}, {}
    for panel in panels:
        name = panel['name']
        path = args.root/'panels'/name
        coords, eigenvalues = read_pca(path)
        info = json.loads((path/'panel.json').read_text())
        ids, ref, changed = match_coordinates(baseline, coords, panel['keep'])
        if int(info['sample_count']) != len(ids) or int(info['pca_snps']) <= 0:
            raise ValueError(f'Incorrect panel counts: {name}')
        comparison, aligned = compare_coordinates(ref, changed)
        separation, samples = centroid_separation(ids, changed, [all_groups[s] for s in ids])
        row = {'panel': name, 'kind': panel['kind'], 'mode': panel['mode'],
               'omitted': ';'.join(panel.get('omitted', [])), 'sample_count': len(ids),
               'pca_snps': int(info['pca_snps']), **comparison, **separation}
        rows.append(row)
        centroid_rows.extend({'panel': name, **sample} for sample in samples)
        raw_eigenvalues[name] = eigenvalues.tolist()
        aligned_by_panel[name] = aligned
        if aligned is not None:
            for i, sample in enumerate(ids):
                aligned_rows.append({'panel': name, 'sample': sample, 'group': all_groups[sample],
                                     'raw_PC1': float(changed[i, 0]), 'raw_PC2': float(changed[i, 1]),
                                     'aligned_PC1': float(aligned[i, 0]), 'aligned_PC2': float(aligned[i, 1]),
                                     'baseline_PC1': float(ref[i, 0]), 'baseline_PC2': float(ref[i, 1])})
    table = pd.DataFrame(rows)
    table, fst_summary = summarize_fst(args.root, table, args.out)
    table.to_csv(args.out/'panel_influence.tsv', sep='\t', index=False)
    pd.DataFrame(centroid_rows).to_csv(args.out/'centroid_distances.tsv', sep='\t', index=False)
    pd.DataFrame(aligned_rows).to_csv(args.out/'aligned_pca.tsv', sep='\t', index=False)
    original = args.root/'baseline_original'
    reproduction = {'status': 'not_supplied'}
    if (original/'pca.eigenvec').exists():
        previous, previous_eigenvalues = read_pca(original)
        _, ref, new = match_coordinates(previous, baseline, lookup['baseline_main']['keep'])
        metrics, _ = compare_coordinates(ref, new)
        reproduction = {'status': 'compared', **metrics,
                        'max_absolute_eigenvalue_difference': float(np.max(np.abs(previous_eigenvalues-baseline_eigenvalues))) if previous_eigenvalues.shape == baseline_eigenvalues.shape else None,
                        'previous_eigenvalues': previous_eigenvalues.tolist(),
                        'current_eigenvalues': baseline_eigenvalues.tolist(),
                        'method': 'Separate baseline rerun vs archived main14, matched by retained sample IDs; sign/rotation/scale invariance checked.'}
    ranges = []
    for (kind, mode), group in table.groupby(['kind', 'mode'], sort=True):
        ranges.append({'kind': kind, 'mode': mode, 'panels': len(group),
                       'distance_r_min': finite_stat(group.PC12_pairwise_distance_r, 'min'),
                       'distance_r_max': finite_stat(group.PC12_pairwise_distance_r, 'max'),
                       'relative_procrustes_residual_max': finite_stat(group.procrustes_relative_residual, 'max'),
                       'min_normalized_centroid_margin': finite_stat(group.min_normalized_centroid_margin, 'min'),
                       'all_own_centroid_nearest_panels': int((group.own_centroid_nearest_count == group.sample_count).sum())})
    summary = {'panels': len(panels), 'ranges': ranges, 'baseline_reproduction': reproduction, 'fst': fst_summary,
               'raw_eigenvalues': raw_eigenvalues,
               'definitions': {
                   'comparison_reference': 'baseline_main restricted to the same retained sample IDs; the baseline itself is not refitted for comparison',
                   'PC12_pairwise_distance_r': 'Pearson correlation between the two sets of all shared-sample Euclidean PC1/2 distances; not a p-value',
                   'procrustes_relative_residual': 'Frobenius norm of fitted-PC1/2 minus baseline-PC1/2 divided by centered baseline Frobenius norm; fit allows translation, orthogonal rotation/reflection and one nonnegative uniform scale',
                   'centroid_margin': 'Distance to nearest other-group full centroid minus distance to own-group centroid excluding the focal sample',
                   'normalized_centroid_margin': 'Centroid margin divided by global PC1/2 root-mean-square radius for that panel',
                   'own_centroid_nearest_count': 'Number of samples whose own leave-one-out group centroid is strictly nearer than every other-group centroid; same fitted samples, not held-out classification',
                   'variance_labels': 'Raw eigenvalues retained including negatives. PC percentages are not labelled explained variance.',
                   'scope': 'Descriptive influence on the 14 observed study samples; 14 deletions and 4 candidate-pair controls in each of two site-selection modes, not bootstrap samples or independent population validation.'}}
    (args.out/'influence_summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    plot_overlays(panels, baseline, aligned_by_panel, all_groups, args.out)
    plot_metrics(panels, table, args.out)
    (args.out/'SUMMARY_COMPLETE').write_text('38 PCA panels summarized; sample IDs/counts and raw eigenvalues checked. FST status: '+fst_summary['status']+'.\n')


if __name__ == '__main__':
    main()
