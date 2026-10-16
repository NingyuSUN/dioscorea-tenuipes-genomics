"""Audit P50 fixed nuclear quartets independently using Bio.Phylo four-point sums.

Copy into P50/workflow before use, or supply --root. Reads original P24 trees.
Branch lengths and support labels are deliberately ignored: every structural
edge has unit length, including an edge whose stored length is exactly zero.
No production parser or descendant-mask code is imported.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import traceback

from Bio import Phylo


R = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path('workspace/wgs/dt_research_p24_20260928')
FOCAL = (
    'Dioscorea_tenuipes_historical',
    'Dioscorea_tokoro',
    'Dioscorea_gracillima',
    'Dioscorea_zingiberensis',
)
CLASSES = ('DT_TOK', 'DT_GRA', 'DT_ZIN', 'unresolved')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def unique_rows(rows, label):
    out = {}
    for row in rows:
        gene = row['gene']
        require(isinstance(gene, str) and gene and Path(gene).name == gene,
                f'{label}: unsafe or invalid gene identifier {gene!r}')
        require(gene not in out, f'{label}: repeated locus {gene}')
        out[gene] = row
    return out


def validate_tips(tree, expected, label):
    terminals = tree.get_terminals()
    names = [leaf.name for leaf in terminals]
    require(len(names) == len(set(names)), f'{label}: duplicate leaf name')
    require(set(names) == set(expected),
            f'{label}: tip mismatch; absent={sorted(set(expected)-set(names))}; '
            f'unexpected={sorted(set(names)-set(expected))}')
    return dict(zip(names, terminals))


def classify_four_point(tree, expected, label):
    """Unique minimum of the three unit-edge distance sums identifies a split."""
    tips = validate_tips(tree, expected, label)
    if not set(FOCAL).issubset(tips):
        return 'missing_taxa', None
    # Bio.Phylo get_path excludes root and includes the terminal node.
    # Common path-prefix length gives the LCA depth without branch lengths.
    paths = [tree.get_path(tips[name]) for name in FOCAL]

    def distance(i, j):
        common = 0
        for left, right in zip(paths[i], paths[j]):
            if left is not right:
                break
            common += 1
        return len(paths[i]) + len(paths[j]) - 2 * common

    sums = (
        distance(0, 1) + distance(2, 3),
        distance(0, 2) + distance(1, 3),
        distance(0, 3) + distance(1, 2),
    )
    ordered = sorted(sums)
    # For an additive tree metric the largest two sums must be equal.
    require(ordered[1] == ordered[2], f'{label}: non-tree four-point sums {sums}')
    if ordered[0] == ordered[2]:
        return 'unresolved', list(sums)
    require(sums.count(ordered[0]) == 1, f'{label}: ambiguous minimum {sums}')
    return CLASSES[sums.index(ordered[0])], list(sums)


def compare_counts(counts, expected_total, label):
    require(isinstance(counts, dict) and set(counts) == set(CLASSES),
            f'{label}: expected exactly the four quartet count keys')
    require(all(type(n) is int and n >= 0 for n in counts.values()),
            f'{label}: counts must be nonnegative integers (not booleans)')
    require(sum(counts.values()) == expected_total,
            f'{label}: counts do not sum to {expected_total}')


def audit(root, source, result_path, manifest_path, original_audit_path, receipt_path):
    require(not receipt_path.exists(), f'Refuse to overwrite existing audit: {receipt_path}')
    started = time.monotonic()
    receipt = {
        'status': 'RUNNING',
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'method': 'Independent Bio.Phylo root-to-tip paths and unit-edge four-point distance sums',
        'structural_zero_length_edges': 'retained',
        'focal_taxa_order': list(FOCAL),
        'sources': {'production_result': str(result_path), 'manifest': str(manifest_path),
                    'original_audit': str(original_audit_path), 'source_root': str(source)},
        'records': [],
    }
    try:
        result = read_json(result_path)
        manifest = read_json(manifest_path)
        source_audit = read_json(original_audit_path)
        production = unique_rows(result['loci'], 'production')
        originals = unique_rows(manifest['loci'], 'original manifest')
        previous = unique_rows(source_audit['records'], 'original audit')
        require(len(originals) == 173, f'Expected 173 original loci, got {len(originals)}')
        require('2661at3193' not in originals, 'Excluded locus 2661 found in original manifest')
        require(set(production) == set(originals) == set(previous), 'Locus sets disagree')
        receipt['input_sha256'] = {str(path): digest(path) for path in
                                  (result_path, manifest_path, original_audit_path)}
        totals = Counter()
        points = Counter()
        pooled = Counter({key: 0 for key in CLASSES})
        for gene, original in originals.items():
            supplied = production[gene]
            expected = original['retained_taxa']
            require(len(expected) == len(set(expected)), f'{gene}: duplicate manifest taxa')
            require(len(supplied['expected_taxa']) == len(set(supplied['expected_taxa'])),
                    f'{gene}: duplicate production expected_taxa')
            require(set(supplied['expected_taxa']) == set(expected),
                    f'{gene}: expected taxa differ from the original manifest')
            missing = sorted(set(FOCAL) - set(expected))
            require(sorted(supplied['missing_taxa']) == missing,
                    f'{gene}: production missing taxa disagree')
            eligible = not missing
            require(type(supplied['eligible']) is bool and supplied['eligible'] == eligible,
                    f'{gene}: eligibility disagrees with original taxon set')
            folder = source / 'gene_trees' / gene
            point_file = folder / 'gene.treefile'
            boot_file = folder / 'gene.ufboot'
            point_hash, boot_hash = digest(point_file), digest(boot_file)
            require(point_hash == previous[gene]['tree_sha256'], f'{gene}: original point tree hash mismatch')
            require(boot_hash == previous[gene]['bootstrap_sha256'], f'{gene}: original bootstrap hash mismatch')
            point = Phylo.read(point_file, 'newick')
            point_split, point_sums = classify_four_point(point, expected, gene + ':point')
            require(point_split == supplied['point_split'], f'{gene}: point split disagrees')
            row = {'gene': gene, 'eligible': eligible, 'missing_taxa': missing,
                   'point_split': point_split, 'point_distance_sums': point_sums,
                   'tree_sha256': point_hash, 'bootstrap_sha256': boot_hash}
            totals['total_loci'] += 1
            points[point_split] += 1
            if not eligible:
                require(supplied['bootstrap_counts'] is None,
                        f'{gene}: missing-taxon locus must have null bootstrap counts')
                row['bootstrap_counts'] = None
                row['bootstrap_trees_parsed'] = 0
                totals['ineligible_loci'] += 1
            else:
                counts = Counter({key: 0 for key in CLASSES})
                n = 0
                for n, tree in enumerate(Phylo.parse(boot_file, 'newick'), 1):
                    require(n <= 1000, f'{gene}: more than 1000 bootstrap trees')
                    split, _ = classify_four_point(tree, expected, f'{gene}:bootstrap:{n}')
                    require(split in CLASSES, f'{gene}: eligible bootstrap lost a focal taxon')
                    counts[split] += 1
                require(n == 1000, f'{gene}: expected 1000 bootstrap trees, found {n}')
                compare_counts(supplied['bootstrap_counts'], 1000, gene + ':production')
                require(dict(counts) == supplied['bootstrap_counts'],
                        f'{gene}: independent bootstrap counts disagree: {dict(counts)}')
                row['bootstrap_counts'] = dict(counts)
                row['bootstrap_trees_parsed'] = n
                totals['eligible_loci'] += 1
                totals['bootstrap_trees_parsed'] += n
                pooled.update(counts)
            receipt['records'].append(row)
        scalar_totals = {key: totals[key] for key in
                         ('total_loci', 'eligible_loci', 'ineligible_loci', 'bootstrap_trees_parsed')}
        require(scalar_totals['total_loci'] == scalar_totals['eligible_loci'] + scalar_totals['ineligible_loci'],
                'Scalar locus conservation failed')
        require(scalar_totals['bootstrap_trees_parsed'] == 1000 * scalar_totals['eligible_loci'],
                'Scalar bootstrap conservation failed')
        require(sum(points.values()) == scalar_totals['total_loci'], 'Point total conservation failed')
        compare_counts(dict(pooled), scalar_totals['bootstrap_trees_parsed'], 'pooled independent counts')
        # The production schema may omit aggregates. If it publishes these
        # explicitly named totals, compare rather than silently trusting them.
        aggregate_checks = []
        for field, computed in (('scalar_totals', scalar_totals),
                                ('point_counts', {key: points[key] for key in CLASSES + ('missing_taxa',)}),
                                ('pooled_bootstrap_counts', dict(pooled))):
            if field in result:
                require(result[field] == computed, f'Production {field} disagrees')
                aggregate_checks.append(field)
        receipt.update(status='PASS_INDEPENDENT_FOUR_POINT_AUDIT', scalar_totals=scalar_totals,
                       point_counts={key: points[key] for key in CLASSES + ('missing_taxa',)},
                       pooled_bootstrap_counts=dict(pooled), production_aggregate_fields_checked=aggregate_checks,
                       comparison_scope='Every original point tree; every bootstrap tree for four-taxon-eligible loci; '
                       'all 173 point/bootstrap source hashes; exact original expected leaf sets',
                       limits=['Ineligible bootstrap corpora were hash-checked but not parsed in this audit',
                               'Pooled replicate counts are computational totals, not independent biological observations',
                               'Agreement establishes computational consistency, not orthology or evolutionary cause'])
    except Exception as exc:
        receipt.update(status='FAIL_INDEPENDENT_FOUR_POINT_AUDIT', error=str(exc),
                       traceback=traceback.format_exc())
        raise
    finally:
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        receipt['elapsed_seconds'] = time.monotonic() - started
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        with receipt_path.open('x') as handle:
            json.dump(receipt, handle, indent=2)
            handle.write('\n')
    print(json.dumps({'status': receipt['status'], 'receipt': str(receipt_path),
                      'scalar_totals': receipt['scalar_totals']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=R)
    parser.add_argument('--source-root', type=Path, default=SOURCE_ROOT)
    parser.add_argument('--result', type=Path)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--original-audit', type=Path)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    audit(args.root, args.source_root,
          args.result or args.root / 'results/NUCLEAR_QUARTETS.json',
          args.manifest or args.source_root / 'metadata/LOCUS_MANIFEST.json',
          args.original_audit or args.source_root / 'validation/ALL_LOCUS_AUDIT.json',
          args.receipt or args.root / 'validation/INDEPENDENT_QUARTET_AUDIT.json')
