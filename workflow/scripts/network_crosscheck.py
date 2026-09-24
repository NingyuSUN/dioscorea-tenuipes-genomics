#!/usr/bin/env python3
"""Glue for the Neighbor-Net cross-implementation check (Snakemake only).

Produces `results/supplement17/cycle.txt` and `network_cross_validation.json`,
which script 49 reads. On the original run's inputs it reproduces the original
files exactly.

  cycle    CYCLE line of script 45's network.nex -> one 1-based index per line,
           the fixed order read by 47b_network_fixed_order.R.
  compare  Match SciPy-NNLS splits (script 45) with phangorn same-order and
           automatic-order splits (47b / 47) by bipartition; report counts and
           the maximum absolute weight difference over the union of splits,
           an absent split weighing 0.
"""
import argparse, json, re
from pathlib import Path
import pandas as pd


def cycle(nexus, out):
    found = re.findall(r'^\s*CYCLE\s+([\d\s]+);', Path(nexus).read_text(), flags=re.M | re.I)
    if len(found) != 1:
        raise ValueError(f'Expected one CYCLE line in {nexus}, found {len(found)}')
    order = found[0].split()
    if sorted(map(int, order)) != list(range(1, len(order) + 1)):
        raise ValueError('CYCLE is not a permutation of 1..n')
    Path(out).write_text(''.join(i + '\n' for i in order))


def load_splits(path, taxa):
    """Canonical key: the side that excludes the first taxon in sorted order."""
    anchor = min(taxa)
    splits = {}
    table = pd.read_csv(path, sep='\t')
    for side, weight in zip(table['side1'], table['weight']):
        part = frozenset(side.split(','))
        if not part <= taxa or not part or part == taxa:
            raise ValueError(f'Invalid split in {path}: {side}')
        key = frozenset(taxa - part) if anchor in part else part
        if key in splits:
            raise ValueError(f'Duplicate split in {path}: {side}')
        splits[key] = float(weight)
    return splits


def compare(a, b):
    diff = max(abs(a.get(k, 0.) - b.get(k, 0.)) for k in a.keys() | b.keys())
    return len(a), len(b), len(a.keys() & b.keys()), diff


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('cycle')
    c.add_argument('--nexus', type=Path, required=True)
    c.add_argument('--out', type=Path, required=True)
    m = sub.add_parser('compare')
    m.add_argument('--distance', type=Path, required=True, help="script 45 distance.tsv (taxon set)")
    m.add_argument('--primary', type=Path, required=True, help="script 45 network_splits.tsv")
    m.add_argument('--same-order', type=Path, required=True, help="47b phangorn_same_order_splits.tsv")
    m.add_argument('--automatic', type=Path, required=True, help="47 phangorn_splits.tsv")
    m.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.cmd == 'cycle':
        cycle(a.nexus, a.out)
        return
    taxa = frozenset(pd.read_csv(a.distance, sep='\t', index_col=0).index)
    primary = load_splits(a.primary, taxa)
    n1, n2, shared, diff = compare(primary, load_splits(a.same_order, taxa))
    if not n1 == n2 == shared:
        raise ValueError(f'Same-order split systems differ: {n1}/{n2}, shared {shared}')
    auto = compare(primary, load_splits(a.automatic, taxa))
    result = {'same_order_primary_splits': n1, 'same_order_phangorn_splits': n2, 'maximum_weight_difference': diff,
              'automatic_order_comparison': dict(zip(['primary_splits', 'phangorn_splits', 'shared_splits', 'maximum_weight_difference'], auto))}
    a.out.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
