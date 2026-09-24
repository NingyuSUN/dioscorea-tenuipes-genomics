import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

MODULE = Path(__file__).resolve().parents[1] / 'scripts' / '54_summarize_influence.py'
if not MODULE.exists():  # Isolated development directory; repository uses scripts/.
    MODULE = Path(__file__).resolve().parents[1] / '54_summarize_influence.py'
spec = importlib.util.spec_from_file_location('influence_pca', MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class PCAInfluenceTests(unittest.TestCase):
    def test_rotation_reflection_scale_translation_invariance(self):
        reference = np.array([[0., 0.], [2., 0.], [0., 1.], [2., 2.]])
        transform = np.array([[0., -1.], [-1., 0.]])
        altered = 3.2 * reference @ transform + np.array([7., -8.])
        result, aligned = mod.compare_coordinates(reference, altered)
        self.assertAlmostEqual(result['PC12_pairwise_distance_r'], 1.)
        self.assertLess(result['procrustes_relative_residual'], 1e-12)
        np.testing.assert_allclose(aligned, reference, atol=1e-12)

    def test_common_samples_reordered_and_subset(self):
        baseline = pd.DataFrame({'PC1': [1., 2., 3., 8.], 'PC2': [4., 7., 8., 1.]},
                                index=['a', 'b', 'c', 'd'])
        other = baseline.loc[['d', 'a', 'c']] * 2 + 4
        ids, reference, altered = mod.match_coordinates(baseline, other, ['c', 'd', 'a'])
        self.assertEqual(ids, ['c', 'd', 'a'])
        np.testing.assert_array_equal(reference, baseline.loc[ids].to_numpy())
        metrics, aligned = mod.compare_coordinates(reference, altered)
        self.assertLess(metrics['procrustes_relative_residual'], 1e-12)
        self.assertRaises(ValueError, mod.match_coordinates, baseline, other, ['a', 'b', 'c'])

    def test_collapsed_or_too_few_coordinates_are_explicit(self):
        result, aligned = mod.compare_coordinates(np.zeros((4, 2)), np.ones((4, 2)))
        self.assertIsNone(result['PC12_pairwise_distance_r'])
        self.assertIsNone(result['procrustes_relative_residual'])
        self.assertIsNone(aligned)
        self.assertEqual(result['comparison_status'], 'zero_spread')
        result, aligned = mod.compare_coordinates(np.eye(2), np.eye(2))
        self.assertEqual(result['comparison_status'], 'fewer_than_three_samples')
        self.assertIsNone(result['PC12_pairwise_distance_r'])

    def test_centroid_own_sample_is_excluded(self):
        points = np.array([[0., 0.], [2., 0.], [5., 0.], [7., 0.]])
        groups = ['A', 'A', 'B', 'B']
        result, rows = mod.centroid_separation(['a', 'b', 'c', 'd'], points, groups)
        self.assertAlmostEqual(rows[0]['own_loo_centroid_distance'], 2.)
        self.assertAlmostEqual(rows[0]['nearest_other_centroid_distance'], 6.)
        self.assertEqual(result['own_centroid_nearest_count'], 4)
        _, changed = mod.centroid_separation(['a', 'b', 'c', 'd'], points * 10 + 3, groups)
        self.assertAlmostEqual(rows[0]['normalized_centroid_margin'], changed[0]['normalized_centroid_margin'])

    def test_singleton_group_and_nonfinite_rejected(self):
        result, rows = mod.centroid_separation(['a', 'b', 'c'], np.array([[0., 0.], [5., 0.], [7., 0.]]), ['A', 'B', 'B'])
        self.assertEqual(result['centroid_evaluable_count'], 2)
        self.assertIsNone(rows[0]['normalized_centroid_margin'])
        self.assertRaises(ValueError, mod.compare_coordinates,
                          np.array([[0., 0.], [1., 2.], [3., np.nan]]), np.eye(3, 2))

    def test_header_pca_and_raw_negative_eigenvalues(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path/'pca.eigenvec').write_text('FID IID PC1 PC2\n0 a 1 2\n0 b 3 4\n0 c 5 6\n')
            (path/'pca.eigenval').write_text('2\n1\n-0.01\n')
            frame, eigenvalues = mod.read_pca(path)
            self.assertEqual(list(frame.index), ['a', 'b', 'c'])
            self.assertEqual(eigenvalues.tolist(), [2., 1., -0.01])

    def test_actual_fst_schema_mode_agrees_and_survives_merge(self):
        panels = [{'panel': 'baseline_main', 'kind': 'baseline', 'mode': 'fixed', 'omitted': ''},
                  {'panel': 'baseline_complete', 'kind': 'baseline', 'mode': 'complete', 'omitted': ''}]
        panels += [{'panel': f'{kind}_{i:02d}_{mode}', 'kind': kind, 'mode': mode, 'omitted': 'example'}
                   for kind, count in [('loo', 14), ('pair', 4)] for i in range(count) for mode in ['fixed', 'refilter']]
        table = pd.DataFrame(panels)
        # Production engine includes mode plus diagnostic counts. Exercise that
        # actual integration boundary, not a reduced hand-written schema.
        data = pd.DataFrame([{'panel': p['panel'], 'pair': pair, 'weighted_fst': value,
                              'numerator': value*10, 'denominator': 10., 'selected_sites': 20,
                              'defined_sites': 18, 'mode': p['mode'], 'undefined_sites': 2,
                              'mean_fst': value, 'negative_sites': 1, 'infinite_sites': 0}
                             for p in panels for pair, value in zip(mod.FST_PAIRS, [.2, .21, .1])])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data.to_csv(root/'fst_summary.tsv', sep='\t', index=False)
            result, summary = mod.summarize_fst(root, table, root)
            self.assertEqual(len(result), 38)
            self.assertEqual(summary['panel_pair_rows'], 114)
            self.assertEqual(summary['both_AIC_pairs_exceed_HCO_IZU_panels'], 38)
            merged = pd.read_csv(root/'fst_influence.tsv', sep='\t')
            self.assertIn('mode', merged.columns)
            self.assertNotIn('mode_x', merged.columns)
            self.assertEqual(merged.groupby('mode').size().to_dict(), {'complete': 3, 'fixed': 57, 'refilter': 54})
            data.loc[0, 'mode'] = 'refilter'
            data.to_csv(root/'fst_summary.tsv', sep='\t', index=False)
            with self.assertRaisesRegex(ValueError, 'mode'):
                mod.summarize_fst(root, table, root)


if __name__ == '__main__':
    unittest.main()
