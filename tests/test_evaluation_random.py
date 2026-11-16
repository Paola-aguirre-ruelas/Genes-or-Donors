"""Donor-level macro-F1, false-positive accounting, random-reference panels, data guard."""
import unittest

import numpy as np

from compact_panels import evaluation, random_reference
from compact_panels.evaluation import score_individual

C = ["NK_16hi", "CD4.Naive", "CD4.CM", "CD8.Naive", "CD8.TE"]


class TestMacroF1(unittest.TestCase):
    def test_non_evaluable_cell_counts_as_false_positive(self):
        yt = np.array([C[0]] * 20 + [C[1]] * 20 + [C[2]] * 20 + [C[3]])
        yp = yt.copy(); yp[-1] = C[0]
        got, terms = score_individual(yt, yp, C[:3], C)
        self.assertAlmostEqual(terms[C[0]], 40 / 41, places=15)
        self.assertAlmostEqual(got, ((40 / 41) + 1 + 1) / 3, places=15)

    def test_class_below_threshold_excluded_from_average(self):
        truth = np.array([C[0]] * 25 + [C[1]] * 25 + [C[2]] * 25 + [C[3]] * 19)
        self.assertEqual(evaluation.evaluable_classes(truth, C, 20), C[:3])
        pred = truth.copy(); pred[-19:] = C[4]                 # all CD8.Naive cells wrong: no effect on C_d terms
        got, terms = score_individual(truth, pred, C[:3], C)
        self.assertEqual(set(terms), set(C[:3]))
        self.assertEqual(got, 1.0)

    def test_fewer_than_three_evaluable_classes_rejected(self):
        truth = np.array([C[0]] * 25 + [C[1]] * 25 + [C[2]] * 5)
        with self.assertRaises(ValueError):
            evaluation.score_donors(truth, np.array(["d"] * len(truth)), {"w0": truth}, ["d"], C, 20, 3)

    def test_donor_equal_mean_differs_from_pooled_cells(self):
        # Donor A: 200 cells perfectly predicted; donor B: 60 cells, poorly predicted.
        ta = np.array([C[0]] * 100 + [C[1]] * 50 + [C[2]] * 50)
        tb = np.array([C[0]] * 20 + [C[1]] * 20 + [C[2]] * 20)
        pb = tb.copy(); pb[:10] = C[1]; pb[20:30] = C[2]; pb[40:50] = C[0]
        truth = np.concatenate([ta, tb]); pred = np.concatenate([ta, pb])
        donors = np.array(["A"] * len(ta) + ["B"] * len(tb))
        rows = evaluation.score_donors(truth, donors, {"w0": pred}, ["A", "B"], C, 20, 3)
        donor_equal = np.mean([r["macro_f1_w0"] for r in rows])
        pooled, _ = score_individual(truth, pred, C[:3], C)
        self.assertAlmostEqual(rows[0]["macro_f1_w0"], 1.0)
        self.assertNotAlmostEqual(donor_equal, pooled, places=3)
        self.assertEqual([r["n_test_cells"] for r in rows], [200, 60])


class TestRandomReference(unittest.TestCase):
    def test_deterministic_bin_matched_panels(self):
        rng = np.random.default_rng(12345)
        means = rng.normal(size=2000)
        ref = rng.choice(2000, size=64, replace=False)
        kw = dict(n_genes=2000, k=64, bins=20, draws=20, seed=20261115, salt=444)
        a, counts = random_reference.random_panels(means, ref, 0, **kw)
        b, _ = random_reference.random_panels(means, ref, 0, **kw)
        c, _ = random_reference.random_panels(means, ref, 1, **kw)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)                              # fold enters the seed
        self.assertEqual(len({tuple(p) for p in a}), 20)
        order = np.lexsort((np.arange(2000), means))
        bins = np.empty(2000, int); bins[order] = np.arange(2000) // 100
        for p in a:
            self.assertEqual(len(set(p)), 64)
            self.assertEqual(np.bincount(bins[p], minlength=20).tolist(), counts)

    def test_invalid_reference_rejected(self):
        with self.assertRaises(ValueError):
            random_reference.random_panels(np.zeros(2000), [1] * 64, 0, n_genes=2000, k=64, bins=20, draws=1, seed=1, salt=1)


class TestDataGuard(unittest.TestCase):
    def test_guard_refuses_forbidden_groups_and_rows(self):
        from compact_panels import data
        from helpers import synthetic_dataset
        info = synthetic_dataset()
        self.assertTrue(data.guard_selftest(info["h5ad"]))
        obs, var, accessed = data.load_metadata(info["h5ad"])
        self.assertEqual({a.split("/")[0] for a in accessed}, {"obs", "var"})
        names = sorted(obs.index[:5])
        with self.assertRaises(data.LeakageGuardError):
            data.read_expression(info["h5ad"], names, obs.index, allowed=set(names[:4]))
        blk = data.read_expression(info["h5ad"], names[::-1], obs.index, allowed=set(names))
        import anndata as ad
        full = ad.read_h5ad(info["h5ad"])
        self.assertTrue(np.array_equal(blk.matrix.toarray(), full[names[::-1]].X.toarray()))  # order preserved


if __name__ == "__main__":
    unittest.main()
