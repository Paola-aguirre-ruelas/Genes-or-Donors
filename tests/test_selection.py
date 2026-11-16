"""Targets, alpha search (termination, fallback, warm starts, resume) and fold conventions."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compact_panels import selection, targets
from compact_panels.selection import SearchStop, advance, initial_state, replay_trace
from helpers import study_cfg, synth_cfg

SEL = {"k": 64, "max_steps": 25, "extra_steps": 25, "fallback_trigger": 8, "lower_frac": 1e-4, "first_alpha_frac": 0.1}
CONV0 = {"search_variant": "fold0_recovery", "w0_lower_frac": 0.01}


def run_nnz(variant, arm, nnz_fn, amax=1.0, steps=60):
    """Drive the state machine with a synthetic monotone nnz(alpha) curve."""
    state = initial_state(variant, arm, amax, SEL, CONV0)
    alphas = []
    for step in range(1, steps + 1):
        a = state["next_alpha"]
        alphas.append(a)
        nz = nnz_fn(a)
        norms = np.zeros(2000)
        norms[:nz] = np.linspace(2, 1, nz) if nz else 0
        state, _, out = advance(state, alpha=a, nz=nz, norms=norms, best_norms=np.zeros(2000), step=step, sel_cfg=SEL)
        if out is not None:
            return alphas, out, state
    return alphas, None, state


class TestTargets(unittest.TestCase):
    def test_blocks_have_unit_total_variance(self):
        rng = np.random.default_rng(0)
        cls = rng.choice(list("ABCDE"), size=500)
        ind = rng.choice([f"I{i}" for i in range(10)], size=500)
        y0 = targets.targets_for_arm("w0", cls, ind, sorted(set(ind)), list("ABCDE"))
        y1 = targets.targets_for_arm("w1", cls, ind, sorted(set(ind)), list("ABCDE"))
        self.assertEqual(y0.shape, (500, 5))
        self.assertEqual(y1.shape, (500, 15))
        self.assertTrue(np.isclose(np.sum(y1[:, :5] ** 2), 500) and np.isclose(np.sum(y1[:, 5:] ** 2), 500))
        self.assertTrue(np.allclose(y0, y1[:, :5]))
        self.assertTrue(np.allclose(y1.mean(axis=0), 0))
        with self.assertRaises(ValueError):
            targets.block(["A", "Z"], ["A", "B"])


class TestHelpers(unittest.TestCase):
    def test_top_k_tie_break_by_lower_index(self):
        norms = np.array([1.0, 3.0, 3.0, 0.5, 3.0])
        self.assertEqual(selection.top_k_by_norm(norms, 2).tolist(), [1, 2])

    def test_alpha_max_zeroes_all_coefficients(self):
        from sklearn.linear_model import MultiTaskLasso
        rng = np.random.default_rng(1)
        X = rng.normal(size=(80, 12)); Y = rng.normal(size=(80, 3))
        am = selection.alpha_max(X, Y)
        self.assertTrue(np.all(MultiTaskLasso(alpha=am * 1.0001).fit(X, Y).coef_ == 0))
        self.assertTrue(np.any(MultiTaskLasso(alpha=am * 0.9).fit(X, Y).coef_ != 0))


class TestSearchStateMachine(unittest.TestCase):
    def test_initial_brackets_per_variant(self):
        s = initial_state("fold0_recovery", "w0", 2.0, SEL, CONV0)
        self.assertEqual((s["lo"], s["hi"], s["next_alpha"], s["skip_halving"]), (0.02, 2.0, 0.2, True))
        s = initial_state("fold0_recovery", "w1", 2.0, SEL, CONV0)
        self.assertEqual((s["lo"], s["next_alpha"], s["skip_halving"]), (2.0e-4, 0.2, False))
        s = initial_state("cv_halving", "w0", 2.0, SEL, CONV0)
        self.assertEqual((s["lo"], s["skip_halving"]), (2.0e-4, False))

    def test_exact_target_terminates(self):
        alphas, out, _ = run_nnz("cv_halving", "w0", lambda a: 64)
        self.assertEqual(len(alphas), 1)
        self.assertTrue(out.exact and out.discarded == 0 and len(out.panel_idx) == 64)

    def test_halving_then_bisection(self):
        f = lambda a: 64 if 0.031 < a < 0.032 else (100 if a <= 0.031 else 30)
        alphas, out, _ = run_nnz("cv_halving", "w1", f)
        self.assertEqual(alphas[:3], [0.1, 0.05, 0.025])                 # halving until > 64
        self.assertAlmostEqual(alphas[3], np.sqrt(0.025 * 0.05))         # then geometric bisection
        self.assertTrue(out.exact)

    def test_fold0_w0_bisects_from_second_fit(self):
        f = lambda a: 64 if 0.031 < a < 0.032 else (100 if a <= 0.031 else 30)
        alphas, out, _ = run_nnz("fold0_recovery", "w0", f)
        self.assertAlmostEqual(alphas[1], np.sqrt(0.01 * 0.1))           # no halving for fold-0 w0
        self.assertTrue(out.exact)

    def test_fallback_at_step_25_when_few_discarded(self):
        # nnz jumps from 67 to 61: exact 64 never reached; best alpha above target has 67 genes.
        f = lambda a: 67 if a <= 0.0494 else 61
        alphas, out, state = run_nnz("cv_halving", "w1", f)
        self.assertEqual(len(alphas), 25)
        self.assertFalse(out.exact)
        self.assertEqual(out.discarded, 3)
        self.assertEqual(out.alpha_final, state["best_alpha"])
        self.assertEqual(out.panel_idx, list(range(64)))                 # 64 largest norms of the best fit

    def test_more_than_8_discarded_continues_to_50(self):
        f = lambda a: 80 if a <= 0.0494 else 50
        alphas, out, _ = run_nnz("cv_halving", "w1", f)
        self.assertEqual(len(alphas), 50)
        self.assertEqual(out.discarded, 16)
        alphas0, out0, _ = run_nnz("fold0_recovery", "w1", f)
        self.assertEqual(len(alphas0), 50)

    def test_no_fit_above_target_differs_by_variant(self):
        f = lambda a: 10                                                  # never above 64
        with self.assertRaises(SearchStop):                               # cannot halve below the lower bracket
            run_nnz("cv_halving", "w1", f, steps=60)
        state = initial_state("fold0_recovery", "w0", 1.0, SEL, CONV0)
        with self.assertRaises(SearchStop):                               # fold-0 stops at fit 25 without fallback
            for step in range(1, 26):
                a = state["next_alpha"]
                state, _, _ = advance(state, alpha=a, nz=10, norms=None, best_norms=None, step=step, sel_cfg=SEL)

    def test_replay_trace_round_trip(self):
        f = lambda a: 64 if 0.031 < a < 0.032 else (100 if a <= 0.031 else 30)
        alphas, out, _ = run_nnz("cv_halving", "w1", f)
        trace = [(a, f(a)) for a in alphas]
        derived, out2, _ = replay_trace("cv_halving", "w1", 1.0, trace, SEL, CONV0)
        self.assertEqual(derived, alphas)
        self.assertTrue(out2.exact)


class TestSelectorWarmStart(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(3)
        self.X = np.asfortranarray(rng.normal(size=(200, 40)))
        B = np.zeros((40, 3)); B[:6] = rng.normal(size=(6, 3))
        self.Y = self.X @ B + 0.1 * rng.normal(size=(200, 3))

    def test_zero_iteration_refit_is_judged_by_gap(self):
        s = selection.Selector(self.X, self.Y, tol=1e-4, max_iter=1000, retry_max_iter=10000, label="t", log=lambda m: None)
        am = selection.alpha_max(self.X, self.Y)
        s.fit_at(0.3 * am, 1, "t")
        s.fit_at(0.3 * am, 2, "t")                                        # identical alpha, warm start
        rec = s.fits[-1]
        self.assertEqual(rec["n_iter"], 0)                                # pre-loop gap check returns immediately
        self.assertTrue(rec["converged"])                                 # accepted because gap <= threshold
        self.assertEqual(s.fits[0]["n_nonzero"], rec["n_nonzero"])

    def test_retry_restores_starting_coefficients(self):
        s = selection.Selector(self.X, self.Y, tol=1e-12, max_iter=2, retry_max_iter=5, label="t", log=lambda m: None)
        am = selection.alpha_max(self.X, self.Y)
        s.fit_at(0.1 * am, 1, "t")
        rec = s.fits[-1]
        self.assertTrue(rec["retried"])
        self.assertFalse(rec["converged"])
        self.assertEqual(rec["first_attempt"]["max_iter"], 2)
        self.assertEqual(rec["max_iter"], 5)
        self.assertEqual(s.est.max_iter, 2)                               # restored to the default budget


class TestResumeEquivalence(unittest.TestCase):
    def test_interrupted_search_equals_uninterrupted(self):
        cfg = synth_cfg()
        rng = np.random.default_rng(5)
        X = np.asfortranarray(rng.normal(size=(300, 100)))
        B = np.zeros((100, 5)); B[:20] = rng.normal(size=(20, 5))
        Y = X @ B + rng.normal(size=(300, 5))
        genes = np.array([f"G{i}" for i in range(100)])
        quiet = lambda m: None
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            selection.run_search(X, Y, arm="w0", fold=1, genes=genes, cfg=cfg, dest=Path(d1), max_new_fits=50,
                                 provenance={"p": 1}, log=quiet)
            while not (Path(d2) / "panel_w0.json").exists():
                selection.run_search(X, Y, arm="w0", fold=1, genes=genes, cfg=cfg, dest=Path(d2), max_new_fits=1,
                                     provenance={"p": 1}, log=quiet)
            p1 = json.loads((Path(d1) / "panel_w0.json").read_text())
            p2 = json.loads((Path(d2) / "panel_w0.json").read_text())
            self.assertEqual(p1["panel_idx"], p2["panel_idx"])
            self.assertEqual(len(p1["panel_idx"]), cfg.k)
            with np.load(Path(d1) / "checkpoint_w0.npz") as a, np.load(Path(d2) / "checkpoint_w0.npz") as b:
                self.assertTrue(np.array_equal(a["coef"], b["coef"]))     # bit-identical warm-start path
                fa = [f["alpha"] for f in json.loads(str(a["meta"].item()))["fits"]]
                fb = [f["alpha"] for f in json.loads(str(b["meta"].item()))["fits"]]
                self.assertEqual(fa, fb)
            with self.assertRaises(Exception):                            # incompatible provenance refused
                selection.run_search(X, Y, arm="w0", fold=1, genes=genes, cfg=cfg, dest=Path(d1), max_new_fits=1,
                                     provenance={"p": 2}, log=quiet)


class TestHistoricalConventions(unittest.TestCase):
    def test_study_configuration_encodes_fold_specific_paths(self):
        cfg = study_cfg()
        c0, c1, c3 = cfg.convention(0), cfg.convention(1), cfg.convention(3)
        self.assertEqual((c0["search_variant"], c0["hvg_convention"], c0["w0_lower_frac"]), ("fold0_recovery", "pilot", 0.01))
        self.assertEqual((c1["search_variant"], c1["hvg_convention"]), ("cv_halving", "cv"))
        self.assertEqual(c3, c1)
        s = cfg["selector"]
        self.assertEqual((s["k"], s["tol"], s["max_iter"], s["retry_max_iter"], s["max_steps"], s["extra_steps"],
                          s["fallback_trigger"], s["first_alpha_frac"], s["lower_frac"]),
                         (64, 1e-4, 1000, 10000, 25, 25, 8, 0.1, 1e-4))
        self.assertEqual((cfg.seed, cfg["subsample"]["selector_cap"], cfg["subsample"]["classifier_cap"]), (20261112, 50, 200))
        self.assertEqual((cfg["random_reference"]["seed"], cfg["random_reference"]["salt"]), (20261115, 444))
        self.assertEqual(cfg["stability"]["g4_seed"], 20261115)


if __name__ == "__main__":
    unittest.main()
