"""End-to-end pipeline on synthetic data: training-only transformations, isolation of
classifier/evaluation samples, write-once scoring, report regeneration."""
import json
import unittest

import numpy as np

from compact_panels import data, preprocessing
from compact_panels.io_utils import IntegrityError, read_json
from compact_panels.pipeline import Study
from compact_panels.reporting import build_results
from compact_panels.validate import run_validation
from helpers import synth_cfg, synthetic_paths


class TestEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = synth_cfg()
        cls.paths = synthetic_paths("e2e")
        cls.st = Study(cls.cfg, cls.paths, log=lambda m: None)
        cls.validation = run_validation(cls.cfg, cls.paths, refs_path=None, acceptance=None, deep=True)
        for k in range(cls.cfg.n_folds):
            while not cls.st.run_panels(k, 50):
                pass
            cls.st.run_classifiers(k)
            cls.st.run_random(k)
            cls.st.run_score(k)
        cls.st.run_g3(0, "w0")
        cls.st.run_g4(50)

    def test_validation_passes_without_fitting(self):
        self.assertEqual(self.validation["status"], "PASS", [c for c in self.validation["checks"] if not c["passed"]])
        self.assertEqual(self.validation["models_fitted"], 0)

    def test_preprocessing_uses_fitting_cells_only(self):
        k = 1
        plan = self.st.plan(k)
        blk = self.st.read_rows(plan.fitting, allowed=set(plan.fitting))
        hv = self.cfg["hvg"]
        pre, F = preprocessing.fit_preprocessing(blk.matrix, blk.var_names, int(hv["n_top_genes"]), hv["flavor"], "cv")
        with np.load(self.st.fold_dir(k) / "preprocessing.npz") as z:
            self.assertTrue(np.array_equal(z["means"], pre.means) and np.array_equal(z["sds"], pre.sds))
        # adding held-out rows changes the statistics, so equality above proves they were not used
        both = self.st.read_rows(plan.fitting + plan.test, allowed=set(plan.fitting) | set(plan.test))
        mu_all = both.matrix[:, pre.candidate_var_idx].toarray().astype(np.float64).mean(axis=0)
        self.assertFalse(np.allclose(mu_all, pre.means, rtol=0, atol=1e-12))

    def test_training_stage_cannot_read_held_out_rows(self):
        plan = self.st.plan(2)
        with self.assertRaises(data.LeakageGuardError):
            self.st.read_rows(plan.test[:3], allowed=set(plan.fitting))

    def test_classifier_and_evaluation_samples_are_separate(self):
        for k in range(self.cfg.n_folds):
            plan = self.st.plan(k)
            self.assertTrue(set(plan.classifier) <= set(plan.fitting))
            self.assertFalse(set(plan.classifier) & set(plan.test))
            self.assertFalse(set(plan.selector) & set(plan.test))
            self.assertFalse(set(plan.withheld) & set(plan.classifier))
            ev = read_json(self.st.fold_dir(k) / "evaluation.json")
            self.assertEqual({r["individual_id"] for r in ev["per_individual"]}, set(plan.test_people))
            self.assertEqual(sum(r["n_test_cells"] for r in ev["per_individual"]), len(plan.test))

    def test_panels_have_k_genes_and_traces(self):
        for k in range(self.cfg.n_folds):
            for arm in ("w0", "w1"):
                p = read_json(self.st.fold_dir(k) / f"panel_{arm}.json")
                self.assertEqual(len(set(p["panel_idx"])), self.cfg.k)
                self.assertTrue(p["all_converged"])
                self.assertEqual(p["variant"], "fold0_recovery" if k == 0 else "cv_halving")

    def test_scoring_is_write_once(self):
        with self.assertRaises(IntegrityError):
            self.st.run_score(0)

    def test_training_refused_after_scoring(self):
        with self.assertRaises(IntegrityError):
            self.st.run_panels(0, 1)

    def test_report_regeneration_is_deterministic(self):
        root = self.paths.run_dir
        a = build_results(root, root / "results_a", self.cfg)
        b = build_results(root, root / "results_b", self.cfg)
        self.assertEqual((root / "results_a" / "SHA256SUMS.txt").read_bytes(), (root / "results_b" / "SHA256SUMS.txt").read_bytes())
        res = json.loads((root / "results_a" / "results.json").read_text())
        s = res["summary"]
        per = [json.loads((root / f"fold{k}" / "evaluation.json").read_text()) for k in range(self.cfg.n_folds)]
        rows = [r for e in per for r in e["per_individual"]]
        self.assertAlmostEqual(s["paired_delta"], float(np.mean([r["paired_delta"] for r in rows])), places=14)
        self.assertEqual(s["n_individuals"], len(rows))
        self.assertIn("g3_fold0_w0", res["stability"]["checks"])
        self.assertIn("not confirmatory", res["stability"]["label"])
        self.assertIn("g4", res["stability"]["checks"])
        table = (root / "results_a" / "table_01_individuals_anonymized.csv").read_text()
        self.assertNotIn("IND", table)                                    # source identifiers anonymized


if __name__ == "__main__":
    unittest.main()


class TestReportRendering(unittest.TestCase):
    def test_template_renders_from_results_only(self):
        from compact_panels.reporting import render_report
        TestEndToEnd.setUpClass() if not hasattr(TestEndToEnd, "st") else None
        root = TestEndToEnd.paths.run_dir
        out = root / "results_render"
        if not out.exists():
            build_results(root, out, TestEndToEnd.cfg)
        results = json.loads((out / "results.json").read_text())
        from helpers import REPO
        text = render_report((REPO / "report" / "report_template.md").read_text(), results, None, "0" * 64, "2026-01-01")
        self.assertNotIn("{{", text)
        self.assertIn(f"{results['summary']['n_individuals']} individuals", text)
