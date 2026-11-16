"""Donor isolation, split determinism, eligibility and Stage B withholding."""
import math
import unittest

import numpy as np
import pandas as pd

from compact_panels import identity, partition
from helpers import synth_cfg, toy_obs


class TestIdentity(unittest.TestCase):
    def test_unmapped_patient_raises(self):
        obs = pd.DataFrame({"patient_id": ["A", "B"]}, index=["c1", "c2"])
        mapping = pd.DataFrame({"h5ad_patient_id": ["A"], "individual_id": ["I1"]})
        with self.assertRaises(KeyError):
            identity.map_individuals(obs, mapping)

    def test_multi_id_individual_collapses(self):
        obs = pd.DataFrame({"patient_id": ["A", "A2", "B"]}, index=["c1", "c2", "c3"])
        mapping = pd.DataFrame({"h5ad_patient_id": ["A", "A2", "B"], "individual_id": ["I1", "I1", "I2"]})
        self.assertEqual(identity.map_individuals(obs, mapping).tolist(), ["I1", "I1", "I2"])

    def test_exclusion_requires_listed_cells_in_pool(self):
        pool = pd.DataFrame({"x": [1, 2]}, index=["a", "b"])
        self.assertEqual(identity.apply_exclusion(pool, ["a"]).index.tolist(), ["b"])
        with self.assertRaises(ValueError):
            identity.apply_exclusion(pool, ["zzz"])

    def test_eligibility_counts_and_threshold(self):
        rows = []
        for ind, counts in {"I1": (20, 20, 20), "I2": (20, 19, 25), "I3": (5, 30, 30)}.items():
            for c, n in zip("ABC", counts):
                rows += [{"individual_id": ind, "author_cell_type": c}] * n
        pool = pd.DataFrame(rows)
        counts, k, E = identity.eligibility(pool, ["A", "B", "C"], 20, 3)
        self.assertEqual(k.to_dict(), {"I1": 3, "I2": 2, "I3": 2})
        self.assertEqual(E, ["I1"])        # exactly 20 cells counts; 19 does not
        counts, k, E = identity.eligibility(pool, ["A", "B", "C"], 20, 2)
        self.assertEqual(E, ["I1", "I2", "I3"])


class TestFolds(unittest.TestCase):
    def test_fold_determinism_and_disjointness(self):
        inds = [f"I{i:02d}" for i in range(24)]
        sites = [("S1", "S2", "S3")[i % 3] for i in range(24)]
        a = partition.make_folds(inds, sites, 4, True, 20261112)
        b = partition.make_folds(inds, sites, 4, True, 20261112)
        tests = [sorted(inds[j] for j in te) for _, te in a]
        self.assertEqual(tests, [sorted(inds[j] for j in te) for _, te in b])
        flat = [x for t in tests for x in t]
        self.assertEqual(sorted(flat), inds)               # each individual in exactly one test fold
        for tr, te in a:
            self.assertFalse(set(tr) & set(te))
        for t in tests:                                    # site stratification: 2 per site per fold
            self.assertEqual(sorted(sites[inds.index(x)] for x in t), ["S1", "S1", "S2", "S2", "S3", "S3"])
        c = partition.make_folds(inds, sites, 4, True, 20261113)
        self.assertNotEqual(tests, [sorted(inds[j] for j in te) for _, te in c])

    def test_unsorted_individuals_rejected(self):
        with self.assertRaises(ValueError):
            partition.make_folds(["B", "A"], ["S", "S"], 2, True, 1)

    def test_canonical_hash_stable(self):
        s1, h1 = partition.canonical_assignment([["b", "a"], ["c"]], 7, 2)
        s2, h2 = partition.canonical_assignment([["a", "b"], ["c"]], 7, 2)
        self.assertEqual(h1, h2)
        self.assertIn('"test_folds":[["a","b"],["c"]]', s1)


class TestStageB(unittest.TestCase):
    def test_withholding_rule_integer_form(self):
        for n in range(0, 100001):
            self.assertEqual(partition.n_withheld(n), (2 * n + 5) // 10)
        for n in range(0, 2001):
            self.assertEqual(partition.n_withheld(n), math.floor(0.2 * n + 0.5))
        self.assertEqual([partition.n_withheld(n) for n in (1, 2, 3, 8, 13, 23)], [0, 0, 1, 2, 3, 5])

    def test_withhold_deterministic_and_within_training(self):
        obs = toy_obs()
        w1, rows = partition.withhold(obs, ["A", "B", "C"], [20261112, 1])
        w2, _ = partition.withhold(obs, ["A", "B", "C"], [20261112, 1])
        w3, _ = partition.withhold(obs, ["A", "B", "C"], [20261112, 2])
        self.assertEqual(w1, w2)
        self.assertNotEqual(w1, w3)
        self.assertTrue(set(w1) <= set(obs.index))
        for r in rows:
            sub = obs[(obs.individual_id == r["individual_id"]) & (obs.author_cell_type == r["class"])]
            self.assertEqual(int(sub.index.isin(w1).sum()), partition.n_withheld(len(sub)))

    def test_matches_independent_reference_including_small_strata(self):
        # Reference written directly from the frozen rule; strata with n = 1, 2 (m = 0) included.
        obs = toy_obs(n_ind=4, per=10)
        extra = pd.DataFrame({"individual_id": ["I00", "I01", "I01"], "author_cell_type": ["D", "D", "D"],
                              "Site": ["S1"] * 3, "patient_id": ["P00", "P01", "P01"], "Status": ["Healthy"] * 3},
                             index=["zz_a", "zz_b", "zz_c"])
        cells = pd.concat([obs, extra])
        classes = ["D", "A", "B", "C"]
        got, _ = partition.withhold(cells, classes, [20261112, 3])
        rng = np.random.default_rng([20261112, 3])
        ref = []
        for d in sorted(cells.individual_id.unique()):
            for c in classes:
                names = sorted(cells.index[(cells.individual_id == d) & (cells.author_cell_type == c)])
                if names:
                    idx = rng.choice(len(names), size=(2 * len(names) + 5) // 10, replace=False)
                    ref += [names[i] for i in idx]
        self.assertEqual(got, sorted(ref))
        self.assertFalse({"zz_a", "zz_b", "zz_c"} & set(got))

    def test_subsample_cap_and_independent_streams(self):
        obs = toy_obs(per=60)
        sel, _ = partition.draw_subsample(obs, ["A", "B", "C"], 50, [7, 0, 1])
        clf, _ = partition.draw_subsample(obs, ["A", "B", "C"], 200, [7, 0, 2])
        self.assertTrue(set(sel) <= set(obs.index) and set(clf) <= set(obs.index))
        counts = obs.loc[sel].groupby(["individual_id", "author_cell_type"]).size()
        self.assertEqual(int(counts.max()), 50)
        self.assertEqual(len(clf), len(obs))               # cap above every stratum size keeps all cells
        self.assertEqual(sel, sorted(sel))


class TestFoldPlanOnSynthetic(unittest.TestCase):
    def test_plan_isolation(self):
        from compact_panels.pipeline import Study
        from helpers import synthetic_paths
        cfg = synth_cfg()
        st = Study(cfg, synthetic_paths("plan_test"))
        co = st.cohort()
        self.assertNotIn("IND9LPS", co.eligible)                       # excluded-status individual removed
        self.assertTrue(all(not i.endswith("07") for i in co.eligible))  # ineligible individuals removed
        asg = st.assignment()
        for k in range(cfg.n_folds):
            p = st.plan(k)
            self.assertFalse(set(p.train_people) & set(p.test_people))
            self.assertFalse(set(p.train) & set(p.test))
            self.assertEqual(set(p.fitting) | set(p.withheld), set(p.train))
            self.assertTrue(set(p.selector) <= set(p.fitting) and set(p.classifier) <= set(p.fitting))
            pid_tr = set(co.obs.loc[p.train].patient_id)
            pid_te = set(co.obs.loc[p.test].patient_id)
            self.assertFalse(pid_tr & pid_te)                          # multi-ID individuals never split
            self.assertEqual(p.seeds["stage_b"], [cfg.seed, k])
            self.assertEqual(p.seeds["selector"], [cfg.seed, k, 1])
            self.assertEqual(p.seeds["classifier"], [cfg.seed, k, 2])
        st2 = Study(cfg, synthetic_paths("plan_test2"))
        self.assertEqual(st2.assignment().canonical_sha256, asg.canonical_sha256)
        self.assertEqual(st2.plan(1).provenance_sha(), st.plan(1).provenance_sha())


if __name__ == "__main__":
    unittest.main()


class TestDerivedInputs(unittest.TestCase):
    def test_mapping_from_sdrf_and_pair_hash(self):
        import tempfile, os
        sdrf = ("Source Name\tCharacteristics[individual]\tExtra\n"
                "BGCV01_CV0001\tCV0001\tx\nBGCV02_CV0002\tCV0002\tx\nnewcastle7\tC-1\tx\nnewcastle7v2\tC-1\tx\nnewcastle7\tC-1\ty\n")
        with tempfile.NamedTemporaryFile("w", suffix=".sdrf.txt", delete=False) as f:
            f.write(sdrf)
        try:
            m = identity.mapping_from_sdrf(f.name, ["CV0001", "newcastle7", "newcastle7v2", "CV0002", "CV0001"])
            self.assertEqual(m.values.tolist(), [["CV0001", "CV0001"], ["CV0002", "CV0002"], ["newcastle7", "C-1"],
                                                 ["newcastle7v2", "C-1"]])
            shuffled = m.sample(frac=1, random_state=0).assign(extra="ignored")
            self.assertEqual(identity.mapping_pairs_sha(m), identity.mapping_pairs_sha(shuffled))
            with self.assertRaises(KeyError):
                identity.mapping_from_sdrf(f.name, ["UNKNOWN"])
        finally:
            os.unlink(f.name)

    def test_qc_rule(self):
        obs = pd.DataFrame({"Status": ["LPS", "Covid", "Covid", "Covid"], "patient_id": ["A", "A", "A", "B"],
                            "Collection_Day": ["D28", "D28", "D0", "D28"]}, index=["c1", "c2", "c3", "c4"])
        self.assertEqual(identity.qc_exclusion_from_rule(obs, {"status_is_not": "LPS", "patient_id": "A",
                                                               "Collection_Day": "D28"}), ["c2"])
