"""Optional: replay the recorded historical alpha traces through the state machine.

Needs read access to the historical artifacts (set CGP_LEGACY_ROOT); skipped otherwise,
so CI runs without any private data. No model is fitted: the recorded (alpha, nnz) pairs
are fed to :func:`replay_trace` and the derived alpha sequence must be bit-identical.
"""
import json
import os
import unittest
from pathlib import Path

import numpy as np

from compact_panels.selection import replay_trace
from helpers import REPO, study_cfg

LEGACY = os.environ.get("CGP_LEGACY_ROOT")


@unittest.skipUnless(LEGACY and Path(LEGACY).is_dir(), "historical artifacts not available (set CGP_LEGACY_ROOT)")
class TestHistoricalAlphaTraces(unittest.TestCase):
    def test_all_eight_traces_replay_exactly(self):
        from compact_panels.reference import Legacy
        cfg = study_cfg()
        legacy = Legacy(Path(LEGACY), REPO / "configs" / "historical_references.toml")
        for k in range(4):
            conv = cfg.convention(k)
            for arm in ("w0", "w1"):
                ck = legacy.checkpoint(k, arm)
                meta, panel = ck["meta"], legacy.panel(k, arm)
                trace = [(f["alpha"], f["n_nonzero"]) for f in meta["fits"]]
                derived, out, state = replay_trace(conv["search_variant"], arm, meta["alpha_max"], trace,
                                                   cfg["selector"], conv)
                with self.subTest(fold=k, arm=arm):
                    self.assertEqual(derived, [t[0] for t in trace])          # bit-identical alpha sequence
                    self.assertIsNotNone(out)
                    self.assertEqual(len(derived), panel["n_fits"])
                    self.assertEqual(out.exact, panel["exact64"])
                    self.assertEqual(out.discarded, panel["discarded"])
                    self.assertEqual(out.alpha_final, panel["alpha_final"])
                    if not out.exact:                                        # fallback = top-64 of stored best norms
                        from compact_panels.selection import top_k_by_norm
                        self.assertEqual(top_k_by_norm(ck["best_norms"], 64).tolist(), sorted(panel["panel_idx"]))


if __name__ == "__main__":
    unittest.main()
