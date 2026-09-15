import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import triadic_ei_matched_controls as audit


class AuditTests(unittest.TestCase):
    def test_models_are_valid_tpms(self):
        for matrix in audit.make_models().values():
            self.assertTrue(np.all(matrix >= -1e-14))
            self.assertTrue(np.allclose(matrix.sum(axis=1), 1.0))

    def test_designated_partition_closure(self):
        truth = np.repeat(np.arange(8), 6)
        models = audit.make_models()
        for name in ("eta_0", "eta_0.02", "eta_0.2", "eta_0.8"):
            self.assertLess(audit.closure_tv(models[name], truth), 1e-12)
        self.assertGreater(audit.closure_tv(models["isospectral_broken_8partition"], truth), 1e-4)

    def test_balanced_partition(self):
        labels = audit.balanced_partition(np.random.default_rng(1))
        self.assertTrue(np.array_equal(np.bincount(labels), np.repeat(6, 8)))

    def test_null_preserves_margins(self):
        rng = np.random.default_rng(2)
        counts = audit.sample_counts(audit.make_models()["eta_0.2"], 100, rng)
        shuffled = audit.shuffle_futures_exact(counts, np.random.default_rng(3))
        self.assertTrue(np.array_equal(shuffled.sum(axis=0), counts.sum(axis=0)))
        self.assertTrue(np.array_equal(shuffled.sum(axis=1), counts.sum(axis=1)))

    def test_small_run_writes_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            result = audit.run(Path(directory), random_partitions=10, repeats=2, seed=7)
            self.assertEqual(len(result["population"]), 5)
            self.assertTrue((Path(directory) / "triadic_ei_matched_controls_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
