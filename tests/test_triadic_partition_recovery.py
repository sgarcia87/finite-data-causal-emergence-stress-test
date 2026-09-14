import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import tempfile
import unittest
from pathlib import Path

import numpy as np
from sklearn.metrics import adjusted_rand_score

from triadic_partition_recovery import (
    discover_partition,
    empirical_tpm,
    make_models,
    population_closure_residual,
    run,
)


class TriadicPartitionRecoveryTests(unittest.TestCase):
    def test_hidden_partition_is_exact_except_broken_control(self):
        truth = np.repeat(np.arange(8), 6)
        models = make_models()
        for name, p in models.items():
            residual = population_closure_residual(p, truth)
            if name == "isospectral_broken_8partition":
                self.assertGreater(residual, 1e-5)
            else:
                self.assertLess(residual, 1e-12)

    def test_noiseless_rank8_case_recovers_hidden_families(self):
        truth = np.repeat(np.arange(8), 6)
        p = make_models()["eta_0"]
        discovered = discover_partition(p, seed=0)
        self.assertEqual(adjusted_rand_score(truth, discovered), 1.0)

    def test_broken_control_is_isospectral(self):
        models = make_models()
        reference = np.linalg.svd(models["eta_0.02"], compute_uv=False)
        broken = np.linalg.svd(models["isospectral_broken_8partition"], compute_uv=False)
        self.assertLess(float(np.max(np.abs(reference - broken))), 1e-12)

    def test_small_run_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run(Path(tmp), samples=[100], repeats=2, random_controls=3, seed=9)
            self.assertEqual(result["independent_repeats"], 2)
            self.assertTrue((Path(tmp) / "triadic_partition_recovery_runs.csv").exists())
            self.assertTrue((Path(tmp) / "triadic_partition_recovery_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
