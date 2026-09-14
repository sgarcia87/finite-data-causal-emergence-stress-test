import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import tempfile
import unittest
from pathlib import Path

import numpy as np

import triadic_exact_hierarchy as hierarchy
import triadic_blind_discovery as discovery


class TriadicBlindDiscoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        triads = hierarchy.enumerate_triads4()
        graph = hierarchy.adjacency(triads)
        cls.laplacian = 6 * np.eye(48, dtype=int) - graph

    def test_population_gap_selects_four_not_hidden_levels(self):
        for tau in (0.1, 0.4, 1.0):
            p = hierarchy.heat_kernel(self.laplacian, tau)
            values, _ = discovery.spectral_values_vectors(p)
            self.assertEqual(discovery.select_k_by_raw_gap(values), 4)

    def test_selection_stays_within_frozen_range(self):
        p = hierarchy.heat_kernel(self.laplacian, 0.4)
        values, _ = discovery.spectral_values_vectors(p)
        self.assertIn(discovery.select_k_by_raw_gap(values), range(2, 13))

    def test_t_partition_is_cube_parity_not_a_coordinate_cut(self):
        triads = hierarchy.enumerate_triads4()
        family, t_class = hierarchy.triadic_partitions(triads)
        parity = ((family >> 0) & 1) ^ ((family >> 1) & 1) ^ ((family >> 2) & 1)
        self.assertTrue(np.array_equal(t_class, parity))

    def test_small_run_writes_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            result = discovery.run(Path(directory), [500], repeats=2, seed=9)
            self.assertEqual(result["hidden_exact_hierarchy"], [48, 8, 2])
            self.assertEqual(len(result["summary"]), 9)
            for name in (
                "triadic_blind_discovery_runs.csv",
                "triadic_blind_discovery_summary.csv",
                "triadic_blind_discovery_summary.json",
            ):
                self.assertTrue((Path(directory) / name).exists())


if __name__ == "__main__":
    unittest.main()
