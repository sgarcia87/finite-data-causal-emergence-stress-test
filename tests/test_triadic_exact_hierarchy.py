import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import tempfile
import unittest
from pathlib import Path

import numpy as np

import triadic_exact_hierarchy as hierarchy


class TriadicExactHierarchyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.triads = hierarchy.enumerate_triads4()
        cls.graph = hierarchy.adjacency(cls.triads)
        cls.laplacian = 6 * np.eye(48, dtype=int) - cls.graph
        cls.family, cls.t_class = hierarchy.triadic_partitions(cls.triads)

    def test_graph_and_partition_sizes(self):
        self.assertEqual(len(self.triads), 48)
        self.assertEqual(int(self.graph.sum() // 2), 144)
        self.assertEqual(np.bincount(self.family).tolist(), [6] * 8)
        self.assertEqual(np.bincount(self.t_class).tolist(), [24, 24])

    def test_exact_nested_quotients(self):
        q8 = hierarchy.exact_generator_quotient(self.laplacian, self.family)
        q2 = hierarchy.exact_generator_quotient(self.laplacian, self.t_class)
        self.assertTrue(all(len(set(self.t_class[self.family == k])) == 1 for k in range(8)))
        self.assertTrue(np.allclose(np.linalg.eigvalsh(q8), [0, 2, 2, 2, 4, 4, 4, 6]))
        self.assertTrue(np.array_equal(q2, [[3, -3], [-3, 3]]))

    def test_full_rank_with_exact_closure(self):
        for tau in (0.1, 0.4, 1.0):
            p = hierarchy.heat_kernel(self.laplacian, tau)
            self.assertEqual(np.linalg.matrix_rank(p, tol=1e-13), 48)
            self.assertLess(hierarchy.quotient_and_closure(p, self.family)[1], 1e-12)
            self.assertLess(hierarchy.quotient_and_closure(p, self.t_class)[1], 1e-12)

    def test_small_run_writes_all_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            result = hierarchy.run(Path(directory), [100], repeats=2, splits=30, seed=7)
            self.assertEqual(result["exact_nested_hierarchy"], [48, 8, 2])
            for name in (
                "triadic_hierarchy_population.csv",
                "triadic_resolved_modes_runs.csv",
                "triadic_resolved_modes_summary.csv",
                "triadic_hierarchy_summary.json",
            ):
                self.assertTrue((Path(directory) / name).exists())


if __name__ == "__main__":
    unittest.main()
