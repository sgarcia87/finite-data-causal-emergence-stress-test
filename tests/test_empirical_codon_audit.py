import tempfile
import unittest
from pathlib import Path

import numpy as np

import empirical_codon_audit as audit


SOURCE = Path("data/empirical_ecm")


class EmpiricalCodonAuditTests(unittest.TestCase):
    def test_source_hashes_and_parse(self):
        for filename, expected in audit.SOURCE_HASHES.items():
            path = SOURCE / filename
            self.assertEqual(audit.sha256(path), expected)
            codons, exchangeability, frequencies = audit.parse_ecm(path)
            self.assertEqual(len(codons), 61)
            self.assertEqual(exchangeability.shape, (61, 61))
            self.assertTrue(np.allclose(exchangeability, exchangeability.T))
            self.assertAlmostEqual(float(frequencies.sum()), 1.0)
            self.assertTrue(all(audit.GENETIC_CODE[c] != "STOP" for c in codons))

    def test_partition_branch(self):
        codons, _, _ = audit.parse_ecm(SOURCE / "ECMrest.dat")
        labels = audit.partitions(codons)
        self.assertEqual({k: int(v.max() + 1) for k, v in labels.items()},
                         {"B16": 16, "AA20": 20, "L23": 23})
        self.assertTrue(audit.refines(labels["L23"], labels["B16"]))
        self.assertTrue(audit.refines(labels["L23"], labels["AA20"]))
        self.assertFalse(audit.refines(labels["B16"], labels["AA20"]))
        self.assertFalse(audit.refines(labels["AA20"], labels["B16"]))

    def test_generator_and_transition_matrices(self):
        for filename in audit.SOURCE_HASHES:
            _, exchangeability, frequencies = audit.parse_ecm(SOURCE / filename)
            q = audit.build_generator(exchangeability, frequencies)
            self.assertLess(float(np.abs(q.sum(axis=1)).max()), 1e-12)
            self.assertAlmostEqual(float(-frequencies @ np.diag(q)), 1.0)
            self.assertLess(float(np.abs(frequencies[:, None] * q - frequencies[None, :] * q.T).max()), 1e-12)
            for horizon in audit.HORIZONS:
                p = audit.transition_matrix(q, horizon)
                self.assertTrue(np.all(p >= 0))
                self.assertTrue(np.allclose(p.sum(axis=1), 1.0))

    def test_small_end_to_end_run(self):
        with tempfile.TemporaryDirectory() as directory:
            result = audit.run(SOURCE, Path(directory), controls=5, repeats=2, seed=7)
            self.assertEqual(len(result["closure"]), 30)
            self.assertEqual(len(result["recovery_summary"]), 60)
            self.assertTrue((Path(directory) / "summary.json").exists())


if __name__ == "__main__":
    unittest.main()
