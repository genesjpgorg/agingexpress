"""Behavior checks for the scientific scoring and input safeguards."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("evaluate_gi", Path(__file__).resolve().parents[1] / "scripts/evaluate_gi.py")
evaluation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluation)


class EvaluationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.rows = []
        requests, targets = [], []
        for index, effect in enumerate((-1, 0, 1), 1):
            gene = f"ENSG{index:011d}"
            targets.append({"ensembl_gene_id": gene, "gene_symbol": f"G{index}", "scaled_log2_change_age20_to75": effect, "adjusted_p_below_0_05": "0"})
            for condition, value in (("young", 10), ("old", (11 * 2 ** effect) - 1)):
                request_id = f"{gene}__primary_{condition}"
                requests.append({"request_id": request_id, "ensembl_gene_id": gene, "variant": "primary", "model": "g0-expression"})
                self.rows.append({"request_id": request_id, "tpm": str(value), "sequence_sha256": "a" * 64, "resolved_ensembl_gene_id": gene, "model": "g0-expression", "run_id": "synthetic_test_only"})
        for filename, rows in (("gi_requests.tsv", requests), ("primary_targets.tsv", targets)):
            with (self.root / filename).open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

    def tearDown(self):
        self.temp.cleanup()

    def test_known_order_and_rounded_zero_direction(self):
        metrics, pairs = evaluation.score(self.root, self.rows, iterations=20)
        result = metrics["primary_all_available_genes"]
        self.assertAlmostEqual(result["spearman_rho"], 1)
        self.assertEqual(result["direction_n_nonzero_reference"], 2)
        self.assertEqual(result["direction_accuracy_nonzero_reference"], 1)
        self.assertEqual(len(pairs), 3)

    def test_changed_sequence_is_rejected(self):
        self.rows[1]["sequence_sha256"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "Different sequences"):
            evaluation.score(self.root, self.rows, iterations=0)

    def test_failed_prediction_is_not_zero(self):
        self.rows[0]["tpm"] = ""
        with self.assertRaisesRegex(ValueError, "Missing 1 primary"):
            evaluation.score(self.root, self.rows, iterations=0)
        metrics, pairs = evaluation.score(self.root, self.rows, allow_incomplete=True, iterations=0)
        self.assertEqual(metrics["primary_pair_count"], 2)
        self.assertEqual(len(metrics["missing_primary_request_ids"]), 1)
        self.assertEqual(len(pairs), 2)

    def test_constant_age_prediction_is_undefined(self):
        for row in self.rows:
            row["tpm"] = "10"
        metrics, _ = evaluation.score(self.root, self.rows, iterations=20)
        self.assertIsNone(metrics["primary_all_available_genes"]["spearman_rho"])
        self.assertEqual(metrics["primary_all_available_genes"]["direction_accuracy_nonzero_reference"], 0)

    def test_invalid_and_duplicate_rows_are_rejected(self):
        for value in ("-1", "nan", "inf"):
            invalid = [dict(r) for r in self.rows]
            invalid[0]["tpm"] = value
            with self.assertRaisesRegex(ValueError, "Invalid TPM"):
                evaluation.score(self.root, invalid, iterations=0)
        with self.assertRaisesRegex(ValueError, "duplicate request_id"):
            evaluation.score(self.root, self.rows + [self.rows[0]], iterations=0)


if __name__ == "__main__":
    unittest.main()
