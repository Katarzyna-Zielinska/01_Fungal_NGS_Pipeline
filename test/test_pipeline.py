import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"


def load_script(name):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


alignment_qc = load_script("alignment_qc")
coverage_summary = load_script("coverage_summary")
depth_summary = load_script("depth_summary")


class AlignmentQCTests(unittest.TestCase):
    def test_extract_metric_with_percentage(self):
        count, percentage = alignment_qc.extract_metric(
            "100 + 2 primary mapped (91.23% : N/A)",
            r"^(\d+) \+ \d+ primary mapped \(([\d.]+)%"
        )
        self.assertEqual(count, 100)
        self.assertAlmostEqual(percentage, 91.23)

    def test_extract_metric_missing_raises(self):
        with self.assertRaises(ValueError):
            alignment_qc.extract_metric(
                "no matching metric",
                r"^(\d+) primary$"
            )


class CoverageSummaryTests(unittest.TestCase):
    def make_coverage_file(self, directory, rows):
        path = Path(directory) / "coverage.tsv"
        pd.DataFrame(rows).to_csv(
            path, sep="\t", header=False, index=False
        )
        return path

    def test_valid_coverage_table(self):
        rows = [
            ["contig1", 0, 100, "gene=locus_tag=GENE1", 10, 80, 100, 0.8],
            ["contig1", 100, 250, "gene=locus_tag=GENE2", 20, 150, 150, 1.0],
        ]

        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_coverage_file(tmp, rows)
            original_expected = coverage_summary.EXPECTED_GENES
            coverage_summary.EXPECTED_GENES = 2
            try:
                df = coverage_summary.load_coverage_file(path, "TEST")
            finally:
                coverage_summary.EXPECTED_GENES = original_expected

        self.assertEqual(len(df), 2)
        self.assertEqual(set(df["gene_id"]), {"GENE1", "GENE2"})

    def test_invalid_coverage_fraction_raises(self):
        rows = [
            ["contig1", 0, 100, "gene=locus_tag=GENE1", 10, 80, 100, 0.5],
        ]

        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_coverage_file(tmp, rows)
            original_expected = coverage_summary.EXPECTED_GENES
            coverage_summary.EXPECTED_GENES = 1
            try:
                with self.assertRaises(ValueError):
                    coverage_summary.load_coverage_file(path, "TEST")
            finally:
                coverage_summary.EXPECTED_GENES = original_expected


class DepthSummaryTests(unittest.TestCase):
    def test_calculate_genome_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_file = Path(tmp) / "genome.txt"
            genome_file.write_text(
                "contig1\t100\ncontig2\t250\n",
                encoding="utf-8"
            )
            contigs, length = depth_summary.calculate_genome_length(
                genome_file
            )

        self.assertEqual(contigs, 2)
        self.assertEqual(length, 350)

    def test_invalid_contig_length_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            genome_file = Path(tmp) / "genome.txt"
            genome_file.write_text("contig1\t0\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                depth_summary.calculate_genome_length(genome_file)


if __name__ == "__main__":
    unittest.main(verbosity=2)
