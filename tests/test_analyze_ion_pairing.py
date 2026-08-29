import tempfile
from pathlib import Path
import sys
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import analyze_ion_pairing as analysis


class IonPairingTests(unittest.TestCase):
    def test_classification_rule_matches_reference_method(self):
        self.assertEqual(analysis.classify_neighbor_count(0), "SIP/SSIP")
        self.assertEqual(analysis.classify_neighbor_count(1), "CIP")
        self.assertEqual(analysis.classify_neighbor_count(2), "AGG")
        self.assertEqual(analysis.classify_neighbor_count(5), "AGG")

    def test_neighbor_population_analysis(self):
        na_atoms = [1, 2]
        cutoffs = [0.5]
        columns = [(0.5, 1), (0.5, 2)]
        rows = [
            [0.0, 0.0, 1.0],
            [1.0, 2.0, 1.0],
        ]
        frame_rows, summaries, histogram_rows = analysis.analyze_neighbor_rows(
            rows, columns, cutoffs, na_atoms, requested_blocks=2
        )
        self.assertEqual(len(frame_rows), 2)
        summary = summaries[0]
        self.assertEqual(summary["na_environments"], 4)
        self.assertAlmostEqual(summary["sip_ssip_fraction"], 0.25)
        self.assertAlmostEqual(summary["cip_fraction"], 0.50)
        self.assertAlmostEqual(summary["agg_fraction"], 0.25)
        self.assertEqual(sum(row[2] for row in histogram_rows), 4)

    def test_first_peak_and_minimum_detection(self):
        combined = []
        for index in range(61):
            radius = index * 0.01
            peak = 8.0 * pow(2.718281828, -((radius - 0.34) / 0.035) ** 2)
            valley = -0.4 * pow(2.718281828, -((radius - 0.45) / 0.025) ** 2)
            combined.append([radius, 1.0 + peak + valley, radius + 0.01, radius])
        result = analysis.find_first_peak_and_minimum(combined)
        self.assertAlmostEqual(result["peak_r_nm"], 0.34, delta=0.02)
        self.assertAlmostEqual(result["minimum_r_nm"], 0.45, delta=0.03)

    def test_xvg_parser_ignores_metadata(self):
        with tempfile.NamedTemporaryFile(suffix=".xvg", delete=False) as temporary:
            path = Path(temporary.name)
        try:
            path.write_text(
                "@ title \"test\"\n# comment\n0.0 1.0\n0.1 2.0\n",
                encoding="utf-8",
            )
            self.assertEqual(
                analysis.parse_xvg(path), [[0.0, 1.0], [0.1, 2.0]]
            )
        finally:
            path.unlink(missing_ok=True)

    def test_cutoff_deduplication(self):
        self.assertEqual(
            analysis.unique_cutoffs([0.5, 0.405, 0.5, 0.4, 0.45]),
            [0.4, 0.405, 0.45, 0.5],
        )

    def test_default_sensitivity_cutoffs_include_ec_dmc_references(self):
        self.assertEqual(
            analysis.DEFAULT_SENSITIVITY_CUTOFFS_NM,
            (0.35, 0.4, 0.405, 0.45, 0.6),
        )

    def test_streaming_neighbor_analysis_matches_in_memory_result(self):
        na_atoms = [1, 2]
        cutoffs = [0.5]
        columns = [(0.5, 1), (0.5, 2)]
        rows = [[0.0, 0.0, 1.0], [1.0, 2.0, 1.0]]
        _, expected_summaries, expected_histogram = analysis.analyze_neighbor_rows(
            rows, columns, cutoffs, na_atoms, requested_blocks=2
        )
        with tempfile.NamedTemporaryFile(suffix=".xvg", delete=False) as source:
            source_path = Path(source.name)
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as output:
            output_path = Path(output.name)
        try:
            source_path.write_text(
                "# neighbor counts\n0.0 0.0 1.0\n1.0 2.0 1.0\n",
                encoding="utf-8",
            )
            summaries, histogram, frames = analysis.analyze_neighbor_file(
                source_path,
                output_path,
                columns,
                cutoffs,
                na_atoms,
                requested_blocks=2,
            )
            self.assertEqual(frames, 2)
            self.assertEqual(histogram, expected_histogram)
            self.assertEqual(summaries, expected_summaries)
            self.assertEqual(len(output_path.read_text().splitlines()), 3)
        finally:
            source_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)

    def test_compact_intermediate_list_never_includes_canonical_csv(self):
        names = {
            path.name for path in analysis.analysis_intermediate_paths(Path("analysis"))
        }
        self.assertIn("na_pf6_neighbor_counts.xvg", names)
        self.assertNotIn("ion_pairing_summary.csv", names)
        self.assertNotIn("rdf_coordination_Na_P.csv", names)


if __name__ == "__main__":
    unittest.main()
