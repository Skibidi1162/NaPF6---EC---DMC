import sys
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import plot_ion_pairing as plotting


class IonPairingPlotTests(unittest.TestCase):
    def test_moving_average_uses_available_neighbors_at_edges(self):
        self.assertEqual(
            plotting.moving_average([0.0, 10.0, 20.0, 30.0, 40.0], 3),
            [5.0, 10.0, 20.0, 30.0, 35.0],
        )

    def test_select_cutoff_rows_uses_exact_requested_cutoff(self):
        rows = [
            {"cutoff_nm": 0.45, "value": 1.0},
            {"cutoff_nm": 0.5, "value": 2.0},
            {"cutoff_nm": 0.5, "value": 3.0},
        ]
        selected = plotting.select_cutoff_rows(rows, 0.5)
        self.assertEqual([row["value"] for row in selected], [2.0, 3.0])

    def test_missing_primary_cutoff_is_actionable(self):
        with self.assertRaisesRegex(plotting.PlotError, "available: 0.4, 0.6"):
            plotting.select_cutoff_rows(
                [{"cutoff_nm": 0.4}, {"cutoff_nm": 0.6}], 0.5
            )

    def test_output_paths_cover_all_plots_and_formats(self):
        paths = plotting.output_paths(Path("analysis"), ["png", "svg"])
        self.assertEqual(len(paths), 2 * len(plotting.PLOT_STEMS))
        self.assertIn(
            Path("analysis/figure_ion_pairing_primary.svg"), paths
        )


if __name__ == "__main__":
    unittest.main()
