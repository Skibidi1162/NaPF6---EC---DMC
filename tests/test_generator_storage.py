import sys
from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import electrolyte_generator_v2 as generator


class GeneratorStorageTests(unittest.TestCase):
    def test_workflow_disables_gromacs_backup_accumulation(self):
        case = generator.from_counts(2, 3, 1, 1.2, "test")
        script = generator.workflow_script(case, "storage_test")
        self.assertIn("export GMX_MAXBACKUP=0", script)

    def test_compaction_preserves_production_and_npt_restart_files(self):
        case = generator.from_counts(2, 3, 1, 1.2, "test")
        script = generator.workflow_script(case, "storage_test")
        compact_block = script.split("COMPACT_FILES=(", 1)[1].split(")", 1)[0]
        self.assertIn('"$NPT.xtc"', compact_block)
        self.assertNotIn('"$NPT.gro"', compact_block)
        self.assertNotIn('"$NPT.cpt"', compact_block)
        self.assertNotIn('"$PROD.tpr"', compact_block)
        self.assertNotIn('"$PROD.xtc"', compact_block)
        self.assertIn("compact-preview", script)

    def test_preparation_and_npt_disable_coordinate_trajectories(self):
        for name in ("prep.mdp", "npt_eq.mdp"):
            text = (PROJECT_ROOT / name).read_text(encoding="utf-8")
            settings = {
                line.split("=", 1)[0].strip(): line.split("=", 1)[1].strip()
                for line in text.splitlines()
                if "=" in line and not line.lstrip().startswith(";")
            }
            self.assertEqual(settings["nstxout-compressed"], "0")


if __name__ == "__main__":
    unittest.main()
