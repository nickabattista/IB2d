from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from generate_case import build_case
from main2d import run_case


class SmokeTests(unittest.TestCase):
    def test_two_step_headless_run(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            build_case(case_dir, nx=32, ny=32, dt=2.5e-5, steps=2)
            states = []

            def collect_state(**state):
                states.append((state["step"], state["time"], state["velocity_x"].shape))

            history = run_case(
                case_dir,
                state_callback=collect_state,
                state_callback_interval=1,
            )

            with history.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertAlmostEqual(float(rows[0]["time"]), 1.25e-5)
            self.assertAlmostEqual(float(rows[-1]["time"]), 3.75e-5)
            self.assertLess(float(rows[-1]["force_balance_linf"]), 1e-9)
            self.assertEqual([state[0] for state in states], [0, 1, 2])
            self.assertEqual([state[2] for state in states], [(32, 32)] * 3)


if __name__ == "__main__":
    unittest.main()
