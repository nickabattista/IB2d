from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from force_history import ForceHistoryWriter


class ForceHistoryTests(unittest.TestCase):
    def test_reaction_and_lagrangian_balance_are_recorded(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "forces.csv"
            force_density = np.array([[2.0, -1.0], [4.0, 3.0]])
            ds = 0.25
            fluid_force = ds * force_density.sum(axis=0)
            target_info = np.array([
                [0.0, 0.1, 0.2, 1.0],
                [1.0, 0.3, 0.4, 1.0],
            ])
            positions = target_info[:, 1:3].copy()

            with ForceHistoryWriter(output, lx=1.0, ly=1.0) as writer:
                writer.record(
                    step=0,
                    time=0.5,
                    body_force=-fluid_force,
                    fluid_force=fluid_force,
                    lagrangian_force_density=force_density,
                    lagrangian_weight=ds,
                    lagrangian_positions=positions,
                    target_info=target_info,
                )

            with output.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertAlmostEqual(float(row["fx_body"]), -1.5)
            self.assertAlmostEqual(float(row["fy_body"]), -0.5)
            self.assertAlmostEqual(float(row["force_balance_linf"]), 0.0)
            self.assertAlmostEqual(float(row["max_target_error"]), 0.0)


if __name__ == "__main__":
    unittest.main()
