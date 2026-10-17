from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parents[1]
CORE = HERE.parents[1] / "IBM_Blackbox"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CORE))

import IBM_Driver
import main2d


class PeriodicForceBalanceTests(unittest.TestCase):
    def test_mean_removal_preserves_force_fluctuations_and_inputs(self):
        force_x = np.array([[4.0, 7.0], [-2.0, 3.0]])
        force_y = np.array([[1.5, -3.5], [8.5, 5.5]])
        original_x = force_x.copy()
        original_y = force_y.copy()

        balanced_x, balanced_y = IBM_Driver._subtract_spatial_mean_force(
            force_x, force_y
        )

        self.assertAlmostEqual(float(np.mean(balanced_x)), 0.0)
        self.assertAlmostEqual(float(np.mean(balanced_y)), 0.0)
        np.testing.assert_allclose(
            balanced_x - balanced_x[0, 0],
            force_x - force_x[0, 0],
        )
        np.testing.assert_allclose(
            balanced_y - balanced_y[0, 0],
            force_y - force_y[0, 0],
        )
        np.testing.assert_array_equal(force_x, original_x)
        np.testing.assert_array_equal(force_y, original_y)

    def test_run_case_forwards_explicit_constraint_to_driver(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            for name in ("input2d", "plate.vertex", "plate.target"):
                (case_dir / name).write_text("placeholder\n", encoding="utf-8")
            (case_dir / "case_metadata.json").write_text(
                json.dumps({"domain": {"lx": 1.0, "ly": 1.0}}),
                encoding="utf-8",
            )

            with mock.patch.object(
                main2d, "please_Initialize_Simulation", return_value=(object(),)
            ), mock.patch.object(main2d.Driver, "main") as driver:
                main2d.run_case(case_dir, balance_periodic_force=True)

            self.assertTrue(driver.call_args.kwargs["balance_periodic_force"])


if __name__ == "__main__":
    unittest.main()
