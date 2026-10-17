from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from generate_case import build_case
from update_Target_Point_Positions import reset_reference, update_Target_Point_Positions


class KinematicsTests(unittest.TestCase):
    def test_absolute_update_does_not_accumulate_drift(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            build_case(case_dir, nx=32, ny=32, dt=2.5e-5, steps=4)
            vertices = np.loadtxt(case_dir / "plate.vertex", skiprows=1)
            count = vertices.shape[0]
            targets = np.column_stack((
                np.arange(count),
                vertices,
                np.full(count, 5.0e7),
            ))

            previous_directory = Path.cwd()
            try:
                os.chdir(case_dir)
                reset_reference()
                first = update_Target_Point_Positions(2.5e-5, 0.0, targets.copy())
                repeated = update_Target_Point_Positions(2.5e-5, 0.0, first.copy())
            finally:
                os.chdir(previous_directory)

            np.testing.assert_allclose(first, repeated, rtol=0.0, atol=1e-15)
            self.assertGreater(first[0, 1], vertices[0, 0])
            np.testing.assert_allclose(first[:, 2], vertices[:, 1])


if __name__ == "__main__":
    unittest.main()
