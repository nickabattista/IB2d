from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
CORE = HERE.parents[1] / "IBM_Blackbox"
sys.path.insert(0, str(CORE))

from Supp import give_Delta_Kernel


class DeltaKernelTests(unittest.TestCase):
    def test_key_four_point_kernel_values(self):
        dx = 0.25
        offsets = np.array([[0.0, dx, 2.0 * dx, 3.0 * dx]])
        values = give_Delta_Kernel(offsets, dx)
        np.testing.assert_allclose(
            values,
            np.array([[0.5 / dx, 0.25 / dx, 0.0, 0.0]]),
            rtol=0.0,
            atol=1e-14,
        )


if __name__ == "__main__":
    unittest.main()
