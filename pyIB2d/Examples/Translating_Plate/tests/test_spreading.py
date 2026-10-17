from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
CORE = HERE.parents[1] / "IBM_Blackbox"
sys.path.insert(0, str(CORE))

from please_Find_Lagrangian_Forces_On_Eulerian_grid import (
    give_Me_Delta_Function_Approximations_For_Force_Calc,
    spread_Lagrangian_Forces_To_Eulerian_Grid,
)


class SpreadingTests(unittest.TestCase):
    def test_local_spread_matches_legacy_dense_product(self):
        nx = 16
        ny = 16
        lx = 1.0
        ly = 1.0
        dx = lx / nx
        dy = ly / ny
        x = np.arange(nx) * dx
        y = np.arange(ny) * dy
        x_lag = np.array([0.01, 0.37, 0.99])
        y_lag = np.array([0.98, 0.52, 0.02])
        fx = np.array([2.0, -1.0, 0.25])
        fy = np.array([-0.5, 3.0, 1.5])
        ds = dx / 2.0
        grid_info = [nx, ny, lx, ly, dx, dy, 4, x_lag.size, ds]

        delta_x, delta_y = give_Me_Delta_Function_Approximations_For_Force_Calc(
            x, y, grid_info, x_lag, y_lag
        )
        dense_fx = ((fx * ds) * delta_y) @ delta_x
        dense_fy = ((fy * ds) * delta_y) @ delta_x
        local_fx, local_fy = spread_Lagrangian_Forces_To_Eulerian_Grid(
            fx, fy, ds, x, y, grid_info, x_lag, y_lag
        )

        np.testing.assert_allclose(local_fx, dense_fx, rtol=2e-15, atol=1e-14)
        np.testing.assert_allclose(local_fy, dense_fy, rtol=2e-15, atol=1e-14)
        self.assertAlmostEqual(dx * dy * local_fx.sum(), ds * fx.sum(), places=13)
        self.assertAlmostEqual(dx * dy * local_fy.sum(), ds * fy.sum(), places=13)


if __name__ == "__main__":
    unittest.main()
