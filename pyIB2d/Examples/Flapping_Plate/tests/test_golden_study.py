from __future__ import annotations

import copy
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from run_re64_golden_study import (
    DEFAULT_CYCLES,
    U_REF,
    _write_fixed_axis_svg,
    design_stage,
    detect_operational_steady_cycle,
    fixed_axis_forces,
    force_oscillation_diagnostic,
    inferred_domain_mean_velocity,
    select_case,
    validate_manifest,
)


class GoldenStudyDefinitionTests(unittest.TestCase):
    def test_stage_manifests_are_deterministic_and_have_expected_matrix(self):
        expected_counts = {
            "stiffness": 4,
            "timestep": 3,
            "box": 3,
            "grid": 3,
            "startup": 1,
            "production": 6,
        }
        for stage, expected_count in expected_counts.items():
            first = design_stage(stage)
            second = design_stage(stage)
            self.assertEqual(first, second)
            self.assertEqual(len(first["cases"]), expected_count)
            self.assertEqual(first["configuration_sha256"], second["configuration_sha256"])
            self.assertEqual(first["baseline_requested"]["cycles"], DEFAULT_CYCLES[stage])
            self.assertTrue(first["coupling_controls"]["periodic_mean_momentum_constraint"])
            for case in first["cases"]:
                self.assertFalse(case["damped_springs"])
                self.assertEqual(case["target_damping"], 0.0)
                self.assertTrue(case["balance_periodic_force"])
                self.assertAlmostEqual(case["ds_over_dx"], 0.5)
                self.assertAlmostEqual(case["steps"] * case["dt"], case["cycles"] * 0.032725)
            validate_manifest(first)

        grid = design_stage("grid")
        self.assertEqual(
            {case["dt"] for case in grid["cases"]},
            {grid["cases"][0]["dt"]},
        )
        production = design_stage("production")
        self.assertEqual(
            sum(bool(case["save_visual_snapshots"]) for case in production["cases"]),
            1,
        )
        self.assertEqual(
            {
                (case["beta_degrees"], case["pitch_amplitude_degrees"])
                for case in production["cases"]
            },
            {(0.0, 30.0), (0.0, 45.0), (30.0, 30.0), (30.0, 45.0),
             (62.8, 30.0), (62.8, 45.0)},
        )

    def test_manifest_tampering_and_bad_index_are_rejected(self):
        manifest = design_stage("startup")
        altered = copy.deepcopy(manifest)
        altered["cases"][0]["dt"] *= 2.0
        with self.assertRaisesRegex(ValueError, "manifest configuration hash"):
            validate_manifest(altered)
        with self.assertRaisesRegex(ValueError, "case index"):
            select_case(manifest, 1)

    def test_balance_can_be_disabled_explicitly(self):
        manifest = design_stage("startup", balance_periodic_force=False)
        self.assertFalse(
            manifest["coupling_controls"]["periodic_mean_momentum_constraint"]
        )
        self.assertFalse(manifest["cases"][0]["balance_periodic_force"])
        self.assertIn("No mean-force balance", manifest["analysis_note"])


class GoldenForceConventionTests(unittest.TestCase):
    def test_fixed_axis_signs_and_normalization_are_exact(self):
        force = fixed_axis_forces([2.0, -4.0], [3.0, -5.0], force_scale=2.0)
        np.testing.assert_allclose(force["lift"], [3.0, -5.0])
        np.testing.assert_allclose(force["drag"], [-2.0, 4.0])
        np.testing.assert_allclose(force["cl"], [1.5, -2.5])
        np.testing.assert_allclose(force["cd"], [-1.0, 2.0])
        self.assertGreater(U_REF, 0.0)

    def test_unbalanced_bulk_drift_integral(self):
        u, v = inferred_domain_mean_velocity(
            [2.0, 2.0], [-1.0, -1.0], dt=0.1, rho=2.0, area=5.0
        )
        np.testing.assert_allclose(u, [-0.02, -0.04])
        np.testing.assert_allclose(v, [0.01, 0.02])

    def test_fixed_axis_svg_uses_task_labels(self):
        curves = {
            "time": np.arange(4, dtype=float),
            "phase": (np.arange(4, dtype=float) + 0.5) / 4.0,
            "lift": np.asarray([1.0, 2.0, -1.0, -2.0]),
            "drag": np.asarray([-2.0, 0.0, 2.0, 0.0]),
            "cl": np.asarray([0.5, 1.0, -0.5, -1.0]),
            "cd": np.asarray([-1.0, 0.0, 1.0, 0.0]),
        }
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = _write_fixed_axis_svg(
                Path(temporary_directory) / "force.svg",
                curves,
                {"case_id": "test", "steps_per_cycle": 4},
                1,
                result_status="operationally_steady",
            )
            text = path.read_text(encoding="utf-8")
            self.assertIn("L = Fy", text)
            self.assertIn("D = -Fx", text)
            self.assertIn("coefficient using Uref", text)

    @staticmethod
    def _cycle_rows() -> list[dict[str, float | int]]:
        return [
            {"cycle": cycle, "mean_cl": 1.0 + 0.001 * cycle,
             "mean_cd": 0.001 * cycle}
            for cycle in range(1, 5)
        ]

    def test_steady_detector_requires_final_two_pairs(self):
        comparisons = [
            {"earlier_cycle": 1, "later_cycle": 2,
             "cl_relative_l2": 0.2, "cd_relative_l2": 0.2},
            {"earlier_cycle": 2, "later_cycle": 3,
             "cl_relative_l2": 0.01, "cd_relative_l2": 0.01},
            {"earlier_cycle": 3, "later_cycle": 4,
             "cl_relative_l2": 0.015, "cd_relative_l2": 0.015},
        ]
        result = detect_operational_steady_cycle(self._cycle_rows(), comparisons)
        self.assertTrue(result["reached"])
        self.assertEqual(result["averaging_cycle"], 4)

        failed = copy.deepcopy(comparisons)
        failed[-1]["cd_relative_l2"] = 0.021
        result = detect_operational_steady_cycle(self._cycle_rows(), failed)
        self.assertFalse(result["reached"])
        self.assertIsNone(result["averaging_cycle"])

    def test_alternating_force_diagnostic_detects_nyquist_component(self):
        values = (-1.0) ** np.arange(100)
        diagnostic = force_oscillation_diagnostic(values)
        self.assertTrue(math.isclose(
            diagnostic["nyquist_amplitude_over_fluctuation_rms"], 1.0
        ))
        self.assertGreater(
            diagnostic["second_difference_rms_over_fluctuation_rms"], 3.9
        )


if __name__ == "__main__":
    unittest.main()
