from __future__ import annotations

import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from analyze_forces import (
    OUTPUT_HEADER,
    analyze_case,
    build_summary,
    compare_cycles,
    derive_coefficients,
    motion_and_force_parameters,
)


def metadata() -> dict:
    return {
        "case_name": "synthetic_flapping_plate",
        "fluid": {"rho": 2.0, "mu": 0.01},
        "geometry": {"nominal_chord": 0.5},
        "motion": {
            "kind": "centered_sinusoidal_stroke_and_pitch",
            "period": 2.0,
            "frequency": 0.5,
            "stroke_midpoint": [0.8, 0.8],
            "stroke_amplitude": 4.0 / math.pi,
            "stroke_angle_degrees": 30.0,
            "mean_angle_of_attack_degrees": 60.0,
            "pitch_amplitude_degrees": 45.0,
            "pitch_phase_radians": 0.0,
            "pivot_fraction": 0.5,
        },
        "force": {
            # Umax = pi*A/T = 2, Urms = sqrt(2).
            "umax_force_scale": 2.0,
            "urms_force_scale": 1.0,
        },
    }


class ForceAnalysisTests(unittest.TestCase):
    def test_parameter_resolution_and_normalizations(self):
        parameters = motion_and_force_parameters(metadata())
        self.assertAlmostEqual(parameters["umax"], 2.0)
        self.assertAlmostEqual(parameters["urms"], math.sqrt(2.0))
        self.assertAlmostEqual(parameters["umax_force_scale"], 2.0)
        self.assertAlmostEqual(parameters["urms_force_scale"], 1.0)

    def test_true_motion_axes_and_reversal_flag(self):
        # At phase 1/4 the centre moves opposite the beta stroke direction.
        beta = math.radians(30.0)
        tangent = -np.asarray([math.cos(beta), math.sin(beta)])
        normal = np.asarray([-tangent[1], tangent[0]])
        force = 2.0 * (-2.0 * tangent + 3.0 * normal)
        history = {
            "step": np.asarray([0.0, 1.0]),
            "time": np.asarray([0.0, 0.5]),
            "fx_body": np.asarray([0.0, force[0]]),
            "fy_body": np.asarray([0.0, force[1]]),
        }
        curves = derive_coefficients(history, metadata())
        self.assertTrue(bool(curves["reversal"][0]))
        self.assertTrue(math.isnan(float(curves["cd_umax"][0])))
        self.assertTrue(math.isnan(float(curves["cl_umax"][0])))
        self.assertFalse(bool(curves["reversal"][1]))
        self.assertAlmostEqual(float(curves["cd_umax"][1]), 2.0)
        self.assertAlmostEqual(float(curves["cl_umax"][1]), 3.0)
        self.assertAlmostEqual(float(curves["cd_urms"][1]), 4.0)
        self.assertAlmostEqual(float(curves["cl_urms"][1]), 6.0)
        self.assertAlmostEqual(float(curves["ch_umax"][1]), force[0] / 2.0)
        self.assertAlmostEqual(float(curves["cv_umax"][1]), force[1] / 2.0)

    def test_midpoint_cells_touching_reversal_are_masked(self):
        times = np.arange(0.1, 2.0, 0.2)
        history = {
            "step": np.arange(times.size, dtype=float),
            "time": times,
            "fx_body": np.ones(times.size),
            "fy_body": np.ones(times.size),
        }
        curves = derive_coefficients(
            history, metadata(), reversal_speed_fraction=0.0
        )
        expected = np.zeros(times.size, dtype=bool)
        expected[[0, 4, 5, 9]] = True
        np.testing.assert_array_equal(curves["reversal"], expected)
        self.assertTrue(np.all(np.isnan(curves["cd_urms"][expected])))
        self.assertTrue(np.all(np.isfinite(curves["cd_urms"][~expected])))

    def test_cycle_comparison_does_not_interpolate_across_reversal_gaps(self):
        samples_per_cycle = 20
        times = (
            np.arange(2 * samples_per_cycle, dtype=float) + 0.5
        ) * (2.0 / samples_per_cycle)
        phase = np.mod(times / 2.0, 1.0)
        history = {
            "step": np.arange(times.size, dtype=float),
            "time": times,
            "fx_body": 1.0 + 0.25 * np.cos(2.0 * math.pi * phase),
            "fy_body": 2.0 + 0.50 * np.sin(2.0 * math.pi * phase),
        }
        curves = derive_coefficients(
            history, metadata(), reversal_speed_fraction=0.0
        )
        comparison = compare_cycles(curves, 1, 2)
        self.assertIsNotNone(comparison)

        for name in ("cd_umax", "cl_umax", "cd_urms", "cl_urms"):
            series = comparison["series"][name]
            self.assertIsNotNone(series)
            self.assertGreater(series["excluded_phase_sample_count"], 0)
            self.assertLess(
                series["comparison_sample_count"], series["phase_sample_count"]
            )
            self.assertLess(series["relative_l2"], 1.0e-12)

        for name in ("ch_umax", "cv_umax", "ch_urms", "cv_urms"):
            series = comparison["series"][name]
            self.assertIsNotNone(series)
            self.assertEqual(series["excluded_phase_sample_count"], 0)
            self.assertLess(series["relative_l2"], 1.0e-12)

    def test_phase_and_one_based_cycle_columns(self):
        history = {
            "step": np.arange(4, dtype=float),
            "time": np.asarray([0.25, 2.25, 4.25, 6.25]),
            "fx_body": np.ones(4),
            "fy_body": np.zeros(4),
        }
        curves = derive_coefficients(history, metadata())
        np.testing.assert_array_equal(curves["cycle"], [1, 2, 3, 4])
        np.testing.assert_allclose(curves["phase"], 0.125)
        np.testing.assert_allclose(curves["time_over_period"], [0.125, 1.125, 2.125, 3.125])

    def test_end_to_end_artifacts_and_cycle_three_four_summary(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            (case_dir / "case_metadata.json").write_text(
                json.dumps(metadata()), encoding="utf-8"
            )
            samples_per_cycle = 16
            times = (
                np.arange(4 * samples_per_cycle, dtype=float) + 0.5
            ) * (2.0 / samples_per_cycle)
            phase = np.mod(times / 2.0, 1.0)
            with (case_dir / "force_history.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                fields = (
                    "step", "time", "fx_body", "fy_body", "force_balance_linf",
                    "max_target_error", "rms_target_error",
                )
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for step, (time, position) in enumerate(zip(times, phase)):
                    writer.writerow({
                        "step": step,
                        "time": time,
                        "fx_body": 2.0 * math.cos(2.0 * math.pi * position),
                        "fy_body": 4.0 + math.sin(2.0 * math.pi * position),
                        "force_balance_linf": 1.0e-12,
                        "max_target_error": 2.0e-4,
                        "rms_target_error": 1.0e-4,
                    })

            result = analyze_case(case_dir)
            output_dir = case_dir / "force_analysis"
            for name in (
                "force_coefficients.csv", "cycle_summary.json",
                "force_coefficients.svg", "cycle_3_4_overlay.svg",
            ):
                self.assertTrue((output_dir / name).is_file(), name)
            self.assertEqual(
                (output_dir / "force_coefficients.svg").read_text(
                    encoding="utf-8"
                ).lstrip()[:5],
                "<?xml",
            )
            with (output_dir / "force_coefficients.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(tuple(reader.fieldnames or ()), OUTPUT_HEADER)
                rows = list(reader)
            self.assertEqual(len(rows), 64)
            self.assertEqual(rows[0]["cycle"], "1")

            summary = result["summary"]
            self.assertEqual(summary["requested_analysis_cycles"], [3, 4])
            self.assertEqual(
                [item["available"] for item in summary["cycle_summaries"]],
                [True, True],
            )
            cycle_three = summary["cycle_summaries"][0]
            self.assertAlmostEqual(
                cycle_three["series"]["fy_body"]["mean"], 4.0, places=12
            )
            comparison = summary["cycle_comparison"]
            self.assertIsNotNone(comparison)
            self.assertLess(
                comparison["series"]["cv_urms"]["relative_l2"], 1.0e-12
            )

    def test_short_pilot_reports_missing_requested_cycles_without_failure(self):
        history = {
            "step": np.arange(8, dtype=float),
            "time": (np.arange(8, dtype=float) + 0.5) * 0.25,
            "fx_body": np.ones(8),
            "fy_body": np.zeros(8),
        }
        curves = derive_coefficients(history, metadata())
        summary = build_summary(curves, metadata(), analysis_cycles=(3, 4))
        self.assertEqual(summary["available_cycles"], [1])
        self.assertFalse(summary["cycle_summaries"][0]["available"])
        self.assertFalse(summary["cycle_summaries"][1]["available"])
        self.assertIsNone(summary["cycle_comparison"])


if __name__ == "__main__":
    unittest.main()
