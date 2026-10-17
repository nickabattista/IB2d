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

from case_config import NOMINAL_CHORD, PEAK_TRANSLATION_SPEED, REYNOLDS_NUMBER
from flapping_visual import FlappingVisualSnapshotWriter, build_visual_artifacts
from generate_case import build_case
from kinematics import (
    FlappingMotion,
    body_coordinates_from_world,
    world_coordinates_from_body,
)
from main2d import run_case


class GeneratedCaseTests(unittest.TestCase):
    def test_grid_marker_and_force_scaling_contract(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = build_case(case_dir, nx=64, ny=64, steps=8)
            dx = metadata["grid"]["dx"]
            ds = metadata["geometry"]["marker_weight"]
            self.assertAlmostEqual(ds / dx, 0.5)
            self.assertAlmostEqual(
                metadata["fluid"]["reynolds_number"], REYNOLDS_NUMBER
            )
            self.assertAlmostEqual(
                metadata["fluid"]["rho"]
                * PEAK_TRANSLATION_SPEED
                * NOMINAL_CHORD
                / metadata["fluid"]["mu"],
                REYNOLDS_NUMBER,
            )
            self.assertAlmostEqual(
                metadata["force"]["umax_force_scale"],
                2.0 * metadata["force"]["urms_force_scale"],
            )
            self.assertEqual(
                metadata["motion"]["angle_convention"],
                "theta=beta-alpha, leading_edge_to_trailing_edge",
            )

            vertex = np.loadtxt(case_dir / "plate.vertex", skiprows=1)
            distances = np.linalg.norm(np.diff(vertex, axis=0), axis=1)
            np.testing.assert_allclose(distances, ds, rtol=1.0e-13, atol=1.0e-15)
            self.assertAlmostEqual(vertex.shape[0] * ds / NOMINAL_CHORD, 1.0)

    def test_real_two_step_solver_and_flapping_snapshots(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = build_case(
                case_dir, nx=32, ny=32, dt=1.25e-6, steps=2
            )
            writer = FlappingVisualSnapshotWriter.from_case(
                case_dir, every_steps=1
            )
            history_path = run_case(
                case_dir,
                state_callback=writer.record,
                state_callback_interval=1,
            )
            with history_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            for row in rows:
                self.assertTrue(
                    all(math.isfinite(float(value)) for value in row.values())
                )

            snapshots = sorted((case_dir / "visual_snapshots").glob("*.npz"))
            self.assertEqual(len(snapshots), 3)
            with np.load(snapshots[-1], allow_pickle=False) as snapshot:
                self.assertEqual(int(snapshot["step"]), 2)
                self.assertAlmostEqual(float(snapshot["target_time"]), 2.5e-6)
                target = np.asarray(snapshot["target_positions"], dtype=float)
            motion = FlappingMotion.from_mapping(metadata["motion"])
            initial = np.loadtxt(case_dir / "plate.vertex", skiprows=1)
            body = body_coordinates_from_world(initial, 0.0, motion)
            # The saved full-step target must be the same rigid body evaluated
            # at the full-step state time, not the solver's midpoint timestamp.
            expected = world_coordinates_from_body(body, 2.5e-6, motion)
            np.testing.assert_allclose(target, expected, atol=1.0e-15)

            artifacts = build_visual_artifacts(
                case_dir, make_mp4=False, max_width=480
            )
            self.assertEqual(
                artifacts["artifact_type"], "flapping_plate_visual_pilot"
            )
            self.assertEqual(artifacts["manifest"], "visual_manifest.json")
            diagnostics = json.loads(
                (case_dir / "visual_pilot" / "oscillation_diagnostics.json")
                .read_text(encoding="utf-8")
            )
            self.assertIn("CH", diagnostics["series"])
            self.assertIn("CV", diagnostics["series"])
            self.assertNotIn("CD", diagnostics["series"])


if __name__ == "__main__":
    unittest.main()
