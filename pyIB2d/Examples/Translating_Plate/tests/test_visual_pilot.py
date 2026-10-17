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

from force_history import HEADER
from case_config import translation_distance
from run_re128_visual_pilot import (
    PILOT_CASE,
    SNAPSHOT_EVERY_STEPS,
    expected_snapshot_steps,
)
from visual_pilot import (
    VisualSnapshotWriter,
    build_visual_artifacts,
    periodic_vorticity,
)


class VisualPilotTests(unittest.TestCase):
    def test_periodic_vorticity_matches_smooth_field(self):
        nx = ny = 64
        lx = ly = 1.0
        dx, dy = lx / nx, ly / ny
        x = np.arange(nx) * dx
        y = np.arange(ny) * dy
        xx, yy = np.meshgrid(x, y)
        u = np.sin(2.0 * math.pi * yy)
        v = np.sin(2.0 * math.pi * xx)
        exact = 2.0 * math.pi * (
            np.cos(2.0 * math.pi * xx) - np.cos(2.0 * math.pi * yy)
        )
        error = np.max(np.abs(periodic_vorticity(u, v, dx, dy) - exact))
        self.assertLess(error, 0.03)

    def test_snapshot_and_dependency_free_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = {
                "case_name": "visual_test",
                "domain": {"lx": 1.0, "ly": 1.0},
                "grid": {"nx": 8, "ny": 8, "dx": 0.125, "dy": 0.125},
                "geometry": {"nominal_chord": 0.1},
                "motion": {"speed": 1.0},
                "fluid": {"rho": 1.0},
                "force": {"coefficient_force_scale": 0.05},
            }
            (case_dir / "case_metadata.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            writer = VisualSnapshotWriter.from_case(case_dir, every_steps=2)
            coordinate = np.arange(8) / 8.0
            xx, yy = np.meshgrid(coordinate, coordinate)
            positions = np.asarray([[0.4, 0.45], [0.5, 0.55]])
            targets = np.asarray([
                [0.0, 0.4, 0.45, 1.0],
                [1.0, 0.5, 0.55, 1.0],
            ])
            for step in range(3):
                writer.record(
                    step=step,
                    time=0.01 * step,
                    velocity_x=np.sin(2.0 * math.pi * yy),
                    velocity_y=np.sin(2.0 * math.pi * xx),
                    lagrangian_positions=positions,
                    target_info=targets,
                    dx=0.125,
                    dy=0.125,
                    lx=1.0,
                    ly=1.0,
                )
            self.assertEqual(
                len(list((case_dir / "visual_snapshots").glob("*.npz"))), 2
            )

            with (case_dir / "force_history.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                csv_writer = csv.DictWriter(handle, fieldnames=HEADER)
                csv_writer.writeheader()
                for step in range(16):
                    row = {name: 0.0 for name in HEADER}
                    row.update({
                        "step": step,
                        "time": (step + 0.5) * 0.001,
                        "fx_body": -math.sin(step / 8.0),
                        "fy_body": math.cos(step / 8.0),
                        "max_target_error": 1.0e-3,
                        "rms_target_error": 5.0e-4,
                    })
                    csv_writer.writerow(row)

            manifest = build_visual_artifacts(
                case_dir, make_mp4=False, max_width=480
            )
            output_dir = case_dir / "visual_pilot"
            first_frame = output_dir / "frames" / "frame_000000.png"
            self.assertEqual(first_frame.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertTrue((output_dir / "animation.html").is_file())
            self.assertTrue((output_dir / "force_and_target_history.svg").is_file())
            self.assertTrue((output_dir / "oscillation_diagnostics.json").is_file())
            self.assertEqual(manifest["frame_count"], 2)

    def test_pilot_snapshot_schedule_includes_both_endpoints(self):
        steps = expected_snapshot_steps(PILOT_CASE, SNAPSHOT_EVERY_STEPS)
        self.assertEqual(steps[0], 0)
        self.assertEqual(steps[-1], PILOT_CASE.steps)
        self.assertEqual(len(steps), 101)

    def test_midpoint_target_is_shifted_to_snapshot_time(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = {
                "domain": {"lx": 1.0, "ly": 1.0},
                "grid": {"nx": 4, "ny": 4, "dx": 0.25, "dy": 0.25},
                "motion": {
                    "kind": "ramped_translation_x",
                    "speed": 24.0,
                    "ramp_rate": 200.0,
                    "stop_time": 0.025,
                },
            }
            (case_dir / "case_metadata.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            writer = VisualSnapshotWriter.from_case(case_dir)
            target_time = 0.005
            state_time = 0.010
            reference_x = 0.2
            target_x = reference_x + translation_distance(target_time)
            path = writer.record(
                step=1,
                time=state_time,
                target_time=target_time,
                velocity_x=np.zeros((4, 4)),
                velocity_y=np.zeros((4, 4)),
                lagrangian_positions=np.asarray([[target_x, 0.4]]),
                target_info=np.asarray([[0.0, target_x, 0.4, 1.0]]),
            )
            self.assertIsNotNone(path)
            with np.load(path, allow_pickle=False) as snapshot:
                self.assertAlmostEqual(float(snapshot["target_time"]), state_time)
                self.assertAlmostEqual(
                    float(snapshot["target_positions"][0, 0]),
                    reference_x + translation_distance(state_time),
                )


if __name__ == "__main__":
    unittest.main()
