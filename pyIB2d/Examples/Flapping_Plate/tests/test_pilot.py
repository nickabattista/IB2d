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

from generate_case import build_case
from run_re64_flapping_pilot import (
    PILOT_CASES,
    PILOT_DYNAMIC_VISCOSITY,
    PILOT_PEAK_SPEED,
    PilotCase,
    _postprocessor_source_fingerprints,
    _snapshot_producer_source_fingerprints,
    _validate_completed_case,
    _expected_snapshot_steps,
    _write_analysis_manifest,
    _write_comparison_csv,
    _write_comparison_svg,
    compare_timestep_cases,
    pilot_manifest,
    run_pilot_case,
)


class PilotDefinitionTests(unittest.TestCase):
    def test_paired_cases_are_cycle_and_space_aligned(self):
        coarse, fine = PILOT_CASES
        self.assertEqual((coarse.steps, fine.steps), (52360, 104720))
        self.assertEqual(fine.dt, 0.5 * coarse.dt)
        self.assertAlmostEqual(coarse.tfinal, fine.tfinal)
        self.assertAlmostEqual(coarse.actual_cycles, 2.0)
        self.assertAlmostEqual(fine.actual_cycles, 2.0)
        self.assertEqual(coarse.snapshot_interval * 2, fine.snapshot_interval)
        for case in PILOT_CASES:
            self.assertEqual((case.nx, case.ny), (512, 512))
            self.assertAlmostEqual(case.ds / case.dx, 0.5)
            self.assertAlmostEqual(case.tether_number, 0.025 / (4 ** case.index))
            self.assertEqual(len(_expected_snapshot_steps(case)), 41)

        manifest = pilot_manifest()
        self.assertEqual(manifest["analysis"]["movie_case_index"], 1)
        self.assertAlmostEqual(
            manifest["physics"]["peak_translation_speed"], PILOT_PEAK_SPEED
        )
        self.assertAlmostEqual(
            manifest["physics"]["dynamic_viscosity"], PILOT_DYNAMIC_VISCOSITY
        )
        self.assertEqual(manifest["motion"]["stroke_angle_degrees"], 62.8)
        self.assertEqual(manifest["motion"]["pitch_amplitude_degrees"], 45.0)
        self.assertEqual(len(manifest["configuration_sha256"]), 64)

    @staticmethod
    def _write_synthetic_case(case_dir: Path, samples_per_cycle: int) -> None:
        period = 0.04
        dt = period / samples_per_cycle
        steps = 2 * samples_per_cycle
        build_case(case_dir, nx=32, ny=32, dt=dt, steps=steps, period=period)
        times = (np.arange(steps, dtype=float) + 0.5) * dt
        phase = 2.0 * math.pi * times / period
        force_x = 3.0 * np.sin(phase) + 0.2 * np.cos(2.0 * phase)
        force_y = -2.0 * np.cos(phase) + 0.1 * np.sin(3.0 * phase)
        with (case_dir / "force_history.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.writer(handle)
            writer.writerow(("step", "time", "fx_body", "fy_body"))
            for index, (time_value, fx, fy) in enumerate(
                zip(times, force_x, force_y)
            ):
                writer.writerow((index, time_value, fx, fy))

    def test_synthetic_cross_timestep_comparison_and_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            coarse, fine = root / "coarse", root / "fine"
            self._write_synthetic_case(coarse, 80)
            self._write_synthetic_case(fine, 160)
            comparison, table = compare_timestep_cases((coarse, fine))
            self.assertEqual(comparison["cycle"], 2)
            self.assertEqual(len(table["phase"]), 1441)
            for field in ("ch_urms", "cv_urms"):
                self.assertLess(
                    comparison["series"][field]["relative_l2_using_fine"],
                    2.0e-3,
                )
            reversal_index = int(np.argmin(np.abs(table["phase"] - 0.5)))
            self.assertTrue(
                math.isnan(float(table["coarse_cd_urms"][reversal_index]))
            )
            self.assertGreater(
                comparison["excluded_common_phase_sample_count"], 0
            )
            self.assertLess(
                comparison["series"]["cd_urms"]["comparison_sample_count"],
                1441,
            )

            csv_path = _write_comparison_csv(root / "comparison.csv", table)
            svg_path = _write_comparison_svg(root / "comparison.svg", table)
            self.assertEqual(len(csv_path.read_text(encoding="utf-8").splitlines()), 1442)
            self.assertIn("timestep comparison", svg_path.read_text(encoding="utf-8"))

    def test_small_run_is_atomic_and_rejects_stale_source_provenance(self):
        case = PilotCase(
            index=99,
            case_id="test_two_step_case",
            dt=1.25e-6,
            steps=2,
            snapshot_interval=1,
            nx=32,
            ny=32,
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            case_dir = run_pilot_case(case, root)
            completion = _validate_completed_case(case_dir, case)
            self.assertEqual(completion["snapshot_steps"], [0, 1, 2])
            self.assertEqual(
                len(completion["snapshot_producer_source_set_sha256"]), 64
            )
            self.assertFalse(any((root / "cases").glob("*.partial-*")))

            metadata_path = case_dir / "case_metadata.json"
            saved = metadata_path.read_text(encoding="utf-8")
            metadata_value = json.loads(saved)
            self.assertEqual(
                metadata_value["pilot"]["snapshot_producer_source_sha256"],
                _snapshot_producer_source_fingerprints(),
            )
            self.assertEqual(
                metadata_value["pilot"]["postprocessor_source_sha256_at_run"],
                _postprocessor_source_fingerprints(),
            )
            metadata_value["pilot"]["snapshot_producer_source_sha256"] = {}
            metadata_path.write_text(
                json.dumps(metadata_value), encoding="utf-8"
            )
            with self.assertRaisesRegex(
                ValueError, "different snapshot-producing source"
            ):
                _validate_completed_case(case_dir, case)
            metadata_path.write_text(saved, encoding="utf-8")

            metadata_value = json.loads(saved)
            metadata_value["pilot"]["numerical_source_sha256"] = {}
            metadata_path.write_text(
                json.dumps(metadata_value), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "different numerical source"):
                _validate_completed_case(case_dir, case)
            metadata_path.write_text(saved, encoding="utf-8")

    def test_analysis_manifest_records_current_postprocessor_provenance(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory).resolve()
            summary_dir = root / "summary"
            artifact_path = root / "cases" / "fine" / "force_history.csv"
            manifest = _write_analysis_manifest(
                summary_dir,
                {
                    "schema_version": 1,
                    "artifact": str(artifact_path),
                },
                root,
            )

            saved = json.loads(
                (summary_dir / "pilot_analysis_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(saved, manifest)
            self.assertEqual(
                saved["postprocessor_source_sha256"],
                _postprocessor_source_fingerprints(),
            )
            self.assertEqual(
                len(saved["postprocessor_source_set_sha256"]), 64
            )
            self.assertEqual(saved["artifact"], "cases/fine/force_history.csv")


if __name__ == "__main__":
    unittest.main()
