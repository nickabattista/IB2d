from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from analyze_re128_study import (
    difference_metrics,
    high_frequency_energy_fraction,
    observed_order,
)
from generate_case import build_case
from force_history import HEADER
from re128_study import (
    FINAL_TIME,
    MU_RE128,
    NOMINAL_CHORD,
    PRODUCTION_CASES,
    REYNOLDS_NUMBER,
    RHO,
    TRANSLATION_SPEED,
    case_record,
    centered_leading_edge,
    generate_study_case,
    run_study_case,
)


class Re128ManifestTests(unittest.TestCase):
    def test_production_families_hold_the_intended_controls(self):
        grid = sorted(
            (case for case in PRODUCTION_CASES if "grid" in case.study_axes),
            key=lambda case: case.cells_per_chord,
        )
        boundary = sorted(
            (case for case in PRODUCTION_CASES if "boundary" in case.study_axes),
            key=lambda case: case.box_chords,
        )

        self.assertEqual([case.cells_per_chord for case in grid], [32.0, 64.0, 128.0])
        self.assertEqual([case.box_chords for case in grid], [16.0, 16.0, 16.0])
        self.assertEqual([case.box_chords for case in boundary], [10.0, 16.0, 24.0, 32.0])
        self.assertTrue(
            all(math.isclose(case.cells_per_chord, 64.0) for case in boundary)
        )
        self.assertTrue(
            all(math.isclose(case.ds, 0.5 * case.dx) for case in PRODUCTION_CASES)
        )
        self.assertEqual([case.dt for case in grid], [2.5e-6] * 3)
        self.assertEqual([case.target_stiffness for case in grid], [1.0e9] * 3)
        self.assertTrue(
            all(
                math.isclose(
                    case.target_stiffness,
                    grid[0].target_stiffness,
                )
                for case in grid
            )
        )
        self.assertEqual(
            [round(case.penalty_coupling_number, 12) for case in grid],
            [0.1, 0.2, 0.4],
        )
        self.assertTrue(
            all(
                math.isclose(case.dt, boundary[0].dt)
                and math.isclose(case.target_stiffness, boundary[0].target_stiffness)
                and math.isclose(case.ds, boundary[0].ds)
                for case in boundary
            )
        )
        self.assertEqual(len({case.case_id for case in PRODUCTION_CASES}), len(PRODUCTION_CASES))
        self.assertAlmostEqual(
            RHO * TRANSLATION_SPEED * NOMINAL_CHORD / MU_RE128,
            REYNOLDS_NUMBER,
        )

    def test_generated_study_case_has_exact_chord_re_and_centering(self):
        case = PRODUCTION_CASES[2]
        with tempfile.TemporaryDirectory() as temporary_directory:
            case_dir = Path(temporary_directory)
            metadata = generate_study_case(case, case_dir)
            loaded = json.loads((case_dir / "case_metadata.json").read_text())

            self.assertAlmostEqual(
                metadata["geometry"]["effective_quadrature_chord"], NOMINAL_CHORD
            )
            self.assertAlmostEqual(
                loaded["fluid"]["nominal_reynolds_number"], REYNOLDS_NUMBER
            )
            self.assertEqual(loaded["geometry"]["marker_count"], 128)
            self.assertEqual(loaded["grid"]["cells_per_nominal_chord"], 64.0)

            x0, y0, placement = centered_leading_edge(case)
            self.assertEqual(loaded["geometry"]["leading_edge"], [x0, y0])
            self.assertAlmostEqual(
                x0, placement["minimum_x_clearance"], places=14
            )
            self.assertAlmostEqual(
                y0, placement["minimum_y_clearance"], places=14
            )

    def test_generator_rejects_non_square_cells(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            with self.assertRaisesRegex(ValueError, "square Eulerian cells"):
                build_case(
                    Path(temporary_directory), nx=32, ny=32, lx=1.0, ly=2.0
                )

    def test_completed_case_must_match_requested_configuration(self):
        case = replace(PRODUCTION_CASES[0], dt=FINAL_TIME)
        with tempfile.TemporaryDirectory() as temporary_directory:
            final_dir = (
                Path(temporary_directory) / "cases" / case.phase / case.case_id
            )
            final_dir.mkdir(parents=True)
            generate_study_case(case, final_dir)
            (final_dir / "COMPLETE.json").write_text(
                json.dumps(case_record(case)), encoding="utf-8"
            )
            values = {name: 0.0 for name in HEADER}
            values["step"] = 0
            values["time"] = 0.5 * case.dt
            (final_dir / "force_history.csv").write_text(
                ",".join(HEADER) + "\n"
                + ",".join(str(values[name]) for name in HEADER) + "\n",
                encoding="utf-8",
            )
            self.assertEqual(
                run_study_case(case, Path(temporary_directory)), final_dir
            )
            with self.assertRaisesRegex(FileExistsError, "different configuration"):
                run_study_case(
                    replace(case, dt=0.5 * case.dt), Path(temporary_directory)
                )


class Re128AnalysisTests(unittest.TestCase):
    @staticmethod
    def _loaded(case, coefficient):
        time = np.linspace(0.0, 1.0, coefficient.size)
        return {"case": case, "history": {"time": time, "cd": coefficient}}

    def test_observed_order_for_quadratic_synthetic_error(self):
        cases = sorted(
            (case for case in PRODUCTION_CASES if "grid" in case.study_axes),
            key=lambda case: case.cells_per_chord,
        )
        time = np.linspace(0.0, 1.0, 51)
        exact = np.sin(2.0 * math.pi * time)
        loaded = [
            self._loaded(case, exact + case.dx * case.dx * np.cos(time))
            for case in cases
        ]
        result = observed_order(*loaded, field="cd")
        self.assertAlmostEqual(result["observed_order"], 2.0, places=10)

    def test_normalized_difference_is_zero_for_identical_signal(self):
        case = PRODUCTION_CASES[0]
        values = np.linspace(-2.0, 3.0, 21)
        reference = self._loaded(case, values)
        candidate = self._loaded(case, values.copy())
        result = difference_metrics(reference, candidate, "cd")
        self.assertEqual(result["relative_l2"], 0.0)
        self.assertEqual(result["relative_linf"], 0.0)

    def test_high_frequency_diagnostic_separates_step_ringing(self):
        sample = np.arange(2048, dtype=float)
        slow = np.sin(2.0 * math.pi * sample / 512.0)
        ringing = slow + 0.5 * np.power(-1.0, sample)
        self.assertLess(high_frequency_energy_fraction(slow), 0.01)
        self.assertGreater(high_frequency_energy_fraction(ringing), 0.10)


if __name__ == "__main__":
    unittest.main()
