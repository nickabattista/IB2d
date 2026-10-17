#!/usr/bin/env python3
"""Manifest-driven Re=128 translating-plate study for one case per process.

The current validated benchmark is a ramped translating plate.  This module
does not invent the still-unspecified heave/pitch law for the eventual
flapping case.  It provides a clean spatial/periodic-box study and a required
penalty-stiffness/timestep calibration that can later be reused by that case.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import socket
import time
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from case_config import (
    ANGLE_DEGREES,
    NOMINAL_CHORD,
    RAMP_RATE,
    RHO,
    SUPPORT,
    TRANSLATION_SPEED,
    TRANSLATION_STOP_TIME,
    translation_distance,
)
from force_history import HEADER
from generate_case import build_case


STUDY_NAME = "re128_translating_plate"
REYNOLDS_NUMBER = 128.0
MU_RE128 = RHO * abs(TRANSLATION_SPEED) * NOMINAL_CHORD / REYNOLDS_NUMBER
FINAL_TIME = 0.025
PRODUCTION_DT = 2.5e-6
PRODUCTION_TARGET_STIFFNESS = 1.0e9
REFERENCE_CELLS_PER_CHORD = 64.0


@dataclass(frozen=True)
class StudyCase:
    index: int
    case_id: str
    phase: str
    study_axes: tuple[str, ...]
    box_chords: float
    nx: int
    ny: int
    dt: float
    target_stiffness: float

    @property
    def lx(self) -> float:
        return self.box_chords * NOMINAL_CHORD

    @property
    def ly(self) -> float:
        return self.lx

    @property
    def dx(self) -> float:
        return self.lx / self.nx

    @property
    def ds(self) -> float:
        return 0.5 * self.dx

    @property
    def cells_per_chord(self) -> float:
        # Since lx=box_chords*chord, this equivalent ratio avoids a second
        # floating-point division through lx and keeps the fixed-resolution
        # periodic-box family's nominal controls identical.
        return self.nx / self.box_chords

    @property
    def steps(self) -> int:
        raw_steps = FINAL_TIME / self.dt
        steps = int(round(raw_steps))
        if steps < 1 or not math.isclose(
            steps * self.dt, FINAL_TIME, rel_tol=1e-12, abs_tol=1e-15
        ):
            raise ValueError(f"{self.case_id}: final time is not divisible by dt")
        return steps

    @property
    def penalty_coupling_number(self) -> float:
        """Explicit tether coupling k*dt^2*ds/(rho*dx*dy)."""

        return (
            self.target_stiffness * self.dt * self.dt * self.ds
            / (RHO * self.dx * self.dx)
        )


def _production_parameters(
    cells_per_chord: float,
    *,
    reference_dt: float = PRODUCTION_DT,
    reference_stiffness: float = PRODUCTION_TARGET_STIFFNESS,
) -> tuple[float, float]:
    """Return the fixed physical and temporal controls for the grid family.

    The IB target force is a Lagrangian force density and spreading supplies
    the quadrature factor ``ds``.  Keeping ``k`` fixed therefore keeps the
    continuum plate compliance fixed as the marker spacing is refined.  A
    common, independently calibrated ``dt`` makes the grid comparison a
    fixed-timestep spatial study.  Its explicit spread/interpolate tether
    number is checked separately and is largest on the finest grid.
    """

    if (
        float(cells_per_chord) <= 0.0
        or reference_dt <= 0.0
        or reference_stiffness <= 0.0
    ):
        raise ValueError("production parameters must be positive")
    return float(reference_dt), float(reference_stiffness)


def _production_case(index: int, case_id: str, axes: tuple[str, ...],
                     box_chords: float, n: int) -> StudyCase:
    cells_per_chord = n / box_chords
    dt, stiffness = _production_parameters(cells_per_chord)
    return StudyCase(
        index, case_id, "production", axes, box_chords, n, n, dt, stiffness
    )


CALIBRATION_CASES = (
    StudyCase(0, "cal_k5e8_dt1p25e6", "calibration", ("stiffness",),
              16.0, 1024, 1024, 1.25e-6, 5.0e8),
    StudyCase(1, "cal_k1e9_dt1p25e6", "calibration", ("stiffness", "timestep"),
              16.0, 1024, 1024, 1.25e-6, 1.0e9),
    StudyCase(2, "cal_k2e9_dt1p25e6", "calibration", ("stiffness",),
              16.0, 1024, 1024, 1.25e-6, 2.0e9),
    StudyCase(3, "cal_k1e9_dt2p5e6", "calibration", ("timestep",),
              16.0, 1024, 1024, 2.5e-6, 1.0e9),
    StudyCase(4, "cal_k1e9_dt5e6", "calibration", ("timestep",),
              16.0, 1024, 1024, 5.0e-6, 1.0e9),
)


PRODUCTION_CASES = (
    _production_case(0, "grid_L16_P32_N512", ("grid",), 16.0, 512),
    _production_case(1, "boundary_L10_P64_N640", ("boundary",), 10.0, 640),
    _production_case(
        2, "shared_L16_P64_N1024", ("grid", "boundary"), 16.0, 1024
    ),
    _production_case(3, "boundary_L24_P64_N1536", ("boundary",), 24.0, 1536),
    _production_case(4, "grid_L16_P128_N2048", ("grid",), 16.0, 2048),
    _production_case(5, "boundary_L32_P64_N2048", ("boundary",), 32.0, 2048),
)


def cases_for_phase(phase: str, *, dt: float | None = None,
                    target_stiffness: float | None = None) -> tuple[StudyCase, ...]:
    if phase == "calibration":
        if dt is not None or target_stiffness is not None:
            raise ValueError("calibration cases have fixed dt and stiffness values")
        return CALIBRATION_CASES
    if phase != "production":
        raise ValueError(f"unknown study phase: {phase}")

    chosen_dt = PRODUCTION_DT if dt is None else float(dt)
    chosen_stiffness = (
        PRODUCTION_TARGET_STIFFNESS
        if target_stiffness is None else float(target_stiffness)
    )
    if chosen_dt <= 0.0 or chosen_stiffness <= 0.0:
        raise ValueError("production dt and target stiffness must be positive")
    selected_cases = []
    for case in PRODUCTION_CASES:
        case_dt, case_stiffness = _production_parameters(
            case.cells_per_chord,
            reference_dt=chosen_dt,
            reference_stiffness=chosen_stiffness,
        )
        selected_cases.append(
            replace(case, dt=case_dt, target_stiffness=case_stiffness)
        )
    return tuple(selected_cases)


def centered_leading_edge(case: StudyCase) -> tuple[float, float, dict]:
    """Center the complete continuous plate sweep in the periodic box."""

    angle = math.radians(ANGLE_DEGREES)
    plate_dx = NOMINAL_CHORD * math.cos(angle)
    plate_dy = NOMINAL_CHORD * math.sin(angle)
    travel = translation_distance(
        TRANSLATION_STOP_TIME,
        speed=TRANSLATION_SPEED,
        ramp_rate=RAMP_RATE,
        stop_time=TRANSLATION_STOP_TIME,
    )

    relative_x = (0.0, plate_dx, travel, travel + plate_dx)
    relative_y = (0.0, plate_dy)
    x_min, x_max = min(relative_x), max(relative_x)
    y_min, y_max = min(relative_y), max(relative_y)
    swept_width = x_max - x_min
    swept_height = y_max - y_min
    if swept_width >= case.lx or swept_height >= case.ly:
        raise ValueError(f"{case.case_id}: swept plate does not fit in the box")

    x0 = 0.5 * (case.lx - swept_width) - x_min
    y0 = 0.5 * (case.ly - swept_height) - y_min
    placement = {
        "strategy": "center_complete_swept_envelope",
        "translation_distance": travel,
        "swept_width": swept_width,
        "swept_height": swept_height,
        "minimum_x_clearance": 0.5 * (case.lx - swept_width),
        "minimum_y_clearance": 0.5 * (case.ly - swept_height),
    }
    return x0, y0, placement


def case_record(case: StudyCase) -> dict:
    record = asdict(case)
    leading_x, leading_y, placement = centered_leading_edge(case)
    record.update({
        "lx": case.lx,
        "ly": case.ly,
        "dx": case.dx,
        "dy": case.dx,
        "ds": case.ds,
        "cells_per_chord": case.cells_per_chord,
        "steps": case.steps,
        "final_time": FINAL_TIME,
        "reynolds_number": REYNOLDS_NUMBER,
        "rho": RHO,
        "mu": MU_RE128,
        "reference_speed": TRANSLATION_SPEED,
        "nominal_chord": NOMINAL_CHORD,
        "angle_degrees": ANGLE_DEGREES,
        "ramp_rate": RAMP_RATE,
        "translation_stop_time": TRANSLATION_STOP_TIME,
        "ib_kernel_support": SUPPORT,
        "leading_edge": [leading_x, leading_y],
        "placement": placement,
        "advective_cfl": abs(TRANSLATION_SPEED) * case.dt / case.dx,
        "viscous_number": (MU_RE128 / RHO) * case.dt / (case.dx * case.dx),
        "target_stiffness_times_ds": case.target_stiffness * case.ds,
        "penalty_coupling_number": case.penalty_coupling_number,
        "refinement_law": {
            "marker_spacing": "ds=dx/2",
            "grid_family_dt": "constant after independent timestep calibration",
            "grid_family_target_stiffness": (
                "constant continuum target-force-density coefficient"
            ),
            "reference_cells_per_chord": REFERENCE_CELLS_PER_CHORD,
        },
    })
    encoded = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    record["configuration_sha256"] = hashlib.sha256(encoded).hexdigest()
    return record


def write_manifest(path: Path, cases: tuple[StudyCase, ...]) -> Path:
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 2,
        "study_name": STUDY_NAME,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cases": [case_record(case) for case in cases],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def generate_study_case(case: StudyCase, case_dir: Path) -> dict:
    x0, y0, placement = centered_leading_edge(case)
    metadata = build_case(
        case_dir,
        nx=case.nx,
        ny=case.ny,
        dt=case.dt,
        steps=case.steps,
        target_stiffness=case.target_stiffness,
        lx=case.lx,
        ly=case.ly,
        rho=RHO,
        mu=MU_RE128,
        nominal_chord=NOMINAL_CHORD,
        leading_edge_x=x0,
        leading_edge_y=y0,
        angle_degrees=ANGLE_DEGREES,
        translation_speed=TRANSLATION_SPEED,
        ramp_rate=RAMP_RATE,
        translation_stop_time=TRANSLATION_STOP_TIME,
    )
    if not math.isclose(
        metadata["geometry"]["effective_quadrature_chord"],
        NOMINAL_CHORD,
        rel_tol=1e-12,
        abs_tol=1e-15,
    ):
        raise ValueError(f"{case.case_id}: grid does not represent the chord exactly")

    metadata["study"] = case_record(case)
    metadata["study"]["placement"] = placement
    metadata["study"]["status"] = "generated"
    metadata_path = Path(case_dir) / "case_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def validate_history(history_path: Path, case: StudyCase) -> dict:
    with Path(history_path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in HEADER if name not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"history is missing columns: {', '.join(missing)}")
        rows = list(reader)
    if len(rows) != case.steps:
        raise ValueError(
            f"{case.case_id}: expected {case.steps} rows, found {len(rows)}"
        )

    values = {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in HEADER
    }
    if not all(np.all(np.isfinite(column)) for column in values.values()):
        raise FloatingPointError(f"{case.case_id}: non-finite history value")
    expected_steps = np.arange(case.steps, dtype=int)
    if not np.array_equal(values["step"].astype(int), expected_steps):
        raise ValueError(f"{case.case_id}: step sequence is incomplete")
    expected_final_time = (case.steps - 0.5) * case.dt
    if not math.isclose(
        values["time"][-1], expected_final_time, rel_tol=1e-12, abs_tol=1e-14
    ):
        raise ValueError(f"{case.case_id}: unexpected final force timestamp")

    return {
        "rows": len(rows),
        "final_force_time": float(values["time"][-1]),
        "max_force_balance_linf": float(np.max(values["force_balance_linf"])),
        "max_target_error": float(np.max(values["max_target_error"])),
        "max_target_error_over_chord": float(
            np.max(values["max_target_error"]) / NOMINAL_CHORD
        ),
        "max_target_error_over_dx": float(
            np.max(values["max_target_error"]) / case.dx
        ),
    }


def run_study_case(case: StudyCase, output_root: Path) -> Path:
    output_root = Path(output_root).resolve()
    phase_root = output_root / "cases" / case.phase
    final_dir = phase_root / case.case_id
    complete_path = final_dir / "COMPLETE.json"
    if complete_path.is_file():
        completed = json.loads(complete_path.read_text(encoding="utf-8"))
        requested_hash = case_record(case)["configuration_sha256"]
        if completed.get("configuration_sha256") != requested_hash:
            raise FileExistsError(
                f"completed case has a different configuration: {final_dir}; "
                "use a new output root or move the old case aside"
            )
        metadata_path = final_dir / "case_metadata.json"
        history_path = final_dir / "force_history.csv"
        if not metadata_path.is_file() or not history_path.is_file():
            raise FileExistsError(
                f"completed marker is missing required artifacts: {final_dir}"
            )
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("study", {}).get("configuration_sha256") != requested_hash:
            raise FileExistsError(
                f"completed metadata has a different configuration: {final_dir}"
            )
        validate_history(history_path, case)
        print(f"Already complete: {final_dir}")
        return final_dir
    if final_dir.exists():
        raise FileExistsError(
            f"incomplete final directory exists: {final_dir}; move it aside before retrying"
        )

    attempt_id = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if array_id is not None:
        attempt_id += f"-{array_id}"
    attempt_dir = phase_root / f"{case.case_id}.partial-{attempt_id}"
    if attempt_dir.exists():
        raise FileExistsError(f"attempt directory already exists: {attempt_dir}")
    attempt_dir.mkdir(parents=True, exist_ok=False)

    started = datetime.now(timezone.utc)
    wall_start = time.perf_counter()
    metadata = generate_study_case(case, attempt_dir)
    # Importing the legacy driver prints capability messages and allocates
    # module-level state, so defer it until a simulation will actually run.
    from main2d import run_case
    history_path = run_case(attempt_dir)
    diagnostics = validate_history(history_path, case)
    elapsed = time.perf_counter() - wall_start

    metadata["study"]["status"] = "complete"
    metadata["study"]["runtime"] = {
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": elapsed,
        "hostname": socket.gethostname(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": array_id,
    }
    metadata["study"]["completion_diagnostics"] = diagnostics
    (attempt_dir / "case_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (attempt_dir / "COMPLETE.json").write_text(
        json.dumps(metadata["study"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    phase_root.mkdir(parents=True, exist_ok=True)
    attempt_dir.replace(final_dir)
    print(f"Completed: {final_dir}")
    return final_dir


def select_case(cases: tuple[StudyCase, ...], *, index: int | None,
                case_id: str | None) -> StudyCase:
    if (index is None) == (case_id is None):
        raise ValueError("select exactly one of --case-index or --case-id")
    if index is not None:
        matches = [case for case in cases if case.index == index]
    else:
        matches = [case for case in cases if case.case_id == case_id]
    if len(matches) != 1:
        requested = index if index is not None else case_id
        raise ValueError(f"case not found: {requested}")
    return matches[0]


def print_cases(cases: tuple[StudyCase, ...]) -> None:
    print(
        "index,case_id,axes,L_over_c,N,cells_per_chord,dx,ds,dt,steps,"
        "target_stiffness,kT_times_ds,penalty_coupling"
    )
    for case in cases:
        print(
            f"{case.index},{case.case_id},{'+'.join(case.study_axes)},"
            f"{case.box_chords:g},{case.nx},{case.cells_per_chord:g},"
            f"{case.dx:.8g},{case.ds:.8g},{case.dt:.8g},{case.steps},"
            f"{case.target_stiffness:.8g},"
            f"{case.target_stiffness * case.ds:.8g},"
            f"{case.penalty_coupling_number:.8g}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("calibration", "production"),
                        default="production")
    parser.add_argument("--output-root", type=Path,
                        default=Path("runs/re128_study"))
    parser.add_argument("--dt", type=float,
                        help="override the common production timestep")
    parser.add_argument("--target-stiffness", type=float,
                        help="override the common production stiffness")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list", action="store_true")
    action.add_argument("--write-manifest", type=Path)
    action.add_argument("--case-index", type=int)
    action.add_argument("--case-id")
    args = parser.parse_args()

    try:
        cases = cases_for_phase(
            args.phase, dt=args.dt, target_stiffness=args.target_stiffness
        )
        if args.list:
            print_cases(cases)
            return
        if args.write_manifest is not None:
            print(write_manifest(args.write_manifest, cases))
            return
        case = select_case(cases, index=args.case_index, case_id=args.case_id)
        run_study_case(case, args.output_root)
    except (ValueError, FileExistsError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
