#!/usr/bin/env python3
"""Run and compare the two-case Re=64 flapping-plate Puma pilot.

The pilot is intentionally narrow: it checks one paper-based tilted motion at
two timesteps before a kinematic sweep is attempted.  ``show`` is read-only,
``run`` executes one Slurm-array member, and ``analyze`` compares the completed
second-cycle histories and renders the fine-timestep movie.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
import platform
import socket
import sys
import time
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from analyze_forces import analyze_case, derive_coefficients, read_force_history
from case_config import (
    MEAN_ANGLE_OF_ATTACK_DEGREES,
    NOMINAL_CHORD,
    PITCH_AMPLITUDE_DEGREES,
    PITCH_PHASE_RADIANS,
    PIVOT_FRACTION,
    RHO,
    STROKE_AMPLITUDE_OVER_CHORD,
    STROKE_ANGLE_DEGREES,
    SUPPORT,
    TARGET_STIFFNESS,
)
from flapping_visual import (
    FlappingVisualSnapshotWriter,
    build_visual_artifacts,
    coefficient_history,
    write_oscillation_diagnostics,
)
from generate_case import build_case


HERE = Path(__file__).resolve().parent
STUDY_ID = "re64_tilted_flapping_timestep_pilot"
REYNOLDS_NUMBER = 64.0
GRID_POINTS = 512
BOX_CHORDS = 16.0
REQUESTED_CYCLES = 2
SNAPSHOTS_PER_CYCLE = 20
ANALYSIS_CYCLE = 2
COMPARISON_SAMPLES = 1441
# This differs from pi*(0.25)/24 by only 7.65e-8 seconds (2.34e-6 relatively),
# but makes each cycle contain exactly 26,180 coarse or 52,360 fine steps.
# build_case derives Umax and mu from this value, so Re remains exactly 64.
PILOT_PERIOD = 0.032725
PILOT_STROKE_AMPLITUDE = 2.5 * NOMINAL_CHORD
PILOT_PEAK_SPEED = math.pi * PILOT_STROKE_AMPLITUDE / PILOT_PERIOD
PILOT_DYNAMIC_VISCOSITY = (
    RHO * PILOT_PEAK_SPEED * NOMINAL_CHORD / REYNOLDS_NUMBER
)


@dataclass(frozen=True)
class PilotCase:
    index: int
    case_id: str
    dt: float
    steps: int
    snapshot_interval: int
    nominal_cycles: int = REQUESTED_CYCLES
    nx: int = GRID_POINTS
    ny: int = GRID_POINTS
    lx: float = BOX_CHORDS * NOMINAL_CHORD
    ly: float = BOX_CHORDS * NOMINAL_CHORD
    target_stiffness: float = TARGET_STIFFNESS

    @property
    def tfinal(self) -> float:
        return self.steps * self.dt

    @property
    def actual_cycles(self) -> float:
        return self.tfinal / PILOT_PERIOD

    @property
    def dx(self) -> float:
        return self.lx / self.nx

    @property
    def ds(self) -> float:
        return 0.5 * self.dx

    @property
    def tether_number(self) -> float:
        """Return kT*dt^2*ds/(rho*dx^2), a stiffness diagnostic."""

        return (
            self.target_stiffness * self.dt * self.dt * self.ds
            / (RHO * self.dx * self.dx)
        )

    def record(self) -> dict[str, object]:
        result = asdict(self)
        result.update({
            "tfinal": self.tfinal,
            "actual_cycles": self.actual_cycles,
            "period": PILOT_PERIOD,
            "dx": self.dx,
            "ds": self.ds,
            "ds_over_dx": self.ds / self.dx,
            "cells_per_chord": NOMINAL_CHORD / self.dx,
            "tether_number": self.tether_number,
        })
        return result


def _make_case(index: int, case_id: str, dt: float) -> PilotCase:
    exact_steps = REQUESTED_CYCLES * PILOT_PERIOD / dt
    steps = int(round(exact_steps))
    if not math.isclose(exact_steps, steps, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError("pilot period must contain an integer number of steps")
    cadence = max(
        1, int(round(PILOT_PERIOD / dt / SNAPSHOTS_PER_CYCLE))
    )
    return PilotCase(index, case_id, dt, steps, cadence)


PILOT_CASES = (
    _make_case(0, "re64_tilted_n512_dt1p25e-6", 1.25e-6),
    _make_case(1, "re64_tilted_n512_dt6p25e-7", 6.25e-7),
)


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _physics_controls() -> dict[str, object]:
    return {
        "reynolds_number": REYNOLDS_NUMBER,
        "nominal_chord": NOMINAL_CHORD,
        "rho": RHO,
        "period": PILOT_PERIOD,
        "peak_translation_speed": PILOT_PEAK_SPEED,
        "dynamic_viscosity": PILOT_DYNAMIC_VISCOSITY,
        "reynolds_definition": "rho*Umax*nominal_chord/mu",
        "requested_cycles": REQUESTED_CYCLES,
    }


def _motion_controls() -> dict[str, object]:
    return {
        "stroke_amplitude": PILOT_STROKE_AMPLITUDE,
        "stroke_amplitude_over_chord": STROKE_AMPLITUDE_OVER_CHORD,
        "stroke_angle_degrees": STROKE_ANGLE_DEGREES,
        "mean_angle_of_attack_degrees": MEAN_ANGLE_OF_ATTACK_DEGREES,
        "pitch_amplitude_degrees": PITCH_AMPLITUDE_DEGREES,
        "pitch_phase_radians": PITCH_PHASE_RADIANS,
        "pivot_fraction": PIVOT_FRACTION,
        "angle_convention": "theta=beta-alpha, leading_edge_to_trailing_edge",
        "evaluation_time": "step_midpoint",
    }


def _case_definition(case: PilotCase) -> dict[str, object]:
    """Return every physical and numerical input that determines a run."""

    return {
        "case": case.record(),
        "physics": _physics_controls(),
        "motion": _motion_controls(),
        "coupling": {
            "delta_support": SUPPORT,
            "lagrangian_spacing_rule": "ds=dx/2",
        },
    }


def _case_configuration_sha256(case: PilotCase) -> str:
    return _sha256_json(_case_definition(case))


def pilot_manifest() -> dict[str, object]:
    """Return the complete planned experiment without touching the filesystem."""

    cases = [case.record() for case in PILOT_CASES]
    manifest = {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "purpose": (
            "Implementation, target-following, and timestep check for the exact "
            "tilted Re=64 flapping case; not a spatial-convergence result"
        ),
        "physics": _physics_controls(),
        "motion": _motion_controls(),
        "fixed_controls": {
            "grid": [GRID_POINTS, GRID_POINTS],
            "box_chords": BOX_CHORDS,
            "target_stiffness": TARGET_STIFFNESS,
            "delta_support": SUPPORT,
            "lagrangian_spacing_rule": "ds=dx/2",
            "snapshots_per_cycle_requested": SNAPSHOTS_PER_CYCLE,
        },
        "cases": cases,
        "analysis": {
            "case_cycle_summaries": [1, 2],
            "cross_timestep_comparison_cycle": ANALYSIS_CYCLE,
            "cross_timestep_phase_samples": COMPARISON_SAMPLES,
            "movie_case_index": 1,
            "scientific_status": "human_review_required",
        },
    }
    manifest["configuration_sha256"] = _sha256_json(manifest)
    return manifest


def _case_for_index(index: int) -> PilotCase:
    for case in PILOT_CASES:
        if case.index == int(index):
            return case
    raise ValueError(f"case index must be one of {[case.index for case in PILOT_CASES]}")


def _fingerprints(paths: Sequence[Path]) -> dict[str, str]:
    project_root = HERE.parents[2]
    return {
        str(path.relative_to(project_root)).replace("\\", "/"):
            hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(paths)
    }


def _numerical_source_fingerprints() -> dict[str, str]:
    candidates = [
        HERE / "case_config.py",
        HERE / "generate_case.py",
        HERE / "kinematics.py",
        HERE / "update_Target_Point_Positions.py",
        HERE / "force_history.py",
        HERE / "main2d.py",
    ]
    candidates.extend(sorted((HERE.parents[1] / "IBM_Blackbox").glob("*.py")))
    return _fingerprints(candidates)


def _snapshot_producer_source_fingerprints() -> dict[str, str]:
    """Hash code that determines the stored snapshot contents."""

    candidates = [
        HERE / "flapping_visual.py",
        HERE.parent / "Translating_Plate" / "visual_pilot.py",
    ]
    return _fingerprints(candidates)


def _postprocessor_source_fingerprints() -> dict[str, str]:
    """Hash code that derives force, comparison, and rendered artifacts."""

    candidates = [
        HERE / "run_re64_flapping_pilot.py",
        HERE / "analyze_forces.py",
        HERE / "flapping_visual.py",
        HERE.parent / "Translating_Plate" / "visual_pilot.py",
    ]
    return _fingerprints(candidates)


def _read_json(path: Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _write_json(path: Path, value: object) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def _portable_paths(value: object, root: Path) -> object:
    """Convert artifact paths under ``root`` to download-safe relative paths."""

    if isinstance(value, dict):
        return {key: _portable_paths(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable_paths(item, root) for item in value]
    if isinstance(value, str):
        candidate = Path(value)
        if candidate.is_absolute():
            try:
                return candidate.relative_to(root).as_posix()
            except ValueError:
                return value
    return value


def _validate_completed_case(case_dir: Path, case: PilotCase) -> dict:
    completion_path = case_dir / "COMPLETE.json"
    if not completion_path.is_file():
        raise FileNotFoundError(f"missing completion marker: {completion_path}")
    completion = _read_json(completion_path)
    expected_hash = _case_configuration_sha256(case)
    if completion.get("case_configuration_sha256") != expected_hash:
        raise ValueError(f"completed case configuration does not match: {case_dir}")
    metadata = _read_json(case_dir / "case_metadata.json")
    pilot = metadata.get("pilot")
    if not isinstance(pilot, Mapping):
        raise ValueError(f"case metadata has no pilot provenance: {case_dir}")
    if pilot.get("case_configuration_sha256") != expected_hash:
        raise ValueError(f"case metadata configuration is stale: {case_dir}")
    current_numerical_sources = _numerical_source_fingerprints()
    if pilot.get("numerical_source_sha256") != current_numerical_sources:
        raise ValueError(
            f"completed case was produced by different numerical source: {case_dir}"
        )
    if completion.get("numerical_source_set_sha256") != _sha256_json(
        current_numerical_sources
    ):
        raise ValueError(f"completion source fingerprint is stale: {case_dir}")
    current_snapshot_sources = _snapshot_producer_source_fingerprints()
    if pilot.get("snapshot_producer_source_sha256") != current_snapshot_sources:
        raise ValueError(
            "completed case was produced by different snapshot-producing "
            f"source: {case_dir}"
        )
    if completion.get("snapshot_producer_source_set_sha256") != _sha256_json(
        current_snapshot_sources
    ):
        raise ValueError(
            f"completion snapshot-source fingerprint is stale: {case_dir}"
        )
    if int(metadata["time"]["steps"]) != case.steps:
        raise ValueError(f"step count does not match pilot definition: {case_dir}")
    history = read_force_history(case_dir / "force_history.csv")
    if history["time"].size != case.steps:
        raise ValueError(
            f"force history has {history['time'].size} rows; expected {case.steps}"
        )
    snapshot_steps = _read_snapshot_steps(case_dir / "visual_snapshots")
    expected_snapshot_steps = _expected_snapshot_steps(case)
    if snapshot_steps != expected_snapshot_steps:
        raise ValueError(
            f"snapshot cadence is incomplete in {case_dir}: "
            f"got {snapshot_steps}, expected {expected_snapshot_steps}"
        )
    return completion


def _expected_snapshot_steps(case: PilotCase) -> list[int]:
    steps = list(range(0, case.steps + 1, case.snapshot_interval))
    if not steps or steps[-1] != case.steps:
        steps.append(case.steps)
    return steps


def _read_snapshot_steps(snapshot_dir: Path) -> list[int]:
    paths = sorted(Path(snapshot_dir).glob("snapshot_step*.npz"))
    values = []
    for path in paths:
        with np.load(path, allow_pickle=False) as snapshot:
            values.append(int(snapshot["step"]))
    return sorted(values)


def run_pilot_case(case: PilotCase, output_root: Path) -> Path:
    """Generate and execute one case, publishing it only after validation."""

    output_root = Path(output_root).resolve()
    cases_root = output_root / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    final_dir = cases_root / case.case_id
    if final_dir.exists():
        _validate_completed_case(final_dir, case)
        print(f"Case already complete and validated: {final_dir}")
        return final_dir

    partial_dir = cases_root / (
        f".{case.case_id}.partial-{os.environ.get('SLURM_JOB_ID', os.getpid())}-"
        f"{uuid.uuid4().hex[:8]}"
    )
    partial_dir.mkdir(parents=False, exist_ok=False)
    started = time.time()
    metadata = build_case(
        partial_dir,
        nx=case.nx,
        ny=case.ny,
        dt=case.dt,
        steps=case.steps,
        target_stiffness=case.target_stiffness,
        lx=case.lx,
        ly=case.ly,
        reynolds_number=REYNOLDS_NUMBER,
        period=PILOT_PERIOD,
        stroke_amplitude_over_chord=STROKE_AMPLITUDE_OVER_CHORD,
        stroke_angle_degrees=STROKE_ANGLE_DEGREES,
        mean_angle_of_attack_degrees=MEAN_ANGLE_OF_ATTACK_DEGREES,
        pitch_amplitude_degrees=PITCH_AMPLITUDE_DEGREES,
        pitch_phase_radians=PITCH_PHASE_RADIANS,
    )
    numerical_sources = _numerical_source_fingerprints()
    snapshot_sources = _snapshot_producer_source_fingerprints()
    metadata["case_name"] = case.case_id
    metadata["pilot"] = {
        "study_id": STUDY_ID,
        "case_index": case.index,
        "case_configuration": _case_definition(case),
        "case_configuration_sha256": _case_configuration_sha256(case),
        "numerical_source_sha256": numerical_sources,
        "snapshot_producer_source_sha256": snapshot_sources,
        "postprocessor_source_sha256_at_run": (
            _postprocessor_source_fingerprints()
        ),
    }
    _write_json(partial_dir / "case_metadata.json", metadata)

    snapshot_writer = FlappingVisualSnapshotWriter.from_case(
        partial_dir, every_steps=case.snapshot_interval
    )
    # Delay the heavy solver import until the run action is actually requested.
    from main2d import run_case

    history_path = run_case(
        partial_dir,
        state_callback=snapshot_writer.record,
        state_callback_interval=case.snapshot_interval,
    )
    history = read_force_history(history_path)
    if history["time"].size != case.steps:
        raise RuntimeError(
            f"solver wrote {history['time'].size} force rows; expected {case.steps}"
        )
    snapshot_steps = _read_snapshot_steps(partial_dir / "visual_snapshots")
    expected_snapshot_steps = _expected_snapshot_steps(case)
    if snapshot_steps != expected_snapshot_steps:
        raise RuntimeError(
            f"snapshot steps are {snapshot_steps}; expected {expected_snapshot_steps}"
        )

    elapsed = time.time() - started
    completion = {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "case_id": case.case_id,
        "case_index": case.index,
        "case_configuration_sha256": _case_configuration_sha256(case),
        "numerical_source_set_sha256": _sha256_json(numerical_sources),
        "snapshot_producer_source_set_sha256": _sha256_json(snapshot_sources),
        "force_row_count": int(history["time"].size),
        "snapshot_count": len(snapshot_steps),
        "snapshot_steps": snapshot_steps,
        "all_core_force_values_finite": bool(
            all(np.all(np.isfinite(history[name])) for name in ("time", "fx_body", "fy_body"))
        ),
        "elapsed_seconds": elapsed,
        "runtime": {
            "hostname": socket.gethostname(),
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        },
    }
    _write_json(partial_dir / "COMPLETE.json", completion)
    if final_dir.exists():
        raise FileExistsError(f"case appeared while this run was active: {final_dir}")
    partial_dir.rename(final_dir)
    print(f"Completed case: {final_dir}")
    return final_dir


COMPARISON_FIELDS = (
    "ch_umax", "cv_umax", "ch_urms", "cv_urms",
    "cd_umax", "cl_umax", "cd_urms", "cl_urms",
)


def _cycle_series(
    curves: Mapping[str, np.ndarray], cycle_number: int, field: str
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.asarray(curves["cycle"], dtype=int) == int(cycle_number)
    phase = np.asarray(curves["phase"], dtype=float)[mask]
    values = np.asarray(curves[field], dtype=float)[mask]
    finite = np.isfinite(phase) & np.isfinite(values)
    if np.count_nonzero(finite) < 2:
        raise ValueError(f"not enough finite cycle-{cycle_number} samples for {field}")
    phase, values = phase[finite], values[finite]
    order = np.argsort(phase)
    return phase[order], values[order]


def compare_timestep_cases(
    case_dirs: Sequence[Path], cycle_number: int = ANALYSIS_CYCLE
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    """Compare coarse and fine histories on a shared phase grid."""

    if len(case_dirs) != 2:
        raise ValueError("exactly two case directories are required")
    curve_sets = []
    metadata_sets = []
    for case_dir in case_dirs:
        metadata = _read_json(Path(case_dir) / "case_metadata.json")
        history = read_force_history(Path(case_dir) / "force_history.csv")
        metadata_sets.append(metadata)
        curve_sets.append(derive_coefficients(history, metadata))

    ranges = []
    for curves in curve_sets:
        mask = np.asarray(curves["cycle"], dtype=int) == cycle_number
        phase = np.asarray(curves["phase"], dtype=float)[mask]
        if phase.size < 2:
            raise ValueError(f"cycle {cycle_number} is unavailable")
        ranges.append((float(np.min(phase)), float(np.max(phase))))
    low = max(value[0] for value in ranges)
    high = min(value[1] for value in ranges)
    if high <= low:
        raise ValueError("completed cases have no common phase interval")
    common_phase = np.linspace(low, high, COMPARISON_SAMPLES)
    reversal_exclusion_half_width = max(
        float(metadata["time"]["dt"]) / float(metadata["motion"]["period"])
        for metadata in metadata_sets
    )
    distance_to_reversal = np.minimum.reduce((
        np.abs(common_phase),
        np.abs(common_phase - 0.5),
        np.abs(common_phase - 1.0),
    ))
    reversal_exclusion = distance_to_reversal <= reversal_exclusion_half_width

    table: dict[str, np.ndarray] = {"phase": common_phase}
    metrics: dict[str, object] = {}
    for field in COMPARISON_FIELDS:
        phase_coarse, coarse = _cycle_series(curve_sets[0], cycle_number, field)
        phase_fine, fine = _cycle_series(curve_sets[1], cycle_number, field)
        coarse_i = np.interp(common_phase, phase_coarse, coarse)
        fine_i = np.interp(common_phase, phase_fine, fine)
        if field.startswith(("cd_", "cl_")):
            coarse_i[reversal_exclusion] = np.nan
            fine_i[reversal_exclusion] = np.nan
        difference = coarse_i - fine_i
        finite = np.isfinite(coarse_i) & np.isfinite(fine_i)
        if np.count_nonzero(finite) < 2:
            raise ValueError(f"not enough non-reversal comparison points for {field}")
        coarse_valid = coarse_i[finite]
        fine_valid = fine_i[finite]
        difference_valid = difference[finite]
        fine_rms = float(np.sqrt(np.mean(fine_valid * fine_valid)))
        rms_difference = float(
            np.sqrt(np.mean(difference_valid * difference_valid))
        )
        table[f"coarse_{field}"] = coarse_i
        table[f"fine_{field}"] = fine_i
        table[f"coarse_minus_fine_{field}"] = difference
        metrics[field] = {
            "rms_difference": rms_difference,
            "relative_l2_using_fine": (
                rms_difference / fine_rms if fine_rms > 0.0 else None
            ),
            "maximum_absolute_difference": float(
                np.max(np.abs(difference_valid))
            ),
            "coarse_mean": float(np.mean(coarse_valid)),
            "fine_mean": float(np.mean(fine_valid)),
            "comparison_sample_count": int(np.count_nonzero(finite)),
        }
    comparison = {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "scientific_status": "human_review_required",
        "interpretation": (
            "Coarse-minus-fine differences after phase interpolation over cycle 2. "
            "Motion-relative CD/CL omit a documented window around each "
            "direction reversal; CH/CV remain defined there. "
            "These diagnostics do not constitute spatial convergence or an "
            "automatic physics-validation pass."
        ),
        "cycle": cycle_number,
        "common_phase_interval": [low, high],
        "phase_sample_count": COMPARISON_SAMPLES,
        "motion_relative_reversal_exclusion_half_width_phase": (
            reversal_exclusion_half_width
        ),
        "excluded_common_phase_sample_count": int(
            np.count_nonzero(reversal_exclusion)
        ),
        "coarse_case": PILOT_CASES[0].record(),
        "fine_case": PILOT_CASES[1].record(),
        "series": metrics,
    }
    return comparison, table


def _write_comparison_csv(path: Path, table: Mapping[str, np.ndarray]) -> Path:
    names = tuple(table)
    rows = zip(*(np.asarray(table[name]) for name in names))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in rows:
            writer.writerow(
                "" if not math.isfinite(float(value))
                else format(float(value), ".17g")
                for value in row
            )
    return path


def _write_comparison_svg(path: Path, table: Mapping[str, np.ndarray]) -> Path:
    phase = np.asarray(table["phase"], dtype=float)
    panels = (("ch_urms", "CH (Urms normalization)"),
              ("cv_urms", "CV (Urms normalization)"))
    elements = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="760" viewBox="0 0 1100 760">',
        '<rect width="100%" height="100%" fill="#fafafa"/>',
        '<style>.title{font:bold 22px sans-serif}.label{font:15px sans-serif}.legend{font:14px sans-serif}</style>',
        '<text x="550" y="32" text-anchor="middle" class="title">Re=64 flapping pilot: cycle-2 timestep comparison</text>',
    ]
    for panel_index, (field, label) in enumerate(panels):
        box = (80.0, 65.0 + panel_index * 330.0, 970.0, 270.0)
        coarse = np.asarray(table[f"coarse_{field}"], dtype=float)
        fine = np.asarray(table[f"fine_{field}"], dtype=float)
        all_values = np.concatenate((coarse, fine))
        # Give both lines in a panel the same range.
        ymin, ymax = float(np.min(all_values)), float(np.max(all_values))
        if math.isclose(ymin, ymax):
            ymin, ymax = ymin - 1.0, ymax + 1.0
        left, top, width, height = box
        def points(values: np.ndarray) -> str:
            px = left + width * (phase - phase[0]) / (phase[-1] - phase[0])
            py = top + height * (ymax - values) / (ymax - ymin)
            return " ".join(f"{x:.2f},{y:.2f}" for x, y in zip(px, py))
        elements.extend((
            f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="white" stroke="#777"/>',
            f'<polyline points="{points(coarse)}" fill="none" stroke="#d95f02" stroke-width="1.5"/>',
            f'<polyline points="{points(fine)}" fill="none" stroke="#1b9e77" stroke-width="1.5"/>',
            f'<text x="{left}" y="{top - 10}" class="label">{html.escape(label)}</text>',
            f'<text x="{left + 8}" y="{top + 20}" class="legend" fill="#d95f02">coarse dt</text>',
            f'<text x="{left + 100}" y="{top + 20}" class="legend" fill="#1b9e77">fine dt</text>',
            f'<text x="{left - 8}" y="{top + 12}" text-anchor="end" class="legend">{ymax:.4g}</text>',
            f'<text x="{left - 8}" y="{top + height}" text-anchor="end" class="legend">{ymin:.4g}</text>',
        ))
    elements.extend((
        '<text x="565" y="730" text-anchor="middle" class="label">cycle phase</text>',
        '</svg>',
    ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path


def _write_analysis_manifest(
    summary_dir: Path, result: Mapping[str, object], output_root: Path
) -> dict[str, object]:
    """Stamp current postprocessor provenance and write a portable manifest."""

    postprocessor_sources = _postprocessor_source_fingerprints()
    manifest = dict(result)
    manifest["postprocessor_source_sha256"] = postprocessor_sources
    manifest["postprocessor_source_set_sha256"] = _sha256_json(
        postprocessor_sources
    )
    portable_manifest = _portable_paths(manifest, Path(output_root).resolve())
    _write_json(
        Path(summary_dir) / "pilot_analysis_manifest.json", portable_manifest
    )
    return portable_manifest


def analyze_pilot(output_root: Path, *, make_mp4: bool = True, fps: float = 12.0) -> dict:
    output_root = Path(output_root).resolve()
    case_dirs = [output_root / "cases" / case.case_id for case in PILOT_CASES]
    completions = [
        _validate_completed_case(case_dir, case)
        for case_dir, case in zip(case_dirs, PILOT_CASES)
    ]
    per_case = {}
    for case, case_dir in zip(PILOT_CASES, case_dirs):
        case_analysis = analyze_case(case_dir, analysis_cycles=(1, 2))
        metadata = _read_json(case_dir / "case_metadata.json")
        history = read_force_history(case_dir / "force_history.csv")
        fixed_frame_curves = coefficient_history(history, metadata)
        ringing_path = case_dir / "force_analysis" / (
            "full_resolution_oscillation_diagnostics.json"
        )
        ringing_report = write_oscillation_diagnostics(
            ringing_path, fixed_frame_curves, metadata
        )
        case_analysis["artifacts"]["full_resolution_oscillation_diagnostics"] = (
            str(ringing_path)
        )
        case_analysis["full_resolution_oscillation_diagnostics"] = ringing_report
        per_case[case.case_id] = case_analysis

    summary_dir = output_root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    comparison, table = compare_timestep_cases(case_dirs)
    comparison_csv = _write_comparison_csv(
        summary_dir / "cycle2_timestep_comparison.csv", table
    )
    comparison_svg = _write_comparison_svg(
        summary_dir / "cycle2_timestep_comparison.svg", table
    )
    comparison_path = _write_json(
        summary_dir / "cycle2_timestep_comparison.json", comparison
    )
    fine_visual_dir = case_dirs[1] / "visual_pilot"
    fine_visual = build_visual_artifacts(
        case_dirs[1],
        output_dir=fine_visual_dir,
        fps=fps,
        make_mp4=make_mp4,
    )
    result = {
        "schema_version": 1,
        "study_id": STUDY_ID,
        "scientific_status": "human_review_required",
        "case_completions": completions,
        "per_case_force_analysis": per_case,
        "cross_timestep_comparison": {
            "json": str(comparison_path),
            "csv": str(comparison_csv),
            "svg": str(comparison_svg),
        },
        "fine_case_visual": {
            "base_directory": str(fine_visual_dir),
            "artifacts": fine_visual,
        },
    }
    return _write_analysis_manifest(summary_dir, result, output_root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("show", help="print the exact planned configuration")
    run_parser = subparsers.add_parser("run", help="run one Slurm-array case")
    run_parser.add_argument("--case-index", type=int, required=True)
    run_parser.add_argument("--output-root", type=Path, required=True)
    analyze_parser = subparsers.add_parser("analyze", help="compare completed cases")
    analyze_parser.add_argument("--output-root", type=Path, required=True)
    analyze_parser.add_argument("--fps", type=float, default=12.0)
    analyze_parser.add_argument("--no-mp4", action="store_true")
    args = parser.parse_args()

    if args.action == "show":
        print(json.dumps(pilot_manifest(), indent=2, sort_keys=True))
    elif args.action == "run":
        run_pilot_case(_case_for_index(args.case_index), args.output_root)
    else:
        result = analyze_pilot(
            args.output_root, make_mp4=not args.no_mp4, fps=args.fps
        )
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
