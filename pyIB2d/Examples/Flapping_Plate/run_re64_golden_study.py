#!/usr/bin/env python3
"""Build, run, and analyze the staged Re=64 flapping-plate study.

The workflow is manifest driven.  ``manifest`` creates an immutable description
of one stage, ``run`` executes one indexed case (suitable for a Slurm array), and
``analyze`` summarizes every completed case while tolerating failed or missing
array members.  Structural dampers are not used.  By default, the fluid update
uses an explicitly documented periodic-box mean-momentum constraint while force
output retains the original immersed-boundary reaction.
"""

from __future__ import annotations

import argparse
import copy
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
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from analyze_forces import analyze_case, read_force_history
from case_config import NOMINAL_CHORD, RHO, SUPPORT
from generate_case import build_case


HERE = Path(__file__).resolve().parent
SCHEMA_VERSION = 1
REYNOLDS_NUMBER = 64.0
PERIOD = 0.032725
STROKE_AMPLITUDE_OVER_CHORD = 2.5
STROKE_AMPLITUDE = STROKE_AMPLITUDE_OVER_CHORD * NOMINAL_CHORD
U_REF = math.pi * STROKE_AMPLITUDE / PERIOD
DYNAMIC_VISCOSITY = RHO * U_REF * NOMINAL_CHORD / REYNOLDS_NUMBER
MEAN_ANGLE_OF_ATTACK_DEGREES = 60.0
PITCH_PHASE_RADIANS = 0.0
BASELINE_BETA_DEGREES = 62.8
BASELINE_PITCH_AMPLITUDE_DEGREES = 45.0

DEFAULT_TARGET_STIFFNESS = 1.0e9
DEFAULT_DT = 6.25e-7
DEFAULT_BOX_CHORDS = 16.0
DEFAULT_CELLS_PER_CHORD = 32
DEFAULT_STIFFNESS_VALUES = (2.5e8, 5.0e8, 1.0e9, 2.0e9)
DEFAULT_DT_VALUES = (2.5e-6, 1.25e-6, 6.25e-7)
DEFAULT_BOX_VALUES = (12.0, 16.0, 24.0)
DEFAULT_GRID_VALUES = (16, 32, 64)
DEFAULT_BETA_VALUES = (0.0, 30.0, 62.8)
DEFAULT_PITCH_VALUES = (30.0, 45.0)
DEFAULT_CYCLES = {
    "stiffness": 1,
    "timestep": 2,
    "box": 4,
    "grid": 2,
    "startup": 8,
    "production": 8,
}
STAGES = tuple(DEFAULT_CYCLES)


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256_json(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _without_key(value: Mapping[str, object], key: str) -> dict[str, object]:
    result = copy.deepcopy(dict(value))
    result.pop(key, None)
    return result


def _number_label(value: float) -> str:
    """Return a compact filesystem-safe label for a finite number."""

    number = float(value)
    if not math.isfinite(number):
        raise ValueError("case labels require finite values")
    label = format(number, ".8g").lower().replace("-", "m").replace("+", "")
    return label.replace(".", "p")


def _positive_values(values: Sequence[float], name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if not result or not all(math.isfinite(value) and value > 0.0 for value in result):
        raise ValueError(f"{name} must contain finite positive values")
    return result


def _cycle_aligned_dt(nominal_dt: float) -> tuple[float, int]:
    """Return a nearby dt with an integer number of steps per cycle."""

    if not math.isfinite(float(nominal_dt)) or nominal_dt <= 0.0:
        raise ValueError("dt must be finite and positive")
    steps_per_cycle = max(1, int(round(PERIOD / float(nominal_dt))))
    return PERIOD / steps_per_cycle, steps_per_cycle


def _case_record(
    *,
    index: int,
    stage: str,
    target_stiffness: float,
    nominal_dt: float,
    box_chords: float,
    cells_per_chord: int,
    cycles: int,
    beta_degrees: float = BASELINE_BETA_DEGREES,
    pitch_amplitude_degrees: float = BASELINE_PITCH_AMPLITUDE_DEGREES,
    balance_periodic_force: bool = True,
    save_visual_snapshots: bool = False,
) -> dict[str, object]:
    if stage not in STAGES:
        raise ValueError(f"unknown stage {stage!r}")
    target_stiffness = _positive_values((target_stiffness,), "target stiffness")[0]
    box_chords = _positive_values((box_chords,), "box size")[0]
    if int(cells_per_chord) != cells_per_chord or cells_per_chord < 4:
        raise ValueError("cells_per_chord must be an integer of at least 4")
    if int(cycles) != cycles or cycles < 1:
        raise ValueError("cycles must be a positive integer")
    cells_per_chord = int(cells_per_chord)
    cycles = int(cycles)
    grid_points_float = box_chords * cells_per_chord
    grid_points = int(round(grid_points_float))
    if not math.isclose(grid_points_float, grid_points, abs_tol=1.0e-10):
        raise ValueError("box_chords*cells_per_chord must be an integer")
    if grid_points < 8 or grid_points % 2:
        raise ValueError("the resulting square grid must have an even size of at least 8")
    dt, steps_per_cycle = _cycle_aligned_dt(nominal_dt)
    lx = box_chords * NOMINAL_CHORD
    dx = lx / grid_points
    ds = 0.5 * dx
    steps = cycles * steps_per_cycle
    case_id = (
        f"{stage}_k{_number_label(target_stiffness)}_dt{_number_label(dt)}_"
        f"l{_number_label(box_chords)}c_p{cells_per_chord}_"
        f"b{_number_label(beta_degrees)}_a{_number_label(pitch_amplitude_degrees)}_"
        f"cy{cycles}"
    )
    record: dict[str, object] = {
        "index": int(index),
        "case_id": case_id,
        "stage": stage,
        "nx": grid_points,
        "ny": grid_points,
        "lx": lx,
        "ly": lx,
        "box_chords": box_chords,
        "cells_per_chord": cells_per_chord,
        "dx": dx,
        "ds": ds,
        "ds_over_dx": ds / dx,
        "nominal_dt": float(nominal_dt),
        "dt": dt,
        "steps_per_cycle": steps_per_cycle,
        "cycles": cycles,
        "steps": steps,
        "tfinal": steps * dt,
        "target_stiffness": target_stiffness,
        "damped_springs": False,
        "target_damping": 0.0,
        "balance_periodic_force": bool(balance_periodic_force),
        "beta_degrees": float(beta_degrees),
        "mean_angle_of_attack_degrees": MEAN_ANGLE_OF_ATTACK_DEGREES,
        "pitch_amplitude_degrees": float(pitch_amplitude_degrees),
        "pitch_phase_radians": PITCH_PHASE_RADIANS,
        "save_visual_snapshots": bool(save_visual_snapshots),
        "visual_frames_per_cycle": 20 if save_visual_snapshots else 0,
        "tether_number": target_stiffness * dt * dt * ds / (RHO * dx * dx),
        "translation_cfl": U_REF * dt / dx,
        "viscous_number": (DYNAMIC_VISCOSITY / RHO) * dt / (dx * dx),
    }
    record["case_configuration_sha256"] = _sha256_json(record)
    return record


def design_stage(
    stage: str,
    *,
    target_stiffness: float = DEFAULT_TARGET_STIFFNESS,
    dt: float = DEFAULT_DT,
    box_chords: float = DEFAULT_BOX_CHORDS,
    cells_per_chord: int = DEFAULT_CELLS_PER_CHORD,
    cycles: int | None = None,
    stiffness_values: Sequence[float] = DEFAULT_STIFFNESS_VALUES,
    dt_values: Sequence[float] = DEFAULT_DT_VALUES,
    box_values: Sequence[float] = DEFAULT_BOX_VALUES,
    grid_values: Sequence[int] = DEFAULT_GRID_VALUES,
    beta_values: Sequence[float] = DEFAULT_BETA_VALUES,
    pitch_values: Sequence[float] = DEFAULT_PITCH_VALUES,
    grid_dt_rule: str = "constant",
    balance_periodic_force: bool = True,
) -> dict[str, object]:
    """Return a deterministic, JSON-compatible manifest for one stage."""

    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}")
    cycles = DEFAULT_CYCLES[stage] if cycles is None else int(cycles)
    baseline = {
        "target_stiffness": float(target_stiffness),
        "dt": float(dt),
        "box_chords": float(box_chords),
        "cells_per_chord": int(cells_per_chord),
        "cycles": cycles,
    }
    specs: list[dict[str, object]] = []
    if stage == "stiffness":
        for stiffness in _positive_values(stiffness_values, "stiffness_values"):
            specs.append({**baseline, "target_stiffness": stiffness})
    elif stage == "timestep":
        for candidate_dt in _positive_values(dt_values, "dt_values"):
            specs.append({**baseline, "dt": candidate_dt})
    elif stage == "box":
        for candidate_box in _positive_values(box_values, "box_values"):
            specs.append({**baseline, "box_chords": candidate_box})
    elif stage == "grid":
        if grid_dt_rule not in {"constant_tether", "constant"}:
            raise ValueError("grid_dt_rule must be 'constant_tether' or 'constant'")
        for candidate_grid in grid_values:
            if int(candidate_grid) != candidate_grid or int(candidate_grid) < 4:
                raise ValueError("grid_values must contain integers of at least 4")
            candidate_grid = int(candidate_grid)
            candidate_dt = float(dt)
            if grid_dt_rule == "constant_tether":
                candidate_dt *= math.sqrt(float(cells_per_chord) / candidate_grid)
            specs.append({
                **baseline,
                "cells_per_chord": candidate_grid,
                "dt": candidate_dt,
            })
    elif stage == "startup":
        specs.append(baseline)
    else:
        for beta in beta_values:
            for pitch in pitch_values:
                if not (math.isfinite(float(beta)) and math.isfinite(float(pitch))):
                    raise ValueError("production angles must be finite")
                specs.append({
                    **baseline,
                    "beta_degrees": float(beta),
                    "pitch_amplitude_degrees": float(pitch),
                    "save_visual_snapshots": (
                        math.isclose(float(beta), BASELINE_BETA_DEGREES)
                        and math.isclose(
                            float(pitch), BASELINE_PITCH_AMPLITUDE_DEGREES
                        )
                    ),
                })
        if specs and not any(bool(spec.get("save_visual_snapshots")) for spec in specs):
            specs[-1]["save_visual_snapshots"] = True

    cases = []
    for index, spec in enumerate(specs):
        cases.append(_case_record(
            index=index,
            stage=stage,
            target_stiffness=float(spec["target_stiffness"]),
            nominal_dt=float(spec["dt"]),
            box_chords=float(spec["box_chords"]),
            cells_per_chord=int(spec["cells_per_chord"]),
            cycles=int(spec["cycles"]),
            beta_degrees=float(spec.get("beta_degrees", BASELINE_BETA_DEGREES)),
            pitch_amplitude_degrees=float(
                spec.get("pitch_amplitude_degrees", BASELINE_PITCH_AMPLITUDE_DEGREES)
            ),
            balance_periodic_force=balance_periodic_force,
            save_visual_snapshots=bool(spec.get("save_visual_snapshots", False)),
        ))

    purpose = {
        "stiffness": "screen target-point penalty stiffness at fixed timestep",
        "timestep": "measure immersed-boundary coupling timestep sensitivity",
        "box": "measure periodic-box-size sensitivity at fixed spatial resolution",
        "grid": "measure Eulerian/Lagrangian mesh sensitivity at fixed box size",
        "startup": "identify an operationally repeatable post-startup force cycle",
        "production": "evaluate the requested beta and pitch-amplitude matrix",
    }[stage]
    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "study_id": f"re64_flapping_golden_{stage}",
        "stage": stage,
        "purpose": purpose,
        "physics": {
            "reynolds_number": REYNOLDS_NUMBER,
            "rho": RHO,
            "nominal_chord": NOMINAL_CHORD,
            "stroke_amplitude": STROKE_AMPLITUDE,
            "stroke_amplitude_over_chord": STROKE_AMPLITUDE_OVER_CHORD,
            "period": PERIOD,
            "u_ref": U_REF,
            "dynamic_viscosity": DYNAMIC_VISCOSITY,
            "reynolds_definition": "rho*U_ref*c/mu",
            "u_ref_definition": "pi*S/T",
        },
        "force_convention": {
            "lift": "L=Fy",
            "drag": "D=-Fx",
            "coefficient_scale": "q=0.5*rho*U_ref^2*c",
            "positive_drag_meaning": "drag; negative CD is thrust",
        },
        "motion_controls": {
            "mean_angle_of_attack_degrees": MEAN_ANGLE_OF_ATTACK_DEGREES,
            "pitch_phase_radians": PITCH_PHASE_RADIANS,
            "pivot": "midpoint",
        },
        "coupling_controls": {
            "delta_support": SUPPORT,
            "lagrangian_spacing": "ds=dx/2",
            "damped_springs": False,
            "target_damping": 0.0,
            "periodic_mean_momentum_constraint": bool(balance_periodic_force),
            "mean_force_balance_method": (
                "subtract spatial mean of Eulerian IB force before the fluid update; "
                "retain original IB reaction for plate-force output"
                if balance_periodic_force else "none"
            ),
        },
        "baseline_requested": baseline,
        "grid_dt_rule": grid_dt_rule if stage == "grid" else "not_applicable",
        "cases": cases,
        "analysis_note": (
            "The periodic mean-momentum constraint removes the spatially uniform IB "
            "force before the fluid update while retaining the original plate load. "
            "Reported cumulative drift is therefore the counterfactual unbalanced "
            "bulk-flow drift (the negative of accumulated balancing impulse divided "
            "by fluid mass), not "
            "the actual domain-mean velocity.  Force-cycle agreement is still an "
            "operational post-startup criterion."
            if balance_periodic_force else
            "No mean-force balance is applied.  A nonzero cycle-mean force in a finite "
            "periodic box changes domain-mean fluid momentum, so exact full-flow "
            "periodicity is not expected; force-cycle agreement is operational only."
        ),
    }
    manifest["configuration_sha256"] = _sha256_json(manifest)
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest: Mapping[str, object]) -> None:
    """Validate hashes and invariants for a study manifest."""

    if int(manifest.get("schema_version", -1)) != SCHEMA_VERSION:
        raise ValueError("unsupported manifest schema")
    stage = manifest.get("stage")
    if stage not in STAGES:
        raise ValueError("manifest has an unknown stage")
    expected = _sha256_json(_without_key(manifest, "configuration_sha256"))
    if manifest.get("configuration_sha256") != expected:
        raise ValueError("manifest configuration hash is invalid")
    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("manifest must contain at least one case")
    identifiers = set()
    for expected_index, case in enumerate(cases):
        if not isinstance(case, Mapping):
            raise ValueError("each case must be a JSON object")
        if int(case.get("index", -1)) != expected_index:
            raise ValueError("case indices must be consecutive and zero based")
        if case.get("stage") != stage:
            raise ValueError("case stage does not match manifest stage")
        identifier = str(case.get("case_id", ""))
        if not identifier or identifier in identifiers:
            raise ValueError("case identifiers must be nonempty and unique")
        identifiers.add(identifier)
        if case.get("damped_springs") is not False or float(case.get("target_damping", -1)) != 0.0:
            raise ValueError("this study must not enable structural or target dampers")
        if not isinstance(case.get("balance_periodic_force"), bool):
            raise ValueError("each case must state whether periodic force is balanced")
        if not math.isclose(float(case.get("ds_over_dx", -1)), 0.5, abs_tol=1.0e-14):
            raise ValueError("each case must use ds=dx/2")
        expected_case_hash = _sha256_json(_without_key(case, "case_configuration_sha256"))
        if case.get("case_configuration_sha256") != expected_case_hash:
            raise ValueError(f"case configuration hash is invalid for {identifier}")


def read_manifest(path: Path) -> dict[str, object]:
    with Path(path).open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("manifest must contain a JSON object")
    validate_manifest(value)
    return value


def write_manifest(path: Path, manifest: Mapping[str, object]) -> Path:
    validate_manifest(manifest)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.write_text(text, encoding="utf-8")
    return path


def select_case(manifest: Mapping[str, object], index: int) -> dict[str, object]:
    validate_manifest(manifest)
    cases = manifest["cases"]
    if int(index) != index or not (0 <= int(index) < len(cases)):
        raise ValueError(f"case index must be between 0 and {len(cases) - 1}")
    return copy.deepcopy(cases[int(index)])


def fixed_axis_forces(
    fx: Sequence[float] | np.ndarray,
    fy: Sequence[float] | np.ndarray,
    force_scale: float,
) -> dict[str, np.ndarray]:
    """Apply the task's lab-frame signs: L=Fy and D=-Fx."""

    fx_array = np.asarray(fx, dtype=float)
    fy_array = np.asarray(fy, dtype=float)
    if fx_array.shape != fy_array.shape or fx_array.ndim != 1:
        raise ValueError("Fx and Fy must be one-dimensional arrays of equal length")
    if not (np.all(np.isfinite(fx_array)) and np.all(np.isfinite(fy_array))):
        raise FloatingPointError("force values must be finite")
    if not math.isfinite(float(force_scale)) or force_scale <= 0.0:
        raise ValueError("force_scale must be finite and positive")
    lift = fy_array.copy()
    drag = -fx_array
    return {
        "fx": fx_array,
        "fy": fy_array,
        "lift": lift,
        "drag": drag,
        "cl": lift / force_scale,
        "cd": drag / force_scale,
    }


def inferred_domain_mean_velocity(
    fx_body: Sequence[float] | np.ndarray,
    fy_body: Sequence[float] | np.ndarray,
    *,
    dt: float,
    rho: float,
    area: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the bulk drift that an *unbalanced* periodic box would acquire.

    When the periodic mean-momentum constraint is active, this is a
    counterfactual diagnostic (and the negative of the cumulative balancing
    impulse divided by fluid mass), not the simulated domain-mean velocity.
    """

    forces = fixed_axis_forces(fx_body, fy_body, 1.0)
    scales = (float(dt), float(rho), float(area))
    if not all(math.isfinite(value) and value > 0.0 for value in scales):
        raise ValueError("dt, rho, and area must be finite and positive")
    # F_body is fluid-on-body, hence the force on the fluid is -F_body.
    denominator = float(rho) * float(area)
    u_mean = -np.cumsum(forces["fx"]) * float(dt) / denominator
    v_mean = -np.cumsum(forces["fy"]) * float(dt) / denominator
    return u_mean, v_mean


def detect_operational_steady_cycle(
    cycle_summaries: Sequence[Mapping[str, object]],
    consecutive_comparisons: Sequence[Mapping[str, object]],
    *,
    waveform_relative_l2_tolerance: float = 0.02,
    relative_mean_tolerance: float = 0.01,
    near_zero_absolute_mean_tolerance: float = 0.01,
) -> dict[str, object]:
    """Require two successive repeatable force-cycle pairs.

    Both CL and CD waveforms must meet the relative-L2 tolerance.  Their cycle
    means must change by no more than one percent, with an absolute 0.01
    allowance when the mean is near zero.  The final cycle is selected for
    averaging only when the final two consecutive cycle pairs both pass.
    """

    tolerances = (
        waveform_relative_l2_tolerance,
        relative_mean_tolerance,
        near_zero_absolute_mean_tolerance,
    )
    if not all(math.isfinite(float(value)) and value >= 0.0 for value in tolerances):
        raise ValueError("steady-cycle tolerances must be finite and non-negative")
    summaries = {int(row["cycle"]): row for row in cycle_summaries}
    evaluations: list[dict[str, object]] = []
    for comparison in consecutive_comparisons:
        earlier_cycle = int(comparison["earlier_cycle"])
        later_cycle = int(comparison["later_cycle"])
        if earlier_cycle not in summaries or later_cycle not in summaries:
            raise ValueError("repeatability comparison refers to an unknown cycle")
        evaluation: dict[str, object] = {
            "earlier_cycle": earlier_cycle,
            "later_cycle": later_cycle,
        }
        pair_passes = True
        for name in ("cl", "cd"):
            waveform_value = comparison.get(f"{name}_relative_l2")
            waveform_passes = (
                waveform_value is not None
                and math.isfinite(float(waveform_value))
                and float(waveform_value) <= waveform_relative_l2_tolerance
            )
            earlier_mean = float(summaries[earlier_cycle][f"mean_{name}"])
            later_mean = float(summaries[later_cycle][f"mean_{name}"])
            absolute_change = abs(later_mean - earlier_mean)
            reference_mean = max(abs(earlier_mean), abs(later_mean))
            mean_threshold = max(
                near_zero_absolute_mean_tolerance,
                relative_mean_tolerance * reference_mean,
            )
            mean_passes = absolute_change <= mean_threshold
            evaluation.update({
                f"{name}_waveform_relative_l2": waveform_value,
                f"{name}_waveform_passes": waveform_passes,
                f"{name}_absolute_mean_change": absolute_change,
                f"{name}_mean_change_threshold": mean_threshold,
                f"{name}_mean_passes": mean_passes,
            })
            pair_passes = pair_passes and waveform_passes and mean_passes
        evaluation["pair_passes"] = pair_passes
        evaluations.append(evaluation)

    averaging_cycle = None
    if len(evaluations) >= 2:
        previous, current = evaluations[-2:]
        if (
            bool(previous["pair_passes"])
            and bool(current["pair_passes"])
            and int(previous["later_cycle"]) == int(current["earlier_cycle"])
        ):
            averaging_cycle = int(current["later_cycle"])
    return {
        "reached": averaging_cycle is not None,
        "averaging_cycle": averaging_cycle,
        "required_successive_passing_pairs": 2,
        "criteria": {
            "cl_and_cd_waveform_relative_l2_max": waveform_relative_l2_tolerance,
            "relative_cycle_mean_change_max": relative_mean_tolerance,
            "near_zero_absolute_cycle_mean_change_max": (
                near_zero_absolute_mean_tolerance
            ),
        },
        "pair_evaluations": evaluations,
    }


def force_oscillation_diagnostic(
    values: Sequence[float] | np.ndarray,
) -> dict[str, float | int | None]:
    """Return threshold-free indicators of grid-scale temporal alternation."""

    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size < 3 or not np.all(np.isfinite(array)):
        raise ValueError("oscillation diagnostic needs at least three finite samples")
    centered = array - np.mean(array)
    fluctuation_rms = float(np.sqrt(np.mean(centered * centered)))
    second_difference_rms = float(np.sqrt(np.mean(np.diff(array, n=2) ** 2)))
    alternating_amplitude = float(
        abs(np.mean(centered * ((-1.0) ** np.arange(array.size))))
    )
    return {
        "sample_count": int(array.size),
        "fluctuation_rms": fluctuation_rms,
        "second_difference_rms": second_difference_rms,
        "second_difference_rms_over_fluctuation_rms": (
            second_difference_rms / fluctuation_rms
            if fluctuation_rms > 0.0 else None
        ),
        "nyquist_alternating_amplitude": alternating_amplitude,
        "nyquist_amplitude_over_fluctuation_rms": (
            alternating_amplitude / fluctuation_rms
            if fluctuation_rms > 0.0 else None
        ),
    }


def _numerical_source_fingerprints() -> dict[str, str]:
    project_python_root = HERE.parents[2]
    candidates = [
        HERE / "case_config.py",
        HERE / "generate_case.py",
        HERE / "kinematics.py",
        HERE / "update_Target_Point_Positions.py",
        HERE / "force_history.py",
        HERE / "main2d.py",
    ]
    candidates.extend(sorted((HERE.parents[1] / "IBM_Blackbox").glob("*.py")))
    return {
        path.relative_to(project_python_root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(candidates)
    }


def _snapshot_source_fingerprints() -> dict[str, str]:
    project_python_root = HERE.parents[2]
    candidates = [
        HERE / "flapping_visual.py",
        HERE.parent / "Translating_Plate" / "visual_pilot.py",
    ]
    return {
        path.relative_to(project_python_root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in candidates
    }


def _ensure_manifest_copy(output_root: Path, manifest: Mapping[str, object]) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    destination = output_root / "study_manifest.json"
    text = json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    try:
        with destination.open("x", encoding="utf-8") as handle:
            handle.write(text)
    except FileExistsError:
        existing = read_manifest(destination)
        if existing["configuration_sha256"] != manifest["configuration_sha256"]:
            raise ValueError(f"{output_root} already belongs to a different manifest")


def _read_json(path: Path) -> dict[str, object]:
    with Path(path).open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def _write_json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def _portable_paths(value: object, root: Path) -> object:
    """Convert artifact paths under root to portable relative strings."""

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


def _validate_completed_case(
    case_dir: Path,
    case: Mapping[str, object],
    manifest: Mapping[str, object],
) -> dict[str, object]:
    completion_path = case_dir / "COMPLETE.json"
    if not completion_path.is_file():
        raise FileNotFoundError(f"missing completion marker: {completion_path}")
    completion = _read_json(completion_path)
    if completion.get("case_configuration_sha256") != case["case_configuration_sha256"]:
        raise ValueError("completed case does not match the selected case definition")
    if completion.get("manifest_configuration_sha256") != manifest["configuration_sha256"]:
        raise ValueError("completed case does not match this study manifest")
    current_sources = _numerical_source_fingerprints()
    if completion.get("numerical_source_set_sha256") != _sha256_json(current_sources):
        raise ValueError("completed case was produced by different numerical source")
    if bool(case["save_visual_snapshots"]):
        current_snapshot_sources = _snapshot_source_fingerprints()
        if completion.get("snapshot_source_set_sha256") != _sha256_json(
            current_snapshot_sources
        ):
            raise ValueError("completed visual case used different snapshot source")
    metadata = _read_json(case_dir / "case_metadata.json")
    provenance = metadata.get("golden_study")
    if not isinstance(provenance, Mapping):
        raise ValueError("case metadata is missing golden-study provenance")
    if provenance.get("case_configuration_sha256") != case["case_configuration_sha256"]:
        raise ValueError("case metadata does not match the selected case definition")
    history = read_force_history(case_dir / "force_history.csv")
    if history["time"].size != int(case["steps"]):
        raise ValueError("force history row count does not match the case definition")
    return completion


def run_case_from_manifest(
    manifest: Mapping[str, object], index: int, output_root: Path
) -> Path:
    """Run one case atomically and leave failures in a diagnostic partial folder."""

    validate_manifest(manifest)
    case = select_case(manifest, index)
    output_root = Path(output_root).resolve()
    _ensure_manifest_copy(output_root, manifest)
    cases_root = output_root / "cases"
    cases_root.mkdir(parents=True, exist_ok=True)
    final_dir = cases_root / str(case["case_id"])
    if final_dir.exists():
        _validate_completed_case(final_dir, case, manifest)
        print(f"Case already complete and validated: {final_dir}")
        return final_dir

    job_token = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    partial_dir = cases_root / (
        f".{case['case_id']}.partial-{job_token}-{uuid.uuid4().hex[:8]}"
    )
    partial_dir.mkdir(parents=False, exist_ok=False)
    started = time.time()
    try:
        metadata = build_case(
            partial_dir,
            nx=int(case["nx"]),
            ny=int(case["ny"]),
            dt=float(case["dt"]),
            steps=int(case["steps"]),
            target_stiffness=float(case["target_stiffness"]),
            lx=float(case["lx"]),
            ly=float(case["ly"]),
            rho=RHO,
            reynolds_number=REYNOLDS_NUMBER,
            nominal_chord=NOMINAL_CHORD,
            period=PERIOD,
            stroke_amplitude_over_chord=STROKE_AMPLITUDE_OVER_CHORD,
            stroke_angle_degrees=float(case["beta_degrees"]),
            mean_angle_of_attack_degrees=MEAN_ANGLE_OF_ATTACK_DEGREES,
            pitch_amplitude_degrees=float(case["pitch_amplitude_degrees"]),
            pitch_phase_radians=PITCH_PHASE_RADIANS,
        )
        sources = _numerical_source_fingerprints()
        snapshot_sources = (
            _snapshot_source_fingerprints()
            if bool(case["save_visual_snapshots"]) else {}
        )
        metadata["case_name"] = case["case_id"]
        metadata["golden_study"] = {
            "study_id": manifest["study_id"],
            "stage": manifest["stage"],
            "case_index": int(case["index"]),
            "case_configuration": case,
            "case_configuration_sha256": case["case_configuration_sha256"],
            "manifest_configuration_sha256": manifest["configuration_sha256"],
            "numerical_source_sha256": sources,
            "snapshot_source_sha256": snapshot_sources,
            "damped_springs": False,
            "target_damping": 0.0,
            "periodic_mean_momentum_constraint": bool(
                case["balance_periodic_force"]
            ),
            "force_output_before_mean_balance": True,
        }
        _write_json(partial_dir / "case_metadata.json", metadata)

        from main2d import run_case
        steps_per_cycle = int(case["steps_per_cycle"])
        state_callback_interval = steps_per_cycle
        snapshot_count_expected = 0
        snapshot_writer = None
        final_cycle_start = int(case["steps"])
        if bool(case["save_visual_snapshots"]):
            from flapping_visual import FlappingVisualSnapshotWriter

            requested_frames = int(case["visual_frames_per_cycle"])
            frame_count = next(
                count for count in range(requested_frames, 0, -1)
                if steps_per_cycle % count == 0
            )
            state_callback_interval = steps_per_cycle // frame_count
            final_cycle_start = int(case["steps"]) - steps_per_cycle
            snapshot_writer = FlappingVisualSnapshotWriter.from_case(
                partial_dir, every_steps=state_callback_interval
            )
            snapshot_count_expected = frame_count + 1

        def progress_and_visual_callback(**payload: object):
            step = int(payload["step"])
            if step > 0 and step % steps_per_cycle == 0:
                mean_u = float(np.mean(np.asarray(payload["velocity_x"], dtype=float)))
                mean_v = float(np.mean(np.asarray(payload["velocity_y"], dtype=float)))
                print(
                    f"PROGRESS case={case['case_id']} "
                    f"cycle={step // steps_per_cycle}/{case['cycles']} "
                    f"step={step}/{case['steps']} mean_u={mean_u:.6e} "
                    f"mean_v={mean_v:.6e} elapsed_s={time.time() - started:.1f}",
                    flush=True,
                )
            if snapshot_writer is not None and step >= final_cycle_start:
                return snapshot_writer.record(**payload)
            return None

        history_path = run_case(
            partial_dir,
            state_callback=progress_and_visual_callback,
            state_callback_interval=state_callback_interval,
            balance_periodic_force=bool(case["balance_periodic_force"]),
        )
        history = read_force_history(history_path)
        core_names = ("time", "fx_body", "fy_body")
        if history["time"].size != int(case["steps"]):
            raise RuntimeError("solver did not write the requested number of force samples")
        if not all(np.all(np.isfinite(history[name])) for name in core_names):
            raise FloatingPointError("solver produced non-finite time or force values")
        snapshot_count = len(list((partial_dir / "visual_snapshots").glob("snapshot_step*.npz")))
        if snapshot_count != snapshot_count_expected:
            raise RuntimeError(
                f"visual case wrote {snapshot_count} snapshots; expected "
                f"{snapshot_count_expected}"
            )
        completion = {
            "schema_version": SCHEMA_VERSION,
            "study_id": manifest["study_id"],
            "stage": manifest["stage"],
            "case_id": case["case_id"],
            "case_index": case["index"],
            "case_configuration_sha256": case["case_configuration_sha256"],
            "manifest_configuration_sha256": manifest["configuration_sha256"],
            "numerical_source_set_sha256": _sha256_json(sources),
            "snapshot_source_set_sha256": (
                _sha256_json(snapshot_sources) if snapshot_sources else None
            ),
            "force_row_count": int(history["time"].size),
            "all_core_force_values_finite": True,
            "periodic_mean_momentum_constraint": bool(
                case["balance_periodic_force"]
            ),
            "snapshot_count": snapshot_count,
            "elapsed_seconds": time.time() - started,
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
    except BaseException as error:
        failure = {
            "case_id": case["case_id"],
            "case_index": case["index"],
            "elapsed_seconds": time.time() - started,
            "error_type": type(error).__name__,
            "error": str(error),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        }
        try:
            _write_json(partial_dir / "FAILED.json", failure)
        finally:
            raise


def _write_fixed_axis_csv(
    path: Path,
    history: Mapping[str, np.ndarray],
    case: Mapping[str, object],
) -> dict[str, np.ndarray]:
    q = 0.5 * RHO * U_REF * U_REF * NOMINAL_CHORD
    force = fixed_axis_forces(history["fx_body"], history["fy_body"], q)
    steps = np.asarray(history["step"], dtype=int)
    times = np.asarray(history["time"], dtype=float)
    steps_per_cycle = int(case["steps_per_cycle"])
    cycles = steps // steps_per_cycle + 1
    phase = ((steps % steps_per_cycle) + 0.5) / steps_per_cycle
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("step", "time", "cycle", "phase", "Fx", "Fy", "L", "D", "CL", "CD"))
        for row in zip(
            steps, times, cycles, phase, force["fx"], force["fy"],
            force["lift"], force["drag"], force["cl"], force["cd"],
        ):
            writer.writerow(
                int(value) if position in (0, 2) else format(float(value), ".17g")
                for position, value in enumerate(row)
            )
    return {**force, "step": steps, "time": times, "cycle": cycles, "phase": phase}


def _write_fixed_axis_svg(
    path: Path,
    curves: Mapping[str, np.ndarray],
    case: Mapping[str, object],
    cycle_number: int,
    *,
    result_status: str,
) -> Path:
    """Plot the task-defined lab-frame forces and coefficients for one cycle."""

    steps_per_cycle = int(case["steps_per_cycle"])
    start = (int(cycle_number) - 1) * steps_per_cycle
    stop = start + steps_per_cycle
    if stop > len(np.asarray(curves["time"])):
        raise ValueError("requested fixed-axis plot cycle is unavailable")
    sample_count = min(2400, steps_per_cycle)
    indices = np.unique(
        np.rint(np.linspace(start, stop - 1, sample_count)).astype(int)
    )
    phase = np.asarray(curves["phase"], dtype=float)[indices]
    panels = (
        (("lift", "L = Fy", "#1b9e77"), ("drag", "D = -Fx", "#d95f02"),
         "force per unit span"),
        (("cl", "CL", "#1b9e77"), ("cd", "CD", "#d95f02"),
         "coefficient using Uref"),
    )
    width, height = 1100.0, 760.0
    left, panel_width, panel_height = 92.0, 960.0, 270.0
    elements = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{int(width)}" height="{int(height)}" viewBox="0 0 {int(width)} {int(height)}">',
        '<rect width="100%" height="100%" fill="#fafafa"/>',
        '<style>.title{font:bold 21px sans-serif}.label{font:15px sans-serif}.tick{font:12px sans-serif}.legend{font:14px sans-serif}</style>',
        f'<text x="550" y="30" text-anchor="middle" class="title">{html.escape(str(case["case_id"]))}: fixed-axis cycle {int(cycle_number)}</text>',
        f'<text x="550" y="52" text-anchor="middle" class="label">{html.escape(result_status)}; q = 0.5 rho Uref^2 c</text>',
    ]
    for panel_index, (first, second, y_label) in enumerate(panels):
        top = 82.0 + panel_index * 325.0
        series = (first, second)
        all_values = np.concatenate([
            np.asarray(curves[name], dtype=float)[indices] for name, _label, _color in series
        ])
        ymin, ymax = float(np.min(all_values)), float(np.max(all_values))
        if math.isclose(ymin, ymax):
            padding = max(1.0, abs(ymin) * 0.1)
        else:
            padding = 0.06 * (ymax - ymin)
        ymin, ymax = ymin - padding, ymax + padding

        def points(values: np.ndarray) -> str:
            x = left + panel_width * phase
            y = top + panel_height * (ymax - values) / (ymax - ymin)
            return " ".join(f"{px:.2f},{py:.2f}" for px, py in zip(x, y))

        elements.extend((
            f'<rect x="{left}" y="{top}" width="{panel_width}" height="{panel_height}" fill="white" stroke="#777"/>',
            f'<text x="{left - 10}" y="{top + 12}" text-anchor="end" class="tick">{ymax:.4g}</text>',
            f'<text x="{left - 10}" y="{top + panel_height}" text-anchor="end" class="tick">{ymin:.4g}</text>',
            f'<text x="20" y="{top + panel_height / 2}" class="label">{html.escape(y_label)}</text>',
        ))
        if ymin <= 0.0 <= ymax:
            zero_y = top + panel_height * ymax / (ymax - ymin)
            elements.append(
                f'<line x1="{left}" y1="{zero_y:.2f}" x2="{left + panel_width}" y2="{zero_y:.2f}" stroke="#bbb"/>'
            )
        for legend_index, (name, label, color) in enumerate(series):
            values = np.asarray(curves[name], dtype=float)[indices]
            elements.extend((
                f'<polyline points="{points(values)}" fill="none" stroke="{color}" stroke-width="1.4"/>',
                f'<text x="{left + 15 + 95 * legend_index}" y="{top + 22}" class="legend" fill="{color}">{html.escape(label)}</text>',
            ))
    elements.extend((
        '<text x="572" y="730" text-anchor="middle" class="label">phase in selected complete flap</text>',
        '</svg>',
    ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path


def summarize_fixed_axis_case(
    history: Mapping[str, np.ndarray],
    case: Mapping[str, object],
) -> dict[str, object]:
    """Return cycle means, operational repeatability, and mean-flow drift."""

    q = 0.5 * RHO * U_REF * U_REF * NOMINAL_CHORD
    curves = fixed_axis_forces(history["fx_body"], history["fy_body"], q)
    steps_per_cycle = int(case["steps_per_cycle"])
    requested_cycles = int(case["cycles"])
    expected_samples = requested_cycles * steps_per_cycle
    if curves["cl"].size != expected_samples:
        raise ValueError("history length is not the requested integer number of cycles")
    dx = float(case["dx"])
    max_error = np.asarray(history.get("max_target_error", np.full(expected_samples, np.nan)))
    rms_error = np.asarray(history.get("rms_target_error", np.full(expected_samples, np.nan)))
    force_balance = np.asarray(
        history.get("force_balance_linf", np.full(expected_samples, np.nan))
    )
    unbalanced_u_mean, unbalanced_v_mean = inferred_domain_mean_velocity(
        history["fx_body"], history["fy_body"],
        dt=float(case["dt"]), rho=RHO, area=float(case["lx"]) * float(case["ly"]),
    )
    cycle_summaries = []
    for cycle_number in range(1, requested_cycles + 1):
        start = (cycle_number - 1) * steps_per_cycle
        stop = cycle_number * steps_per_cycle
        row: dict[str, object] = {"cycle": cycle_number, "sample_count": steps_per_cycle}
        for name in ("fx", "fy", "lift", "drag", "cl", "cd"):
            values = curves[name][start:stop]
            row[f"mean_{name}"] = float(np.mean(values))
            row[f"rms_{name}"] = float(np.sqrt(np.mean(values * values)))
            row[f"min_{name}"] = float(np.min(values))
            row[f"max_{name}"] = float(np.max(values))
        finite_max = max_error[start:stop][np.isfinite(max_error[start:stop])]
        finite_rms = rms_error[start:stop][np.isfinite(rms_error[start:stop])]
        finite_balance = force_balance[start:stop][
            np.isfinite(force_balance[start:stop])
        ]
        row["max_target_error_over_dx"] = (
            float(np.max(finite_max) / dx) if finite_max.size else None
        )
        row["rms_target_error_over_dx"] = (
            float(np.sqrt(np.mean(finite_rms * finite_rms)) / dx)
            if finite_rms.size else None
        )
        row["max_force_balance_linf"] = (
            float(np.max(finite_balance)) if finite_balance.size else None
        )
        row["max_force_balance_linf_over_q"] = (
            float(np.max(finite_balance) / q) if finite_balance.size else None
        )
        row["unbalanced_bulk_u_at_cycle_end"] = float(
            unbalanced_u_mean[stop - 1]
        )
        row["unbalanced_bulk_v_at_cycle_end"] = float(
            unbalanced_v_mean[stop - 1]
        )
        cycle_summaries.append(row)

    comparisons = []
    for second_cycle in range(2, requested_cycles + 1):
        first_start = (second_cycle - 2) * steps_per_cycle
        second_start = (second_cycle - 1) * steps_per_cycle
        comparison: dict[str, object] = {
            "earlier_cycle": second_cycle - 1,
            "later_cycle": second_cycle,
        }
        for name in ("cl", "cd"):
            earlier = curves[name][first_start:first_start + steps_per_cycle]
            later = curves[name][second_start:second_start + steps_per_cycle]
            rms_difference = float(np.sqrt(np.mean((later - earlier) ** 2)))
            later_rms = float(np.sqrt(np.mean(later ** 2)))
            comparison[f"{name}_relative_l2"] = (
                rms_difference / later_rms if later_rms > 0.0 else None
            )
            comparison[f"{name}_cycle_mean_change"] = float(
                np.mean(later) - np.mean(earlier)
            )
        comparisons.append(comparison)

    mean_force_nonzero = bool(
        abs(cycle_summaries[-1]["mean_fx"]) > 1.0e-12
        or abs(cycle_summaries[-1]["mean_fy"]) > 1.0e-12
    )
    balance_active = bool(case.get("balance_periodic_force", False))
    drift_speed = math.hypot(unbalanced_u_mean[-1], unbalanced_v_mean[-1])
    steady = detect_operational_steady_cycle(cycle_summaries, comparisons)
    last_start = (requested_cycles - 1) * steps_per_cycle
    oscillation = {
        name: force_oscillation_diagnostic(
            curves[name][last_start:last_start + steps_per_cycle]
        )
        for name in ("cl", "cd")
    }
    reported_average = None
    if steady["reached"]:
        selected = cycle_summaries[int(steady["averaging_cycle"]) - 1]
        reported_average = {
            "cycle": selected["cycle"],
            "mean_fx": selected["mean_fx"],
            "mean_fy": selected["mean_fy"],
            "mean_lift": selected["mean_lift"],
            "mean_drag": selected["mean_drag"],
            "mean_cl": selected["mean_cl"],
            "mean_cd": selected["mean_cd"],
        }
    return {
        "force_definition": {
            "lift": "L=Fy",
            "drag": "D=-Fx",
            "coefficient_scale": q,
            "coefficient_scale_definition": "0.5*rho*U_ref^2*c",
        },
        "cycle_summaries": cycle_summaries,
        "consecutive_cycle_force_repeatability": comparisons,
        "operational_steady_cycle": steady,
        "reported_post_startup_average": reported_average,
        "reported_result_status": (
            "operationally_steady" if steady["reached"] else "transient_not_reported"
        ),
        "last_cycle_high_frequency_diagnostic": {
            "interpretation": (
                "Threshold-free comparison aid. Large second-difference or Nyquist "
                "ratios relative to smaller-dt/stiffness cases can indicate coupled "
                "target-force oscillation; this is not an automatic pass/fail test."
            ),
            **oscillation,
        },
        "periodicity_status": "operational_force_repeatability_only",
        "periodicity_warning": (
            "The periodic mean-momentum constraint was active.  The original IB "
            "reaction is reported as plate force, while the spatially uniform part of "
            "the force applied to the fluid was removed.  Exact wake periodicity is "
            "not asserted; consecutive-force-cycle agreement is an operational "
            "post-startup criterion."
            if balance_active else
            "No mean-force balance was applied.  A nonzero cycle-mean force in a "
            "finite periodic box changes domain-mean fluid momentum, so exact full-flow "
            "periodicity is not possible; consecutive-force-cycle agreement is only an "
            "operational post-startup criterion."
        ),
        "last_cycle_mean_force_is_nonzero": mean_force_nonzero,
        "bulk_momentum_diagnostic": {
            "method": "Umean(t)=-integral(F_body)dt/(rho*Lx*Ly)",
            "periodic_mean_momentum_constraint": balance_active,
            "interpretation": (
                "counterfactual unbalanced bulk-flow drift; equal to the negative "
                "of cumulative balancing impulse divided by fluid mass; not actual "
                "mean velocity"
                if balance_active else "inferred simulated domain-mean velocity drift"
            ),
            "actual_mean_velocity_expectation": (
                "held at its initially zero value up to solver roundoff"
                if balance_active else "not constrained"
            ),
            "final_unbalanced_bulk_u": float(unbalanced_u_mean[-1]),
            "final_unbalanced_bulk_v": float(unbalanced_v_mean[-1]),
            "final_unbalanced_bulk_speed_over_u_ref": float(drift_speed / U_REF),
        },
    }


def _partial_directories(output_root: Path) -> list[str]:
    cases_root = output_root / "cases"
    if not cases_root.is_dir():
        return []
    return sorted(path.name for path in cases_root.glob(".*.partial-*"))


def _periodic_phase_profile(
    values: Sequence[float] | np.ndarray,
    case: Mapping[str, object],
    cycle_number: int,
    *,
    samples: int = 721,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate one complete, uniformly sampled cycle to a shared phase grid."""

    values_array = np.asarray(values, dtype=float)
    steps_per_cycle = int(case["steps_per_cycle"])
    start = (int(cycle_number) - 1) * steps_per_cycle
    stop = start + steps_per_cycle
    cycle_values = values_array[start:stop]
    if cycle_values.size != steps_per_cycle or not np.all(np.isfinite(cycle_values)):
        raise ValueError("requested phase profile is incomplete or non-finite")
    source_phase = (np.arange(steps_per_cycle, dtype=float) + 0.5) / steps_per_cycle
    phase = (np.arange(int(samples), dtype=float) + 0.5) / int(samples)
    extended_phase = np.concatenate((source_phase - 1.0, source_phase, source_phase + 1.0))
    extended_values = np.tile(cycle_values, 3)
    return phase, np.interp(phase, extended_phase, extended_values)


def _relative_l2(values: np.ndarray, reference: np.ndarray) -> float | None:
    difference_rms = float(np.sqrt(np.mean((values - reference) ** 2)))
    reference_rms = float(np.sqrt(np.mean(reference ** 2)))
    return difference_rms / reference_rms if reference_rms > 0.0 else None


def _calibration_comparisons(
    completed: Sequence[tuple[Mapping[str, object], Mapping[str, np.ndarray], Mapping[str, object]]],
    manifest: Mapping[str, object],
) -> list[dict[str, object]]:
    """Compare adjacent calibration cases to their next refinement/reference."""

    if manifest["stage"] not in {"stiffness", "timestep", "box", "grid"}:
        return []
    by_id = {str(case["case_id"]): (case, history, summary) for case, history, summary in completed}
    comparisons: list[dict[str, object]] = []
    planned = manifest["cases"]
    for earlier_case, reference_case in zip(planned, planned[1:]):
        earlier_id = str(earlier_case["case_id"])
        reference_id = str(reference_case["case_id"])
        row: dict[str, object] = {
            "case_id": earlier_id,
            "reference_case_id": reference_id,
            "reference_meaning": {
                "stiffness": "next_higher_stiffness",
                "timestep": "next_finer_timestep",
                "box": "next_larger_box",
                "grid": "next_finer_grid",
            }[str(manifest["stage"])],
        }
        if earlier_id not in by_id or reference_id not in by_id:
            row["status"] = "unavailable_missing_case"
            comparisons.append(row)
            continue
        case, history, summary = by_id[earlier_id]
        ref_case, ref_history, ref_summary = by_id[reference_id]
        case_cycle = int(case["cycles"])
        ref_cycle = int(ref_case["cycles"])
        q = 0.5 * RHO * U_REF * U_REF * NOMINAL_CHORD
        force = fixed_axis_forces(history["fx_body"], history["fy_body"], q)
        ref_force = fixed_axis_forces(
            ref_history["fx_body"], ref_history["fy_body"], q
        )
        row.update({
            "status": "complete",
            "case_cycle": case_cycle,
            "reference_cycle": ref_cycle,
            "case_operationally_steady": bool(
                summary["operational_steady_cycle"]["reached"]
            ),
            "reference_operationally_steady": bool(
                ref_summary["operational_steady_cycle"]["reached"]
            ),
        })
        for name in ("cl", "cd"):
            _, profile = _periodic_phase_profile(force[name], case, case_cycle)
            _, reference_profile = _periodic_phase_profile(
                ref_force[name], ref_case, ref_cycle
            )
            row.update({
                f"{name}_relative_l2_vs_reference": _relative_l2(
                    profile, reference_profile
                ),
                f"{name}_case_mean": float(np.mean(profile)),
                f"{name}_reference_mean": float(np.mean(reference_profile)),
                f"{name}_case_min": float(np.min(profile)),
                f"{name}_case_max": float(np.max(profile)),
                f"{name}_reference_min": float(np.min(reference_profile)),
                f"{name}_reference_max": float(np.max(reference_profile)),
            })
        comparisons.append(row)
    return comparisons


def _rotational_invariance_report(
    completed: Sequence[tuple[Mapping[str, object], Mapping[str, np.ndarray], Mapping[str, object]]],
) -> list[dict[str, object]]:
    """Compare production forces after rotating each beta case to stroke axes."""

    records: dict[tuple[float, float], tuple[Mapping[str, object], Mapping[str, np.ndarray], Mapping[str, object]]] = {}
    for case, history, summary in completed:
        records[(float(case["pitch_amplitude_degrees"]), float(case["beta_degrees"]))] = (
            case, history, summary
        )
    q = 0.5 * RHO * U_REF * U_REF * NOMINAL_CHORD
    report: list[dict[str, object]] = []
    for pitch in DEFAULT_PITCH_VALUES:
        reference_key = (float(pitch), 0.0)
        if reference_key not in records:
            continue
        ref_case, ref_history, ref_summary = records[reference_key]
        ref_cycle = int(
            ref_summary["operational_steady_cycle"]["averaging_cycle"]
            or ref_case["cycles"]
        )
        ref_fx = np.asarray(ref_history["fx_body"], dtype=float) / q
        ref_fy = np.asarray(ref_history["fy_body"], dtype=float) / q
        _, ref_along = _periodic_phase_profile(ref_fx, ref_case, ref_cycle)
        _, ref_normal = _periodic_phase_profile(ref_fy, ref_case, ref_cycle)
        for beta in DEFAULT_BETA_VALUES[1:]:
            key = (float(pitch), float(beta))
            row: dict[str, object] = {
                "pitch_amplitude_degrees": float(pitch),
                "beta_degrees": float(beta),
                "reference_beta_degrees": 0.0,
            }
            if key not in records:
                row["status"] = "unavailable_missing_case"
                report.append(row)
                continue
            case, history, summary = records[key]
            cycle = int(
                summary["operational_steady_cycle"]["averaging_cycle"]
                or case["cycles"]
            )
            beta_radians = math.radians(float(beta))
            fx = np.asarray(history["fx_body"], dtype=float) / q
            fy = np.asarray(history["fy_body"], dtype=float) / q
            along = math.cos(beta_radians) * fx + math.sin(beta_radians) * fy
            normal = -math.sin(beta_radians) * fx + math.cos(beta_radians) * fy
            _, along_profile = _periodic_phase_profile(along, case, cycle)
            _, normal_profile = _periodic_phase_profile(normal, case, cycle)
            row.update({
                "status": "complete",
                "case_cycle": cycle,
                "reference_cycle": ref_cycle,
                "case_operationally_steady": bool(
                    summary["operational_steady_cycle"]["reached"]
                ),
                "reference_operationally_steady": bool(
                    ref_summary["operational_steady_cycle"]["reached"]
                ),
                "stroke_axis_force_relative_l2": _relative_l2(
                    along_profile, ref_along
                ),
                "stroke_normal_force_relative_l2": _relative_l2(
                    normal_profile, ref_normal
                ),
            })
            report.append(row)
    return report


def analyze_stage(
    manifest: Mapping[str, object], output_root: Path, *, make_mp4: bool = True
) -> dict[str, object]:
    """Analyze completed cases and report missing/failed cases without aborting."""

    validate_manifest(manifest)
    output_root = Path(output_root).resolve()
    _ensure_manifest_copy(output_root, manifest)
    summary_root = output_root / "summary"
    summary_root.mkdir(parents=True, exist_ok=True)
    statuses: list[dict[str, object]] = []
    table_rows: list[dict[str, object]] = []
    completed_data: list[
        tuple[Mapping[str, object], Mapping[str, np.ndarray], Mapping[str, object]]
    ] = []
    for case in manifest["cases"]:
        case_dir = output_root / "cases" / str(case["case_id"])
        status: dict[str, object] = {
            "index": case["index"],
            "case_id": case["case_id"],
        }
        if not case_dir.is_dir():
            status["status"] = "missing"
            statuses.append(status)
            continue
        try:
            completion = _validate_completed_case(case_dir, case, manifest)
            history = read_force_history(case_dir / "force_history.csv")
            analysis_dir = case_dir / "force_analysis"
            fixed_path = analysis_dir / "fixed_axis_forces.csv"
            fixed_curves = _write_fixed_axis_csv(fixed_path, history, case)
            fixed_summary = summarize_fixed_axis_case(history, case)
            fixed_summary_path = _write_json(
                analysis_dir / "fixed_axis_cycle_summary.json", fixed_summary
            )
            plot_cycle = int(
                fixed_summary["operational_steady_cycle"]["averaging_cycle"]
                or case["cycles"]
            )
            fixed_plot_path = _write_fixed_axis_svg(
                analysis_dir / "fixed_axis_force_cycle.svg",
                fixed_curves,
                case,
                plot_cycle,
                result_status=str(fixed_summary["reported_result_status"]),
            )
            try:
                legacy = analyze_case(
                    case_dir,
                    output_dir=analysis_dir / "supplemental_motion_relative",
                    analysis_cycles=tuple(range(1, int(case["cycles"]) + 1)),
                )
            except (FloatingPointError, KeyError, OSError, RuntimeError, ValueError) as error:
                legacy = {
                    "status": "optional_supplemental_analysis_failed",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            visual_artifacts = None
            if bool(case["save_visual_snapshots"]):
                from flapping_visual import build_visual_artifacts

                try:
                    visual_artifacts = build_visual_artifacts(
                        case_dir,
                        output_dir=case_dir / "representative_flow_animation",
                        fps=10.0,
                        make_mp4=make_mp4,
                    )
                except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
                    visual_artifacts = {
                        "status": "optional_visual_render_failed",
                        "error_type": type(error).__name__,
                        "error": str(error),
                    }
            status.update({
                "status": "complete",
                "completion": completion,
                "fixed_axis_history": str(fixed_path.relative_to(output_root).as_posix()),
                "fixed_axis_summary": str(fixed_summary_path.relative_to(output_root).as_posix()),
                "fixed_axis_force_plot": str(
                    fixed_plot_path.relative_to(output_root).as_posix()
                ),
                "supplemental_analysis": legacy,
                "representative_flow_animation": visual_artifacts,
            })
            last = fixed_summary["cycle_summaries"][-1]
            steady = fixed_summary["operational_steady_cycle"]
            reported = fixed_summary["reported_post_startup_average"]
            oscillation = fixed_summary["last_cycle_high_frequency_diagnostic"]
            prior_comparison = fixed_summary["consecutive_cycle_force_repeatability"]
            latest_repeatability = prior_comparison[-1] if prior_comparison else {}
            table_rows.append({
                "index": case["index"],
                "case_id": case["case_id"],
                "target_stiffness": case["target_stiffness"],
                "dt": case["dt"],
                "box_chords": case["box_chords"],
                "cells_per_chord": case["cells_per_chord"],
                "beta_degrees": case["beta_degrees"],
                "pitch_amplitude_degrees": case["pitch_amplitude_degrees"],
                "last_cycle": last["cycle"],
                "last_cycle_mean_cl": last["mean_cl"],
                "last_cycle_mean_cd": last["mean_cd"],
                "last_cycle_rms_cl": last["rms_cl"],
                "last_cycle_rms_cd": last["rms_cd"],
                "last_cycle_min_cl": last["min_cl"],
                "last_cycle_max_cl": last["max_cl"],
                "last_cycle_min_cd": last["min_cd"],
                "last_cycle_max_cd": last["max_cd"],
                "operationally_steady": steady["reached"],
                "averaging_cycle": steady["averaging_cycle"],
                "reported_mean_cl": reported["mean_cl"] if reported else None,
                "reported_mean_cd": reported["mean_cd"] if reported else None,
                "reported_mean_lift": reported["mean_lift"] if reported else None,
                "reported_mean_drag": reported["mean_drag"] if reported else None,
                "max_target_error_over_dx": last["max_target_error_over_dx"],
                "rms_target_error_over_dx": last["rms_target_error_over_dx"],
                "max_force_balance_linf": last["max_force_balance_linf"],
                "max_force_balance_linf_over_q": last["max_force_balance_linf_over_q"],
                "cl_relative_l2_previous_cycle": latest_repeatability.get("cl_relative_l2"),
                "cd_relative_l2_previous_cycle": latest_repeatability.get("cd_relative_l2"),
                "cl_second_difference_ratio": oscillation["cl"]["second_difference_rms_over_fluctuation_rms"],
                "cd_second_difference_ratio": oscillation["cd"]["second_difference_rms_over_fluctuation_rms"],
                "cl_nyquist_ratio": oscillation["cl"]["nyquist_amplitude_over_fluctuation_rms"],
                "cd_nyquist_ratio": oscillation["cd"]["nyquist_amplitude_over_fluctuation_rms"],
                "unbalanced_bulk_drift_over_u_ref": fixed_summary["bulk_momentum_diagnostic"]["final_unbalanced_bulk_speed_over_u_ref"],
            })
            completed_data.append((case, history, fixed_summary))
        except (
            FileNotFoundError, FloatingPointError, KeyError, OSError,
            RuntimeError, ValueError,
        ) as error:
            status.update({
                "status": "invalid_or_analysis_failed",
                "error_type": type(error).__name__,
                "error": str(error),
            })
        statuses.append(status)

    table_path = summary_root / "case_summary.csv"
    fieldnames = (
        "index", "case_id", "target_stiffness", "dt", "box_chords",
        "cells_per_chord", "beta_degrees", "pitch_amplitude_degrees",
        "last_cycle", "last_cycle_mean_cl", "last_cycle_mean_cd",
        "last_cycle_rms_cl", "last_cycle_rms_cd", "last_cycle_min_cl",
        "last_cycle_max_cl", "last_cycle_min_cd", "last_cycle_max_cd",
        "operationally_steady", "averaging_cycle", "reported_mean_cl",
        "reported_mean_cd", "reported_mean_lift", "reported_mean_drag",
        "max_target_error_over_dx", "rms_target_error_over_dx",
        "max_force_balance_linf", "max_force_balance_linf_over_q",
        "cl_relative_l2_previous_cycle", "cd_relative_l2_previous_cycle",
        "cl_second_difference_ratio", "cd_second_difference_ratio",
        "cl_nyquist_ratio", "cd_nyquist_ratio",
        "unbalanced_bulk_drift_over_u_ref",
    )
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in table_rows:
            writer.writerow(row)
    completed_count = sum(row["status"] == "complete" for row in statuses)
    calibration_comparisons = _calibration_comparisons(completed_data, manifest)
    comparison_path = summary_root / "cross_case_comparison.csv"
    comparison_fields = sorted({key for row in calibration_comparisons for key in row})
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        if comparison_fields:
            writer = csv.DictWriter(handle, fieldnames=comparison_fields)
            writer.writeheader()
            writer.writerows(calibration_comparisons)
    rotational_report = (
        _rotational_invariance_report(completed_data)
        if manifest["stage"] == "production" else []
    )
    result: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "study_id": manifest["study_id"],
        "stage": manifest["stage"],
        "manifest_configuration_sha256": manifest["configuration_sha256"],
        "planned_case_count": len(manifest["cases"]),
        "completed_case_count": completed_count,
        "partial": completed_count != len(manifest["cases"]),
        "case_statuses": statuses,
        "leftover_partial_directories": _partial_directories(output_root),
        "case_summary_csv": table_path.relative_to(output_root).as_posix(),
        "cross_case_comparisons": calibration_comparisons,
        "cross_case_comparison_csv": comparison_path.relative_to(output_root).as_posix(),
        "production_rotational_invariance": rotational_report,
        "periodicity_warning": manifest["analysis_note"],
        "selection_status": "human_review_required",
    }
    portable_result = _portable_paths(result, output_root)
    _write_json(summary_root / "stage_analysis.json", portable_result)
    return portable_result


def _parse_float_list(text: str) -> tuple[float, ...]:
    try:
        return tuple(float(part.strip()) for part in text.split(",") if part.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _parse_int_list(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part.strip()) for part in text.split(",") if part.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    return values


def _add_design_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--target-stiffness", type=float, default=DEFAULT_TARGET_STIFFNESS)
    parser.add_argument("--dt", type=float, default=DEFAULT_DT)
    parser.add_argument("--box-chords", type=float, default=DEFAULT_BOX_CHORDS)
    parser.add_argument("--cells-per-chord", type=int, default=DEFAULT_CELLS_PER_CHORD)
    parser.add_argument("--cycles", type=int)
    parser.add_argument("--stiffness-values", type=_parse_float_list, default=DEFAULT_STIFFNESS_VALUES)
    parser.add_argument("--dt-values", type=_parse_float_list, default=DEFAULT_DT_VALUES)
    parser.add_argument("--box-values", type=_parse_float_list, default=DEFAULT_BOX_VALUES)
    parser.add_argument("--grid-values", type=_parse_int_list, default=DEFAULT_GRID_VALUES)
    parser.add_argument("--beta-values", type=_parse_float_list, default=DEFAULT_BETA_VALUES)
    parser.add_argument("--pitch-values", type=_parse_float_list, default=DEFAULT_PITCH_VALUES)
    parser.add_argument("--grid-dt-rule", choices=("constant_tether", "constant"), default="constant")
    parser.add_argument(
        "--no-balance-periodic-force",
        dest="balance_periodic_force",
        action="store_false",
        help="disable the periodic-box mean-momentum constraint for a sensitivity run",
    )
    parser.set_defaults(balance_periodic_force=True)


def _manifest_from_args(args: argparse.Namespace) -> dict[str, object]:
    return design_stage(
        args.stage,
        target_stiffness=args.target_stiffness,
        dt=args.dt,
        box_chords=args.box_chords,
        cells_per_chord=args.cells_per_chord,
        cycles=args.cycles,
        stiffness_values=args.stiffness_values,
        dt_values=args.dt_values,
        box_values=args.box_values,
        grid_values=args.grid_values,
        beta_values=args.beta_values,
        pitch_values=args.pitch_values,
        grid_dt_rule=args.grid_dt_rule,
        balance_periodic_force=args.balance_periodic_force,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    manifest_parser = subparsers.add_parser("manifest", help="write a deterministic stage manifest")
    _add_design_arguments(manifest_parser)
    manifest_parser.add_argument("--output", type=Path, required=True)
    show_parser = subparsers.add_parser("show", help="print and validate a manifest")
    show_parser.add_argument("--manifest", type=Path, required=True)
    count_parser = subparsers.add_parser("count", help="print the manifest case count")
    count_parser.add_argument("--manifest", type=Path, required=True)
    run_parser = subparsers.add_parser("run", help="run one manifest case")
    run_parser.add_argument("--manifest", type=Path, required=True)
    run_parser.add_argument("--case-index", type=int, required=True)
    run_parser.add_argument("--output-root", type=Path, required=True)
    analyze_parser = subparsers.add_parser("analyze", help="analyze all available cases")
    analyze_parser.add_argument("--manifest", type=Path, required=True)
    analyze_parser.add_argument("--output-root", type=Path, required=True)
    analyze_parser.add_argument(
        "--no-mp4", action="store_true",
        help="still build the representative HTML animation but skip optional MP4",
    )
    args = parser.parse_args()

    try:
        if args.action == "manifest":
            manifest = _manifest_from_args(args)
            output = write_manifest(args.output, manifest)
            print(output)
        elif args.action == "show":
            print(json.dumps(read_manifest(args.manifest), indent=2, sort_keys=True))
        elif args.action == "count":
            print(len(read_manifest(args.manifest)["cases"]))
        elif args.action == "run":
            run_case_from_manifest(
                read_manifest(args.manifest), args.case_index, args.output_root
            )
        else:
            result = analyze_stage(
                read_manifest(args.manifest), args.output_root,
                make_mp4=not args.no_mp4,
            )
            print(json.dumps(result, indent=2, sort_keys=True))
    except (FileNotFoundError, FloatingPointError, KeyError, OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
