#!/usr/bin/env python3
"""Summarize Re=128 calibration, grid-convergence, and box-size runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from force_history import HEADER
from re128_study import (
    FINAL_TIME,
    NOMINAL_CHORD,
    PRODUCTION_DT,
    PRODUCTION_TARGET_STIFFNESS,
    REYNOLDS_NUMBER,
    StudyCase,
    case_record,
    cases_for_phase,
)


def read_history(path: Path) -> dict[str, np.ndarray]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in HEADER if name not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} has no data rows")
    data = {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in HEADER
    }
    if not all(np.all(np.isfinite(column)) for column in data.values()):
        raise ValueError(f"{path} contains non-finite values")
    return data


def load_case(output_root: Path, case: StudyCase) -> dict:
    case_dir = Path(output_root) / "cases" / case.phase / case.case_id
    if not (case_dir / "COMPLETE.json").is_file():
        raise FileNotFoundError(f"case is not complete: {case_dir}")
    metadata = json.loads((case_dir / "case_metadata.json").read_text(encoding="utf-8"))
    requested_hash = case_record(case)["configuration_sha256"]
    completed_hash = metadata.get("study", {}).get("configuration_sha256")
    if completed_hash != requested_hash:
        raise ValueError(
            f"case configuration does not match the requested manifest: {case_dir}"
        )
    history = read_history(case_dir / "force_history.csv")
    force_scale = float(metadata["force"]["coefficient_force_scale"])
    speed = float(metadata["motion"]["speed"])
    chord = float(metadata["geometry"]["nominal_chord"])
    history["tau"] = history["time"] * abs(speed) / chord
    history["cd"] = -history["fx_body"] / force_scale
    history["cl"] = history["fy_body"] / force_scale
    return {"case": case, "directory": str(case_dir.resolve()),
            "metadata": metadata, "history": history}


def signal_statistics(values: np.ndarray, tau: np.ndarray,
                      midpoint_spacing: float) -> dict:
    return {
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "mean": float(np.mean(values)),
        "rms": float(np.sqrt(np.mean(values * values))),
        "impulse_over_convective_time": float(midpoint_spacing * np.sum(values)),
    }


def high_frequency_energy_fraction(values: np.ndarray,
                                   cutoff_cycles_per_step: float = 0.1) -> float:
    """Return detrended spectral energy above a step-scale frequency cutoff.

    This is a numerical-ringing diagnostic, not a physical-frequency filter.
    Cubic detrending prevents the start-up transient from dominating the
    low-frequency denominator.
    """

    values = np.asarray(values, dtype=float)
    if values.size < 8:
        return 0.0
    coordinate = np.linspace(-1.0, 1.0, values.size)
    degree = min(3, values.size - 1)
    trend = np.polyval(np.polyfit(coordinate, values, degree), coordinate)
    spectrum = np.fft.rfft(values - trend)
    energy = np.abs(spectrum) ** 2
    frequencies = np.fft.rfftfreq(values.size, d=1.0)
    nonzero = frequencies > 0.0
    denominator = float(np.sum(energy[nonzero]))
    if denominator <= np.finfo(float).tiny:
        return 0.0
    return float(
        np.sum(energy[frequencies >= cutoff_cycles_per_step]) / denominator
    )


def case_summary(loaded: dict) -> dict:
    case = loaded["case"]
    metadata = loaded["metadata"]
    history = loaded["history"]
    max_target_error = float(np.max(history["max_target_error"]))
    max_rms_target_error = float(np.max(history["rms_target_error"]))
    convective_dt = (
        abs(float(metadata["motion"]["speed"])) * case.dt / NOMINAL_CHORD
    )
    return {
        "case_id": case.case_id,
        "phase": case.phase,
        "study_axes": list(case.study_axes),
        "nx": case.nx,
        "ny": case.ny,
        "box_chords": case.box_chords,
        "cells_per_chord": case.cells_per_chord,
        "dx": case.dx,
        "ds": case.ds,
        "dt": case.dt,
        "advective_cfl": abs(float(metadata["motion"]["speed"])) * case.dt / case.dx,
        "viscous_number": (
            float(metadata["fluid"]["mu"]) / float(metadata["fluid"]["rho"])
        ) * case.dt / (case.dx * case.dx),
        "steps": int(history["step"].size),
        "target_stiffness": case.target_stiffness,
        "target_stiffness_times_ds": case.target_stiffness * case.ds,
        "penalty_coupling_number": case.penalty_coupling_number,
        "nominal_reynolds_number": float(
            metadata["fluid"]["nominal_reynolds_number"]
        ),
        "effective_chord": float(
            metadata["geometry"]["effective_quadrature_chord"]
        ),
        "max_target_error": max_target_error,
        "max_target_error_over_chord": max_target_error / NOMINAL_CHORD,
        "max_target_error_over_dx": max_target_error / case.dx,
        "max_rms_target_error": max_rms_target_error,
        "max_rms_target_error_over_chord": max_rms_target_error / NOMINAL_CHORD,
        "max_rms_target_error_over_dx": max_rms_target_error / case.dx,
        "max_force_balance_linf": float(np.max(history["force_balance_linf"])),
        "high_frequency_energy_fraction": {
            "drag_coefficient": high_frequency_energy_fraction(history["cd"]),
            "lift_coefficient": high_frequency_energy_fraction(history["cl"]),
            "cutoff_cycles_per_step": 0.1,
            "detrend_degree": 3,
        },
        "drag_coefficient": signal_statistics(
            history["cd"], history["tau"], convective_dt
        ),
        "lift_coefficient": signal_statistics(
            history["cl"], history["tau"], convective_dt
        ),
        "directory": loaded["directory"],
    }


def common_reference_samples(reference: dict, candidate: dict,
                             field: str) -> tuple[np.ndarray, np.ndarray]:
    ref_time = reference["history"]["time"]
    candidate_time = candidate["history"]["time"]
    lower = max(ref_time[0], candidate_time[0])
    upper = min(ref_time[-1], candidate_time[-1])
    mask = (ref_time >= lower) & (ref_time <= upper)
    if np.count_nonzero(mask) < 2:
        raise ValueError("histories have no usable common time interval")
    reference_values = reference["history"][field][mask]
    candidate_values = np.interp(
        ref_time[mask], candidate_time, candidate["history"][field]
    )
    return reference_values, candidate_values


def difference_metrics(reference: dict, candidate: dict, field: str) -> dict:
    ref, cand = common_reference_samples(reference, candidate, field)
    difference = cand - ref
    l2 = float(np.sqrt(np.mean(difference * difference)))
    linf = float(np.max(np.abs(difference)))
    reference_rms = max(float(np.sqrt(np.mean(ref * ref))), np.finfo(float).tiny)
    reference_linf = max(float(np.max(np.abs(ref))), np.finfo(float).tiny)
    return {
        "l2": l2,
        "linf": linf,
        "relative_l2": l2 / reference_rms,
        "relative_linf": linf / reference_linf,
        "samples": int(ref.size),
    }


def compare_cases(reference: dict, candidate: dict) -> dict:
    return {
        "reference": reference["case"].case_id,
        "candidate": candidate["case"].case_id,
        "drag_coefficient": difference_metrics(reference, candidate, "cd"),
        "lift_coefficient": difference_metrics(reference, candidate, "cl"),
    }


def absolute_difference_norm(first: dict, second: dict, field: str) -> float:
    reference_values, candidate_values = common_reference_samples(first, second, field)
    difference = candidate_values - reference_values
    return float(np.sqrt(np.mean(difference * difference)))


def observed_order(coarse: dict, medium: dict, fine: dict, field: str) -> dict:
    coarse_medium = absolute_difference_norm(medium, coarse, field)
    medium_fine = absolute_difference_norm(fine, medium, field)
    refinement_ratio = coarse["case"].dx / medium["case"].dx
    next_ratio = medium["case"].dx / fine["case"].dx
    if not math.isclose(refinement_ratio, next_ratio, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError("observed order requires a constant refinement ratio")
    order = None
    if coarse_medium > 0.0 and medium_fine > 0.0:
        order = math.log(coarse_medium / medium_fine) / math.log(refinement_ratio)
    return {
        "coarse_medium_l2_difference": coarse_medium,
        "medium_fine_l2_difference": medium_fine,
        "refinement_ratio": refinement_ratio,
        "observed_order": order,
    }


def calibration_report(
    loaded: list[dict],
    relative_l2_tolerance: float,
    target_error_tolerance: float,
    target_error_dx_tolerance: float,
    rms_target_error_dx_tolerance: float,
    high_frequency_energy_tolerance: float,
    high_frequency_reference_multiplier: float,
    penalty_coupling_tolerance: float,
) -> dict:
    by_id = {item["case"].case_id: item for item in loaded}
    stiffness_reference = by_id["cal_k2e9_dt1p25e6"]
    timestep_reference = by_id["cal_k1e9_dt1p25e6"]

    stiffness = [
        compare_cases(stiffness_reference, by_id[case_id])
        for case_id in ("cal_k5e8_dt1p25e6", "cal_k1e9_dt1p25e6")
    ]
    timestep = [
        compare_cases(timestep_reference, by_id[case_id])
        for case_id in ("cal_k1e9_dt2p5e6", "cal_k1e9_dt5e6")
    ]
    nearest_stiffness = stiffness[-1]
    proposed_timestep = timestep[0]
    proposed_case = by_id["cal_k1e9_dt2p5e6"]
    proposed_summary = case_summary(proposed_case)
    timestep_reference_summary = case_summary(timestep_reference)
    force_converged = all(
        nearest_stiffness[field]["relative_l2"] <= relative_l2_tolerance
        for field in ("drag_coefficient", "lift_coefficient")
    )
    timestep_converged = all(
        proposed_timestep[field]["relative_l2"] <= relative_l2_tolerance
        for field in ("drag_coefficient", "lift_coefficient")
    )
    target_rigid_enough = all((
        proposed_summary["max_target_error_over_chord"] <= target_error_tolerance,
        proposed_summary["max_target_error_over_dx"] <= target_error_dx_tolerance,
        proposed_summary["max_rms_target_error_over_dx"]
        <= rms_target_error_dx_tolerance,
    ))
    proposed_hf = proposed_summary["high_frequency_energy_fraction"]
    reference_hf = timestep_reference_summary["high_frequency_energy_fraction"]
    ringing_acceptable = all(
        proposed_hf[field] <= high_frequency_energy_tolerance
        and proposed_hf[field] <= max(
            high_frequency_reference_multiplier * reference_hf[field], 1.0e-12
        )
        for field in ("drag_coefficient", "lift_coefficient")
    )
    penalty_coupling_acceptable = (
        proposed_summary["penalty_coupling_number"] <= penalty_coupling_tolerance
    )
    return {
        "relative_l2_tolerance": relative_l2_tolerance,
        "target_error_over_chord_tolerance": target_error_tolerance,
        "target_error_over_dx_tolerance": target_error_dx_tolerance,
        "rms_target_error_over_dx_tolerance": rms_target_error_dx_tolerance,
        "high_frequency_energy_tolerance": high_frequency_energy_tolerance,
        "high_frequency_reference_multiplier": high_frequency_reference_multiplier,
        "penalty_coupling_tolerance": penalty_coupling_tolerance,
        "stiffness_comparisons": stiffness,
        "timestep_comparisons": timestep,
        "nearest_stiffness_pair_converged": force_converged,
        "proposed_dt_2p5e6_converged": timestep_converged,
        "proposed_case_target_error_acceptable": target_rigid_enough,
        "proposed_case_ringing_acceptable": ringing_acceptable,
        "proposed_case_penalty_coupling_acceptable": penalty_coupling_acceptable,
        "proposed_case_diagnostics": proposed_summary,
        "small_dt_reference_diagnostics": timestep_reference_summary,
        "production_ready": bool(
            force_converged
            and timestep_converged
            and target_rigid_enough
            and ringing_acceptable
            and penalty_coupling_acceptable
        ),
        "provisional_production_parameters": {
            "reference_cells_per_chord": 64.0,
            "reference_target_stiffness": PRODUCTION_TARGET_STIFFNESS,
            "reference_dt": PRODUCTION_DT,
            "grid_family_scaling": {
                "marker_spacing": "ds=dx/2",
                "target_stiffness": "constant",
                "dt": "constant after independent timestep calibration",
            },
        },
    }


def production_report(loaded: list[dict]) -> dict:
    grid = sorted(
        (item for item in loaded if "grid" in item["case"].study_axes),
        key=lambda item: item["case"].cells_per_chord,
    )
    boundary = sorted(
        (item for item in loaded if "boundary" in item["case"].study_axes),
        key=lambda item: item["case"].box_chords,
    )
    if len(grid) != 3 or len(boundary) != 4:
        raise ValueError("production manifest does not contain the expected families")

    finest_grid = grid[-1]
    largest_box = boundary[-1]
    return {
        "grid": {
            "interpretation": (
                "fixed-parameter spatial path: ds=dx/2 with constant dt and "
                "constant continuum target stiffness"
            ),
            "reference": finest_grid["case"].case_id,
            "comparisons": [
                compare_cases(finest_grid, candidate) for candidate in grid[:-1]
            ],
            "observed_order": {
                "drag_coefficient": observed_order(*grid, field="cd"),
                "lift_coefficient": observed_order(*grid, field="cl"),
            },
        },
        "boundary": {
            "reference": largest_box["case"].case_id,
            "comparisons": [
                compare_cases(largest_box, candidate) for candidate in boundary[:-1]
            ],
        },
    }


def write_case_csv(path: Path, summaries: list[dict]) -> None:
    fieldnames = (
        "case_id", "phase", "study_axes", "nx", "ny", "box_chords",
        "cells_per_chord", "dx", "ds", "dt", "steps", "target_stiffness",
        "target_stiffness_times_ds", "penalty_coupling_number",
        "advective_cfl", "viscous_number",
        "nominal_reynolds_number", "effective_chord", "max_target_error",
        "max_target_error_over_chord", "max_target_error_over_dx",
        "max_rms_target_error", "max_rms_target_error_over_chord",
        "max_rms_target_error_over_dx", "cd_high_frequency_energy_fraction",
        "cl_high_frequency_energy_fraction",
        "max_force_balance_linf", "cd_min", "cd_max", "cd_mean", "cd_rms",
        "cd_impulse", "cl_min", "cl_max", "cl_mean", "cl_rms", "cl_impulse",
        "directory",
    )
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for summary in summaries:
            cd = summary["drag_coefficient"]
            cl = summary["lift_coefficient"]
            writer.writerow({
                **{name: summary[name] for name in fieldnames if name in summary},
                "study_axes": "+".join(summary["study_axes"]),
                "cd_min": cd["minimum"], "cd_max": cd["maximum"],
                "cd_mean": cd["mean"], "cd_rms": cd["rms"],
                "cd_impulse": cd["impulse_over_convective_time"],
                "cl_min": cl["minimum"], "cl_max": cl["maximum"],
                "cl_mean": cl["mean"], "cl_rms": cl["rms"],
                "cl_impulse": cl["impulse_over_convective_time"],
                "cd_high_frequency_energy_fraction": summary[
                    "high_frequency_energy_fraction"
                ]["drag_coefficient"],
                "cl_high_frequency_energy_fraction": summary[
                    "high_frequency_energy_fraction"
                ]["lift_coefficient"],
            })


def analyze(output_root: Path, phase: str, *, dt: float | None,
            target_stiffness: float | None, relative_l2_tolerance: float,
            target_error_tolerance: float,
            target_error_dx_tolerance: float = 0.5,
            rms_target_error_dx_tolerance: float = 0.25,
            high_frequency_energy_tolerance: float = 0.05,
            high_frequency_reference_multiplier: float = 3.0,
            penalty_coupling_tolerance: float = 1.0) -> dict:
    cases = cases_for_phase(phase, dt=dt, target_stiffness=target_stiffness)
    loaded = [load_case(output_root, case) for case in cases]
    summaries = [case_summary(item) for item in loaded]
    exact_re = all(
        math.isclose(
            summary["nominal_reynolds_number"], REYNOLDS_NUMBER,
            rel_tol=1e-12, abs_tol=1e-12,
        )
        for summary in summaries
    )
    exact_chord = all(
        math.isclose(
            summary["effective_chord"], NOMINAL_CHORD,
            rel_tol=1e-12, abs_tol=1e-15,
        )
        for summary in summaries
    )
    comparisons = (
        calibration_report(
            loaded,
            relative_l2_tolerance,
            target_error_tolerance,
            target_error_dx_tolerance,
            rms_target_error_dx_tolerance,
            high_frequency_energy_tolerance,
            high_frequency_reference_multiplier,
            penalty_coupling_tolerance,
        )
        if phase == "calibration" else production_report(loaded)
    )
    return {
        "schema_version": 3,
        "phase": phase,
        "final_time": FINAL_TIME,
        "all_cases_exact_re128": exact_re,
        "all_cases_exact_chord": exact_chord,
        "case_summaries": summaries,
        "comparisons": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("calibration", "production"), required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--dt", type=float,
                        help="production dt used when the cases were submitted")
    parser.add_argument("--target-stiffness", type=float,
                        help="production stiffness used when the cases were submitted")
    parser.add_argument("--relative-l2-tolerance", type=float, default=0.02)
    parser.add_argument("--target-error-over-chord-tolerance", type=float, default=0.01)
    parser.add_argument("--target-error-over-dx-tolerance", type=float, default=0.5)
    parser.add_argument("--rms-target-error-over-dx-tolerance", type=float, default=0.25)
    parser.add_argument("--high-frequency-energy-tolerance", type=float, default=0.05)
    parser.add_argument("--high-frequency-reference-multiplier", type=float, default=3.0)
    parser.add_argument("--penalty-coupling-tolerance", type=float, default=1.0)
    args = parser.parse_args()

    try:
        report = analyze(
            args.output_root,
            args.phase,
            dt=args.dt,
            target_stiffness=args.target_stiffness,
            relative_l2_tolerance=args.relative_l2_tolerance,
            target_error_tolerance=args.target_error_over_chord_tolerance,
            target_error_dx_tolerance=args.target_error_over_dx_tolerance,
            rms_target_error_dx_tolerance=args.rms_target_error_over_dx_tolerance,
            high_frequency_energy_tolerance=args.high_frequency_energy_tolerance,
            high_frequency_reference_multiplier=(
                args.high_frequency_reference_multiplier
            ),
            penalty_coupling_tolerance=args.penalty_coupling_tolerance,
        )
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))

    summary_dir = args.output_root.resolve() / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    json_path = summary_dir / f"{args.phase}_report.json"
    csv_path = summary_dir / f"{args.phase}_cases.csv"
    json_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_case_csv(csv_path, report["case_summaries"])
    print(json.dumps(report["comparisons"], indent=2, sort_keys=True))
    print(f"JSON report: {json_path}")
    print(f"Case table: {csv_path}")


if __name__ == "__main__":
    main()
