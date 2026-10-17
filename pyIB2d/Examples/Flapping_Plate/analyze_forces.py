#!/usr/bin/env python3
"""Post-process force histories for a prescribed flapping plate.

The solver records the dimensional force exerted by the fluid on the plate in
laboratory coordinates.  This module keeps those raw forces and derives two
distinct families of coefficients:

* ``CH`` and ``CV`` are signed horizontal and vertical force coefficients.
* ``CD`` and ``CL`` are resolved relative to the instantaneous translational
  velocity of the plate centre.  ``CD`` is positive opposite the direction of
  plate motion and ``CL`` is positive 90 degrees counter-clockwise from it.

Both families are written with the peak-speed and RMS-speed normalizations.
Motion-relative coefficients are deliberately left blank at stroke reversal,
where the direction of translation is undefined, and in midpoint sample cells
that touch one.  ``reversal`` is set to 1 for all such samples.

Only NumPy and the standard library are required.  Typical use is::

    python analyze_forces.py runs/tilted_pilot

The output directory defaults to ``CASE_DIR/force_analysis``.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


OUTPUT_HEADER = (
    "step",
    "time",
    "time_over_period",
    "phase",
    "phase_radians",
    "cycle",
    "fx_body",
    "fy_body",
    "center_velocity_x",
    "center_velocity_y",
    "center_speed",
    "reversal",
    "reversal_cell_touches",
    "ch_umax",
    "cv_umax",
    "ch_urms",
    "cv_urms",
    "cd_umax",
    "cl_umax",
    "cd_urms",
    "cl_urms",
    "force_balance_linf",
    "max_target_error",
    "rms_target_error",
)

SUMMARY_FIELDS = (
    "fx_body",
    "fy_body",
    "ch_umax",
    "cv_umax",
    "ch_urms",
    "cv_urms",
    "cd_umax",
    "cl_umax",
    "cd_urms",
    "cl_urms",
    "force_balance_linf",
    "max_target_error",
    "rms_target_error",
)


def _read_metadata(source: Path | Mapping[str, object]) -> dict:
    if isinstance(source, Mapping):
        return dict(source)
    path = Path(source)
    if path.is_dir():
        path = path / "case_metadata.json"
    with path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError(f"metadata in {path} must be a JSON object")
    return metadata


def read_force_history(path: Path) -> dict[str, np.ndarray]:
    """Read and validate the solver's dimensional force history."""

    path = Path(path)
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"step", "time", "fx_body", "fy_body"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = required.difference(reader.fieldnames or ())
            raise ValueError(
                f"force history {path} lacks required columns: "
                + ", ".join(sorted(missing))
            )
        rows = list(reader)
        fieldnames = tuple(reader.fieldnames)
    if not rows:
        raise ValueError(f"force history {path} contains no samples")

    result: dict[str, np.ndarray] = {}
    for name in fieldnames:
        try:
            result[name] = np.asarray([float(row[name]) for row in rows], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(f"column {name!r} in {path} is not numeric") from error
    for name in required:
        if not np.all(np.isfinite(result[name])):
            raise FloatingPointError(f"column {name!r} in {path} is not finite")
    if np.any(np.diff(result["time"]) <= 0.0):
        raise ValueError("force-history times must be strictly increasing")
    if np.any(np.diff(result["step"]) <= 0.0):
        raise ValueError("force-history steps must be strictly increasing")
    return result


def _positive(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def motion_and_force_parameters(metadata: Mapping[str, object]) -> dict[str, float]:
    """Resolve the kinematic and normalization values used by this analysis."""

    motion = metadata.get("motion")
    fluid = metadata.get("fluid")
    geometry = metadata.get("geometry")
    force = metadata.get("force", {})
    if not isinstance(motion, Mapping):
        raise ValueError("metadata.motion must be an object")
    if not isinstance(fluid, Mapping):
        raise ValueError("metadata.fluid must be an object")
    if not isinstance(geometry, Mapping):
        raise ValueError("metadata.geometry must be an object")
    if not isinstance(force, Mapping):
        raise ValueError("metadata.force must be an object")

    if "period" in motion:
        period = _positive(motion["period"], "motion.period")
    elif "frequency" in motion:
        period = 1.0 / _positive(motion["frequency"], "motion.frequency")
    else:
        raise ValueError("metadata.motion needs period or frequency")
    frequency = 1.0 / period
    if "frequency" in motion and not math.isclose(
        float(motion["frequency"]), frequency, rel_tol=1.0e-10, abs_tol=0.0
    ):
        raise ValueError("motion.period and motion.frequency are inconsistent")

    amplitude_value = motion.get("stroke_amplitude", motion.get("amplitude"))
    if amplitude_value is None:
        raise ValueError("metadata.motion needs stroke_amplitude")
    amplitude = _positive(amplitude_value, "motion.stroke_amplitude")
    beta_degrees = float(
        motion.get("stroke_angle_degrees", motion.get("beta_degrees", 0.0))
    )
    if not math.isfinite(beta_degrees):
        raise ValueError("motion.stroke_angle_degrees must be finite")

    rho = _positive(fluid["rho"], "fluid.rho")
    chord = _positive(geometry["nominal_chord"], "geometry.nominal_chord")
    derived_umax = math.pi * amplitude / period
    umax = _positive(
        motion.get(
            "peak_translation_speed",
            motion.get("umax", motion.get("maximum_translation_speed", derived_umax)),
        ),
        "motion.peak_translation_speed",
    )
    if not math.isclose(umax, derived_umax, rel_tol=1.0e-10, abs_tol=1.0e-14):
        raise ValueError(
            "motion.peak_translation_speed is inconsistent with "
            "pi*motion.stroke_amplitude/motion.period"
        )
    urms = umax / math.sqrt(2.0)
    if "rms_translation_speed" in motion and not math.isclose(
        _positive(motion["rms_translation_speed"], "motion.rms_translation_speed"),
        urms,
        rel_tol=1.0e-10,
        abs_tol=1.0e-14,
    ):
        raise ValueError(
            "motion.rms_translation_speed is inconsistent with Umax/sqrt(2)"
        )
    derived_umax_scale = 0.5 * rho * umax * umax * chord
    derived_urms_scale = 0.5 * rho * urms * urms * chord
    umax_scale = _positive(
        force.get("umax_force_scale", derived_umax_scale),
        "force.umax_force_scale",
    )
    urms_scale = _positive(
        force.get("urms_force_scale", derived_urms_scale),
        "force.urms_force_scale",
    )
    for name, supplied, expected in (
        ("force.umax_force_scale", umax_scale, derived_umax_scale),
        ("force.urms_force_scale", urms_scale, derived_urms_scale),
    ):
        if not math.isclose(supplied, expected, rel_tol=1.0e-10, abs_tol=1.0e-14):
            raise ValueError(f"{name} is inconsistent with rho, speed, and chord")

    time_origin = float(motion.get("time_origin", 0.0))
    if not math.isfinite(time_origin):
        raise ValueError("motion.time_origin must be finite")
    return {
        "period": period,
        "frequency": frequency,
        "amplitude": amplitude,
        "beta_degrees": beta_degrees,
        "beta_radians": math.radians(beta_degrees),
        "time_origin": time_origin,
        "rho": rho,
        "chord": chord,
        "umax": umax,
        "urms": urms,
        "umax_force_scale": umax_scale,
        "urms_force_scale": urms_scale,
    }


def derive_coefficients(
    history: Mapping[str, np.ndarray],
    metadata: Mapping[str, object],
    *,
    reversal_speed_fraction: float = 1.0e-6,
) -> dict[str, np.ndarray]:
    """Return raw, fixed-frame, and motion-relative force time series."""

    if not math.isfinite(reversal_speed_fraction) or reversal_speed_fraction < 0.0:
        raise ValueError("reversal_speed_fraction must be finite and non-negative")
    parameters = motion_and_force_parameters(metadata)
    time = np.asarray(history["time"], dtype=float)
    step = np.asarray(history["step"], dtype=float)
    fx = np.asarray(history["fx_body"], dtype=float)
    fy = np.asarray(history["fy_body"], dtype=float)
    if not (time.ndim == step.ndim == fx.ndim == fy.ndim == 1):
        raise ValueError("force-history arrays must be one-dimensional")
    if len({time.size, step.size, fx.size, fy.size}) != 1 or time.size == 0:
        raise ValueError("force-history arrays must be nonempty and have equal length")
    if not all(np.all(np.isfinite(values)) for values in (time, step, fx, fy)):
        raise FloatingPointError("force-history core arrays must be finite")

    elapsed = time - parameters["time_origin"]
    time_over_period = elapsed / parameters["period"]
    phase = np.mod(time_over_period, 1.0)
    phase_radians = 2.0 * math.pi * phase
    cycle = np.floor(time_over_period).astype(np.int64) + 1

    stroke_x = math.cos(parameters["beta_radians"])
    stroke_y = math.sin(parameters["beta_radians"])
    signed_speed = -parameters["umax"] * np.sin(phase_radians)
    velocity_x = signed_speed * stroke_x
    velocity_y = signed_speed * stroke_y
    speed = np.abs(signed_speed)
    speed_reversal = speed <= reversal_speed_fraction * parameters["umax"]

    # Force samples are evaluated at step midpoints.  With an even number of
    # steps per half-stroke, no sample occurs at the exact reversal: the two
    # neighboring midpoint control volumes touch it instead.  Flag those cells
    # explicitly so the discontinuous motion-relative axes are never plotted
    # or compared as though their direction were defined through reversal.
    if time.size > 1:
        sample_interval = float(np.median(np.diff(time)))
        reversal_half_width_phase = 0.5 * sample_interval / parameters["period"]
    else:
        reversal_half_width_phase = 0.0
    distance_to_reversal = np.minimum.reduce((
        phase,
        np.abs(phase - 0.5),
        1.0 - phase,
    ))
    cell_touches_reversal = distance_to_reversal <= (
        reversal_half_width_phase * (1.0 + 1.0e-9) + 1.0e-15
    )
    reversal = speed_reversal | cell_touches_reversal

    # The velocity unit vector has the sign of the scalar stroke velocity.
    # Fill reversal rows with NaN before division so their direction cannot be
    # accidentally interpreted as physically meaningful.
    direction_sign = np.sign(signed_speed)
    direction_sign = np.where(reversal, np.nan, direction_sign)
    tangent_x = direction_sign * stroke_x
    tangent_y = direction_sign * stroke_y
    normal_x = -tangent_y
    normal_y = tangent_x
    drag_force = -(fx * tangent_x + fy * tangent_y)
    lift_force = fx * normal_x + fy * normal_y

    def optional(name: str) -> np.ndarray:
        if name not in history:
            return np.full(time.shape, np.nan, dtype=float)
        values = np.asarray(history[name], dtype=float)
        if values.shape != time.shape:
            raise ValueError(f"force-history column {name!r} has the wrong length")
        return values

    qmax = parameters["umax_force_scale"]
    qrms = parameters["urms_force_scale"]
    return {
        "step": step.astype(np.int64),
        "time": time,
        "time_over_period": time_over_period,
        "phase": phase,
        "phase_radians": phase_radians,
        "cycle": cycle,
        "fx_body": fx,
        "fy_body": fy,
        "center_velocity_x": velocity_x,
        "center_velocity_y": velocity_y,
        "center_speed": speed,
        "reversal": reversal,
        "reversal_cell_touches": cell_touches_reversal,
        "ch_umax": fx / qmax,
        "cv_umax": fy / qmax,
        "ch_urms": fx / qrms,
        "cv_urms": fy / qrms,
        "cd_umax": drag_force / qmax,
        "cl_umax": lift_force / qmax,
        "cd_urms": drag_force / qrms,
        "cl_urms": lift_force / qrms,
        "force_balance_linf": optional("force_balance_linf"),
        "max_target_error": optional("max_target_error"),
        "rms_target_error": optional("rms_target_error"),
    }


def _format_csv(value: object) -> str:
    if isinstance(value, (bool, np.bool_)):
        return "1" if bool(value) else "0"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    number = float(value)
    return "" if not math.isfinite(number) else format(number, ".17g")


def write_coefficient_history(path: Path, curves: Mapping[str, np.ndarray]) -> Path:
    """Write the complete derived history, using blanks for undefined values."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sample_count = len(np.asarray(curves["time"]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(OUTPUT_HEADER)
        for index in range(sample_count):
            writer.writerow([
                _format_csv(np.asarray(curves[name])[index]) for name in OUTPUT_HEADER
            ])
    return path


def _sample_widths(time: np.ndarray, start: float, end: float) -> np.ndarray:
    """Return clipped midpoint-quadrature widths for one cycle."""

    if time.size == 1:
        return np.asarray([end - start], dtype=float)
    interior = 0.5 * (time[:-1] + time[1:])
    edges = np.concatenate(([start], interior, [end]))
    edges = np.clip(edges, start, end)
    return np.maximum(np.diff(edges), 0.0)


def _series_statistics(
    values: np.ndarray, weights: np.ndarray, period: float
) -> dict[str, object]:
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(finite):
        return {
            "sample_count": 0,
            "coverage_fraction": 0.0,
            "mean": None,
            "rms": None,
            "minimum": None,
            "maximum": None,
            "peak_to_peak": None,
            "time_integral": None,
        }
    selected = values[finite]
    selected_weights = weights[finite]
    duration = float(np.sum(selected_weights))
    integral = float(np.sum(selected_weights * selected))
    mean = integral / duration
    return {
        "sample_count": int(selected.size),
        "coverage_fraction": duration / period,
        "mean": mean,
        "rms": float(np.sqrt(np.sum(selected_weights * selected * selected) / duration)),
        "minimum": float(np.min(selected)),
        "maximum": float(np.max(selected)),
        "peak_to_peak": float(np.ptp(selected)),
        "time_integral": integral,
    }


def summarize_cycle(
    curves: Mapping[str, np.ndarray], cycle_number: int, parameters: Mapping[str, float]
) -> dict[str, object]:
    cycle_values = np.asarray(curves["cycle"], dtype=int)
    selected = cycle_values == int(cycle_number)
    if not np.any(selected):
        return {
            "cycle": int(cycle_number),
            "available": False,
            "sample_count": 0,
        }
    time = np.asarray(curves["time"], dtype=float)[selected]
    start = parameters["time_origin"] + (int(cycle_number) - 1) * parameters["period"]
    end = start + parameters["period"]
    weights = _sample_widths(time, start, end)
    reversal = np.asarray(curves["reversal"], dtype=bool)[selected]
    return {
        "cycle": int(cycle_number),
        "available": True,
        "sample_count": int(time.size),
        "time_start": float(time[0]),
        "time_end": float(time[-1]),
        "reversal_sample_count": int(np.count_nonzero(reversal)),
        "reversal_sample_fraction": float(np.mean(reversal)),
        "series": {
            name: _series_statistics(
                np.asarray(curves[name], dtype=float)[selected], weights, parameters["period"]
            )
            for name in SUMMARY_FIELDS
        },
    }


def compare_cycles(
    curves: Mapping[str, np.ndarray], first_cycle: int, second_cycle: int
) -> dict[str, object] | None:
    """Compare phase-aligned coefficient waveforms for two available cycles."""

    cycle = np.asarray(curves["cycle"], dtype=int)
    first_mask = cycle == int(first_cycle)
    second_mask = cycle == int(second_cycle)
    if not np.any(first_mask) or not np.any(second_mask):
        return None
    fields = (
        "ch_umax", "cv_umax", "ch_urms", "cv_urms",
        "cd_umax", "cl_umax", "cd_urms", "cl_urms",
    )
    comparisons: dict[str, object] = {}
    for name in fields:
        phase_a = np.asarray(curves["phase"])[first_mask]
        phase_b = np.asarray(curves["phase"])[second_mask]
        values_a = np.asarray(curves[name], dtype=float)[first_mask]
        values_b = np.asarray(curves[name], dtype=float)[second_mask]
        finite_a = np.isfinite(phase_a) & np.isfinite(values_a)
        finite_b = np.isfinite(phase_b) & np.isfinite(values_b)
        if np.count_nonzero(finite_a) < 2 or np.count_nonzero(finite_b) < 2:
            comparisons[name] = None
            continue
        low = max(float(np.min(phase_a[finite_a])), float(np.min(phase_b[finite_b])))
        high = min(float(np.max(phase_a[finite_a])), float(np.max(phase_b[finite_b])))
        if high <= low:
            comparisons[name] = None
            continue
        common_phase = np.linspace(low, high, 721)
        first = _interpolate_finite_segments(phase_a, values_a, common_phase)
        second = _interpolate_finite_segments(phase_b, values_b, common_phase)
        comparable = np.isfinite(first) & np.isfinite(second)
        if np.count_nonzero(comparable) < 2:
            comparisons[name] = None
            continue
        first_valid = first[comparable]
        second_valid = second[comparable]
        difference = second_valid - first_valid
        denominator = float(np.sqrt(np.mean(first_valid * first_valid)))
        rms_difference = float(np.sqrt(np.mean(difference * difference)))
        comparisons[name] = {
            "phase_sample_count": int(common_phase.size),
            "comparison_sample_count": int(np.count_nonzero(comparable)),
            "excluded_phase_sample_count": int(
                common_phase.size - np.count_nonzero(comparable)
            ),
            "rms_difference": rms_difference,
            "relative_l2": (
                rms_difference / denominator if denominator > 0.0 else None
            ),
            "mean_difference": float(
                np.mean(second_valid) - np.mean(first_valid)
            ),
            "maximum_absolute_difference": float(np.max(np.abs(difference))),
        }
    return {
        "first_cycle": int(first_cycle),
        "second_cycle": int(second_cycle),
        "series": comparisons,
    }


def _interpolate_finite_segments(
    source_x: np.ndarray,
    source_y: np.ndarray,
    query_x: np.ndarray,
) -> np.ndarray:
    """Interpolate within finite runs without bridging intentional gaps.

    Motion-relative force histories contain NaNs in midpoint cells touching a
    stroke reversal.  Interpolating all finite values at once would reconnect
    the two half-strokes across that undefined change of basis.  This helper
    treats each contiguous finite run as an independent interpolation segment.
    """

    source_x = np.asarray(source_x, dtype=float)
    source_y = np.asarray(source_y, dtype=float)
    query_x = np.asarray(query_x, dtype=float)
    if source_x.ndim != 1 or source_y.ndim != 1 or query_x.ndim != 1:
        raise ValueError("interpolation arrays must be one-dimensional")
    if source_x.size != source_y.size:
        raise ValueError("source interpolation arrays must have equal length")

    result = np.full(query_x.shape, np.nan, dtype=float)
    finite = np.isfinite(source_x) & np.isfinite(source_y)
    padded = np.concatenate(([False], finite, [False]))
    transitions = np.flatnonzero(padded[1:] != padded[:-1])
    for start, stop in transitions.reshape(-1, 2):
        if stop - start < 2:
            continue
        segment_x = source_x[start:stop]
        segment_y = source_y[start:stop]
        if np.any(np.diff(segment_x) <= 0.0):
            raise ValueError("finite interpolation segments must increase")
        selected = (query_x >= segment_x[0]) & (query_x <= segment_x[-1])
        result[selected] = np.interp(
            query_x[selected], segment_x, segment_y
        )
    return result


def build_summary(
    curves: Mapping[str, np.ndarray],
    metadata: Mapping[str, object],
    analysis_cycles: Sequence[int] = (3, 4),
    *,
    reversal_speed_fraction: float = 1.0e-6,
) -> dict[str, object]:
    parameters = motion_and_force_parameters(metadata)
    cycles = tuple(dict.fromkeys(int(value) for value in analysis_cycles))
    if not cycles or any(value < 1 for value in cycles):
        raise ValueError("analysis cycles must be positive integers")
    cycle_summaries = [summarize_cycle(curves, value, parameters) for value in cycles]
    comparison = (
        compare_cycles(curves, cycles[0], cycles[1]) if len(cycles) >= 2 else None
    )
    available = sorted(set(int(value) for value in np.asarray(curves["cycle"])))
    return {
        "schema_version": 1,
        "case_name": metadata.get("case_name", "flapping_plate"),
        "interpretation": (
            "Forces are fluid-on-body forces per unit span. CH and CV are signed "
            "laboratory x/y coefficients. CD is positive opposite plate-centre "
            "translation and CL is positive to the left of that translation. "
            "Motion-relative values are undefined at flagged stroke reversals; "
            "the flag also covers midpoint sample cells that touch a reversal."
        ),
        "normalizations": {
            "umax": {
                "speed": parameters["umax"],
                "force_scale": parameters["umax_force_scale"],
                "definition": "0.5*rho*Umax^2*chord",
            },
            "urms": {
                "speed": parameters["urms"],
                "force_scale": parameters["urms_force_scale"],
                "definition": "0.5*rho*Urms^2*chord, Urms=Umax/sqrt(2)",
            },
        },
        "kinematics": {
            "period": parameters["period"],
            "frequency": parameters["frequency"],
            "stroke_amplitude": parameters["amplitude"],
            "stroke_angle_degrees": parameters["beta_degrees"],
            "velocity_definition": (
                "-Umax*sin(2*pi*phase)*(cos(beta), sin(beta))"
            ),
            "reversal_speed_fraction": reversal_speed_fraction,
        },
        "available_cycles": available,
        "requested_analysis_cycles": list(cycles),
        "cycle_summaries": cycle_summaries,
        "cycle_comparison": comparison,
    }


def _plot_range(series: Iterable[np.ndarray]) -> tuple[float, float]:
    finite_arrays = []
    for values in series:
        selected = np.asarray(values, dtype=float)
        selected = selected[np.isfinite(selected)]
        if selected.size:
            finite_arrays.append(selected)
    if not finite_arrays:
        return -1.0, 1.0
    lower = min(0.0, min(float(np.min(values)) for values in finite_arrays))
    upper = max(0.0, max(float(np.max(values)) for values in finite_arrays))
    if math.isclose(lower, upper):
        padding = max(abs(lower), 1.0) * 0.05
    else:
        padding = 0.05 * (upper - lower)
    return lower - padding, upper + padding


def _svg_segments(
    x: np.ndarray,
    y: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    box: tuple[float, float, float, float],
) -> list[str]:
    left, top, width, height = box
    xmin, xmax = x_range
    ymin, ymax = y_range
    valid = np.isfinite(x) & np.isfinite(y)
    segments: list[str] = []
    start = None
    for index, is_valid in enumerate(np.append(valid, False)):
        if is_valid and start is None:
            start = index
        elif not is_valid and start is not None:
            chosen_x = x[start:index]
            chosen_y = y[start:index]
            px = left + (chosen_x - xmin) / (xmax - xmin) * width
            py = top + (ymax - chosen_y) / (ymax - ymin) * height
            segments.append(" ".join(
                ("M" if offset == 0 else "L") + f"{a:.3f},{b:.3f}"
                for offset, (a, b) in enumerate(zip(px, py))
            ))
            start = None
    return segments


def _svg_panel(
    elements: list[str],
    x: np.ndarray,
    series: Sequence[tuple[str, np.ndarray, str]],
    box: tuple[float, float, float, float],
    ylabel: str,
) -> None:
    left, top, width, height = box
    xmin, xmax = float(np.min(x)), float(np.max(x))
    if math.isclose(xmin, xmax):
        xmin, xmax = xmin - 0.5, xmax + 0.5
    ymin, ymax = _plot_range(values for _name, values, _color in series)
    elements.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="#fff" stroke="#222"/>'
    )
    for fraction in np.linspace(0.0, 1.0, 5):
        px = left + fraction * width
        py = top + fraction * height
        xv = xmin + fraction * (xmax - xmin)
        yv = ymax - fraction * (ymax - ymin)
        elements.extend((
            f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" '
            f'y2="{top + height}" stroke="#e4e4e4"/>',
            f'<line x1="{left}" y1="{py:.2f}" x2="{left + width}" '
            f'y2="{py:.2f}" stroke="#e4e4e4"/>',
            f'<text x="{px:.2f}" y="{top + height + 19}" text-anchor="middle" '
            f'class="tick">{xv:.3g}</text>',
            f'<text x="{left - 9}" y="{py + 4:.2f}" text-anchor="end" '
            f'class="tick">{yv:.3g}</text>',
        ))
    for name, values, color in series:
        for path_data in _svg_segments(
            np.asarray(x), np.asarray(values), (xmin, xmax), (ymin, ymax), box
        ):
            elements.append(
                f'<path d="{path_data}" fill="none" stroke="{color}" '
                'stroke-width="1.7" stroke-linejoin="round"/>'
            )
    elements.append(
        f'<text x="22" y="{top + height / 2}" text-anchor="middle" '
        f'transform="rotate(-90 22 {top + height / 2})" class="axis">'
        f'{html.escape(ylabel)}</text>'
    )
    for index, (name, _values, color) in enumerate(series):
        x0 = left + 12 + index * 230
        elements.extend((
            f'<line x1="{x0}" y1="{top + 18}" x2="{x0 + 27}" '
            f'y2="{top + 18}" stroke="{color}" stroke-width="3"/>',
            f'<text x="{x0 + 34}" y="{top + 23}" class="legend">'
            f'{html.escape(name)}</text>',
        ))


def write_force_svg(
    path: Path, curves: Mapping[str, np.ndarray], metadata: Mapping[str, object]
) -> Path:
    """Write raw and coefficient histories against nondimensional time."""

    x = np.asarray(curves["time_over_period"], dtype=float)
    width, height = 1120, 1020
    elements = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<style>.title{font:bold 21px sans-serif}.axis{font:15px sans-serif}'
        '.tick{font:11px sans-serif;fill:#333}.legend{font:13px sans-serif}</style>',
        '<rect width="100%" height="100%" fill="#fafafa"/>',
        f'<text x="560" y="31" text-anchor="middle" class="title">'
        f'{html.escape(str(metadata.get("case_name", "flapping plate")))}: '
        'force history</text>',
    ]
    _svg_panel(
        elements, x,
        (("Fx body", curves["fx_body"], "#1f77b4"),
         ("Fy body", curves["fy_body"], "#d62728")),
        (85.0, 52.0, 995.0, 245.0), "raw force per unit span",
    )
    _svg_panel(
        elements, x,
        (("CH (Urms)", curves["ch_urms"], "#1f77b4"),
         ("CV (Urms)", curves["cv_urms"], "#d62728")),
        (85.0, 370.0, 995.0, 245.0), "fixed-frame coefficient",
    )
    _svg_panel(
        elements, x,
        (("CD (Urms)", curves["cd_urms"], "#2ca02c"),
         ("CL (Urms)", curves["cl_urms"], "#9467bd")),
        (85.0, 688.0, 995.0, 245.0), "motion-relative coefficient",
    )
    elements.extend((
        '<text x="582" y="983" text-anchor="middle" class="axis">time / period</text>',
        '<text x="582" y="1005" text-anchor="middle" class="tick">'
        'Motion-relative curves contain intentional gaps at stroke reversal.</text>',
        '</svg>',
    ))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path


def write_cycle_overlay_svg(
    path: Path,
    curves: Mapping[str, np.ndarray],
    metadata: Mapping[str, object],
    analysis_cycles: Sequence[int],
) -> Path:
    """Overlay selected cycles against phase for a visual periodicity check."""

    palette = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    cycle_values = np.asarray(curves["cycle"], dtype=int)
    selections = [
        (int(number), cycle_values == int(number)) for number in analysis_cycles
        if np.any(cycle_values == int(number))
    ]
    path = Path(path)
    if not selections:
        return path

    width, height = 1120, 710
    elements = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<style>.title{font:bold 21px sans-serif}.axis{font:15px sans-serif}'
        '.tick{font:11px sans-serif;fill:#333}.legend{font:13px sans-serif}</style>',
        '<rect width="100%" height="100%" fill="#fafafa"/>',
        f'<text x="560" y="31" text-anchor="middle" class="title">'
        f'{html.escape(str(metadata.get("case_name", "flapping plate")))}: '
        'cycle overlay</text>',
    ]
    fixed_series = []
    relative_series = []
    for index, (number, mask) in enumerate(selections):
        color_a = palette[(2 * index) % len(palette)]
        color_b = palette[(2 * index + 1) % len(palette)]
        fixed_series.extend((
            (f"C{number} CH", np.where(mask, curves["ch_urms"], np.nan), color_a),
            (f"C{number} CV", np.where(mask, curves["cv_urms"], np.nan), color_b),
        ))
        relative_series.extend((
            (f"C{number} CD", np.where(mask, curves["cd_urms"], np.nan), color_a),
            (f"C{number} CL", np.where(mask, curves["cl_urms"], np.nan), color_b),
        ))
    phase = np.asarray(curves["phase"], dtype=float)
    _svg_panel(elements, phase, tuple(fixed_series), (85.0, 52.0, 995.0, 245.0),
               "fixed-frame coefficient")
    _svg_panel(elements, phase, tuple(relative_series), (85.0, 370.0, 995.0, 245.0),
               "motion-relative coefficient")
    elements.extend((
        '<text x="582" y="674" text-anchor="middle" class="axis">cycle phase</text>',
        '</svg>',
    ))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path


def analyze_case(
    case_dir: Path,
    output_dir: Path | None = None,
    *,
    analysis_cycles: Sequence[int] = (3, 4),
    reversal_speed_fraction: float = 1.0e-6,
) -> dict[str, object]:
    """Create the coefficient table, summaries, and dependency-light plots."""

    case_dir = Path(case_dir).resolve()
    output_dir = (
        Path(output_dir).resolve() if output_dir is not None
        else case_dir / "force_analysis"
    )
    metadata = _read_metadata(case_dir)
    history = read_force_history(case_dir / "force_history.csv")
    curves = derive_coefficients(
        history, metadata, reversal_speed_fraction=reversal_speed_fraction
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficient_path = write_coefficient_history(
        output_dir / "force_coefficients.csv", curves
    )
    summary = build_summary(
        curves, metadata, analysis_cycles,
        reversal_speed_fraction=reversal_speed_fraction,
    )
    summary_path = output_dir / "cycle_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    force_plot_path = write_force_svg(
        output_dir / "force_coefficients.svg", curves, metadata
    )
    cycle_label = "_".join(str(int(value)) for value in analysis_cycles)
    cycle_plot_path = write_cycle_overlay_svg(
        output_dir / f"cycle_{cycle_label}_overlay.svg",
        curves,
        metadata,
        analysis_cycles,
    )
    artifacts = {
        "coefficient_history": str(coefficient_path),
        "cycle_summary": str(summary_path),
        "force_plot": str(force_plot_path),
        "cycle_overlay": str(cycle_plot_path) if cycle_plot_path.is_file() else None,
    }
    return {"summary": summary, "artifacts": artifacts}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--analysis-cycles", type=int, nargs="+", default=(3, 4),
        help="one or more one-based cycle numbers to summarize (default: 3 4)",
    )
    parser.add_argument(
        "--reversal-speed-fraction", type=float, default=1.0e-6,
        help="flag motion-relative coefficients below this fraction of Umax",
    )
    args = parser.parse_args()
    try:
        result = analyze_case(
            args.case_dir,
            args.output_dir,
            analysis_cycles=args.analysis_cycles,
            reversal_speed_fraction=args.reversal_speed_fraction,
        )
    except (FileNotFoundError, FloatingPointError, KeyError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result["artifacts"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
