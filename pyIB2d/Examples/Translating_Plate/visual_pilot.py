#!/usr/bin/env python3
"""Create dependency-light visual diagnostics for a translating-plate run.

Only NumPy and the Python standard library are required.  The intended solver
integration is a state callback of the form::

    snapshots = VisualSnapshotWriter.from_case(
        case_dir, every_steps=20,
    )
    snapshots.record(
        step=cter,
        time=current_time,
        target_time=current_time,
        velocity_x=U,
        velocity_y=V,
        lagrangian_positions=np.column_stack((xLag, yLag)),
        target_info=target_info,
        dx=dx,
        dy=dy,
        lx=Lx,
        ly=Ly,
    )

``record`` computes vorticity with periodic centered differences and writes a
compressed ``.npz`` file.  After the run, :func:`build_visual_artifacts`
renders all snapshots with one shared vorticity scale, creates a browser-
playable HTML animation, plots force and target error histories as SVG, and
writes descriptive oscillation metrics as JSON.  No diagnostic is assigned an
automatic pass/fail threshold.
"""

from __future__ import annotations

import argparse
import binascii
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
from typing import Iterable, Mapping, Sequence
import zlib

import numpy as np

try:
    from case_config import translation_distance
except ImportError:
    # The flapping example deliberately reuses this dependency-light renderer.
    # In that process its own ``case_config`` is already loaded; its subclass
    # handles target-time alignment analytically and never needs this helper.
    translation_distance = None


SNAPSHOT_SCHEMA_VERSION = 1
ARTIFACT_SCHEMA_VERSION = 1


def periodic_vorticity(u: np.ndarray, v: np.ndarray,
                       dx: float, dy: float) -> np.ndarray:
    """Return ``dv/dx - du/dy`` using periodic centered differences.

    Solver velocity arrays are expected in ``(ny, nx)`` order.  The returned
    array has the same order and is not transposed.
    """

    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if u.ndim != 2 or v.ndim != 2 or u.shape != v.shape:
        raise ValueError("u and v must be same-shaped two-dimensional arrays")
    if u.shape[0] < 2 or u.shape[1] < 2:
        raise ValueError("velocity arrays must contain at least 2x2 cells")
    if not np.isfinite(dx) or not np.isfinite(dy) or dx <= 0.0 or dy <= 0.0:
        raise ValueError("dx and dy must be finite and positive")
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
        raise FloatingPointError("velocity arrays contain non-finite values")

    dvdx = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2.0 * dx)
    dudy = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2.0 * dy)
    return dvdx - dudy


def _read_metadata(source: Path | Mapping[str, object]) -> dict:
    if isinstance(source, Mapping):
        return dict(source)
    path = Path(source)
    if path.is_dir():
        path = path / "case_metadata.json"
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


class VisualSnapshotWriter:
    """A testable state-callback that saves periodic-vorticity snapshots."""

    def __init__(self, snapshot_dir: Path, metadata: Mapping[str, object],
                 every_steps: int = 1):
        if int(every_steps) != every_steps or every_steps <= 0:
            raise ValueError("every_steps must be a positive integer")
        self.snapshot_dir = Path(snapshot_dir).resolve()
        self.metadata = dict(metadata)
        domain = self.metadata["domain"]
        grid = self.metadata["grid"]
        self.lx = float(domain["lx"])
        self.ly = float(domain["ly"])
        self.nx = int(grid["nx"])
        self.ny = int(grid["ny"])
        self.dx = float(grid.get("dx", self.lx / self.nx))
        self.dy = float(grid.get("dy", self.ly / self.ny))
        self.every_steps = int(every_steps)
        time_metadata = self.metadata.get("time", {})
        self.total_steps = (
            int(time_metadata["steps"])
            if isinstance(time_metadata, Mapping) and "steps" in time_metadata
            else None
        )
        if min(self.lx, self.ly, self.dx, self.dy) <= 0.0:
            raise ValueError("domain lengths and cell widths must be positive")
        self.snapshot_dir.mkdir(parents=True, exist_ok=True)

    @classmethod
    def from_case(cls, case_dir: Path, every_steps: int = 1,
                  snapshot_dir: Path | None = None) -> "VisualSnapshotWriter":
        case_dir = Path(case_dir).resolve()
        metadata = _read_metadata(case_dir)
        destination = (
            Path(snapshot_dir).resolve()
            if snapshot_dir is not None else case_dir / "visual_snapshots"
        )
        return cls(destination, metadata, every_steps=every_steps)

    def should_record(self, step: int, final: bool = False) -> bool:
        return bool(final or int(step) % self.every_steps == 0)

    def _target_at_state_time(
        self, target_positions: np.ndarray, target_time: float, state_time: float
    ) -> tuple[np.ndarray, float]:
        """Align this case's midpoint target with the saved full-step state."""

        positions = np.asarray(target_positions, dtype=float).copy()
        if math.isclose(target_time, state_time, rel_tol=0.0, abs_tol=1.0e-15):
            return positions, target_time
        motion = self.metadata.get("motion", {})
        if not isinstance(motion, Mapping) or motion.get("kind") != (
            "ramped_translation_x"
        ):
            # Retain the honest source timestamp for an unknown motion law.
            return positions, target_time
        required = ("speed", "ramp_rate", "stop_time")
        if any(name not in motion for name in required):
            raise ValueError("ramped-translation metadata is incomplete")
        if translation_distance is None:
            raise RuntimeError("translating-case displacement helper is unavailable")
        displacement = lambda value: translation_distance(
            value,
            speed=float(motion["speed"]),
            ramp_rate=float(motion["ramp_rate"]),
            stop_time=float(motion["stop_time"]),
        )
        positions[:, 0] += displacement(state_time) - displacement(target_time)
        return positions, state_time

    def record(self, *, step: int, time: float, target_time: float | None = None,
               velocity_x: np.ndarray | None = None,
               velocity_y: np.ndarray | None = None,
               lagrangian_positions: np.ndarray, target_info: np.ndarray,
               dx: float | None = None, dy: float | None = None,
               lx: float | None = None, ly: float | None = None,
               final: bool = False, **_unused: object) -> Path | None:
        """Save one state when its step is on the configured cadence.

        Extra keyword arguments are accepted so this can coexist with a richer
        ``IBM_Driver`` callback payload.  ``target_info`` follows the existing
        IB2d convention: columns are marker id, target x, target y, stiffness.
        """

        step = int(step)
        time = float(time)
        target_time = time if target_time is None else float(target_time)
        if not self.should_record(step, final=final):
            return None
        if not np.isfinite(time) or not np.isfinite(target_time):
            raise FloatingPointError("snapshot or target time is not finite")

        callback_geometry = {
            "dx": (dx, self.dx), "dy": (dy, self.dy),
            "lx": (lx, self.lx), "ly": (ly, self.ly),
        }
        for name, (supplied, expected) in callback_geometry.items():
            if supplied is not None and not math.isclose(
                    float(supplied), expected, rel_tol=1.0e-12, abs_tol=0.0):
                raise ValueError(
                    f"callback {name}={supplied} disagrees with metadata {expected}"
                )
        if velocity_x is None or velocity_y is None:
            raise TypeError("velocity_x and velocity_y are required")

        u = np.asarray(velocity_x, dtype=float)
        v = np.asarray(velocity_y, dtype=float)
        if u.shape != (self.ny, self.nx) or v.shape != (self.ny, self.nx):
            raise ValueError(
                f"velocity shape must be {(self.ny, self.nx)}, got "
                f"{u.shape} and {v.shape}"
            )
        positions = np.asarray(lagrangian_positions, dtype=float)
        targets = np.asarray(target_info, dtype=float)
        if positions.ndim != 2 or positions.shape[1] != 2:
            raise ValueError("lagrangian_positions must have shape (n, 2)")
        if targets.ndim != 2 or targets.shape[1] < 3:
            raise ValueError("target_info must have at least three columns")
        ids = targets[:, 0].astype(np.int64)
        if ids.size and (ids.min() < 0 or ids.max() >= positions.shape[0]):
            raise IndexError("target marker id is outside lagrangian_positions")
        actual = positions[ids, :]
        target_positions, target_time = self._target_at_state_time(
            targets[:, 1:3], target_time, time
        )
        if not np.all(np.isfinite(actual)) or not np.all(
                np.isfinite(target_positions)):
            raise FloatingPointError("plate coordinates contain non-finite values")

        omega = periodic_vorticity(u, v, self.dx, self.dy).astype(
            np.float32, copy=False
        )
        final_path = self.snapshot_dir / f"snapshot_step{step:09d}.npz"
        temporary_path = final_path.with_suffix(".npz.tmp")
        with temporary_path.open("wb") as handle:
            np.savez_compressed(
                handle,
                schema_version=np.asarray(SNAPSHOT_SCHEMA_VERSION, dtype=np.int64),
                step=np.asarray(step, dtype=np.int64),
                time=np.asarray(time, dtype=np.float64),
                target_time=np.asarray(target_time, dtype=np.float64),
                lx=np.asarray(self.lx, dtype=np.float64),
                ly=np.asarray(self.ly, dtype=np.float64),
                dx=np.asarray(self.dx, dtype=np.float64),
                dy=np.asarray(self.dy, dtype=np.float64),
                vorticity=omega,
                actual_positions=np.asarray(actual, dtype=np.float64),
                target_positions=np.asarray(target_positions, dtype=np.float64),
                target_ids=ids,
            )
        os.replace(temporary_path, final_path)
        if (
            step == 0
            or step % (10 * self.every_steps) == 0
            or (self.total_steps is not None and step == self.total_steps)
        ):
            print(f"Saved visual snapshot at step {step}", flush=True)
        return final_path

    __call__ = record


def _load_snapshot(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as stored:
        required = {
            "schema_version", "step", "time", "target_time", "lx", "ly", "dx", "dy",
            "vorticity", "actual_positions", "target_positions", "target_ids",
        }
        missing = required.difference(stored.files)
        if missing:
            raise ValueError(f"{path} is missing: {', '.join(sorted(missing))}")
        result = {name: np.asarray(stored[name]) for name in required}
    if int(result["schema_version"]) != SNAPSHOT_SCHEMA_VERSION:
        raise ValueError(f"unsupported snapshot schema in {path}")
    omega = result["vorticity"]
    if omega.ndim != 2 or not np.all(np.isfinite(omega)):
        raise ValueError(f"invalid vorticity in {path}")
    return result


def _png_chunk(kind: bytes, payload: bytes) -> bytes:
    checksum = binascii.crc32(kind)
    checksum = binascii.crc32(payload, checksum) & 0xFFFFFFFF
    return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", checksum)


def write_rgb_png(path: Path, pixels: np.ndarray) -> Path:
    """Write an 8-bit RGB array as PNG without Pillow or matplotlib."""

    pixels = np.asarray(pixels)
    if pixels.ndim != 3 or pixels.shape[2] != 3:
        raise ValueError("pixels must have shape (height, width, 3)")
    if pixels.dtype != np.uint8:
        pixels = np.clip(pixels, 0, 255).astype(np.uint8)
    height, width, _ = pixels.shape
    if width <= 0 or height <= 0:
        raise ValueError("PNG dimensions must be positive")
    scanlines = b"".join(b"\x00" + row.tobytes() for row in pixels)
    header = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    encoded = (
        b"\x89PNG\r\n\x1a\n"
        + _png_chunk(b"IHDR", header)
        + _png_chunk(b"IDAT", zlib.compress(scanlines, level=6))
        + _png_chunk(b"IEND", b"")
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(encoded)
    return path


def vorticity_rgb(vorticity: np.ndarray, limit: float) -> np.ndarray:
    """Map signed vorticity to a deterministic blue-white-red palette."""

    values = np.asarray(vorticity, dtype=float)
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("vorticity color limit must be finite and positive")
    scaled = np.clip(values / float(limit), -1.0, 1.0)
    blue = np.asarray([33.0, 102.0, 172.0])
    white = np.asarray([247.0, 247.0, 247.0])
    red = np.asarray([178.0, 24.0, 43.0])
    negative_weight = np.clip(-scaled, 0.0, 1.0)[..., None]
    positive_weight = np.clip(scaled, 0.0, 1.0)[..., None]
    rgb = white + negative_weight * (blue - white) + positive_weight * (red - white)
    return np.rint(rgb).astype(np.uint8)


def _draw_disc(image: np.ndarray, x: int, y: int, radius: int,
               color: Sequence[int]) -> None:
    height, width, _ = image.shape
    radius = max(int(radius), 0)
    x0, x1 = max(0, x - radius), min(width - 1, x + radius)
    y0, y1 = max(0, y - radius), min(height - 1, y + radius)
    radius_squared = radius * radius
    for py in range(y0, y1 + 1):
        for px in range(x0, x1 + 1):
            if (px - x) ** 2 + (py - y) ** 2 <= radius_squared:
                image[py, px, :] = color


def _draw_line(image: np.ndarray, start: tuple[int, int], end: tuple[int, int],
               color: Sequence[int], width: int) -> None:
    x0, y0 = start
    x1, y1 = end
    dx = abs(x1 - x0)
    sx = 1 if x0 < x1 else -1
    dy = -abs(y1 - y0)
    sy = 1 if y0 < y1 else -1
    error = dx + dy
    radius = max(0, int(width) // 2)
    while True:
        _draw_disc(image, x0, y0, radius, color)
        if x0 == x1 and y0 == y1:
            return
        doubled = 2 * error
        if doubled >= dy:
            error += dy
            x0 += sx
        if doubled <= dx:
            error += dx
            y0 += sy


def _to_pixel(point: Sequence[float], lx: float, ly: float,
              width: int, height: int) -> tuple[int, int]:
    x = float(point[0]) / lx * (width - 1)
    y = (1.0 - float(point[1]) / ly) * (height - 1)
    return int(round(x)), int(round(y))


def _draw_periodic_polyline(image: np.ndarray, points: np.ndarray,
                            lx: float, ly: float, color: Sequence[int],
                            width: int) -> None:
    points = np.asarray(points, dtype=float)
    if points.shape[0] == 0:
        return
    wrapped = np.column_stack((points[:, 0] % lx, points[:, 1] % ly))
    height, pixel_width, _ = image.shape
    if wrapped.shape[0] == 1:
        _draw_disc(
            image, *_to_pixel(wrapped[0], lx, ly, pixel_width, height),
            max(1, width), color,
        )
        return
    for point_a, point_b_wrapped in zip(wrapped[:-1], wrapped[1:]):
        delta = point_b_wrapped - point_a
        delta -= np.asarray([lx, ly]) * np.round(delta / np.asarray([lx, ly]))
        point_b = point_a + delta
        for offset_x in (-lx, 0.0, lx):
            for offset_y in (-ly, 0.0, ly):
                offset = np.asarray([offset_x, offset_y])
                _draw_line(
                    image,
                    _to_pixel(point_a + offset, lx, ly, pixel_width, height),
                    _to_pixel(point_b + offset, lx, ly, pixel_width, height),
                    color,
                    width,
                )


def render_snapshot(snapshot: Mapping[str, np.ndarray], output_path: Path,
                    vorticity_limit: float, max_width: int = 960,
                    min_width: int = 480) -> Path:
    """Render one snapshot with actual and target plate overlays."""

    omega = np.asarray(snapshot["vorticity"], dtype=float)
    ny, nx = omega.shape
    lx = float(snapshot["lx"])
    ly = float(snapshot["ly"])
    if max_width <= 0 or min_width <= 0 or min_width > max_width:
        raise ValueError("render widths must satisfy 0 < min_width <= max_width")
    output_width = max(min_width, min(max_width, nx))
    output_height = max(1, int(round(output_width * ly / lx)))
    output_height = min(max_width, output_height)
    source_x = np.rint(np.linspace(0, nx - 1, output_width)).astype(int)
    source_y = np.rint(np.linspace(ny - 1, 0, output_height)).astype(int)
    sampled = omega[source_y[:, None], source_x[None, :]]
    image = vorticity_rgb(sampled, vorticity_limit)

    targets = np.asarray(snapshot["target_positions"], dtype=float)
    actual = np.asarray(snapshot["actual_positions"], dtype=float)
    _draw_periodic_polyline(image, targets, lx, ly, (250, 210, 32), 5)
    _draw_periodic_polyline(image, actual, lx, ly, (15, 15, 15), 7)
    _draw_periodic_polyline(image, actual, lx, ly, (250, 250, 250), 3)
    return write_rgb_png(output_path, image)


def _read_force_history(path: Path) -> dict[str, np.ndarray]:
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "time", "fx_body", "fy_body", "max_target_error",
            "rms_target_error",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"force history {path} lacks required columns")
        rows = list(reader)
    if not rows:
        raise ValueError(f"force history {path} contains no samples")
    result = {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in reader.fieldnames
    }
    if not all(np.all(np.isfinite(values)) for values in result.values()):
        raise FloatingPointError(f"force history {path} contains non-finite values")
    if np.any(np.diff(result["time"]) <= 0.0):
        raise ValueError("force-history times must be strictly increasing")
    return result


def coefficient_history(history: Mapping[str, np.ndarray],
                        metadata: Mapping[str, object]) -> dict[str, np.ndarray]:
    chord = float(metadata["geometry"]["nominal_chord"])
    dx = float(metadata["grid"]["dx"])
    speed = abs(float(metadata["motion"]["speed"]))
    scale = float(metadata["force"].get(
        "coefficient_force_scale",
        0.5 * float(metadata["fluid"]["rho"]) * speed * speed * chord,
    ))
    if scale <= 0.0 or chord <= 0.0 or dx <= 0.0:
        raise ValueError("force scale, chord, and dx must be positive")
    time = np.asarray(history["time"], dtype=float)
    force_metadata = metadata.get("force", {})
    component_mode = (
        force_metadata.get("visual_component_mode", "translation_x_drag_lift")
        if isinstance(force_metadata, Mapping)
        else "translation_x_drag_lift"
    )
    if component_mode == "global_horizontal_vertical":
        first_component = np.asarray(history["fx_body"], dtype=float) / scale
        second_component = np.asarray(history["fy_body"], dtype=float) / scale
    elif component_mode == "translation_x_drag_lift":
        first_component = -np.asarray(history["fx_body"], dtype=float) / scale
        second_component = np.asarray(history["fy_body"], dtype=float) / scale
    else:
        raise ValueError(f"unknown visual force-component mode: {component_mode}")
    return {
        "time": time,
        "tau": time * speed / chord,
        # The historical key names are retained for artifact-schema
        # compatibility.  Labels in the SVG identify the actual requested
        # components (CD/CL for translation or CH/CV for flapping).
        "CD": first_component,
        "CL": second_component,
        "max_target_error_over_chord": (
            np.asarray(history["max_target_error"], dtype=float) / chord
        ),
        "rms_target_error_over_chord": (
            np.asarray(history["rms_target_error"], dtype=float) / chord
        ),
        "max_target_error_over_dx": (
            np.asarray(history["max_target_error"], dtype=float) / dx
        ),
        "rms_target_error_over_dx": (
            np.asarray(history["rms_target_error"], dtype=float) / dx
        ),
    }


def _linear_detrend(values: np.ndarray) -> np.ndarray:
    if values.size < 2:
        return values - np.mean(values)
    coordinate = np.linspace(-1.0, 1.0, values.size)
    design = np.column_stack((coordinate, np.ones(values.size)))
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def _polynomial_detrend(values: np.ndarray, degree: int) -> np.ndarray:
    coordinate = np.linspace(-1.0, 1.0, values.size)
    fitted_degree = min(int(degree), values.size - 1)
    return values - np.polyval(
        np.polyfit(coordinate, values, fitted_degree), coordinate
    )


def series_oscillation_metrics(time: np.ndarray,
                               values: np.ndarray) -> dict[str, object]:
    """Return descriptive difference and spectral metrics, without a verdict."""

    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if time.ndim != 1 or values.ndim != 1 or time.size != values.size:
        raise ValueError("time and values must be same-length one-dimensional arrays")
    if values.size == 0:
        raise ValueError("cannot diagnose an empty series")
    differences = np.diff(values)
    second_differences = np.diff(values, n=2)
    centered = values - np.mean(values)
    standard_deviation = float(np.sqrt(np.mean(centered * centered)))
    scale = max(standard_deviation, np.finfo(float).eps)

    nonzero_signs = np.sign(differences[np.abs(differences) > 0.0])
    turning_points = int(np.count_nonzero(np.diff(nonzero_signs))) \
        if nonzero_signs.size > 1 else 0

    spectral: dict[str, object]
    if values.size >= 4 and time[-1] > time[0]:
        uniform_time = np.linspace(time[0], time[-1], values.size)
        uniform_values = np.interp(uniform_time, time, values)
        detrended = _linear_detrend(uniform_values)
        sample_interval = float(uniform_time[1] - uniform_time[0])
        frequencies = np.fft.rfftfreq(values.size, d=sample_interval)
        power = np.abs(np.fft.rfft(detrended)) ** 2
        if power.size:
            power[0] = 0.0
        total_power = float(np.sum(power))
        nyquist = float(0.5 / sample_interval)

        def band_fraction(start_fraction: float) -> float:
            if total_power == 0.0:
                return 0.0
            return float(np.sum(power[frequencies >= start_fraction * nyquist]) /
                         total_power)

        peak_index = int(np.argmax(power)) if total_power > 0.0 else 0
        centroid = (
            float(np.sum(frequencies * power) / total_power)
            if total_power > 0.0 else 0.0
        )
        spectral = {
            "method": "linear detrend, uniform interpolation, real FFT",
            "nyquist_frequency": nyquist,
            "peak_frequency": float(frequencies[peak_index]),
            "spectral_centroid": centroid,
            "spectral_centroid_over_nyquist": (
                centroid / nyquist if nyquist > 0.0 else 0.0
            ),
            "upper_half_nyquist_power_fraction": band_fraction(0.5),
            "upper_quarter_nyquist_power_fraction": band_fraction(0.75),
        }
        step_power = np.abs(
            np.fft.rfft(_polynomial_detrend(values, degree=3))
        ) ** 2
        step_frequencies = np.fft.rfftfreq(values.size, d=1.0)
        if step_power.size:
            step_power[0] = 0.0
        step_total_power = float(np.sum(step_power))
        spectral.update({
            "step_scale_method": "cubic detrend, real FFT",
            "high_frequency_cutoff_cycles_per_step": 0.1,
            "energy_fraction_at_or_above_0p1_cycles_per_step": (
                float(np.sum(step_power[step_frequencies >= 0.1]) /
                      step_total_power)
                if step_total_power > 0.0 else 0.0
            ),
        })
    else:
        spectral = {
            "method": "insufficient samples for FFT",
            "nyquist_frequency": None,
            "peak_frequency": None,
            "spectral_centroid": None,
            "spectral_centroid_over_nyquist": None,
            "upper_half_nyquist_power_fraction": None,
            "upper_quarter_nyquist_power_fraction": None,
            "step_scale_method": "insufficient samples for FFT",
            "high_frequency_cutoff_cycles_per_step": 0.1,
            "energy_fraction_at_or_above_0p1_cycles_per_step": None,
        }

    return {
        "sample_count": int(values.size),
        "mean": float(np.mean(values)),
        "rms": float(np.sqrt(np.mean(values * values))),
        "standard_deviation": standard_deviation,
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "peak_to_peak": float(np.ptp(values)),
        "first_difference_rms": (
            float(np.sqrt(np.mean(differences * differences)))
            if differences.size else 0.0
        ),
        "first_difference_rms_over_standard_deviation": (
            float(np.sqrt(np.mean(differences * differences)) / scale)
            if differences.size else 0.0
        ),
        "second_difference_rms": (
            float(np.sqrt(np.mean(second_differences * second_differences)))
            if second_differences.size else 0.0
        ),
        "total_variation": float(np.sum(np.abs(differences))),
        "turning_point_count": turning_points,
        "turning_point_fraction": (
            turning_points / max(1, values.size - 2)
        ),
        "spectral": spectral,
    }


def write_oscillation_diagnostics(path: Path, curves: Mapping[str, np.ndarray],
                                  metadata: Mapping[str, object]) -> dict:
    time = np.asarray(curves["time"], dtype=float)
    intervals = np.diff(time)
    force_metadata = metadata.get("force", {})
    component_mode = (
        force_metadata.get("visual_component_mode", "translation_x_drag_lift")
        if isinstance(force_metadata, Mapping)
        else "translation_x_drag_lift"
    )
    if component_mode == "global_horizontal_vertical":
        force_series = (("CH", "CD"), ("CV", "CL"))
    else:
        force_series = (("CD", "CD"), ("CL", "CL"))
    diagnostic_series = {
        output_name: series_oscillation_metrics(
            time, np.asarray(curves[curve_name], dtype=float)
        )
        for output_name, curve_name in force_series
    }
    diagnostic_series.update({
        name: series_oscillation_metrics(
            time, np.asarray(curves[name], dtype=float)
        )
        for name in (
            "max_target_error_over_chord", "rms_target_error_over_chord",
            "max_target_error_over_dx", "rms_target_error_over_dx",
        )
    })
    report = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "diagnostic_type": "descriptive_high_frequency_and_oscillation_metrics",
        "interpretation": (
            "No automatic pass/fail thresholds are applied. Compare these metrics "
            "between time-step and stiffness trials and inspect the plotted traces."
        ),
        "case_name": metadata.get("case_name", "unknown"),
        "force_component_mode": component_mode,
        "sampling": {
            "sample_count": int(time.size),
            "time_start": float(time[0]),
            "time_end": float(time[-1]),
            "minimum_interval": float(np.min(intervals)) if intervals.size else None,
            "median_interval": float(np.median(intervals)) if intervals.size else None,
            "maximum_interval": float(np.max(intervals)) if intervals.size else None,
        },
        "series": diagnostic_series,
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def _plot_range(series: Iterable[np.ndarray], include_zero: bool = True) -> tuple[float, float]:
    arrays = [np.asarray(values, dtype=float) for values in series]
    lower = min(float(np.min(values)) for values in arrays)
    upper = max(float(np.max(values)) for values in arrays)
    if include_zero:
        lower = min(lower, 0.0)
        upper = max(upper, 0.0)
    if lower == upper:
        padding = max(abs(lower), 1.0) * 0.05
    else:
        padding = 0.05 * (upper - lower)
    return lower - padding, upper + padding


def _svg_path(x_values: np.ndarray, y_values: np.ndarray,
              x_range: tuple[float, float], y_range: tuple[float, float],
              box: tuple[float, float, float, float]) -> str:
    left, top, width, height = box
    xmin, xmax = x_range
    ymin, ymax = y_range
    # Preserve every sample.  Decimating by evenly spaced indices can alias an
    # alternating-step penalty oscillation into an apparently smooth trace,
    # which is exactly the failure mode this pilot is meant to expose.
    px = left + (x_values - xmin) / (xmax - xmin) * width
    py = top + (ymax - y_values) / (ymax - ymin) * height
    return " ".join(
        ("M" if index == 0 else "L") + f"{x:.3f},{y:.3f}"
        for index, (x, y) in enumerate(zip(px, py))
    )


def _svg_panel(elements: list[str], x: np.ndarray,
               series: Sequence[tuple[str, np.ndarray, str]],
               box: tuple[float, float, float, float], ylabel: str) -> None:
    left, top, width, height = box
    x_range = (float(x[0]), float(x[-1]))
    if x_range[0] == x_range[1]:
        x_range = (x_range[0] - 0.5, x_range[1] + 0.5)
    y_range = _plot_range([values for _, values, _ in series])
    elements.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" '
        'fill="#ffffff" stroke="#222" stroke-width="1"/>'
    )
    for fraction in np.linspace(0.0, 1.0, 6):
        px = left + fraction * width
        py = top + fraction * height
        x_value = x_range[0] + fraction * (x_range[1] - x_range[0])
        y_value = y_range[1] - fraction * (y_range[1] - y_range[0])
        elements.append(
            f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" '
            f'y2="{top + height}" stroke="#e5e5e5"/>'
        )
        elements.append(
            f'<line x1="{left}" y1="{py:.2f}" x2="{left + width}" '
            f'y2="{py:.2f}" stroke="#e5e5e5"/>'
        )
        elements.append(
            f'<text x="{px:.2f}" y="{top + height + 22}" text-anchor="middle" '
            f'class="tick">{x_value:.3g}</text>'
        )
        elements.append(
            f'<text x="{left - 10}" y="{py + 4:.2f}" text-anchor="end" '
            f'class="tick">{y_value:.3g}</text>'
        )
    for name, values, color in series:
        path_data = _svg_path(x, values, x_range, y_range, box)
        elements.append(
            f'<path d="{path_data}" fill="none" stroke="{color}" '
            'stroke-width="2" stroke-linejoin="round"/>'
        )
    elements.append(
        f'<text x="22" y="{top + height / 2}" text-anchor="middle" '
        f'transform="rotate(-90 22 {top + height / 2})" class="axis">'
        f'{html.escape(ylabel)}</text>'
    )
    legend_x = left + 12
    for index, (name, _values, color) in enumerate(series):
        x0 = legend_x + index * 235
        elements.append(
            f'<line x1="{x0}" y1="{top + 20}" x2="{x0 + 30}" '
            f'y2="{top + 20}" stroke="{color}" stroke-width="3"/>'
        )
        elements.append(
            f'<text x="{x0 + 38}" y="{top + 25}" class="legend">'
            f'{html.escape(name)}</text>'
        )


def write_force_svg(path: Path, curves: Mapping[str, np.ndarray],
                    metadata: Mapping[str, object]) -> Path:
    """Write force coefficients and normalized target errors as plain SVG."""

    time = np.asarray(curves["time"], dtype=float)
    width, height = 1100, 760
    elements = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}">',
        '<style>.title{font:bold 22px sans-serif}.axis{font:16px sans-serif}'
        '.tick{font:12px sans-serif;fill:#333}.legend{font:14px sans-serif}</style>',
        '<rect width="100%" height="100%" fill="#fafafa"/>',
        f'<text x="550" y="32" text-anchor="middle" class="title">'
        f'{html.escape(str(metadata.get("case_name", "plate")))}: force and target history</text>',
    ]
    force_metadata = metadata.get("force", {})
    labels = (
        force_metadata.get("visual_component_labels")
        if isinstance(force_metadata, Mapping)
        else None
    )
    if not isinstance(labels, (list, tuple)) or len(labels) != 2:
        labels = ("CD = -Fx / scale", "CL = Fy / scale")
    _svg_panel(
        elements,
        time,
        ((str(labels[0]), curves["CD"], "#1f77b4"),
         (str(labels[1]), curves["CL"], "#d62728")),
        (82.0, 55.0, 980.0, 275.0),
        "force coefficient",
    )
    _svg_panel(
        elements,
        time,
        (("maximum target error / dx",
          curves["max_target_error_over_dx"], "#9467bd"),
         ("RMS target error / dx",
          curves["rms_target_error_over_dx"], "#2ca02c"),
         ("0.5 dx max review",
          np.full(time.shape, 0.5), "#666666"),
         ("0.25 dx RMS review",
          np.full(time.shape, 0.25), "#999999")),
        (82.0, 405.0, 980.0, 275.0),
        "target error / dx",
    )
    elements.extend([
        '<text x="572" y="730" text-anchor="middle" class="axis">time</text>',
        '</svg>',
    ])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")
    return path


def _write_animation_html(path: Path, frame_paths: Sequence[Path],
                          frame_records: Sequence[dict], color_limit: float,
                          force_svg: Path, mp4_path: Path | None,
                          fps: float) -> Path:
    output_dir = Path(path).parent
    relative_frames = [frame.relative_to(output_dir).as_posix() for frame in frame_paths]
    force_relative = force_svg.relative_to(output_dir).as_posix()
    mp4_relative = mp4_path.relative_to(output_dir).as_posix() if mp4_path else None
    video = (
        f'<video controls loop preload="metadata" src="{html.escape(mp4_relative)}"></video>'
        if mp4_relative else ""
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Plate visual pilot</title>
<style>
body {{ margin:0; font:16px system-ui,sans-serif; background:#18202b; color:#edf2f7; }}
main {{ max-width:1100px; margin:auto; padding:20px; }}
.card {{ background:#242f3e; border-radius:10px; padding:16px; margin-bottom:18px; }}
#frame {{ display:block; width:100%; height:auto; image-rendering:auto; background:white; }}
.controls {{ display:flex; gap:12px; align-items:center; flex-wrap:wrap; margin-top:12px; }}
input[type=range] {{ flex:1; min-width:260px; }}
.bar {{ height:14px; border-radius:7px; background:linear-gradient(90deg,rgb(33,102,172),rgb(247,247,247),rgb(178,24,43)); }}
.labels {{ display:flex; justify-content:space-between; font:13px ui-monospace,monospace; }}
.legend {{ display:flex; gap:18px; margin:8px 0; }}
.swatch {{ display:inline-block; width:24px; height:5px; vertical-align:middle; margin-right:6px; }}
object,video {{ width:100%; background:white; }}
button,select {{ font:inherit; padding:6px 10px; }}
</style>
</head>
<body><main>
<h1>Plate visual pilot</h1>
<section class="card">
<img id="frame" alt="Vorticity and plate snapshot">
<div class="legend"><span><i class="swatch" style="background:#fff;border:1px solid #111"></i>actual plate</span><span><i class="swatch" style="background:rgb(250,210,32)"></i>target plate</span></div>
<div class="bar"></div><div class="labels"><span>-{color_limit:.6g}</span><span>vorticity</span><span>+{color_limit:.6g}</span></div>
<div class="controls"><button id="play">Play</button><input id="slider" type="range" min="0" max="{max(0, len(relative_frames)-1)}" value="0"><label>speed <select id="speed"><option value="2">0.5x</option><option value="1" selected>1x</option><option value="0.5">2x</option><option value="0.25">4x</option></select></label><output id="caption"></output></div>
</section>
{video}
<section class="card"><object type="image/svg+xml" data="{html.escape(force_relative)}">Force plot</object></section>
</main>
<script>
const frames={json.dumps(relative_frames)};
const records={json.dumps(list(frame_records))};
const image=document.getElementById('frame'), slider=document.getElementById('slider');
const caption=document.getElementById('caption'), play=document.getElementById('play');
const speed=document.getElementById('speed'); let timer=null;
function show(i) {{ i=Math.max(0,Math.min(frames.length-1,Number(i))); slider.value=i; image.src=frames[i]; const r=records[i]; const stateTime=Number(r.time), targetTime=Number(r.target_time); const targetLabel=Math.abs(stateTime-targetTime)>1e-14 ? `, target t=${{targetTime.toPrecision(6)}}` : ''; caption.textContent=`step ${{r.step}}, t=${{stateTime.toPrecision(6)}}${{targetLabel}}`; }}
function stop() {{ if(timer!==null) clearTimeout(timer); timer=null; play.textContent='Play'; }}
function tick() {{ let i=Number(slider.value)+1; if(i>=frames.length)i=0; show(i); timer=setTimeout(tick,1000/{float(fps):.12g}*Number(speed.value)); }}
play.onclick=()=>{{ if(timer===null){{play.textContent='Pause';tick();}}else stop(); }};
slider.oninput=()=>{{stop();show(slider.value);}}; if(frames.length)show(0);
</script></body></html>
"""
    Path(path).write_text(document, encoding="utf-8")
    return Path(path)


def _try_make_mp4(frame_dir: Path, frame_count: int, output_path: Path,
                  fps: float) -> tuple[Path | None, str]:
    executable = shutil.which("ffmpeg")
    if executable is None:
        return None, "ffmpeg not found; HTML animation remains available"
    command = [
        executable, "-y", "-loglevel", "error", "-framerate", str(float(fps)),
        "-i", str(frame_dir / "frame_%06d.png"), "-frames:v", str(frame_count),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
        str(output_path),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        message = completed.stderr.strip() or f"ffmpeg exited {completed.returncode}"
        return None, message
    return output_path, "created"


def build_visual_artifacts(case_dir: Path, snapshot_dir: Path | None = None,
                           output_dir: Path | None = None, fps: float = 12.0,
                           make_mp4: bool = True,
                           vorticity_limit: float | None = None,
                           max_width: int = 960) -> dict:
    """Build the complete visual-pilot package for one completed case."""

    case_dir = Path(case_dir).resolve()
    snapshot_dir = (
        Path(snapshot_dir).resolve()
        if snapshot_dir is not None else case_dir / "visual_snapshots"
    )
    output_dir = (
        Path(output_dir).resolve()
        if output_dir is not None else case_dir / "visual_pilot"
    )
    if not np.isfinite(fps) or fps <= 0.0:
        raise ValueError("fps must be finite and positive")
    metadata = _read_metadata(case_dir)
    snapshot_paths = sorted(snapshot_dir.glob("snapshot_step*.npz"))
    if not snapshot_paths:
        raise FileNotFoundError(f"no snapshots found in {snapshot_dir}")
    snapshots = [_load_snapshot(path) for path in snapshot_paths]
    snapshots.sort(key=lambda item: (int(item["step"]), float(item["time"])))

    if vorticity_limit is None:
        color_scale_policy = (
            "one symmetric global scale for every frame; maximum of each "
            "frame's 99.5th absolute-vorticity percentile"
        )
        # A single grid-point extremum can flatten every other movie frame.
        # Use one robust, symmetric scale while retaining the same scale for
        # the complete movie.  Fall back to the absolute maximum for a field
        # whose 99.5th percentile is identically zero.
        vorticity_limit = max(
            float(np.quantile(np.abs(snapshot["vorticity"]), 0.995))
            for snapshot in snapshots
        )
        if vorticity_limit == 0.0:
            vorticity_limit = max(
                float(np.max(np.abs(snapshot["vorticity"])))
                for snapshot in snapshots
            )
        if vorticity_limit == 0.0:
            vorticity_limit = 1.0
    else:
        color_scale_policy = "user-supplied symmetric global scale"
    if not np.isfinite(vorticity_limit) or vorticity_limit <= 0.0:
        raise ValueError("vorticity_limit must be finite and positive")

    output_dir.mkdir(parents=True, exist_ok=True)
    frame_dir = output_dir / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []
    frame_records = []
    for index, snapshot in enumerate(snapshots):
        frame_path = frame_dir / f"frame_{index:06d}.png"
        render_snapshot(
            snapshot, frame_path, float(vorticity_limit), max_width=max_width
        )
        frame_paths.append(frame_path)
        frame_records.append({
            "frame": index,
            "step": int(snapshot["step"]),
            "time": float(snapshot["time"]),
            "target_time": float(snapshot["target_time"]),
            "source": snapshot_paths[index].name,
        })

    force_history = _read_force_history(case_dir / "force_history.csv")
    curves = coefficient_history(force_history, metadata)
    force_svg = write_force_svg(output_dir / "force_and_target_history.svg", curves, metadata)
    diagnostics_path = output_dir / "oscillation_diagnostics.json"
    diagnostics = write_oscillation_diagnostics(diagnostics_path, curves, metadata)

    mp4_path = None
    mp4_status = "disabled"
    if make_mp4:
        mp4_path, mp4_status = _try_make_mp4(
            frame_dir, len(frame_paths), output_dir / "plate_vorticity.mp4", fps
        )
    html_path = _write_animation_html(
        output_dir / "animation.html", frame_paths, frame_records,
        float(vorticity_limit), force_svg, mp4_path, fps,
    )
    manifest = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "artifact_type": "translating_plate_visual_pilot",
        "render_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "case_directory": str(case_dir),
        "snapshot_directory": str(snapshot_dir),
        "frame_count": len(frame_paths),
        "vorticity_color_limits": [-float(vorticity_limit), float(vorticity_limit)],
        "vorticity_color_scale_policy": color_scale_policy,
        "html_animation": str(html_path),
        "force_plot": str(force_svg),
        "oscillation_diagnostics": str(diagnostics_path),
        "mp4": str(mp4_path) if mp4_path else None,
        "mp4_status": mp4_status,
        "diagnostic_interpretation": diagnostics["interpretation"],
        "frames": frame_records,
    }
    manifest_path = output_dir / "visual_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest["manifest"] = str(manifest_path)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--snapshot-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--fps", type=float, default=12.0)
    parser.add_argument("--no-mp4", action="store_true")
    parser.add_argument("--vorticity-limit", type=float)
    parser.add_argument("--max-width", type=int, default=960)
    args = parser.parse_args()
    manifest = build_visual_artifacts(
        args.case_dir,
        snapshot_dir=args.snapshot_dir,
        output_dir=args.output_dir,
        fps=args.fps,
        make_mp4=not args.no_mp4,
        vorticity_limit=args.vorticity_limit,
        max_width=args.max_width,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
