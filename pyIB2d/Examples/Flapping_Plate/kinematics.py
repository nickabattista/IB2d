"""Analytical rigid-body kinematics for the prescribed flapping plate.

The plate center follows a sinusoid along a stroke plane and the plate
rotates about its midchord.  ``theta`` is the directed leading-edge to
trailing-edge angle in laboratory coordinates.  The model convention is

``alpha = alpha_mean + alpha_amplitude*sin(2*pi*t/T + pitch_phase)``
``theta = beta - alpha``.

All public functions are absolute in time.  The IB target update therefore
does not accumulate roundoff or integration drift over many cycles.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class FlappingMotion:
    """Parameters for the centered stroke and midchord pitch motion."""

    stroke_midpoint: tuple[float, float]
    stroke_amplitude: float
    stroke_angle_degrees: float
    mean_angle_of_attack_degrees: float
    pitch_amplitude_degrees: float
    period: float
    pitch_phase_radians: float = 0.0
    pivot_fraction: float = 0.5

    def __post_init__(self) -> None:
        values = (
            *self.stroke_midpoint,
            self.stroke_amplitude,
            self.stroke_angle_degrees,
            self.mean_angle_of_attack_degrees,
            self.pitch_amplitude_degrees,
            self.period,
            self.pitch_phase_radians,
            self.pivot_fraction,
        )
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("all flapping-motion parameters must be finite")
        if self.stroke_amplitude <= 0.0 or self.period <= 0.0:
            raise ValueError("stroke_amplitude and period must be positive")
        if not math.isclose(self.pivot_fraction, 0.5, abs_tol=1.0e-14):
            raise ValueError("this prescribed rigid plate must pivot at midchord")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> "FlappingMotion":
        """Create a motion from the JSON-compatible metadata schema."""

        midpoint = values["stroke_midpoint"]
        if not isinstance(midpoint, Sequence) or len(midpoint) != 2:
            raise ValueError("motion.stroke_midpoint must contain x and y")
        return cls(
            stroke_midpoint=(float(midpoint[0]), float(midpoint[1])),
            stroke_amplitude=float(values["stroke_amplitude"]),
            stroke_angle_degrees=float(values["stroke_angle_degrees"]),
            mean_angle_of_attack_degrees=float(
                values["mean_angle_of_attack_degrees"]
            ),
            pitch_amplitude_degrees=float(values["pitch_amplitude_degrees"]),
            period=float(values["period"]),
            pitch_phase_radians=float(values.get("pitch_phase_radians", 0.0)),
            pivot_fraction=float(values.get("pivot_fraction", 0.5)),
        )


@dataclass(frozen=True)
class PlatePose:
    """Midchord position and directed chord angle at one instant."""

    center_x: float
    center_y: float
    alpha_radians: float
    theta_radians: float

    @property
    def center(self) -> np.ndarray:
        return np.asarray((self.center_x, self.center_y), dtype=float)


@dataclass(frozen=True)
class PlateVelocity:
    """Translational and angular velocity of the prescribed plate."""

    center_velocity_x: float
    center_velocity_y: float
    alpha_rate: float
    theta_rate: float

    @property
    def center_velocity(self) -> np.ndarray:
        return np.asarray(
            (self.center_velocity_x, self.center_velocity_y), dtype=float
        )


def angular_frequency(period: float) -> float:
    """Return ``2*pi/period`` after validating the period."""

    period = float(period)
    if not math.isfinite(period) or period <= 0.0:
        raise ValueError("period must be finite and positive")
    return 2.0 * math.pi / period


def peak_translation_speed(stroke_amplitude: float, period: float) -> float:
    """Return the maximum center speed, ``pi*A0/T``."""

    stroke_amplitude = float(stroke_amplitude)
    period = float(period)
    if (
        not math.isfinite(stroke_amplitude)
        or not math.isfinite(period)
        or stroke_amplitude <= 0.0
        or period <= 0.0
    ):
        raise ValueError("stroke_amplitude and period must be finite and positive")
    return math.pi * stroke_amplitude / period


def period_from_peak_speed(stroke_amplitude: float, peak_speed: float) -> float:
    """Return the period that realizes a prescribed peak center speed."""

    stroke_amplitude = float(stroke_amplitude)
    peak_speed = float(peak_speed)
    if (
        not math.isfinite(stroke_amplitude)
        or not math.isfinite(peak_speed)
        or stroke_amplitude <= 0.0
        or peak_speed <= 0.0
    ):
        raise ValueError("stroke_amplitude and peak_speed must be positive")
    return math.pi * stroke_amplitude / peak_speed


def pose_at_time(time: float, motion: FlappingMotion) -> PlatePose:
    """Evaluate the paper-based plate pose at absolute time ``time``."""

    time = float(time)
    if not math.isfinite(time):
        raise ValueError("time must be finite")

    omega = angular_frequency(motion.period)
    phase = omega * time
    beta = math.radians(motion.stroke_angle_degrees)
    alpha_mean = math.radians(motion.mean_angle_of_attack_degrees)
    alpha_amplitude = math.radians(motion.pitch_amplitude_degrees)

    stroke_coordinate = 0.5 * motion.stroke_amplitude * math.cos(phase)
    center_x = motion.stroke_midpoint[0] + stroke_coordinate * math.cos(beta)
    center_y = motion.stroke_midpoint[1] + stroke_coordinate * math.sin(beta)
    alpha = alpha_mean + alpha_amplitude * math.sin(
        phase + motion.pitch_phase_radians
    )
    theta = beta - alpha
    return PlatePose(center_x, center_y, alpha, theta)


def velocity_at_time(time: float, motion: FlappingMotion) -> PlateVelocity:
    """Evaluate center velocity and directed-chord angular velocity."""

    time = float(time)
    if not math.isfinite(time):
        raise ValueError("time must be finite")

    omega = angular_frequency(motion.period)
    phase = omega * time
    beta = math.radians(motion.stroke_angle_degrees)
    alpha_amplitude = math.radians(motion.pitch_amplitude_degrees)

    stroke_speed = (
        -0.5 * motion.stroke_amplitude * omega * math.sin(phase)
    )
    alpha_rate = alpha_amplitude * omega * math.cos(
        phase + motion.pitch_phase_radians
    )
    return PlateVelocity(
        stroke_speed * math.cos(beta),
        stroke_speed * math.sin(beta),
        alpha_rate,
        -alpha_rate,
    )


def _rotation(angle: float) -> np.ndarray:
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return np.asarray(((cosine, -sine), (sine, cosine)), dtype=float)


def body_coordinates_from_world(
    world_coordinates: np.ndarray,
    time: float,
    motion: FlappingMotion,
) -> np.ndarray:
    """Express world coordinates in the plate frame at ``time``."""

    coordinates = np.asarray(world_coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("world_coordinates must have shape (n, 2)")
    pose = pose_at_time(time, motion)
    return (coordinates - pose.center) @ _rotation(pose.theta_radians)


def world_coordinates_from_body(
    body_coordinates: np.ndarray,
    time: float,
    motion: FlappingMotion,
) -> np.ndarray:
    """Transform plate-frame coordinates to world coordinates at ``time``."""

    coordinates = np.asarray(body_coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("body_coordinates must have shape (n, 2)")
    pose = pose_at_time(time, motion)
    return pose.center + coordinates @ _rotation(pose.theta_radians).T


def point_velocities(
    body_coordinates: np.ndarray,
    time: float,
    motion: FlappingMotion,
) -> np.ndarray:
    """Return the laboratory velocity of each plate-frame point."""

    coordinates = np.asarray(body_coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("body_coordinates must have shape (n, 2)")
    pose = pose_at_time(time, motion)
    velocity = velocity_at_time(time, motion)
    world_offsets = coordinates @ _rotation(pose.theta_radians).T
    rotational = velocity.theta_rate * np.column_stack(
        (-world_offsets[:, 1], world_offsets[:, 0])
    )
    return velocity.center_velocity + rotational
