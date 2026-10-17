"""Authoritative defaults for the fixed-Re=64 tilted flapping plate."""

from __future__ import annotations

import math

from kinematics import period_from_peak_speed


CASE_NAME = "flapping_plate_tilted_re64"

LX = 1.6
LY = 1.6
RHO = 10.0
SUPPORT = 4

REYNOLDS_NUMBER = 64.0
NOMINAL_CHORD = 0.1
PEAK_TRANSLATION_SPEED = 24.0
STROKE_AMPLITUDE_OVER_CHORD = 2.5
STROKE_AMPLITUDE = STROKE_AMPLITUDE_OVER_CHORD * NOMINAL_CHORD
PERIOD = period_from_peak_speed(STROKE_AMPLITUDE, PEAK_TRANSLATION_SPEED)
MU = RHO * PEAK_TRANSLATION_SPEED * NOMINAL_CHORD / REYNOLDS_NUMBER

# Exact tilted/dragonfly-like endpoint from the paper-based study plan.
STROKE_ANGLE_DEGREES = 62.8
MEAN_ANGLE_OF_ATTACK_DEGREES = 60.0
PITCH_AMPLITUDE_DEGREES = 45.0
PITCH_PHASE_RADIANS = 0.0
PIVOT_FRACTION = 0.5
STROKE_MIDPOINT = (0.5 * LX, 0.5 * LY)

TARGET_STIFFNESS = 1.0e9

# These defaults deliberately make only a tiny smoke case.  Production and
# timestep-calibration jobs must supply their resolution, dt, and step count.
DEFAULT_NX = 32
DEFAULT_NY = 32
DEFAULT_DT = 1.25e-6
DEFAULT_STEPS = 8


def reynolds_number(
    rho: float,
    peak_speed: float,
    chord: float,
    dynamic_viscosity: float,
) -> float:
    """Return ``rho*Umax*c/mu`` for the study's Reynolds convention."""

    values = (rho, peak_speed, chord, dynamic_viscosity)
    if not all(math.isfinite(float(value)) and value > 0.0 for value in values):
        raise ValueError("Reynolds-number inputs must be finite and positive")
    return float(rho) * float(peak_speed) * float(chord) / float(dynamic_viscosity)


def viscosity_for_reynolds(
    rho: float,
    peak_speed: float,
    chord: float,
    reynolds: float,
) -> float:
    """Return dynamic viscosity for a desired peak-speed Reynolds number."""

    values = (rho, peak_speed, chord, reynolds)
    if not all(math.isfinite(float(value)) and value > 0.0 for value in values):
        raise ValueError("viscosity-scaling inputs must be finite and positive")
    return float(rho) * float(peak_speed) * float(chord) / float(reynolds)


def marker_spacing(lx: float, nx: int) -> float:
    """Return the required Lagrangian spacing ``ds=dx/2``."""

    if not math.isfinite(float(lx)) or lx <= 0.0:
        raise ValueError("lx must be finite and positive")
    if int(nx) != nx or nx <= 0:
        raise ValueError("nx must be a positive integer")
    return 0.5 * float(lx) / int(nx)


def default_motion_metadata() -> dict:
    """Return JSON-compatible motion parameters for the tilted case."""

    return {
        "kind": "centered_sinusoidal_stroke_and_pitch",
        "stroke_midpoint": list(STROKE_MIDPOINT),
        "stroke_amplitude": STROKE_AMPLITUDE,
        "stroke_amplitude_over_chord": STROKE_AMPLITUDE_OVER_CHORD,
        "stroke_angle_degrees": STROKE_ANGLE_DEGREES,
        "mean_angle_of_attack_degrees": MEAN_ANGLE_OF_ATTACK_DEGREES,
        "pitch_amplitude_degrees": PITCH_AMPLITUDE_DEGREES,
        "pitch_phase_radians": PITCH_PHASE_RADIANS,
        "period": PERIOD,
        "pivot_fraction": PIVOT_FRACTION,
        "angle_convention": "theta=beta-alpha, leading_edge_to_trailing_edge",
        "evaluation_time": "step_midpoint",
    }
