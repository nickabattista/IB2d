"""Authoritative parameters and motion law for the translating-plate baseline."""

from __future__ import annotations

import math


CASE_NAME = "translating_plate"

LX = 1.0
LY = 1.0
RHO = 10.0
MU = 0.187
SUPPORT = 4

NOMINAL_CHORD = 0.1
LEADING_EDGE_X = 0.2
LEADING_EDGE_Y = 0.45
ANGLE_DEGREES = 45.0
TARGET_STIFFNESS = 5.0e7

TRANSLATION_SPEED = 24.0
RAMP_RATE = 200.0
TRANSLATION_STOP_TIME = 0.025

DEFAULT_NX = 32
DEFAULT_NY = 32
DEFAULT_DT = 2.5e-5
DEFAULT_STEPS = 8


def translation_distance(time: float, speed: float = TRANSLATION_SPEED,
                         ramp_rate: float = RAMP_RATE,
                         stop_time: float = TRANSLATION_STOP_TIME) -> float:
    """Return the absolute ramped translation used by the regression case.

    Velocity is ``speed*tanh(ramp_rate*time)`` until ``stop_time`` and zero
    afterward.  Integrating the law analytically avoids accumulated coordinate
    drift and makes MATLAB/Python comparisons deterministic.
    """

    capped_time = min(max(float(time), 0.0), float(stop_time))
    return float(speed) / float(ramp_rate) * math.log(
        math.cosh(float(ramp_rate) * capped_time)
    )
