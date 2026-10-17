"""Absolute midpoint target update for the prescribed flapping plate."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from kinematics import (
    FlappingMotion,
    body_coordinates_from_world,
    world_coordinates_from_body,
)


_body_coordinates = None
_reference_ids = None
_motion = None


def reset_reference() -> None:
    """Clear cached plate-frame coordinates between simulations or tests."""

    global _body_coordinates, _reference_ids, _motion
    _body_coordinates = None
    _reference_ids = None
    _motion = None


def _motion_parameters() -> FlappingMotion:
    with Path("case_metadata.json").open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    return FlappingMotion.from_mapping(metadata["motion"])


def update_Target_Point_Positions(dt, current_time, target_info):
    """Return target coordinates at the force-evaluation midpoint.

    ``plate.vertex`` is expected to contain the prescribed pose at ``t=0``.
    On the first call those coordinates are converted to the rigid plate frame.
    All later positions are reconstructed from that immutable reference, so a
    repeated call at the same time gives bitwise-equivalent target coordinates
    instead of advancing the plate twice.
    """

    global _body_coordinates, _reference_ids, _motion

    target_info = np.asarray(target_info)
    if target_info.ndim != 2 or target_info.shape[1] < 3:
        raise ValueError("target_info must have at least id, x, and y columns")
    if not np.issubdtype(target_info.dtype, np.floating):
        raise TypeError("target_info must use a floating-point dtype")

    dt = float(dt)
    current_time = float(current_time)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")
    if not np.isfinite(current_time):
        raise ValueError("current_time must be finite")

    ids = target_info[:, 0].astype(int)
    if (
        _body_coordinates is None
        or _reference_ids is None
        or not np.array_equal(ids, _reference_ids)
    ):
        _motion = _motion_parameters()
        _body_coordinates = body_coordinates_from_world(
            target_info[:, 1:3], 0.0, _motion
        )
        _reference_ids = ids.copy()

    evaluation_time = current_time + 0.5 * dt
    target_info[:, 1:3] = world_coordinates_from_body(
        _body_coordinates, evaluation_time, _motion
    )
    return target_info
