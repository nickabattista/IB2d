"""Absolute target motion for the translating-plate regression."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from case_config import translation_distance


_reference_coordinates = None
_reference_ids = None


def reset_reference() -> None:
    """Clear cached coordinates, primarily for repeated runs and tests."""

    global _reference_coordinates, _reference_ids
    _reference_coordinates = None
    _reference_ids = None


def _motion_parameters() -> dict:
    with Path("case_metadata.json").open("r", encoding="utf-8") as handle:
        return json.load(handle)["motion"]


def update_Target_Point_Positions(dt, current_time, target_info):
    """Return target coordinates evaluated at the force midpoint."""

    global _reference_coordinates, _reference_ids

    ids = target_info[:, 0].astype(int)
    if (
        _reference_coordinates is None
        or _reference_ids is None
        or not np.array_equal(ids, _reference_ids)
    ):
        _reference_coordinates = target_info[:, 1:3].copy()
        _reference_ids = ids.copy()

    motion = _motion_parameters()
    evaluation_time = float(current_time) + 0.5 * float(dt)
    displacement = translation_distance(
        evaluation_time,
        speed=motion["speed"],
        ramp_rate=motion["ramp_rate"],
        stop_time=motion["stop_time"],
    )

    target_info[:, 1] = _reference_coordinates[:, 0] + displacement
    target_info[:, 2] = _reference_coordinates[:, 1]
    return target_info
