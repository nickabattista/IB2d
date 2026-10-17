"""Flapping-aware wrapper around the dependency-free visual pipeline."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

from kinematics import (
    FlappingMotion,
    body_coordinates_from_world,
    world_coordinates_from_body,
)


HERE = Path(__file__).resolve().parent
TRANSLATING_CASE = HERE.parent / "Translating_Plate"
if str(TRANSLATING_CASE) not in sys.path:
    sys.path.append(str(TRANSLATING_CASE))

from visual_pilot import (  # noqa: E402
    VisualSnapshotWriter as _BaseVisualSnapshotWriter,
    build_visual_artifacts as _build_visual_artifacts,
    coefficient_history,
    write_oscillation_diagnostics,
)


class FlappingVisualSnapshotWriter(_BaseVisualSnapshotWriter):
    """Align midpoint targets with full-step snapshots analytically."""

    def _target_at_state_time(
        self, target_positions: np.ndarray, target_time: float, state_time: float
    ) -> tuple[np.ndarray, float]:
        positions = np.asarray(target_positions, dtype=float).copy()
        if abs(target_time - state_time) <= 1.0e-15:
            return positions, target_time
        motion_values = self.metadata.get("motion", {})
        if motion_values.get("kind") != "centered_sinusoidal_stroke_and_pitch":
            return super()._target_at_state_time(
                positions, target_time, state_time
            )
        motion = FlappingMotion.from_mapping(motion_values)
        body = body_coordinates_from_world(positions, target_time, motion)
        return world_coordinates_from_body(body, state_time, motion), state_time


def build_visual_artifacts(*args, **kwargs) -> dict:
    """Build shared artifacts while identifying them as a flapping pilot."""

    manifest = _build_visual_artifacts(*args, **kwargs)
    manifest["artifact_type"] = "flapping_plate_visual_pilot"
    manifest_path = Path(manifest["manifest"])
    base = manifest_path.parent
    for name in (
        "case_directory", "snapshot_directory", "html_animation", "force_plot",
        "oscillation_diagnostics", "mp4", "manifest",
    ):
        value = manifest.get(name)
        if value is not None and Path(value).is_absolute():
            manifest[name] = os.path.relpath(Path(value), base).replace("\\", "/")
    stored = dict(manifest)
    stored.pop("manifest", None)
    manifest_path.write_text(
        json.dumps(stored, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest

