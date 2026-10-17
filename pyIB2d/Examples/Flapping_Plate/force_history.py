"""Dependency-free signed force-history writer for the flapping plate."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


HEADER = (
    "step",
    "time",
    "fx_body",
    "fy_body",
    "fx_fluid",
    "fy_fluid",
    "fx_body_lagrangian",
    "fy_body_lagrangian",
    "force_balance_linf",
    "max_target_error",
    "rms_target_error",
)


class ForceHistoryWriter:
    """Write IB reaction forces and target-following diagnostics each step.

    Rows are buffered between periodic flushes.  Flushing every solver step is
    especially expensive on a shared network filesystem and does not make an
    interrupted case usable: a valid case still requires its completion
    marker and the expected number of rows.
    """

    def __init__(self, path: Path, lx: float, ly: float,
                 flush_interval: int = 256):
        self.path = Path(path).resolve()
        self.domain_lengths = np.asarray((float(lx), float(ly)))
        if int(flush_interval) != flush_interval or flush_interval < 1:
            raise ValueError("flush_interval must be a positive integer")
        self.flush_interval = int(flush_interval)
        self._handle = None
        self._writer = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = self.path.open("w", newline="", encoding="utf-8")
        self._writer = csv.writer(self._handle)
        self._writer.writerow(HEADER)
        self._handle.flush()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._handle is not None:
            self._handle.close()
        self._handle = None
        self._writer = None
        return False

    @staticmethod
    def _format(value):
        if isinstance(value, (int, np.integer)):
            return str(int(value))
        return format(float(value), ".17g")

    def record(self, *, step, time, body_force, fluid_force,
               lagrangian_force_density, lagrangian_weight,
               lagrangian_positions, target_info):
        if self._writer is None or self._handle is None:
            raise RuntimeError("ForceHistoryWriter must be used as a context manager")

        body_force = np.asarray(body_force, dtype=float)
        fluid_force = np.asarray(fluid_force, dtype=float)
        lagrangian_force_density = np.asarray(
            lagrangian_force_density, dtype=float
        )
        lagrangian_positions = np.asarray(lagrangian_positions, dtype=float)
        target_info = np.asarray(target_info, dtype=float)

        fluid_force_lagrangian = float(lagrangian_weight) * np.sum(
            lagrangian_force_density, axis=0, dtype=np.float64
        )
        body_force_lagrangian = -fluid_force_lagrangian
        force_balance_linf = np.max(
            np.abs(fluid_force - fluid_force_lagrangian)
        )

        ids = target_info[:, 0].astype(int)
        target_error = lagrangian_positions[ids, :] - target_info[:, 1:3]
        target_error -= self.domain_lengths * np.round(
            target_error / self.domain_lengths
        )
        target_error_norm = np.linalg.norm(target_error, axis=1)
        max_target_error = (
            np.max(target_error_norm) if target_error_norm.size else 0.0
        )
        rms_target_error = (
            np.sqrt(np.mean(target_error_norm * target_error_norm))
            if target_error_norm.size else 0.0
        )

        values = (
            step,
            time,
            body_force[0],
            body_force[1],
            fluid_force[0],
            fluid_force[1],
            body_force_lagrangian[0],
            body_force_lagrangian[1],
            force_balance_linf,
            max_target_error,
            rms_target_error,
        )
        if not np.all(np.isfinite(np.asarray(values[1:], dtype=float))):
            raise FloatingPointError("non-finite value in force history")
        self._writer.writerow([self._format(value) for value in values])
        if int(step) % self.flush_interval == 0:
            self._handle.flush()
