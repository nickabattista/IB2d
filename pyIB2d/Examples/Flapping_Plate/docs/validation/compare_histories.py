#!/usr/bin/env python3
"""Compare matched MATLAB and Python IB2d force histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


REQUIRED_COLUMNS = (
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


def read_history(path: Path) -> dict[str, np.ndarray]:
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [name for name in REQUIRED_COLUMNS if name not in (reader.fieldnames or ())]
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    result = {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in REQUIRED_COLUMNS
    }
    result["step"] = result["step"].astype(np.int64)
    if not all(np.all(np.isfinite(values)) for values in result.values()):
        raise ValueError(f"{path} contains non-finite values")
    return result


def history_diagnostics(history: dict[str, np.ndarray]) -> dict:
    expected_steps = np.arange(history["step"].size)
    body_fluid_residual = np.column_stack((
        history["fx_body"] + history["fx_fluid"],
        history["fy_body"] + history["fy_fluid"],
    ))
    return {
        "rows": int(history["step"].size),
        "steps_contiguous_from_zero": bool(np.array_equal(history["step"], expected_steps)),
        "time_strictly_increasing": bool(np.all(np.diff(history["time"]) > 0.0)),
        "max_body_plus_fluid_residual": float(np.max(np.abs(body_fluid_residual))),
        "max_internal_force_balance": float(np.max(history["force_balance_linf"])),
        "max_target_error": float(np.max(history["max_target_error"])),
    }


def field_comparison(reference: np.ndarray, candidate: np.ndarray,
                     atol: float) -> dict:
    difference = candidate - reference
    reference_scale = max(float(np.max(np.abs(reference))), atol)
    denominator = np.maximum(np.abs(reference), atol)
    return {
        "max_abs": float(np.max(np.abs(difference))),
        "rms_abs": float(np.sqrt(np.mean(difference * difference))),
        "max_relative_pointwise": float(np.max(np.abs(difference) / denominator)),
        "normalized_linf": float(np.max(np.abs(difference)) / reference_scale),
    }


def compare(python_path: Path, matlab_path: Path, rtol: float, atol: float) -> dict:
    python_history = read_history(python_path)
    matlab_history = read_history(matlab_path)

    if not np.array_equal(python_history["step"], matlab_history["step"]):
        raise ValueError("MATLAB and Python step columns differ")
    if not np.allclose(
        python_history["time"], matlab_history["time"], rtol=0.0, atol=1e-14
    ):
        raise ValueError("MATLAB and Python force timestamps differ")

    fields = {}
    passed = True
    for field in ("fx_body", "fy_body", "max_target_error"):
        metrics = field_comparison(matlab_history[field], python_history[field], atol)
        metrics["within_tolerance"] = bool(np.allclose(
            python_history[field], matlab_history[field], rtol=rtol, atol=atol
        ))
        fields[field] = metrics
        passed = passed and metrics["within_tolerance"]

    return {
        "schema_version": 1,
        "python_history": str(Path(python_path).resolve()),
        "matlab_history": str(Path(matlab_path).resolve()),
        "rtol": rtol,
        "atol": atol,
        "rows": int(python_history["step"].size),
        "python_diagnostics": history_diagnostics(python_history),
        "matlab_diagnostics": history_diagnostics(matlab_history),
        "field_comparisons": fields,
        "within_requested_tolerance": bool(passed),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_history", type=Path)
    parser.add_argument("matlab_history", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--enforce", action="store_true",
        help="return a failing exit status when the requested tolerance is exceeded",
    )
    args = parser.parse_args()

    report = compare(args.python_history, args.matlab_history, args.rtol, args.atol)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if args.enforce and not report["within_requested_tolerance"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
