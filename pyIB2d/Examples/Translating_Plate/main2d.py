#!/usr/bin/env python3
"""Generate and run the headless translating-plate IB2d regression case."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PYTHON_ROOT = HERE.parents[1]
CORE = PYTHON_ROOT / "IBM_Blackbox"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CORE))

import IBM_Driver as Driver
from force_history import ForceHistoryWriter
from generate_case import build_case
from please_Initialize_Simulation import please_Initialize_Simulation
from update_Target_Point_Positions import reset_reference


def run_case(case_dir: Path, *, state_callback=None,
             state_callback_interval=None) -> Path:
    """Run an already generated case and return its force-history path."""

    case_dir = Path(case_dir).resolve()
    required = ("input2d", "plate.vertex", "plate.target", "case_metadata.json")
    missing = [name for name in required if not (case_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"case is incomplete; missing: {', '.join(missing)}")

    with (case_dir / "case_metadata.json").open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    previous_directory = Path.cwd()
    history_path = case_dir / "force_history.csv"
    try:
        os.chdir(case_dir)
        reset_reference()
        parameters = please_Initialize_Simulation()
        with ForceHistoryWriter(
            history_path,
            lx=metadata["domain"]["lx"],
            ly=metadata["domain"]["ly"],
        ) as force_writer:
            Driver.main(
                *parameters,
                force_callback=force_writer.record,
                state_callback=state_callback,
                state_callback_interval=state_callback_interval,
            )
    finally:
        os.chdir(previous_directory)
    return history_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=HERE / "runs" / "smoke")
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--nx", type=int, default=32)
    parser.add_argument("--ny", type=int, default=32)
    parser.add_argument("--dt", type=float, default=2.5e-5)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--target-stiffness", type=float, default=5.0e7)
    args = parser.parse_args()

    input_path = args.case_dir / "input2d"
    if args.regenerate or not input_path.is_file():
        build_case(
            args.case_dir, nx=args.nx, ny=args.ny, dt=args.dt, steps=args.steps,
            target_stiffness=args.target_stiffness,
        )

    result = run_case(args.case_dir)
    print(f"Force history: {result}")


if __name__ == "__main__":
    main()
