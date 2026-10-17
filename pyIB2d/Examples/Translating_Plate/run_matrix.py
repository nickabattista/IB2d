#!/usr/bin/env python3
"""Generate and run a sequence of Python translating-plate cases."""

from __future__ import annotations

import argparse
from pathlib import Path

from generate_case import build_case
from main2d import run_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[128, 256])
    parser.add_argument("--dt", type=float, default=2.5e-5)
    parser.add_argument("--final-time", type=float, default=0.025)
    parser.add_argument("--output-root", type=Path, default=Path("runs/matrix"))
    parser.add_argument("--target-stiffness", type=float, default=5.0e7)
    args = parser.parse_args()

    steps = round(args.final_time / args.dt)
    if steps < 1 or abs(steps * args.dt - args.final_time) > 1e-12 * max(1.0, args.final_time):
        raise ValueError("final-time must be a positive integer multiple of dt")

    for resolution in args.resolutions:
        case_dir = args.output_root / f"nx{resolution}"
        build_case(
            case_dir, resolution, resolution, args.dt, steps,
            target_stiffness=args.target_stiffness,
        )
        result = run_case(case_dir)
        print(f"nx={resolution}: {result}")


if __name__ == "__main__":
    main()
