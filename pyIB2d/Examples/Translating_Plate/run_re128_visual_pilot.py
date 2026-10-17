#!/usr/bin/env python3
"""Run and render the first Re=128 visual/numerical checkpoint.

This is intentionally a single 1024-square translating-plate case, not the
production convergence matrix.  It exercises a conservative 64-cells/chord
stiffness/timestep reference, records 101 compact vorticity/plate snapshots,
and builds force, target-error, and oscillation diagnostics for human review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import socket
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from re128_study import (
    StudyCase,
    case_record,
    generate_study_case,
    validate_history,
)
from visual_pilot import VisualSnapshotWriter, build_visual_artifacts


PILOT_CASE = StudyCase(
    index=0,
    case_id="visual_L16_P64_N1024",
    phase="visual_pilot",
    study_axes=("visual", "stiffness_timestep_check"),
    box_chords=16.0,
    nx=1024,
    ny=1024,
    dt=1.25e-6,
    target_stiffness=1.0e9,
)
SNAPSHOT_EVERY_STEPS = 200


def numerical_source_sha256() -> str:
    """Fingerprint source that can change the numerical trajectory.

    A completed directory must not be silently reused after a code change just
    because its numerical parameters are unchanged.  Rendering-only code is
    intentionally excluded so plots can be rebuilt without rerunning the
    solver.
    """

    python_root = Path(__file__).resolve().parents[2]
    case_root = Path(__file__).resolve().parent
    case_sources = (
        "case_config.py",
        "force_history.py",
        "generate_case.py",
        "main2d.py",
        "re128_study.py",
        "update_Target_Point_Positions.py",
    )
    files = sorted(
        [
            path for path in (python_root / "IBM_Blackbox").glob("*.py")
            if path.name != "test.py"
        ]
        + [case_root / name for name in case_sources],
        key=lambda path: path.relative_to(python_root).as_posix(),
    )
    digest = hashlib.sha256()
    for path in files:
        relative = path.relative_to(python_root).as_posix().encode("utf-8")
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def expected_snapshot_steps(case: StudyCase, every_steps: int) -> list[int]:
    if every_steps < 1:
        raise ValueError("snapshot interval must be positive")
    steps = [0]
    steps.extend(range(every_steps, case.steps + 1, every_steps))
    if steps[-1] != case.steps:
        steps.append(case.steps)
    return steps


def pilot_record(case: StudyCase = PILOT_CASE,
                 every_steps: int = SNAPSHOT_EVERY_STEPS) -> dict:
    record = {
        "schema_version": 1,
        "purpose": "pre-production visual and force-history checkpoint",
        "case": case_record(case),
        "numerical_source_sha256": numerical_source_sha256(),
        "snapshot_every_steps": int(every_steps),
        "expected_snapshot_steps": expected_snapshot_steps(case, every_steps),
        "review_required": True,
        "review_items": [
            "actual plate follows target without visible deformation or phase lag",
            "vorticity evolution is physically plausible and remains in the box",
            "CD and CL histories do not show step-scale alternating ringing",
            "maximum and RMS target errors remain acceptable",
        ],
    }
    encoded = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    record["visual_configuration_sha256"] = hashlib.sha256(encoded).hexdigest()
    return record


def _snapshot_steps(case_dir: Path) -> list[int]:
    steps = []
    for path in sorted((case_dir / "visual_snapshots").glob("snapshot_step*.npz")):
        with np.load(path, allow_pickle=False) as stored:
            steps.append(int(stored["step"]))
    return steps


def validate_completed_pilot(case_dir: Path, case: StudyCase,
                             every_steps: int) -> dict:
    case_dir = Path(case_dir)
    required = (
        case_dir / "COMPLETE.json",
        case_dir / "case_metadata.json",
        case_dir / "force_history.csv",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError("visual pilot is missing: " + ", ".join(missing))

    complete = json.loads(required[0].read_text(encoding="utf-8"))
    requested = pilot_record(case, every_steps)
    if complete.get("visual_configuration_sha256") != requested[
        "visual_configuration_sha256"
    ]:
        raise ValueError("visual pilot configuration does not match this request")
    history = validate_history(required[2], case)
    actual_steps = _snapshot_steps(case_dir)
    expected_steps = expected_snapshot_steps(case, every_steps)
    if actual_steps != expected_steps:
        raise ValueError(
            f"snapshot sequence is incomplete: expected {expected_steps}, "
            f"found {actual_steps}"
        )
    return {
        "history": history,
        "snapshot_count": len(actual_steps),
        "first_snapshot_step": actual_steps[0],
        "last_snapshot_step": actual_steps[-1],
    }


def run_visual_pilot(output_root: Path, case: StudyCase = PILOT_CASE,
                     every_steps: int = SNAPSHOT_EVERY_STEPS) -> Path:
    output_root = Path(output_root).resolve()
    final_dir = output_root / case.case_id
    complete_path = final_dir / "COMPLETE.json"
    if complete_path.is_file():
        validate_completed_pilot(final_dir, case, every_steps)
        print(f"Already complete: {final_dir}")
        return final_dir
    if final_dir.exists():
        raise FileExistsError(
            f"incomplete visual-pilot directory exists: {final_dir}; move it "
            "aside before retrying"
        )

    output_root.mkdir(parents=True, exist_ok=True)
    attempt_id = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    attempt_dir = output_root / f"{case.case_id}.partial-{attempt_id}"
    if attempt_dir.exists():
        raise FileExistsError(f"attempt directory already exists: {attempt_dir}")
    attempt_dir.mkdir(parents=True, exist_ok=False)

    started = datetime.now(timezone.utc)
    wall_start = time.perf_counter()
    metadata = generate_study_case(case, attempt_dir)
    visual_configuration = pilot_record(case, every_steps)
    metadata["visual_pilot"] = visual_configuration
    (attempt_dir / "case_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    snapshots = VisualSnapshotWriter.from_case(
        attempt_dir, every_steps=every_steps
    )
    # Importing the legacy driver prints capability messages, so defer it
    # until the numerical run rather than polluting ``show`` JSON output.
    from main2d import run_case
    history_path = run_case(
        attempt_dir,
        state_callback=snapshots.record,
        state_callback_interval=every_steps,
    )
    history_diagnostics = validate_history(history_path, case)
    snapshot_steps = _snapshot_steps(attempt_dir)
    expected_steps = expected_snapshot_steps(case, every_steps)
    if snapshot_steps != expected_steps:
        raise ValueError(
            f"snapshot sequence is incomplete: expected {expected_steps}, "
            f"found {snapshot_steps}"
        )

    visual_configuration.update({
        "status": "complete_pending_human_review",
        "configuration_sha256": case_record(case)["configuration_sha256"],
        "runtime": {
            "started_utc": started.isoformat(),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - wall_start,
            "hostname": socket.gethostname(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        },
        "history_diagnostics": history_diagnostics,
        "snapshot_count": len(snapshot_steps),
    })
    metadata["visual_pilot"] = visual_configuration
    (attempt_dir / "case_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (attempt_dir / "COMPLETE.json").write_text(
        json.dumps(visual_configuration, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    attempt_dir.replace(final_dir)
    print(f"Completed numerical pilot: {final_dir}")
    return final_dir


def render_visual_pilot(output_root: Path, *, make_mp4: bool = True,
                        fps: float = 12.0) -> dict:
    case_dir = Path(output_root).resolve() / PILOT_CASE.case_id
    validate_completed_pilot(
        case_dir, PILOT_CASE, SNAPSHOT_EVERY_STEPS
    )
    manifest = build_visual_artifacts(
        case_dir, fps=fps, make_mp4=make_mp4
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("show", "run", "render", "all"))
    parser.add_argument(
        "--output-root", type=Path, default=Path("runs/re128_visual_pilot")
    )
    parser.add_argument("--no-mp4", action="store_true")
    parser.add_argument("--fps", type=float, default=12.0)
    args = parser.parse_args()

    try:
        if args.action == "show":
            print(json.dumps(pilot_record(), indent=2, sort_keys=True))
            return
        if args.action in ("run", "all"):
            run_visual_pilot(args.output_root)
        if args.action in ("render", "all"):
            render_visual_pilot(
                args.output_root, make_mp4=not args.no_mp4, fps=args.fps
            )
    except (FileExistsError, FileNotFoundError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
