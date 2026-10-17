"""Generate a self-contained prescribed flapping-plate IB2d case."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from case_config import (
    DEFAULT_DT,
    DEFAULT_NX,
    DEFAULT_NY,
    DEFAULT_STEPS,
    LX,
    LY,
    MEAN_ANGLE_OF_ATTACK_DEGREES,
    NOMINAL_CHORD,
    PEAK_TRANSLATION_SPEED,
    PITCH_AMPLITUDE_DEGREES,
    PITCH_PHASE_RADIANS,
    REYNOLDS_NUMBER,
    RHO,
    STROKE_AMPLITUDE_OVER_CHORD,
    STROKE_ANGLE_DEGREES,
    SUPPORT,
    TARGET_STIFFNESS,
)
from kinematics import (
    FlappingMotion,
    peak_translation_speed,
    period_from_peak_speed,
    pose_at_time,
    world_coordinates_from_body,
)


def _input_text(nx: int, ny: int, dt: float, steps: int, *,
                lx: float, ly: float, rho: float, mu: float) -> str:
    tfinal = steps * dt
    return f"""% Generated Re=64 prescribed flapping-plate case.
Fluid_Parameters {{
mu = {mu:.17g}
rho = {rho:.17g}
}}

Temporal_Information {{
Tfinal = {tfinal:.17g}
dt = {dt:.17g}
}}

Grid_Parameters {{
Nx = {nx}
Ny = {ny}
Lx = {lx:.17g}
Ly = {ly:.17g}
supp = {SUPPORT}
}}

Lag_Structure_Info {{
springs = 0
update_springs = 0
target_pts = 1
update_target = 1
beams = 0
update_beams = 0
nonInvariant_beams = 0
update_nonInv_beams = 0
FV_LT_muscle = 0
3_element_muscle = 0
arb_ext_force = 0
tracers = 0
mass_pts = 0
gravity = 0
x_gravity_vec_comp = 0
y_gravity_vec_comp = 0
porous_media = 0
concentration = 0
electro_phys = 0
damped_springs = 0
update_damp_springs = 0
boussinesq = 0
expansion_coeff = 0
user_force_model = 0
poroelastic = 0
brinkman = 0
}}

Output_Info {{
print_dump = {steps + 1}
plot_Matlab = 0
plot_LagPts = 0
plot_Velocity = 0
plot_Vorticity = 0
plot_MagVelocity = 0
plot_Pressure = 0
save_Vorticity = 0
save_Pressure = 0
save_uVec = 0
save_uMag = 0
save_uX = 0
save_uY = 0
save_fMag = 0
save_fX = 0
save_fY = 0
save_hier = 0
}}

Lag_Name {{
string_name = "plate"
}}
"""


def build_case(
    output_dir: Path,
    nx: int = DEFAULT_NX,
    ny: int = DEFAULT_NY,
    dt: float = DEFAULT_DT,
    steps: int = DEFAULT_STEPS,
    target_stiffness: float = TARGET_STIFFNESS,
    *,
    lx: float = LX,
    ly: float = LY,
    rho: float = RHO,
    reynolds_number: float = REYNOLDS_NUMBER,
    nominal_chord: float = NOMINAL_CHORD,
    peak_speed: float = PEAK_TRANSLATION_SPEED,
    period: float | None = None,
    stroke_amplitude_over_chord: float = STROKE_AMPLITUDE_OVER_CHORD,
    stroke_angle_degrees: float = STROKE_ANGLE_DEGREES,
    mean_angle_of_attack_degrees: float = MEAN_ANGLE_OF_ATTACK_DEGREES,
    pitch_amplitude_degrees: float = PITCH_AMPLITUDE_DEGREES,
    pitch_phase_radians: float = PITCH_PHASE_RADIANS,
    stroke_midpoint: tuple[float, float] | None = None,
) -> dict:
    """Write solver inputs and return complete, auditable case metadata."""

    output_dir = Path(output_dir).resolve()
    if nx < 8 or ny < 8 or nx % 2 or ny % 2:
        raise ValueError("nx and ny must be even integers of at least 8")
    values = (
        dt, target_stiffness, lx, ly, rho, reynolds_number, nominal_chord,
        peak_speed, stroke_amplitude_over_chord,
    )
    if steps < 1 or not all(math.isfinite(float(v)) and v > 0.0 for v in values):
        raise ValueError("steps and all physical/numerical scales must be positive")

    dx, dy = lx / nx, ly / ny
    if not math.isclose(dx, dy, rel_tol=1.0e-12, abs_tol=1.0e-15):
        raise ValueError("this solver path requires square Eulerian cells")
    ds = 0.5 * dx
    marker_count = max(2, int(round(nominal_chord / ds)))
    effective_chord = marker_count * ds

    amplitude = stroke_amplitude_over_chord * nominal_chord
    if period is None:
        period = period_from_peak_speed(amplitude, peak_speed)
    else:
        period = float(period)
        peak_speed = peak_translation_speed(amplitude, period)
    mu = rho * peak_speed * nominal_chord / reynolds_number
    midpoint = (
        (0.5 * lx, 0.5 * ly)
        if stroke_midpoint is None
        else (float(stroke_midpoint[0]), float(stroke_midpoint[1]))
    )
    motion = FlappingMotion(
        stroke_midpoint=midpoint,
        stroke_amplitude=amplitude,
        stroke_angle_degrees=stroke_angle_degrees,
        mean_angle_of_attack_degrees=mean_angle_of_attack_degrees,
        pitch_amplitude_degrees=pitch_amplitude_degrees,
        period=period,
        pitch_phase_radians=pitch_phase_radians,
        pivot_fraction=0.5,
    )

    beta = math.radians(stroke_angle_degrees)
    half_x = 0.5 * amplitude * abs(math.cos(beta)) + 0.5 * effective_chord
    half_y = 0.5 * amplitude * abs(math.sin(beta)) + 0.5 * effective_chord
    clearances = (
        midpoint[0] - half_x, lx - midpoint[0] - half_x,
        midpoint[1] - half_y, ly - midpoint[1] - half_y,
    )
    if min(clearances) <= 0.0:
        raise ValueError("the conservative swept plate envelope does not fit")

    output_dir.mkdir(parents=True, exist_ok=True)
    body_coordinates = np.column_stack((
        (np.arange(marker_count, dtype=float) + 0.5) * ds
        - 0.5 * effective_chord,
        np.zeros(marker_count, dtype=float),
    ))
    coordinates = world_coordinates_from_body(body_coordinates, 0.0, motion)
    vertex_lines = [str(marker_count)]
    vertex_lines.extend(f"{x:.17e} {y:.17e}" for x, y in coordinates)
    (output_dir / "plate.vertex").write_text(
        "\n".join(vertex_lines) + "\n", encoding="utf-8"
    )
    target_lines = [str(marker_count)]
    target_lines.extend(
        f"{marker_id:d} {target_stiffness:.17e}"
        for marker_id in range(marker_count)
    )
    (output_dir / "plate.target").write_text(
        "\n".join(target_lines) + "\n", encoding="utf-8"
    )
    (output_dir / "input2d").write_text(
        _input_text(nx, ny, dt, steps, lx=lx, ly=ly, rho=rho, mu=mu),
        encoding="utf-8",
    )

    urms = peak_speed / math.sqrt(2.0)
    q_urms = 0.5 * rho * urms * urms * nominal_chord
    q_umax = 0.5 * rho * peak_speed * peak_speed * nominal_chord
    initial_pose = pose_at_time(0.0, motion)
    steps_per_cycle = period / dt
    motion_metadata = {
        "kind": "centered_sinusoidal_stroke_and_pitch",
        "stroke_midpoint": list(midpoint),
        "stroke_amplitude": amplitude,
        "stroke_amplitude_over_chord": stroke_amplitude_over_chord,
        "stroke_angle_degrees": stroke_angle_degrees,
        "mean_angle_of_attack_degrees": mean_angle_of_attack_degrees,
        "pitch_amplitude_degrees": pitch_amplitude_degrees,
        "pitch_phase_radians": pitch_phase_radians,
        "period": period,
        "frequency": 1.0 / period,
        "peak_translation_speed": peak_speed,
        "rms_translation_speed": urms,
        "speed": peak_speed,
        "pivot_fraction": 0.5,
        "angle_convention": "theta=beta-alpha, leading_edge_to_trailing_edge",
        "evaluation_time": "step_midpoint",
    }
    metadata = {
        "schema_version": 2,
        "case_name": "flapping_plate_re64",
        "fluid": {
            "rho": rho,
            "mu": mu,
            "reynolds_number": reynolds_number,
            "effective_reynolds_number": (
                rho * peak_speed * effective_chord / mu
            ),
            "reynolds_definition": "rho*Umax*nominal_chord/mu",
        },
        "domain": {"lx": lx, "ly": ly, "periodic_x": True, "periodic_y": True},
        "grid": {
            "nx": nx, "ny": ny, "dx": dx, "dy": dy, "support": SUPPORT,
            "cells_per_nominal_chord": nominal_chord / dx,
        },
        "time": {
            "dt": dt, "steps": steps, "tfinal": steps * dt,
            "period": period, "steps_per_cycle": steps_per_cycle,
            "cycles": steps * dt / period,
        },
        "geometry": {
            "representation": "open_line_equal_weight_marker_centers",
            "nominal_chord": nominal_chord,
            "effective_quadrature_chord": effective_chord,
            "marker_count": marker_count,
            "marker_weight": ds,
            "target_stiffness": target_stiffness,
            "target_ids_zero_based": True,
            "pivot": "midchord",
            "body_coordinate_first_marker": body_coordinates[0].tolist(),
            "body_coordinate_last_marker": body_coordinates[-1].tolist(),
            "initial_center": initial_pose.center.tolist(),
            "initial_theta_radians": initial_pose.theta_radians,
            "conservative_minimum_boundary_clearance": min(clearances),
        },
        "motion": motion_metadata,
        "force": {
            "reported_quantity": "fluid_on_body_reaction_per_unit_span",
            "laboratory_axes": ["horizontal_x", "vertical_y"],
            "integration": "negative_eulerian_integral_of_ib_force_density",
            "urms_force_scale": q_urms,
            "umax_force_scale": q_umax,
            "coefficient_force_scale": q_urms,
            "primary_normalization": "0.5*rho*Urms^2*nominal_chord",
            "secondary_normalization": "0.5*rho*Umax^2*nominal_chord",
            "reversal_policy": (
                "motion-relative CD and CL are NaN at zero center speed and "
                "in midpoint sample cells touching a stroke reversal"
            ),
            "visual_component_mode": "global_horizontal_vertical",
            "visual_component_labels": ["CH = Fx / q_rms", "CV = Fy / q_rms"],
        },
    }
    (output_dir / "case_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/smoke"))
    parser.add_argument("--nx", type=int, default=DEFAULT_NX)
    parser.add_argument("--ny", type=int, default=DEFAULT_NY)
    parser.add_argument("--dt", type=float, default=DEFAULT_DT)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    parser.add_argument("--target-stiffness", type=float, default=TARGET_STIFFNESS)
    parser.add_argument("--lx", type=float, default=LX)
    parser.add_argument("--ly", type=float, default=LY)
    parser.add_argument("--re", type=float, default=REYNOLDS_NUMBER)
    parser.add_argument("--beta", type=float, default=STROKE_ANGLE_DEGREES)
    parser.add_argument("--alpha0", type=float, default=MEAN_ANGLE_OF_ATTACK_DEGREES)
    parser.add_argument("--pitch-amplitude", type=float, default=PITCH_AMPLITUDE_DEGREES)
    args = parser.parse_args()
    metadata = build_case(
        args.output_dir, nx=args.nx, ny=args.ny, dt=args.dt, steps=args.steps,
        target_stiffness=args.target_stiffness, lx=args.lx, ly=args.ly,
        reynolds_number=args.re, stroke_angle_degrees=args.beta,
        mean_angle_of_attack_degrees=args.alpha0,
        pitch_amplitude_degrees=args.pitch_amplitude,
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
