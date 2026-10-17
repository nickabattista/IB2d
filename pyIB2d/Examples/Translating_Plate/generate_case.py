"""Generate a self-contained translating-plate IB2d run directory."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from case_config import (
    ANGLE_DEGREES,
    CASE_NAME,
    DEFAULT_DT,
    DEFAULT_NX,
    DEFAULT_NY,
    DEFAULT_STEPS,
    LEADING_EDGE_X,
    LEADING_EDGE_Y,
    LX,
    LY,
    MU,
    NOMINAL_CHORD,
    RAMP_RATE,
    RHO,
    SUPPORT,
    TARGET_STIFFNESS,
    TRANSLATION_SPEED,
    TRANSLATION_STOP_TIME,
)


def _input_text(nx: int, ny: int, dt: float, steps: int, *,
                lx: float = LX, ly: float = LY, rho: float = RHO,
                mu: float = MU) -> str:
    tfinal = steps * dt
    return f"""% Generated translating-plate regression case.
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
string_name = \"plate\"
}}
"""


def build_case(output_dir: Path, nx: int = DEFAULT_NX, ny: int = DEFAULT_NY,
               dt: float = DEFAULT_DT, steps: int = DEFAULT_STEPS,
               target_stiffness: float = TARGET_STIFFNESS, *,
               lx: float = LX, ly: float = LY, rho: float = RHO,
               mu: float = MU, nominal_chord: float = NOMINAL_CHORD,
               leading_edge_x: float = LEADING_EDGE_X,
               leading_edge_y: float = LEADING_EDGE_Y,
               angle_degrees: float = ANGLE_DEGREES,
               translation_speed: float = TRANSLATION_SPEED,
               ramp_rate: float = RAMP_RATE,
               translation_stop_time: float = TRANSLATION_STOP_TIME) -> dict:
    """Write solver inputs and return the metadata dictionary."""

    output_dir = Path(output_dir).resolve()
    if nx < 8 or ny < 8 or nx % 2 or ny % 2:
        raise ValueError("nx and ny must be even integers of at least 8")
    if dt <= 0.0 or steps < 1:
        raise ValueError("dt must be positive and steps must be at least one")
    if target_stiffness <= 0.0:
        raise ValueError("target_stiffness must be positive")
    if min(lx, ly, rho, mu, nominal_chord) <= 0.0:
        raise ValueError("domain, fluid, and chord parameters must be positive")
    if ramp_rate <= 0.0 or translation_stop_time < 0.0:
        raise ValueError("ramp_rate must be positive and stop_time non-negative")

    output_dir.mkdir(parents=True, exist_ok=True)

    dx = lx / nx
    dy = ly / ny
    if not math.isclose(dx, dy, rel_tol=1e-12, abs_tol=1e-15):
        raise ValueError(
            "this solver path requires square Eulerian cells (lx/nx == ly/ny)"
        )

    # IBM_Driver uses ds = Lx/(2*Nx).  Enforcing dx == dy above keeps the
    # generated geometry consistent with that solver-side quadrature weight.
    ds = 0.5 * dx
    marker_count = max(2, int(round(nominal_chord / ds)))
    effective_chord = marker_count * ds

    # The markers represent equal-width cells centered along the open plate.
    # Consequently marker_count*ds, rather than endpoint span, is the length
    # used by the legacy constant-ds force quadrature.
    arc = (np.arange(marker_count, dtype=float) + 0.5) * ds
    angle = math.radians(angle_degrees)
    coordinates = np.column_stack((
        leading_edge_x + arc * math.cos(angle),
        leading_edge_y + arc * math.sin(angle),
    ))

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
        _input_text(
            nx=nx, ny=ny, dt=dt, steps=steps,
            lx=lx, ly=ly, rho=rho, mu=mu,
        ),
        encoding="utf-8",
    )

    nominal_reynolds = rho * abs(translation_speed) * nominal_chord / mu
    effective_reynolds = rho * abs(translation_speed) * effective_chord / mu
    force_scale = 0.5 * rho * translation_speed * translation_speed * nominal_chord

    metadata = {
        "schema_version": 1,
        "case_name": CASE_NAME,
        "fluid": {
            "rho": rho,
            "mu": mu,
            "nominal_reynolds_number": nominal_reynolds,
            "effective_reynolds_number": effective_reynolds,
            "reynolds_definition": "rho*abs(translation_speed)*chord/mu",
        },
        "domain": {"lx": lx, "ly": ly, "periodic_x": True, "periodic_y": True},
        "grid": {
            "nx": nx,
            "ny": ny,
            "dx": dx,
            "dy": dy,
            "support": SUPPORT,
            "cells_per_nominal_chord": nominal_chord / dx,
        },
        "time": {"dt": dt, "steps": steps, "tfinal": steps * dt},
        "geometry": {
            "representation": "open_line_equal_weight_marker_centers",
            "nominal_chord": nominal_chord,
            "effective_quadrature_chord": effective_chord,
            "marker_count": marker_count,
            "marker_weight": ds,
            "leading_edge": [leading_edge_x, leading_edge_y],
            "angle_degrees": angle_degrees,
            "target_stiffness": target_stiffness,
            "target_ids_zero_based": True,
        },
        "motion": {
            "kind": "ramped_translation_x",
            "speed": translation_speed,
            "ramp_rate": ramp_rate,
            "stop_time": translation_stop_time,
            "evaluation_time": "step_midpoint",
        },
        "force": {
            "reported_quantity": "fluid_on_body_reaction_per_unit_span",
            "laboratory_axes": ["x", "y"],
            "integration": "negative_eulerian_integral_of_ib_force_density",
            "coefficient_force_scale": force_scale,
            "drag_convention": "drag=-fx_body for translation in +x",
            "lift_convention": "lift=fy_body in laboratory +y",
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
    parser.add_argument("--rho", type=float, default=RHO)
    parser.add_argument("--mu", type=float, default=MU)
    parser.add_argument("--nominal-chord", type=float, default=NOMINAL_CHORD)
    parser.add_argument("--leading-edge-x", type=float, default=LEADING_EDGE_X)
    parser.add_argument("--leading-edge-y", type=float, default=LEADING_EDGE_Y)
    parser.add_argument("--angle-degrees", type=float, default=ANGLE_DEGREES)
    parser.add_argument("--translation-speed", type=float, default=TRANSLATION_SPEED)
    parser.add_argument("--ramp-rate", type=float, default=RAMP_RATE)
    parser.add_argument("--translation-stop-time", type=float, default=TRANSLATION_STOP_TIME)
    args = parser.parse_args()
    metadata = build_case(
        args.output_dir, args.nx, args.ny, args.dt, args.steps,
        target_stiffness=args.target_stiffness,
        lx=args.lx,
        ly=args.ly,
        rho=args.rho,
        mu=args.mu,
        nominal_chord=args.nominal_chord,
        leading_edge_x=args.leading_edge_x,
        leading_edge_y=args.leading_edge_y,
        angle_degrees=args.angle_degrees,
        translation_speed=args.translation_speed,
        ramp_rate=args.ramp_rate,
        translation_stop_time=args.translation_stop_time,
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
