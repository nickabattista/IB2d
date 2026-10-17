# Prescribed flapping plate (headless pyIB2d example)

This example runs a prescribed two-dimensional flapping plate with the pyIB2d
solver in `../../IBM_Blackbox`, headlessly (no plotting) and with NumPy as the
only third-party dependency. The staged study holds `Re=64` fixed and uses
`S/c=2.5`, `alpha_bar=60 degrees`, `phi=0`, and a midchord pivot. It calibrates
the numerical setup on the conservative `beta=62.8 degrees`,
`alpha_a=45 degrees` motion before running the six requested combinations of
`beta={0,30,62.8}` and `alpha_a={30,45}` degrees. The plate uses target points
with `ds=dx/2`; target stiffness is treated as a continuum force-density
coefficient and is not scaled with marker count.

The companion `../Translating_Plate` example is the older regression and
provenance case (a ramped translation, with a MATLAB/Python parity history);
its Re=128 study is documented in `docs/validation/RE128_STUDY.md`.

Like the rest of IB2d, everything here is licensed under GPLv3 (see the
repository `license` file). The example was written in 2026 in Laura Miller's
group; the changes it needed in `IBM_Blackbox` are listed at the end.

## Quick start

From `pyIB2d/Examples/Flapping_Plate`:

```text
python main2d.py --regenerate --nx 32 --steps 8
python -m unittest discover -s tests -v
python run_re64_flapping_pilot.py show
python run_re64_golden_study.py manifest --stage stiffness --output stiffness.json
```

Unlike the other pyIB2d examples, `main2d.py` does not read a hand-written
`input2d`: `generate_case.py` writes `input2d`, `plate.vertex`, `plate.target`
and `case_metadata.json` into the case directory from the command-line
parameters, with `Tfinal = steps * dt`. The smoke run above writes
`force_history.csv` below `runs/smoke` (`runs/` is generated output and is
git-ignored). Its values are dimensional fluid-on-body `Fx` and `Fy` per unit
span plus force-balance and target-following diagnostics. `analyze_forces.py`
derives fixed-frame `CH`/`CV` and motion-relative `CD`/`CL`, using both
peak-speed and RMS-speed normalizations. Motion-relative coefficients are
undefined at stroke reversal; the exact reversal and midpoint sample cells
touching it are deliberately left blank.

## Files

| File | Purpose |
| --- | --- |
| `case_config.py`, `kinematics.py` | physical parameters and the analytical stroke and pitch motion |
| `generate_case.py` | writes a complete case directory (`input2d`, vertex and target files, metadata) |
| `update_Target_Point_Positions.py` | the IB2d hook that moves the target points (absolute positions, evaluated at the force midpoint) |
| `force_history.py` | `ForceHistoryWriter`, the `force_callback` that records the body load every step |
| `main2d.py` | generates (optionally) and runs one case headlessly |
| `analyze_forces.py`, `flapping_visual.py` | post-processing: coefficients, cycle statistics, dependency-free SVG/HTML rendering |
| `run_re64_flapping_pilot.py` | the two-timestep pilot (manifest, run, analyze) |
| `run_re64_golden_study.py` | the staged golden study: stiffness, timestep, box, grid, startup, production |
| `tests/` | unit tests, including a real two-step IBM solve and the periodic force-balance helper |
| `docs/validation/` | study plans and review criteria (`RE64_GOLDEN_STUDY.md`, `RE64_FLAPPING_PILOT.md`, `RE128_STUDY.md`, `README.md`) and `compare_histories.py` |
| `docs/hpc/` | Slurm scripts and submit wrappers for UArizona Puma; start with `README_PUMA.md` |
| `docs/VALIDATION_STATUS.md` | what has and has not been validated (state 2026-10-14) |

## Puma golden study

The production workflow is deliberately staged: stiffness, timestep, box,
grid, startup, and finally the six kinematic cases. Each stage uses immutable
manifests and a dependent analysis job, and later stages should not be launched
until the earlier evidence has been reviewed. On Puma, with the repository
cloned under your allocation, set the allocation-specific values and submit the
first gate from the repository root:

```bash
PI_ACCOUNT=your_pi_account
STUDY_ROOT=/xdisk/your_pi_account/your_netid/plate_re64_golden
cd /xdisk/your_pi_account/your_netid/IB2d
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" stiffness "$STUDY_ROOT" \
  --stiffness-values 2.5e8,5e8,1e9,2e9 \
  --dt 6.25e-7 --box-chords 16 --cells-per-chord 32 --cycles 1
```

The standard study does not use IB2d damped springs: they connect marker pairs
and are not target-point dashpots. It instead selects the smallest converged
target stiffness and largest converged coupling timestep. The periodic fluid
uses an explicitly recorded uniform balancing force to keep its spatially
averaged momentum fixed; the saved plate load remains the original IB
reaction. Exact matrices, decision gates, follow-on commands, resource
estimates, and this modeling convention are in
`docs/validation/RE64_GOLDEN_STUDY.md`.

The older two-timestep flapping pilot remains available in
`docs/validation/RE64_FLAPPING_PILOT.md`. It is an implementation check, not a
spatial-convergence or periodic-image study.

Current evidence and remaining work are recorded in
`docs/VALIDATION_STATUS.md`. The research papers are not included.

The manifests record SHA-256 fingerprints of the numerical source files keyed
by their path relative to the repository root (for example
`pyIB2d/IBM_Blackbox/IBM_Driver.py`), so manifests produced in the original
standalone project layout (`python/...`) do not match manifests produced here.

## Changes to IBM_Blackbox

The example depends on changes to `pyIB2d/IBM_Blackbox` made from 2026-10-10
through 2026-10-14 in Laura Miller's group; they are on the same branch, and
each modified file carries a MODIFICATION NOTICE comment. Without the new
optional arguments the solver behaves exactly as before.

- `IBM_Driver.py`: `main()` accepts optional `force_callback`, `state_callback`
  (with `state_callback_interval`) and `balance_periodic_force` arguments. The
  force callback receives the signed body load every step, evaluated before any
  external, Boussinesq or Brinkman forcing is added; the state callback
  receives read-only velocity and position snapshots at a chosen cadence (used
  by the visual pilots). `balance_periodic_force=True` subtracts the spatial
  mean of the Eulerian force before the periodic fluid update, a uniform
  counter-force that holds the box-mean momentum fixed; recorded loads remain
  the original immersed-boundary reaction. The time loop now runs exactly
  `round(Tfinal/dt)` steps and sets `current_time = step*dt` instead of
  accumulating round-off; `Tfinal` must be an integer multiple of `dt`
  (`ValueError` otherwise). `please_Plot_Results` is imported only when
  plotting is requested, so matplotlib is not needed headlessly. The optional
  compiled writer (`write` extension, `c_write.c`, `setup.py`, `wrtFile.pyx`)
  is unchanged.
- `Supp.py`: the Peskin four-point delta kernel is a vectorized NumPy
  expression; Numba is no longer imported anywhere in `IBM_Blackbox`.
- `please_Find_Lagrangian_Forces_On_Eulerian_grid.py`: force spreading
  accumulates each marker's local four-by-four stencil (`np.add.at`) instead of
  forming dense `Ny x Nb` and `Nb x Nx` delta matrices, so it scales as
  O(supp^2 Nb) instead of O(Nx Ny Nb). `../Translating_Plate/tests/test_spreading.py`
  checks agreement with the dense product. A non-finite Lagrangian state raises
  `FloatingPointError` with a hint to reduce `dt` or the stiffness.
- `give_Me_Coagulation_Force_Densities.py`: the bond check
  `if np.where(checkVec == j)` was always true; it is now `np.any(...)`.
- `please_Read_input2d_File.py`: `FID.close()` and `FID_temp.close()` are
  called (previously the method objects were only referenced).
