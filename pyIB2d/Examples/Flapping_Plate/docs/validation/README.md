# Validation map

The prescribed translating plate was the first cross-language regression.
MATLAB and Python used identical physical inputs, marker coordinates, midpoint
target motion, timestep, and output schema; target IDs were deliberately
one-based in MATLAB and zero-based in Python. That comparison is retained
below as solver evidence, but it is not the final physics case.

The current gate is the Python Re=64 tilted flapping pilot. Its exact controls,
two-timestep Puma workflow, force definitions, review criteria, and download
commands are in `RE64_FLAPPING_PILOT.md`. The earlier translating-grid plan is
retained in `RE128_STUDY.md`.

## Stage A: two-step smoke run

Python, from `pyIB2d/Examples/Translating_Plate`:

```text
python main2d.py --case-dir runs/hpc_smoke --regenerate --nx 32 --ny 32 --steps 2
```

MATLAB, from `matlab/Examples/Translating_Plate`:

```matlab
source_dir = pwd;
case_dir = fullfile(source_dir,'runs','hpc_smoke');
generate_case(case_dir,32,32,2.5e-5,2);
main2d(case_dir);
```

Compare from the workspace root:

```text
python pyIB2d/Examples/Flapping_Plate/docs/validation/compare_histories.py \
  pyIB2d/Examples/Translating_Plate/runs/hpc_smoke/force_history.csv \
  matlab/Examples/Translating_Plate/runs/hpc_smoke/force_history.csv \
  --output pyIB2d/Examples/Flapping_Plate/docs/validation/smoke_comparison.json
```

Do not enforce a tolerance on the first cross-language run.  Inspect the
history and report first, then set a defensible tolerance based on the actual
solver-version differences.

## Stage B: reference history

After the smoke comparison is understood, run `nx=128`, `dt=2.5e-5`, and
`steps=1000` in both languages.  This reaches `t=0.025`, the end of the ramped
translation.  Then repeat at 256.  A 512 run is useful as an offline spatial
reference; reserve 1024 for MATLAB unless the 512 Python timing is acceptable.

The comparison gate checks the complete force histories.  It also reports:

- force-on-body plus force-on-fluid residual;
- Eulerian versus Lagrangian force-integration residual;
- maximum target-following error;
- absolute and normalized MATLAB/Python differences.

The translating case passes only when the sign convention and full histories
are understood.  Matching the legacy `mean(F_Lag)` quantity is not a goal.
