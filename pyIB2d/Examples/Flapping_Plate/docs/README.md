# Flapping-plate study documents

- `validation/` holds the study plans and review criteria: `RE64_GOLDEN_STUDY.md`
  (the staged Re=64 golden study), `RE64_FLAPPING_PILOT.md` (the older
  two-timestep pilot), `RE128_STUDY.md` (the translating-plate Re=128 study),
  `README.md` (the MATLAB/Python parity stages) and `compare_histories.py`, which
  compares two `force_history.csv` files.
- `hpc/` holds the Slurm scripts and submit wrappers for UArizona Puma; start
  with `hpc/README_PUMA.md`.
- `VALIDATION_STATUS.md` records what has and has not been validated.

These files were written for a standalone project whose Python tree was called
`python/` and whose `hpc/` and `validation/` folders sat at the project root.
In this repository the paths have been changed to `pyIB2d/Examples/...` and
`pyIB2d/Examples/Flapping_Plate/docs/{hpc,validation}`, and the submit wrappers
derive the repository root from their own location, so run them from the
repository root (the directory containing `pyIB2d/`).

Account and path components are placeholders: `your_pi_account` is the Slurm
account reported by `va`, `your_netid` the user name, and `C:\path\to\...` a
local download folder.

The paired MATLAB/Python scripts (`hpc/puma_smoke.slurm`, `hpc/puma_matrix.slurm`
and stages A/B of `validation/README.md`) refer to a MATLAB `Translating_Plate`
example under a `matlab/Examples/` tree that is not part of this branch; they are
kept as the record of the language-parity check reported in
`VALIDATION_STATUS.md` and will not run here without that example.
