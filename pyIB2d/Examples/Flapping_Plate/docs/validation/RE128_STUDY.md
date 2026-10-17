# Re = 128 convergence and periodic-box study

## Scope

This study applies to the currently validated **ramped translating plate**.
The final heave/pitch flapping law has not yet been specified, so these runs
must not be described as flapping-cycle results.  The infrastructure is
manifest-driven and can be reused after the flapping motion, reference
velocity, pivot, and averaging window are fixed.

IB2d uses periodic boundaries in both directions.  The "boundary" study is
therefore a **periodic-box/image sensitivity study**, not a solid-wall
clearance study.

The workflow is deliberately staged.  First run one visual pilot and inspect
the plate, vorticity, target following, and unsmoothed force histories.  Do not
submit the calibration or production arrays until that pilot has been reviewed.

## Controlled physical parameters

- Nominal chord: `c = 0.1`
- Density: `rho = 10`
- Translation reference speed: `U = 24`
- Exact nominal Reynolds number: `Re = rho*U*c/mu = 128`
- Dynamic viscosity used by the study: `mu = 0.1875`
- Plate angle: 45 degrees
- Active translation ends at `t = 0.025`, or `tU/c = 6`
- The complete swept envelope is centered independently in each square box.
- Lift and drag normalization per unit span is
  `0.5*rho*U^2*c = 288`.
- Positive drag is opposite the imposed `+x` translation:
  `CD = -Fx_body/288`; `CL = Fy_body/288`.

Every selected grid represents `c=0.1` exactly.  This avoids the small
resolution-dependent chord and Reynolds-number changes in the earlier
128/256 language-parity matrix.

The Lagrangian spacing is always

```text
ds = dx/2
```

The target routine computes the Lagrangian force density `f=kT*(Y-X)` and the
spreading routine applies `f*ds`.  Consequently, the whole-plate force is the
quadrature sum `sum(f*ds)`: refining the marker spacing already adds more
terms.  The same physical continuum plate therefore uses

```text
kT = constant
dt = constant after a separate timestep calibration
```

Scaling `kT` as `1/ds` would make the plate increasingly stiff with grid
refinement.  Constant `kT` and `dt` make the grid family a fixed-parameter
spatial comparison.  The relevant explicit stability diagnostic is the
spread/interpolate tether number

```text
chi = kT*dt^2*ds/(rho*dx^2)
```

not the advective CFL.  With `ds=dx/2`, `chi` grows as `1/dx` at fixed `dt` and
is therefore checked on the finest grid.  We require `chi <= 1`; the
conservative visual pilot gives `chi=0.05` and the finest planned production
grid gives `chi=0.4`.

This choice is consistent with the shipped IB2d examples.  The controlled
Jellyfish, VIV, tethered-jellyfish, and pendulum resolution families keep raw
`kT` constant as marker count increases.  The closest original insect-wing
case uses `rho=10`, `dx=1/1024`, `ds=dx/2`, `kT=5e7`, and `dt=2.5e-5`, for
`chi=1.6`.  Our pilot timestep is twenty times smaller and its `chi` is 32
times smaller while using the stiffer target needed for low plate lag.

## Required Stage 0: visual pilot

The first run is one complete 1024-square case in the 16-chord periodic box:

| Quantity | Value |
| --- | ---: |
| Cells per chord | 64 |
| Eulerian grid | 1024 x 1024 |
| `dx` | 0.0015625 |
| `ds` | 0.00078125 = `dx/2` |
| Target stiffness | 1.0e9 |
| `dt` | 1.25e-6 |
| Steps | 20000 |
| Tether number `chi` | 0.05 |
| Snapshot cadence | every 200 steps, plus initial/final |

The run saves 101 compressed vorticity/plate snapshots and produces:

- a browser-playable `animation.html` with actual and prescribed plates;
- consistent-color PNG frames and, when `ffmpeg` is available, an MP4;
- `force_and_target_history.svg` with raw `CD`, `CL`, maximum target error,
  and RMS target error normalized by `dx` versus time;
- `oscillation_diagnostics.json`, including a cubic-detrended spectral energy
  fraction above 0.1 cycles per step.

This deliberately conservative timestep matches the small-`dt` calibration
reference; the later gate decides whether `2.5e-6` is also acceptable for
production.  This stage has no automatic visual verdict.  Check that the actual plate
tracks the target, the wake is plausible and clear of its periodic images,
and the force curves do not exhibit alternating step-scale ringing.

## Required Stage 1: stiffness/timestep calibration

Target points impose the plate through penalty springs.  Existing results
show that language parity is excellent but do not establish a rigid-plate
limit.  Production is therefore gated by five 1024-square runs in a 16-chord
box:

| Index | Target stiffness | dt | Purpose |
| ---: | ---: | ---: | --- |
| 0 | 5.0e8 | 1.25e-6 | stiffness sequence |
| 1 | 1.0e9 | 1.25e-6 | small-timestep reference + stiffness sequence |
| 2 | 2.0e9 | 1.25e-6 | stiff-limit reference |
| 3 | 1.0e9 | 2.5e-6 | proposed production timestep |
| 4 | 1.0e9 | 5.0e-6 | larger-timestep ringing diagnostic |

The analysis compares complete `CD(t)` and `CL(t)` histories.  Its default
gate requires the nearest stiffness pair and proposed timestep to differ by no
more than 2% in relative L2.  The **actual proposed pair** (`kT=1e9`,
`dt=2.5e-6`) must also satisfy maximum lag below 1% chord and 0.5 `dx`, RMS lag
below 0.25 `dx`, `chi<=1`, and a high-frequency energy gate.  The ringing gate
requires less than 5% energy at or above 0.1 cycles per step and no more than
three times the `dt/2` reference.  These thresholds are documented numerical
acceptance criteria, not universal physical constants.

## Stage 2: spatial and periodic-box study

The provisional common values are `kT=1e9` and `dt=2.5e-6`.  Do not launch
this stage until the visual pilot is acceptable and
the calibration report says `production_ready: true`, or a different pair has
been intentionally reviewed.

| Index | Case | L/c | Grid | Pts/chord | `dt` | `kT` | `chi` | Axis |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | `grid_L16_P32_N512` | 16 | 512² | 32 | 2.5e-6 | 1.0e9 | 0.1 | grid |
| 1 | `boundary_L10_P64_N640` | 10 | 640² | 64 | 2.5e-6 | 1.0e9 | 0.2 | box |
| 2 | `shared_L16_P64_N1024` | 16 | 1024² | 64 | 2.5e-6 | 1.0e9 | 0.2 | grid + box |
| 3 | `boundary_L24_P64_N1536` | 24 | 1536² | 64 | 2.5e-6 | 1.0e9 | 0.2 | box |
| 4 | `grid_L16_P128_N2048` | 16 | 2048² | 128 | 2.5e-6 | 1.0e9 | 0.4 | grid |
| 5 | `boundary_L32_P64_N2048` | 32 | 2048² | 64 | 2.5e-6 | 1.0e9 | 0.2 | box |

The grid family holds `L/c=16`, `dt`, and `kT` fixed while refining `dx` and
maintaining `ds=dx/2`.  This preserves the same continuum plate compliance and
separates spatial error from a changing timestep or penalty law.  The box
family holds `dx`, `ds`, marker count, timestep, stiffness, geometry, and motion
fixed while increasing `L/c`.  Index 2 is shared.  The maximum size is 2048,
not 2056, because 2048 is FFT-friendly and gives the exact refinement ratios
needed here.

The report includes:

- full-history relative L2 and Linf differences for `CD` and `CL`;
- observed grid order from the three systematically refined histories;
- minimum, maximum, mean, RMS, and nondimensional impulse for each force;
- maximum target lag divided by both chord and `dx`;
- immersed-boundary force-balance residuals;
- exact-Re and exact-chord checks.

## Puma commands

Work from the fixed project and result locations:

```bash
cd /xdisk/your_pi_account/your_netid/IB2d
va
```

Replace `YOUR_PI_GROUP` with the account printed by `va`.  Start only the
visual pilot:

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_visual.sh YOUR_PI_GROUP
```

Monitor it with the printed job ID:

```bash
squeue --me
sacct -j VISUAL_JOB_ID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
tail -f /xdisk/your_pi_account/your_netid/plate_re128_visual_pilot/logs/re128_visual_VISUAL_JOB_ID.out
```

When it finishes, package the compact snapshots and rendered diagnostics:

```bash
cd /xdisk/your_pi_account/your_netid
tar -czf plate_re128_visual_pilot_results.tar.gz \
  --exclude='*.partial-*' \
  plate_re128_visual_pilot
sha256sum plate_re128_visual_pilot_results.tar.gz
```

Download in Windows PowerShell.  The destination deliberately has no trailing
backslash inside the quotes:

```powershell
scp your_netid@filexfer.hpc.arizona.edu:/xdisk/your_pi_account/your_netid/plate_re128_visual_pilot_results.tar.gz `
  "C:\path\to\puma_results"
```

Extract it and open
`plate_re128_visual_pilot/visual_L16_P64_N1024/visual_pilot/animation.html`;
the force plot is embedded below the movie.  Only after we review those
artifacts should you submit the
calibration:

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_study.sh YOUR_PI_GROUP calibration
```

The command prints the array and dependent analysis job IDs.  Monitor them:

```bash
squeue --me
sacct -j ARRAY_JOB_ID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
tail -f /xdisk/your_pi_account/your_netid/plate_re128_study/logs/re128_ARRAY_JOB_ID_0.out
```

After the dependent analysis finishes:

```bash
python3 -m json.tool \
  /xdisk/your_pi_account/your_netid/plate_re128_study/summary/calibration_report.json \
  | less
```

If `production_ready` is true, submit production:

```bash
cd /xdisk/your_pi_account/your_netid/IB2d
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_study.sh YOUR_PI_GROUP production
```

If it is false, stop before the expensive jobs and download the calibration
report and case CSVs for review.  To intentionally use a reviewed alternative
production pair, set the values and explicitly acknowledge the failed default
gate, for example:

```bash
export PRODUCTION_DT=2.5e-6
export PRODUCTION_TARGET_STIFFNESS=1e9
export ALLOW_UNCALIBRATED=1
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_study.sh YOUR_PI_GROUP production
```

These environment values are common to every production grid; the driver does
not rescale `dt` or `kT` with marker count.  The
production runs are split into low- and high-resolution arrays, with a
dependent analysis job.  They use the standard partition because the solver
has no checkpoint/restart support and long Windfall jobs can be preempted.

## Results and download

Results are kept outside the source tree:

```text
/xdisk/your_pi_account/your_netid/plate_re128_visual_pilot/
  visual_L16_P64_N1024/
    visual_snapshots/
    visual_pilot/
  logs/

/xdisk/your_pi_account/your_netid/plate_re128_study/
  cases/calibration/
  cases/production/
  logs/
  summary/
  manifest_calibration.json
  manifest_production.json
```

Each case is first written to a unique `.partial-JOBID-TASKID` directory.  A
validated run receives `COMPLETE.json` and is renamed to its final case ID.
There is no numerical restart: an interrupted partial case must be rerun from
time zero.

Package only the compact histories, metadata, reports, manifests, and logs:

```bash
cd /xdisk/your_pi_account/your_netid
tar -czf plate_re128_study_results.tar.gz \
  --exclude='*.partial-*' \
  plate_re128_study
sha256sum plate_re128_study_results.tar.gz
```

Download from Windows PowerShell:

```powershell
scp your_netid@filexfer.hpc.arizona.edu:/xdisk/your_pi_account/your_netid/plate_re128_study_results.tar.gz `
  "C:\path\to\puma_results"
```
