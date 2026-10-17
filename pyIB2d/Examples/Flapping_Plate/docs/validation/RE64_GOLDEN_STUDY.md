# Re=64 flapping-plate golden study

## Purpose

This runbook selects a defensible numerical setup before running the six UAV
kinematic cases.  It separates the questions in this order:

1. target-point stiffness;
2. immersed-boundary coupling timestep;
3. periodic-box size;
4. Eulerian and Lagrangian resolution;
5. removal of startup transients; and
6. the six production kinematics.

Do not launch all stages at once.  Review each stage and pass its selected
values to the next stage.  The representative calibration motion is
`beta=62.8 degrees`, `alpha_a=45 degrees`.  Every case uses `Re=64`,
`c=0.1`, `S/c=2.5`, `T=0.032725`, `rho=10`, `alpha_bar=60 degrees`, and
`phi=0`.  The viscosity is computed from these values rather than set
independently.

For reference, these inputs give `S=0.25`,
`U_ref=pi*S/T=23.9999439`, and `mu=rho*U_ref*c/Re=0.3749991`
(rounded here; the manifest stores the computed values).

The Eulerian grid is square and uniform.  Lagrangian points always use
`ds=dx/2`.  The target stiffness is a continuum force-density penalty, so it
is held constant when the number of Lagrangian points changes.  Multiplication
by `ds` during force spreading preserves the whole-plate force integral.

## Why the standard study has no dampers

Use target-point springs without structural dampers for the primary study.
The existing IB2d `damped_springs` model connects pairs of Lagrangian points;
it is not a dashpot between a marker and its prescribed target.  Turning it on
would add a different structural model and its force would enter the measured
plate load.  The viscous fluid already supplies physical dissipation.

Select the smallest spring stiffness that tracks the prescribed motion and
then select a timestep that removes coupling-scale force oscillations.  If no
spring-only stiffness/timestep pair meets the gates below, stop and revisit
the forcing formulation.  Do not hide the failure by adding an uncalibrated
damper.

The useful coupling diagnostic recorded in every manifest is

```text
tether_number = k * dt^2 * ds / (rho * dx^2)
              = k * dt^2 / (2 * rho * dx), because ds=dx/2.
```

It helps compare cases, but it is not by itself a stability theorem.  Force
and tracking convergence remain the deciding evidence.

## Periodic mean-force convention

A plate with nonzero mean lift applies an equal and opposite nonzero mean
force to the fluid.  In a fully periodic box, that force drives the spatially
uniform velocity mode.  Without a momentum sink, the domain-mean fluid
velocity continually changes and an exact time-periodic full-flow state need
not exist.

By default, the golden-study runner subtracts the spatial mean of the spread
immersed-boundary force before the fluid update.  This is a uniform balancing
body force that holds the initially zero domain-mean velocity fixed, up to
roundoff.  The force saved for the plate is the original, uncorrected
immersed-boundary reaction.  Thus:

- `Fx`, `Fy`, lift, and drag remain the raw fluid-on-plate load;
- the simulated periodic fluid includes a mean-momentum constraint;
- the reported unbalanced bulk-flow drift is counterfactual.  It is the drift
  that would have occurred without the balancing force, or equivalently the
  accumulated balancing impulse divided by fluid mass.

This convention must be disclosed with the final results.  It is not a wall
boundary condition and should not be described as an infinite-domain
simulation.  The box-size study is what demonstrates that the local plate
force is insensitive to the periodic copies and to the uniform balancing
force.  The manifest option `--no-balance-periodic-force` is available only
for a separately labeled sensitivity calculation; it should not be mixed into
the calibrated production set.

## Default calibration matrices

All dimensions below are in chord-based terms.  `P` is Eulerian cells per
chord and the square grid size is `N=(L/c)P`.

| Stage | Values varied | Other starting values | Cycles | Purpose |
| --- | --- | --- | ---: | --- |
| stiffness | `k={2.5e8,5e8,1e9,2e9}` | `dt=6.25e-7`, `L/c=16`, `P=32` (`N=512`) | 1 | Find a nearly rigid but not excessively stiff target penalty. |
| timestep | `dt={2.5e-6,1.25e-6,6.25e-7}` | selected `k`, `L/c=16`, `P=32` | 2 | Test the coupled force-spreading/interpolation stability, not only CFL. |
| box | `L/c={12,16,24}` | selected `k,dt`, `P=32` (`N=384,512,768`) | 4 | Test periodic-image, wake-return, and uniform-balance sensitivity. |
| grid | `P={16,32,64}` | selected `k,dt,L/c`; `ds=dx/2` | 3 | Test spatial convergence.  Keep the selected `dt` constant. |
| startup | one selected setup | representative motion | 8 | Find a repeatable post-startup force cycle. |
| production | `beta={0,30,62.8}` crossed with `alpha_a={30,45}` degrees | all selected numerical values | 8 | Produce the requested six force histories and cycle averages. |

The matrix is deliberately staged.  The `P=64` grid case can be expensive,
but only one representative motion needs it.  Full production should use
`P=32` only if its force history passes comparison with `P=64`.

## Decision gates

The analysis writes threshold-free values as well as the following operational
steady-cycle decision.  Selection is intentionally a human-reviewed decision;
the scripts do not silently certify a setup.

### Basic validity for every case

- The case has `COMPLETE.json`, all requested force samples are finite, and
  no partial directory is being mistaken for a completed result.
- Maximum target error is at most `0.10 dx` and RMS target error is at most
  `0.05 dx` over the compared cycle.
- The force curve has no increasing step-to-step alternating component.  Check
  the Nyquist and second-difference diagnostics together with the force plot;
  a timestep refinement should reduce, not reveal, grid-scale ringing.
- The recorded force-balance residual is small relative to the plate-force
  scale.

### Calibration agreement

For both `CL` and `CD`, compare complete phase-aligned force curves, not only
their averages.  A candidate passes against the next more conservative case
when:

- relative L2 difference of each force-coefficient waveform is at most `2%`;
- the absolute change in each cycle mean is no larger than
  `max(0.01, 0.01*max(abs(mean_1),abs(mean_2)))`; and
- both cases pass the tracking and oscillation checks above.

Choose the smallest passing stiffness, the largest passing timestep, the
smallest passing box, and the coarsest passing mesh.  In the stiffness stage,
agreement should be checked on the stiffness plateau rather than assuming the
largest `k` is best.  An excessively stiff penalty can introduce coupling
oscillations.  If the last two entries of any matrix do not agree, extend that
matrix before proceeding.

For the box stage, also inspect the representative vorticity field.  Reject a
box if a coherent wake structure crosses a periodic boundary and returns to
the plate region during the averaging interval.  Force convergence is the
quantitative gate.

The grid stage normally uses `--grid-dt-rule constant`, which isolates spatial
resolution.  Because the tether number grows as `dx` decreases, first make
sure the selected timestep remains stable at `P=64`.  If it does not, reduce
the common timestep and repeat the timestep and grid checks.  Do not compare
different grids at different timesteps and call the result a pure spatial
convergence test.

### Startup and averaging cycle

The runner declares an operationally repeatable cycle only after two
successive cycle pairs pass all of these checks for both `CL` and `CD`:

- phase-aligned waveform relative L2 difference at most `2%`; and
- cycle-mean change no larger than the `1%`/absolute-`0.01` rule above.

The later cycle of the second passing pair is the averaging cycle.  If the
eight-cycle startup run does not reach this state, rerun a new startup study
with more cycles.  Production averages must not be reported as steady values
for cases that fail this check; extend those cases or label them transient.

Because the fluid is initially quiescent and isotropic, changing `beta`
primarily rotates the same prescribed motion relative to the fixed UAV axes.
The production analysis therefore also rotates the loads to stroke-aligned
coordinates and compares the three `beta` values.  Agreement is a useful
implementation check; fixed-axis lift and drag will differ because their
projections differ.

## Puma instructions

### 1. Set paths once

Replace the generic account and path components below with the values for the
allocation being used.  Keep the study output separate from the source tree.

```bash
PI_ACCOUNT=your_pi_account
PROJECT_ROOT=/xdisk/your_pi_account/your_netid/IB2d
STUDY_ROOT=/xdisk/your_pi_account/your_netid/plate_re64_golden

cd "$PROJECT_ROOT"
```

The submit wrapper creates one immutable output area per stage:
`$STUDY_ROOT/stiffness`, `$STUDY_ROOT/timestep`, and so on.  It runs at most
two array members concurrently and submits a dependent analysis job even when
one member fails, so the partial report identifies what needs to be rerun.

### 2. Screen stiffness

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" stiffness "$STUDY_ROOT" \
  --stiffness-values 2.5e8,5e8,1e9,2e9 \
  --dt 6.25e-7 --box-chords 16 --cells-per-chord 32 --cycles 1
```

Wait for the array and dependent analysis IDs printed by the command.  Review:

```bash
less "$STUDY_ROOT/stiffness/summary/stage_analysis.json"
column -s, -t < "$STUDY_ROOT/stiffness/summary/case_summary.csv" | less -S
column -s, -t < "$STUDY_ROOT/stiffness/summary/cross_case_comparison.csv" | less -S
```

Set the selected stiffness only after applying the gates above:

```bash
K_SELECTED=5e8    # replace with the reviewed result
```

### 3. Refine the timestep

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" timestep "$STUDY_ROOT" \
  --target-stiffness "$K_SELECTED" \
  --dt-values 2.5e-6,1.25e-6,6.25e-7 \
  --box-chords 16 --cells-per-chord 32 --cycles 2
```

Review the three files under `$STUDY_ROOT/timestep/summary`, then set:

```bash
DT_SELECTED=6.25e-7    # replace with the reviewed result
```

If the two finest timesteps do not agree, create a new study root and add a
finer value rather than overwriting the immutable manifest.

### 4. Check the periodic box

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" box "$STUDY_ROOT" \
  --target-stiffness "$K_SELECTED" --dt "$DT_SELECTED" \
  --box-values 12,16,24 --cells-per-chord 32 --cycles 4
```

Review `$STUDY_ROOT/box/summary`, then set:

```bash
BOX_SELECTED=16    # replace with the reviewed result, in chord lengths
```

### 5. Check mesh resolution

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" grid "$STUDY_ROOT" \
  --target-stiffness "$K_SELECTED" --dt "$DT_SELECTED" \
  --box-chords "$BOX_SELECTED" \
  --grid-values 16,32,64 --grid-dt-rule constant --cycles 3
```

Review `$STUDY_ROOT/grid/summary`, then set:

```bash
P_SELECTED=32    # replace with the reviewed cells-per-chord result
```

### 6. Establish the averaging cycle

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" startup "$STUDY_ROOT" \
  --target-stiffness "$K_SELECTED" --dt "$DT_SELECTED" \
  --box-chords "$BOX_SELECTED" --cells-per-chord "$P_SELECTED" \
  --cycles 8
```

Check `operational_steady_cycle` in the fixed-axis summary and the startup
stage table.  If it does not pass, use a new study root and increase `--cycles`.

### 7. Run the six production cases

```bash
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh \
  "$PI_ACCOUNT" production "$STUDY_ROOT" \
  --target-stiffness "$K_SELECTED" --dt "$DT_SELECTED" \
  --box-chords "$BOX_SELECTED" --cells-per-chord "$P_SELECTED" \
  --cycles 8 \
  --beta-values 0,30,62.8 --pitch-values 30,45
```

The production analysis writes the summary table to
`$STUDY_ROOT/production/summary/case_summary.csv`.  Each completed case has
`force_analysis/fixed_axis_forces.csv` with `Fx`, `Fy`, `L`, `D`, `CL`, and
`CD`, and `fixed_axis_cycle_summary.json` with the averaging decision and
cycle statistics.  The `beta=62.8`, `alpha_a=45` case also receives the
representative flow-animation artifacts.

### Monitoring and safe resubmission

The submit command prints exact job IDs and copy-ready `squeue`, `sacct`, and
`tail` commands.  The general forms are:

```bash
squeue -j ARRAY_JOB_ID,ANALYSIS_JOB_ID -o "%.18i %.2t %.10M %.10l %R"
sacct -j ARRAY_JOB_ID,ANALYSIS_JOB_ID \
  --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

After an array has stopped, repeating the identical submit command is safe.
The wrapper reuses the matching immutable manifest and the runner skips cases
whose `COMPLETE.json`, manifest hash, source hash, and output length all
validate.  It does not silently reuse a partial or mismatched case.  Do not
resubmit the same stage while its first array is still active.

## Runtime and resource planning

The wrapper requests the following defaults, all with array concurrency two:

| Stage | Memory per case | Walltime per case |
| --- | ---: | ---: |
| stiffness, timestep | 8 GB | 24 h |
| box, grid | 12 GB | 72 h |
| startup, production | 8 GB | 72 h |
| dependent analysis | 8 GB | 2 h |

A measured Puma run at `N=1024` used about `0.62 GB` maximum resident memory,
so these memory requests leave substantial margin without needlessly asking
for 32 GB.  Runtime is the limiting resource.  At `dt=6.25e-7`, one flap is
52,360 steps.  A prior `N=1024`, 20,000-step run took about 5 h 40 min.  Using
rough FFT-grid scaling only for planning gives these approximate compute times:

- one-cycle `N=512` calibration case: about 3--4 h;
- four-cycle largest box case (`N=768`): about 25--35 h;
- three-cycle finest grid case (`N=1024`): about 40--50 h;
- eight-cycle `N=512` startup or production case: about 25--35 h.

With six production cases and two running at once, an `N=512` production
matrix is roughly three 25--35 h waves, plus queue time.  These are estimates,
not guarantees; use `sacct` from each completed stage to update the next
request.  Eight-cycle `N=1024` production cases are likely to exceed the
default 72 h and are much less practical.  If `P=32` does not pass the `P=64`
comparison, do not launch production blindly.  Reassess the study size and set
`GOLDEN_WALLTIME` to a partition-permitted value only after using the measured
fine-grid rate.

Optional resource overrides do not require editing the scripts:

```bash
export GOLDEN_WALLTIME=96:00:00
export GOLDEN_MEMORY=12gb
export GOLDEN_ARRAY_CONCURRENCY=2
```

Confirm that the selected partition permits the requested walltime before
submission.

