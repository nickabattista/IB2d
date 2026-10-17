# Running the paired validation on UArizona Puma

These scripts run the matched translating-plate regression in Python and
MATLAB and then compare their full force histories.  Submit them from the
repository root (the directory containing `pyIB2d/`) so that `SLURM_SUBMIT_DIR` points at the working copy.  For a
Puma account `your_netid` in the PI allocation `your_pi_account`, the working copy is:

```text
/xdisk/your_pi_account/your_netid/IB2d
```

Upload and unpack the archive in that allocation:

```text
scp IB2d.tar.gz your_netid@filexfer.hpc.arizona.edu:/xdisk/your_pi_account/your_netid/
cd /xdisk/your_pi_account/your_netid
tar -xzf IB2d.tar.gz
cd /xdisk/your_pi_account/your_netid/IB2d
```

The smoke script uses the preemptible `windfall` partition and does not need a
PI account.  It runs only two timesteps at 32 by 32.  Inspect its comparison
before launching the matrix.

```bash
sbatch pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_smoke.slurm
```

The matrix script runs 128 by 128 and 256 by 256 through `t=0.025`.  Supply a
PI group reported by `va` when submitting to the `standard` partition:

```bash
va
sbatch --account=YOUR_PI_GROUP --partition=standard pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_matrix.slurm
```

The module defaults are `python/3.11/3.11.4` and Puma's default `matlab`
module.  For a reproducible run, inspect the versions on an interactive
compute node and override MATLAB with the exact installed name without editing
the script:

```bash
sbatch --export=ALL,MATLAB_MODULE=matlab/EXACT_VERSION pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_smoke.slurm
```

Monitor and inspect jobs with:

```bash
squeue --me
tail -f puma_smoke_JOBID.out
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

If a current MATLAB release requests MathWorks authentication, start a Puma
interactive compute session, load the same MATLAB module, launch `matlab`, and
complete the one-time university sign-in.  Do not launch MATLAB or the solver
on a login node.

## Re = 128 convergence study

The completed 128/256 matrix above is a MATLAB/Python language-parity test, not
a spatial-convergence result.  The staged Python-only Re=128 study, including
its visual pilot, stiffness/timestep gate, `ds=dx/2` refinement with fixed
continuum stiffness and timestep, Slurm arrays, periodic-box controls,
analysis, and download commands, is
documented in `../validation/RE128_STUDY.md`.

The next job is the single visual pilot—not the large arrays:

```bash
cd /xdisk/your_pi_account/your_netid/IB2d
va
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_visual.sh YOUR_PI_GROUP
```

## Re = 64 flapping pilot

The prescribed flapping case has its own two-timestep pilot.  It uses the
exact tilted kinematics (`beta=62.8` degrees, `alpha0=60` degrees), runs two
cycles at `512^2`, and creates the movie and force diagnostics only after both
timestep runs finish successfully.  It does not reuse the older Re=128
translating-plate result directory.

```bash
cd /xdisk/your_pi_account/your_netid/IB2d
bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_flapping_pilot.sh your_pi_account
```

The submit command prints the exact array and analysis job IDs together with
ready-to-paste monitoring commands.  Results are written beneath
`/xdisk/your_pi_account/your_netid/plate_re64_flapping_pilot`.  Full parameters, review
criteria, packaging, and PowerShell download commands are in
`../validation/RE64_FLAPPING_PILOT.md`.
