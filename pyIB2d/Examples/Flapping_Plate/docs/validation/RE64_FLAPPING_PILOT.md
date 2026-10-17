# Re = 64 exact-tilted flapping pilot

## Purpose

This is the first prescribed flapping-plate checkpoint.  It is deliberately a
two-case timestep comparison, not the production dragonfly-to-Drosophila
parameter sweep and not a spatial-convergence result.  Both cases run the same
two-cycle exact-tilted motion; only the timestep changes.  The finer case also
provides the movie used for visual review.

## Fixed physical and numerical controls

| Quantity | Value |
| --- | ---: |
| Reynolds number | 64 |
| Nominal chord, `c` | 0.1 |
| Density, `rho` | 10 |
| Peak translation speed, `Umax` | 23.9999439 (nominally 24) |
| Dynamic viscosity, `mu=rho*Umax*c/Re` | 0.374999123 |
| Stroke amplitude, `A0/c` | 2.5 |
| Stroke-plane angle, `beta` | 62.8 degrees |
| Mean angle of attack, `alpha0` | 60 degrees |
| Pitch amplitude, `B` | 45 degrees |
| Pitch phase | 0 radians |
| Pivot | midchord |
| Period, `pi*A0/Umax` | 0.032725 (cycle-aligned) |
| Simulated duration | 2 periods |
| Periodic box, `L/c` | 16 |
| Eulerian grid | 512 by 512 |
| Cells per chord | 32 |
| `dx` | 0.003125 |
| Lagrangian spacing | `ds=dx/2=0.0015625` |
| Target stiffness | `1e9` in both cases |

The 0.032725 period differs from the nominal-24 value by only `2.34e-6`
relatively, while giving exactly 26,180 coarse or 52,360 fine steps per cycle.
The generator recomputes `Umax` and `mu`, keeping `Re=64` exact. The Python
driver and each completed case record the timestep, step count, complete motion
definition, coefficient normalization, numerical-source fingerprint, and
runtime environment. The manifest—not this rounded table—is authoritative.

## Timestep pair

| Array index | Nominal timestep | Tether number `chi` | Role |
| ---: | ---: | ---: | --- |
| 0 | `1.25e-6` | 0.025 | candidate pilot timestep |
| 1 | `6.25e-7` | 0.00625 | fine-timestep reference and movie |

Here `chi=kT*dt^2*ds/(rho*dx^2)`.  The small advective CFL is not the
acceptance criterion; the comparison targets the force-spreading and
velocity-interpolation stiffness identified for this problem.

## Submit on Puma

Upload the Python-only archive from Windows PowerShell:

```powershell
scp "C:\path\to\IB2d_re64_pilot.tar.gz" `
  your_netid@filexfer.hpc.arizona.edu:/xdisk/your_pi_account/your_netid/
```

Then unpack it into a fresh source directory. This avoids retaining stale
MATLAB or optional C/Cython files from an older extraction:

```bash
cd /xdisk/your_pi_account/your_netid
sha256sum IB2d_re64_pilot.tar.gz
if [[ -e re64_flapping_source ]]; then
  echo "re64_flapping_source already exists; choose a fresh directory name"
else
  mkdir re64_flapping_source
  tar -xzf IB2d_re64_pilot.tar.gz -C re64_flapping_source
  cd /xdisk/your_pi_account/your_netid/re64_flapping_source/IB2d
  bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_flapping_pilot.sh your_pi_account
fi
```

Compare that checksum with `Get-FileHash -Algorithm SHA256` on Windows before
extracting. The archive contains no MATLAB tree, paper PDF, generated run
directory, bytecode cache, or optional C/Cython writer.

The wrapper submits array indices `0` and `1` concurrently on `standard`, then
an `afterok` analysis/render job.  It prints the exact job IDs and copy-ready
`squeue`, `sacct`, and `tail` commands.  The optional argument is the PI account
reported by `va`; the script's default is the placeholder `your_pi_account`, so
either pass the real account or edit that default once.

If a Python module name changes on Puma, override it without editing files:

```bash
PYTHON_MODULE=python/EXACT_VERSION \
  bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_flapping_pilot.sh your_pi_account
```

## Review gate

Do not launch the six-case parameter sweep until all of the following have
been reviewed:

1. The actual plate follows the prescribed plate through both reversals, with
   no visible deformation, marker crossing, or phase lag.
2. The vorticity movie is plausible and does not visibly interact with a
   periodic image during these two cycles.
3. Raw fixed-frame `Fx` and `Fy` and derived `CH`, `CV`, `CD`, and `CL` do not
   contain alternating, step-scale ringing.
4. The two timestep histories agree after phase interpolation, excluding the
   documented motion-relative reversal windows.
5. Target-error and immersed-boundary force-balance diagnostics stay finite;
   inspect them against the nonbinding review guides in the generated plot.
6. Motion-relative `CD` and `CL` midpoint cells touching each reversal are
   flagged and omitted rather than assigned an arbitrary direction.
7. The full-resolution high-frequency diagnostics for both timesteps show no
   unacceptable alternating-step ringing. No numerical threshold is imposed;
   this remains a documented human judgment.

The pilot is a timestep and implementation check only.  Passing it does not
replace the later grid-convergence and periodic-box/image study.

## Generated review artifacts

Each completed case contains its dimensional `force_history.csv`, audited
`case_metadata.json`, `COMPLETE.json`, snapshots, and a `force_analysis/`
directory. The analysis directory contains the full coefficient history,
cycle-1/2 summary, fixed-frame and motion-relative force plot, and a cycle
overlay. The pilot-level `summary/` directory contains:

- `cycle2_timestep_comparison.csv`: phase-aligned coarse, fine, and difference
  histories for `CH`, `CV`, `CD`, and `CL` under both normalizations.
- `cycle2_timestep_comparison.json`: RMS, relative-L2, maximum, and mean
  coarse/fine differences.
- `cycle2_timestep_comparison.svg`: visual comparison of the primary
  RMS-normalized horizontal and vertical coefficients.
- `pilot_analysis_manifest.json`: paths and provenance for every artifact.

Each case's `force_analysis/` directory also contains
`full_resolution_oscillation_diagnostics.json`. It is computed from every
force sample (not the 1,441-point phase table) and labels the fixed-frame
components as `CH` and `CV`.

The fine-timestep case also contains `visual_pilot/animation.html`, PNG frames,
the force/target-error diagnostic, and an MP4 when `ffmpeg` is available. No
script assigns an automatic physics pass: the comparison is explicitly marked
`human_review_required`.

## Package and download

After the dependent analysis job finishes successfully:

```bash
cd /xdisk/your_pi_account/your_netid
tar -czf plate_re64_flapping_pilot_results.tar.gz \
  --exclude='*.partial-*' \
  plate_re64_flapping_pilot
sha256sum plate_re64_flapping_pilot_results.tar.gz
```

Download from Windows PowerShell.  Keep the destination directory outside the
closing quote, with no quoted trailing backslash:

```powershell
New-Item -ItemType Directory -Force "C:\path\to\puma_results"
scp your_netid@filexfer.hpc.arizona.edu:/xdisk/your_pi_account/your_netid/plate_re64_flapping_pilot_results.tar.gz `
  "C:\path\to\puma_results"
```

