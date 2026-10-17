# Validation status

## Completed locally

- The original IB2d archive remains unchanged.
- A separate Python and MATLAB validation workspace has been created.
- Python now runs headlessly with NumPy as its only third-party dependency.
- Target motion is absolute and evaluated at the force midpoint.
- Signed body reaction is logged before external fluid forcing.
- Eulerian and Lagrangian force integrals are recorded independently.
- Python force spreading uses a local four-point stencil.  A regression test
  verifies agreement with the legacy dense product.
- All seventeen translating-case tests pass, including a real two-step solver
  run and the dependency-free visual artifact pipeline.
- The copied Python tree and comparison utilities compile cleanly.

The prescribed flapping implementation is also complete locally:

- The `Re=64` tilted motion uses analytical stroke and pitch, a midchord pivot,
  absolute midpoint target updates, and no accumulated kinematic drift.
- Generated cases enforce `ds=dx/2`, retain `kT=1e9`, and derive
  `mu=rho*Umax*c/Re` from the recorded peak-speed convention.
- Signed dimensional `Fx` and `Fy`, force balance, and target errors are saved
  every step. Post-processing reports `CH`, `CV`, `CD`, and `CL` with both
  `Umax` and `Urms` normalizations and explicit reversal handling.
- All twenty-one flapping-case tests pass. They include analytical kinematics,
  a real two-step IBM solve with flapping-aware snapshots, reversal-safe
  force comparisons, source-provenance guards, and a synthetic coarse/fine
  phase-comparison test.
- The exact Puma manifest contains 52,360 coarse-timestep steps and 104,720
  fine-timestep steps. The cycle-aligned period is `0.032725`, so both runs end
  at `t=0.06545` after exactly two cycles while retaining exact `Re=64`.

The staged golden-study workflow is now implemented locally as well:

- deterministic manifests cover stiffness, timestep, box-size, grid,
  startup, and six-case production stages;
- task-facing output uses `L=Fy`, `D=-Fx`, and the specified `U_ref`
  normalization rather than the supplemental motion-relative convention;
- the optional periodic mean-momentum constraint is explicit in every case
  manifest and leaves the saved unbalanced plate reaction unchanged;
- completed cases receive fixed-axis histories, force plots, cycle summaries,
  tracking and oscillation diagnostics, and cross-case comparisons;
- steady averages are emitted only after two successive cycle pairs meet the
  documented operational repeatability rule; and
- one production case saves the final-cycle flow snapshots used for an HTML
  animation and optional MP4.

All thirty-one flapping tests pass, including a real two-step solve, the
periodic-force-balance helper, golden manifest invariants, fixed force signs,
steady-cycle selection, and prompt-facing plot labels. A separate inexpensive
end-to-end workflow smoke test also completed locally. That smoke result used a
deliberately under-resolved setup and is not scientific evidence.

Local Python timings on this machine:

| Case | Steps | Wall time |
| --- | ---: | ---: |
| 128 x 128 | 20 | 0.35 s |
| 256 x 256 | 20 | 0.91 s |
| 128 x 128 | 1000 | 5.46 s |
| 256 x 256 | 1000 | 38.02 s |

For the 1000-step legacy-stiffness case (`k_T=5e7`), force conservation was
near machine precision.  Maximum target lag was approximately `0.49 dx` at
128 and `0.89 dx` at 256.  This demonstrates that language parity and penalty
rigidity must be tested separately.  Increasing stiffness without reducing the
timestep can become unstable; stiffness and timestep will be refined together.

## Completed on Puma

Puma job `24205603` ran the 128 and 256 cases through `t=0.025` in both
languages.  The complete MATLAB/Python histories agree to about `1e-11`
normalized or better, so the independent implementations now provide a strong
language-parity check for the translating benchmark.

That result is not a grid-convergence result.  The 128 and 256 cases used
slightly different effective chords, the nominal parameters give
`Re=128.342...`, and the target penalty lag is not negligible.  Python force
histories also change substantially between 128 and 256.

The later Re=128 translating visual pilot also finished on Puma and was
downloaded; visual inspection was reported as good. That observation is a
useful implementation check, but it is not evidence for the new periodic
flapping motion or for spatial convergence.

## Next numerical gate

Submit only the stiffness stage documented in
`validation/RE64_GOLDEN_STUDY.md`, review its target tracking and force plateau,
and then pass the selected value to the timestep stage. Grid, box, startup,
and production remain unrun. Do not describe the local smoke test or generated
manifests as completed scientific simulations.
