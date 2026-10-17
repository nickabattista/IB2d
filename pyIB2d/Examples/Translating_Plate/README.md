# Prescribed translating plate (headless pyIB2d example)

A rigid plate driven through the periodic box by a smoothly ramped prescribed
translation, held with target points.  This is the older regression and
provenance case for the `Flapping_Plate` example: it was used for the
MATLAB/Python language-parity check and for the Re=128 study described in
`../Flapping_Plate/docs/validation/RE128_STUDY.md`.

From this folder:

```text
python main2d.py --regenerate --nx 32 --steps 8
python -m unittest discover -s tests -v
```

`main2d.py` generates the case files (`input2d`, `plate.vertex`, `plate.target`,
`case_metadata.json`) into `runs/smoke` with `generate_case.py` and runs the
solver headlessly, recording the body load every step in `force_history.csv`
through the driver's `force_callback`.  `runs/` is generated output and is
git-ignored.  `tests/test_spreading.py` is the regression test for the local
four-point force spreading in `IBM_Blackbox`.  See `../Flapping_Plate/README.md`
for the shared background and the list of changes to `IBM_Blackbox`.
