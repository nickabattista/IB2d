#!/bin/bash
#SBATCH --job-name=plate-re64-flap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-1%2
#SBATCH --mem=8gb
#SBATCH --time=24:00:00

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
PILOT_ROOT="${FLAPPING_PILOT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re64_flapping_pilot}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
RUNNER="${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/run_re64_flapping_pilot.py"

if [[ ! -f "${RUNNER}" ]]; then
    echo "Missing flapping-pilot driver: ${RUNNER}" >&2
    exit 2
fi
if [[ "${SLURM_ARRAY_TASK_ID}" != "0" && "${SLURM_ARRAY_TASK_ID}" != "1" ]]; then
    echo "SLURM_ARRAY_TASK_ID must be 0 or 1" >&2
    exit 2
fi

mkdir -p "${PILOT_ROOT}/cases" "${PILOT_ROOT}/summary"
cd "${PROJECT_ROOT}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg
unset PYTHONPATH

module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found." >&2
    module avail python
    exit 3
fi

echo "Re=64 exact-tilted flapping timestep pilot"
echo "Case index: ${SLURM_ARRAY_TASK_ID}"
echo "Project root: ${PROJECT_ROOT}"
echo "Pilot root: ${PILOT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version
python3 -c "import numpy; print('NumPy', numpy.__version__)"

/usr/bin/time -v python3 -u "${RUNNER}" run \
    --case-index "${SLURM_ARRAY_TASK_ID}" \
    --output-root "${PILOT_ROOT}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Case root: ${PILOT_ROOT}/cases"

