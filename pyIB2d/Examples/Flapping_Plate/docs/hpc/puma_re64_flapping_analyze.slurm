#!/bin/bash
#SBATCH --job-name=plate-re64-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
PILOT_ROOT="${FLAPPING_PILOT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re64_flapping_pilot}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
RUNNER="${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/run_re64_flapping_pilot.py"

if [[ ! -f "${RUNNER}" ]]; then
    echo "Missing flapping-pilot driver: ${RUNNER}" >&2
    exit 2
fi

mkdir -p "${PILOT_ROOT}/summary"
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

echo "Analyzing Re=64 exact-tilted flapping timestep pair"
echo "Project root: ${PROJECT_ROOT}"
echo "Pilot root: ${PILOT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version

/usr/bin/time -v python3 -u "${RUNNER}" analyze \
    --output-root "${PILOT_ROOT}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Review configuration: ${PILOT_ROOT}/pilot_configuration.json"
echo "Review artifacts: ${PILOT_ROOT}/summary"

