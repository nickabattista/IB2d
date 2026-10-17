#!/bin/bash
#SBATCH --job-name=plate-re128-visual
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=24:00:00

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
PILOT_ROOT="${VISUAL_PILOT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re128_visual_pilot}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
PILOT_SCRIPT="${PROJECT_ROOT}/pyIB2d/Examples/Translating_Plate/run_re128_visual_pilot.py"

if [[ ! -f "${PILOT_SCRIPT}" ]]; then
    echo "Missing visual-pilot driver: ${PILOT_SCRIPT}" >&2
    exit 2
fi

mkdir -p "${PILOT_ROOT}" "${PILOT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
unset PYTHONPATH

module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found." >&2
    module avail python
    exit 3
fi

echo "Re=128 visual pilot"
echo "Project root: ${PROJECT_ROOT}"
echo "Pilot root: ${PILOT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version
python3 -c "import numpy; print('NumPy', numpy.__version__)"

/usr/bin/time -v python3 -u "${PILOT_SCRIPT}" all --output-root "${PILOT_ROOT}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Open after download: ${PILOT_ROOT}/visual_L16_P64_N1024/visual_pilot/animation.html"

