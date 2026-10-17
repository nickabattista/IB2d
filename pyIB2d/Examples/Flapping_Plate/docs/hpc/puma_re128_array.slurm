#!/bin/bash
#SBATCH --job-name=plate-re128
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=12:00:00
#SBATCH --array=0-4%2

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
STUDY_ROOT="${STUDY_OUTPUT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re128_study}"
STUDY_PHASE="${STUDY_PHASE:-calibration}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
PRODUCTION_DT="${PRODUCTION_DT:-2.5e-6}" # common fixed-timestep value
PRODUCTION_TARGET_STIFFNESS="${PRODUCTION_TARGET_STIFFNESS:-1e9}" # common value

CASE_SCRIPT="${PROJECT_ROOT}/pyIB2d/Examples/Translating_Plate/re128_study.py"
if [[ ! -f "${CASE_SCRIPT}" ]]; then
    echo "Missing study driver: ${CASE_SCRIPT}" >&2
    exit 2
fi
if [[ "${STUDY_PHASE}" != "calibration" && "${STUDY_PHASE}" != "production" ]]; then
    echo "STUDY_PHASE must be calibration or production" >&2
    exit 2
fi

mkdir -p "${STUDY_ROOT}/cases/${STUDY_PHASE}" "${STUDY_ROOT}/summary"
cd "${PROJECT_ROOT}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH

module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found." >&2
    module avail python
    exit 3
fi

echo "Study phase: ${STUDY_PHASE}"
echo "Case index: ${SLURM_ARRAY_TASK_ID}"
echo "Project root: ${PROJECT_ROOT}"
echo "Study root: ${STUDY_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version
python3 -c "import numpy; print('NumPy', numpy.__version__)"

arguments=(
    --phase "${STUDY_PHASE}"
    --case-index "${SLURM_ARRAY_TASK_ID}"
    --output-root "${STUDY_ROOT}"
)
if [[ "${STUDY_PHASE}" == "production" ]]; then
    arguments+=(
        --dt "${PRODUCTION_DT}"
        --target-stiffness "${PRODUCTION_TARGET_STIFFNESS}"
    )
fi

/usr/bin/time -v python3 "${CASE_SCRIPT}" "${arguments[@]}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Case output: ${STUDY_ROOT}/cases/${STUDY_PHASE}"
