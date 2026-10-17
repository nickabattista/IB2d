#!/bin/bash
#SBATCH --job-name=plate-re128-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
STUDY_ROOT="${STUDY_OUTPUT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re128_study}"
STUDY_PHASE="${STUDY_PHASE:-calibration}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
PRODUCTION_DT="${PRODUCTION_DT:-2.5e-6}"
PRODUCTION_TARGET_STIFFNESS="${PRODUCTION_TARGET_STIFFNESS:-1e9}"

cd "${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH

module purge
module load "${PYTHON_MODULE}"

arguments=(
    --phase "${STUDY_PHASE}"
    --output-root "${STUDY_ROOT}"
)
if [[ "${STUDY_PHASE}" == "production" ]]; then
    arguments+=(
        --dt "${PRODUCTION_DT}"
        --target-stiffness "${PRODUCTION_TARGET_STIFFNESS}"
    )
fi

python3 pyIB2d/Examples/Translating_Plate/analyze_re128_study.py "${arguments[@]}"
echo "Reports written under ${STUDY_ROOT}/summary"
