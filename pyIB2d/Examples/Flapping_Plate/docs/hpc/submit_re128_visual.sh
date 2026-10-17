#!/bin/bash
# Submit the single visual/force-history checkpoint before large study arrays.

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_visual.sh PI_ACCOUNT" >&2
    exit 2
fi

PI_ACCOUNT="$1"
PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
PILOT_ROOT="${VISUAL_PILOT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re128_visual_pilot}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"

if [[ ! -f "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_visual.slurm" ]]; then
    echo "Project not found at ${PROJECT_ROOT}" >&2
    exit 2
fi

mkdir -p "${PILOT_ROOT}/logs"
module purge
module load "${PYTHON_MODULE}"
python3 "${PROJECT_ROOT}/pyIB2d/Examples/Translating_Plate/run_re128_visual_pilot.py" \
    show > "${PILOT_ROOT}/pilot_configuration.json"

common_export="ALL,IB2D_PROJECT_ROOT=${PROJECT_ROOT},VISUAL_PILOT_ROOT=${PILOT_ROOT},PYTHON_MODULE=${PYTHON_MODULE}"
raw_job_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --output="${PILOT_ROOT}/logs/re128_visual_%j.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_visual.slurm")
job_id="${raw_job_id%%;*}"

echo "Visual pilot job: ${job_id}"
echo "Log: ${PILOT_ROOT}/logs/re128_visual_${job_id}.out"

