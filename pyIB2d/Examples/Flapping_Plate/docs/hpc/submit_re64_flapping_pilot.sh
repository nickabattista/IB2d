#!/bin/bash
# Submit the exact tilted Re=64 flapping pilot and its dependent analysis.

set -euo pipefail

if [[ $# -gt 1 ]]; then
    echo "Usage: bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_flapping_pilot.sh [PI_ACCOUNT]" >&2
    exit 2
fi

PI_ACCOUNT="${1:-your_pi_account}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../../../.." && pwd)"
PROJECT_ROOT="${IB2D_PROJECT_ROOT:-${DEFAULT_PROJECT_ROOT}}"
PILOT_ROOT="${FLAPPING_PILOT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re64_flapping_pilot}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
RUNNER="${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/run_re64_flapping_pilot.py"

if [[ ! -f "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_flapping_pilot.slurm" ]]; then
    echo "Project not found at ${PROJECT_ROOT}" >&2
    exit 2
fi
if [[ ! -f "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_flapping_analyze.slurm" ]]; then
    echo "Missing dependent analysis job under ${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc" >&2
    exit 2
fi
if [[ ! -f "${RUNNER}" ]]; then
    echo "Missing flapping-pilot driver: ${RUNNER}" >&2
    exit 2
fi

mkdir -p "${PILOT_ROOT}/logs" "${PILOT_ROOT}/summary"

# Record the exact requested cases before submitting expensive work.  This also
# fails early when the selected Python module cannot import the pilot driver.
module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found." >&2
    module avail python
    exit 3
fi
python3 "${RUNNER}" show > "${PILOT_ROOT}/pilot_configuration.json"

common_export="ALL,IB2D_PROJECT_ROOT=${PROJECT_ROOT},FLAPPING_PILOT_ROOT=${PILOT_ROOT},PYTHON_MODULE=${PYTHON_MODULE}"
raw_array_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --array=0-1%2 --mem=8gb --time=24:00:00 \
    --output="${PILOT_ROOT}/logs/re64_flapping_%A_%a.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_flapping_pilot.slurm")
array_id="${raw_array_id%%;*}"

raw_analysis_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --dependency="afterok:${array_id}" \
    --mem=4gb --time=01:00:00 \
    --output="${PILOT_ROOT}/logs/re64_flapping_analysis_%j.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_flapping_analyze.slurm")
analysis_id="${raw_analysis_id%%;*}"

echo "Re=64 timestep-pair array: ${array_id}"
echo "Dependent analysis/render job: ${analysis_id}"
echo "Array logs: ${PILOT_ROOT}/logs/re64_flapping_${array_id}_{0,1}.out"
echo "Analysis log: ${PILOT_ROOT}/logs/re64_flapping_analysis_${analysis_id}.out"
echo
echo "Copy-ready monitoring commands:"
echo "squeue -j ${array_id},${analysis_id} -o \"%.18i %.2t %.10M %.10l %R\""
echo "sacct -j ${array_id},${analysis_id} --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode"
echo "tail -f ${PILOT_ROOT}/logs/re64_flapping_${array_id}_1.out"

