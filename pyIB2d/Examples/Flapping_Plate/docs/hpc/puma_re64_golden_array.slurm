#!/bin/bash
#SBATCH --job-name=plate-re64-golden
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-0%2
#SBATCH --mem=12gb
#SBATCH --time=24:00:00

set -euo pipefail

: "${IB2D_PROJECT_ROOT:?IB2D_PROJECT_ROOT was not exported by the submit script}"
: "${GOLDEN_OUTPUT_ROOT:?GOLDEN_OUTPUT_ROOT was not exported by the submit script}"
: "${GOLDEN_MANIFEST:?GOLDEN_MANIFEST was not exported by the submit script}"
: "${GOLDEN_STAGE:?GOLDEN_STAGE was not exported by the submit script}"
: "${SLURM_ARRAY_TASK_ID:?This script must be run as a Slurm array job}"

PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
RUNNER="${IB2D_PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/run_re64_golden_study.py"

case "${GOLDEN_STAGE}" in
    stiffness|timestep|box|grid|startup|production) ;;
    *)
        echo "Unknown GOLDEN_STAGE: ${GOLDEN_STAGE}" >&2
        exit 2
        ;;
esac

if [[ ! -f "${RUNNER}" ]]; then
    echo "Missing golden-study driver: ${RUNNER}" >&2
    exit 2
fi
if [[ ! -f "${GOLDEN_MANIFEST}" ]]; then
    echo "Missing immutable study manifest: ${GOLDEN_MANIFEST}" >&2
    exit 2
fi

mkdir -p "${GOLDEN_OUTPUT_ROOT}/logs" "${GOLDEN_OUTPUT_ROOT}/summary"
cd "${IB2D_PROJECT_ROOT}"

# The solver is single-process.  Prevent linked numerical libraries from
# silently starting extra threads on a shared node.
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

case_count=$(python3 "${RUNNER}" count --manifest "${GOLDEN_MANIFEST}")
if ! [[ "${case_count}" =~ ^[1-9][0-9]*$ ]]; then
    echo "The runner returned an invalid case count: ${case_count}" >&2
    exit 3
fi
if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= case_count )); then
    echo "Case index ${SLURM_ARRAY_TASK_ID} is outside 0..$((case_count - 1))." >&2
    exit 3
fi

echo "Re=64 flapping-plate golden study"
echo "Stage: ${GOLDEN_STAGE}"
echo "Case index: ${SLURM_ARRAY_TASK_ID} of $((case_count - 1))"
echo "Manifest: ${GOLDEN_MANIFEST}"
echo "Project root: ${IB2D_PROJECT_ROOT}"
echo "Output root: ${GOLDEN_OUTPUT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version
python3 -c "import numpy; print('NumPy', numpy.__version__)"

# The Python driver owns case-level COMPLETE/FAILED markers and treats a
# resubmission of an already complete case as a safe no-op.
/usr/bin/time -v python3 -u "${RUNNER}" run \
    --manifest "${GOLDEN_MANIFEST}" \
    --case-index "${SLURM_ARRAY_TASK_ID}" \
    --output-root "${GOLDEN_OUTPUT_ROOT}"

echo "Finished: $(date --iso-8601=seconds)"

