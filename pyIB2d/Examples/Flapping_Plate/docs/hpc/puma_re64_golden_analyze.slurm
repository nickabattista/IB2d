#!/bin/bash
#SBATCH --job-name=plate-re64-golden-analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --time=02:00:00

set -euo pipefail

: "${IB2D_PROJECT_ROOT:?IB2D_PROJECT_ROOT was not exported by the submit script}"
: "${GOLDEN_OUTPUT_ROOT:?GOLDEN_OUTPUT_ROOT was not exported by the submit script}"
: "${GOLDEN_MANIFEST:?GOLDEN_MANIFEST was not exported by the submit script}"
: "${GOLDEN_STAGE:?GOLDEN_STAGE was not exported by the submit script}"

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

mkdir -p "${GOLDEN_OUTPUT_ROOT}/summary"
cd "${IB2D_PROJECT_ROOT}"

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

echo "Analyzing Re=64 flapping-plate golden study"
echo "Stage: ${GOLDEN_STAGE}"
echo "Manifest: ${GOLDEN_MANIFEST}"
echo "Project root: ${IB2D_PROJECT_ROOT}"
echo "Output root: ${GOLDEN_OUTPUT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
python3 --version

# This job is submitted with afterany.  The analyzer therefore runs even if a
# case failed, records which cases are incomplete, and must not label a partial
# stage as complete.
/usr/bin/time -v python3 -u "${RUNNER}" analyze \
    --manifest "${GOLDEN_MANIFEST}" \
    --output-root "${GOLDEN_OUTPUT_ROOT}"

echo "Finished: $(date --iso-8601=seconds)"
echo "Review stage artifacts under: ${GOLDEN_OUTPUT_ROOT}/summary"

