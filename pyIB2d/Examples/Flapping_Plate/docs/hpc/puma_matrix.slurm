#!/bin/bash
#SBATCH --job-name=ib2d-matrix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=puma_matrix_%j.out

set -euo pipefail

PROJECT_ROOT="${IB2D_PROJECT_ROOT:-${SLURM_SUBMIT_DIR}}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
MATLAB_MODULE="${MATLAB_MODULE:-matlab}"

if [[ ! -f "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/validation/compare_histories.py" ]]; then
    echo "Submit this job from the IB2d project root." >&2
    exit 2
fi

cd "${PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PYTHONNOUSERSITE=1
unset PYTHONPATH

echo "Project root: ${PROJECT_ROOT}"
echo "Node: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"

module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found. Available modules:" >&2
    module avail python
    exit 3
fi
module list
python3 --version
python3 -c "import numpy; print('NumPy', numpy.__version__)"
python3 pyIB2d/Examples/Translating_Plate/run_matrix.py \
    --resolutions 128 256 \
    --dt 2.5e-5 \
    --final-time 0.025 \
    --target-stiffness 5e7 \
    --output-root pyIB2d/Examples/Translating_Plate/runs/matrix

module purge
if ! module load "${MATLAB_MODULE}"; then
    echo "MATLAB module ${MATLAB_MODULE} was not found. Available modules:" >&2
    module avail matlab
    exit 4
fi
module list
MATLAB_TMP_ROOT="${SLURM_JOBTMP:-/tmp}"
export MATLAB_PREFDIR="${MATLAB_TMP_ROOT}/matlab-prefs-${SLURM_JOB_ID}"
mkdir -p "${MATLAB_PREFDIR}"
matlab -singleCompThread -batch "source_dir=fullfile('${PROJECT_ROOT}','matlab','Examples','Translating_Plate'); addpath(source_dir); output_root=fullfile(source_dir,'runs','matrix'); run_matrix([128 256],2.5e-5,0.025,output_root,5e7);"

module purge
module load "${PYTHON_MODULE}"
for resolution in 128 256; do
    python3 pyIB2d/Examples/Flapping_Plate/docs/validation/compare_histories.py \
        "pyIB2d/Examples/Translating_Plate/runs/matrix/nx${resolution}/force_history.csv" \
        "matlab/Examples/Translating_Plate/runs/matrix/nx${resolution}/force_history.csv" \
        --output "pyIB2d/Examples/Flapping_Plate/docs/validation/matrix_nx${resolution}_comparison.json"
done

echo "Finished: $(date --iso-8601=seconds)"
echo "Reports: ${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/validation/matrix_nx*_comparison.json"
