#!/bin/bash
# Submit the staged Re=128 study and its dependent analysis job.

set -euo pipefail

if [[ $# -lt 2 || $# -gt 2 ]]; then
    echo "Usage: bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re128_study.sh PI_ACCOUNT calibration|production" >&2
    exit 2
fi

PI_ACCOUNT="$1"
STUDY_PHASE="$2"
PROJECT_ROOT="${IB2D_PROJECT_ROOT:-/xdisk/your_pi_account/your_netid/IB2d}"
STUDY_ROOT="${STUDY_OUTPUT_ROOT:-/xdisk/your_pi_account/your_netid/plate_re128_study}"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
PRODUCTION_DT="${PRODUCTION_DT:-2.5e-6}"
PRODUCTION_TARGET_STIFFNESS="${PRODUCTION_TARGET_STIFFNESS:-1e9}"

if [[ "${STUDY_PHASE}" != "calibration" && "${STUDY_PHASE}" != "production" ]]; then
    echo "Phase must be calibration or production." >&2
    exit 2
fi
if [[ ! -f "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_array.slurm" ]]; then
    echo "Project not found at ${PROJECT_ROOT}" >&2
    exit 2
fi

mkdir -p "${STUDY_ROOT}/logs" "${STUDY_ROOT}/summary"

module purge
module load "${PYTHON_MODULE}"
manifest_args=(
    --phase "${STUDY_PHASE}"
    --write-manifest "${STUDY_ROOT}/manifest_${STUDY_PHASE}.json"
)
if [[ "${STUDY_PHASE}" == "production" ]]; then
    manifest_args+=(
        --dt "${PRODUCTION_DT}"
        --target-stiffness "${PRODUCTION_TARGET_STIFFNESS}"
    )
fi
python3 "${PROJECT_ROOT}/pyIB2d/Examples/Translating_Plate/re128_study.py" \
    "${manifest_args[@]}"

common_export="ALL,IB2D_PROJECT_ROOT=${PROJECT_ROOT},STUDY_OUTPUT_ROOT=${STUDY_ROOT},STUDY_PHASE=${STUDY_PHASE},PYTHON_MODULE=${PYTHON_MODULE},PRODUCTION_DT=${PRODUCTION_DT},PRODUCTION_TARGET_STIFFNESS=${PRODUCTION_TARGET_STIFFNESS}"

if [[ "${STUDY_PHASE}" == "calibration" ]]; then
    raw_job_id=$(sbatch --parsable \
        --account="${PI_ACCOUNT}" --partition=standard \
        --array=0-4%2 --mem=12gb --time=24:00:00 \
        --output="${STUDY_ROOT}/logs/re128_%A_%a.out" \
        --export="${common_export}" \
        "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_array.slurm")
    job_id="${raw_job_id%%;*}"
    raw_analysis_id=$(sbatch --parsable \
        --account="${PI_ACCOUNT}" --partition=standard \
        --dependency="afterok:${job_id}" \
        --output="${STUDY_ROOT}/logs/re128_analysis_%j.out" \
        --export="${common_export}" \
        "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_analyze.slurm")
    analysis_id="${raw_analysis_id%%;*}"
    echo "Calibration array: ${job_id}"
    echo "Dependent calibration analysis: ${analysis_id}"
    exit 0
fi

calibration_report="${STUDY_ROOT}/summary/calibration_report.json"
if [[ ! -f "${calibration_report}" ]]; then
    echo "Run and analyze the calibration phase before production." >&2
    exit 4
fi
if [[ "${ALLOW_UNCALIBRATED:-0}" != "1" ]]; then
    if ! python3 -c "import json,sys; report=json.load(open(sys.argv[1])); sys.exit(0 if report['comparisons']['production_ready'] else 1)" "${calibration_report}"; then
        echo "Calibration did not pass the default force/rigidity gates." >&2
        echo "Review ${calibration_report}; set ALLOW_UNCALIBRATED=1 only for an intentional diagnostic run." >&2
        exit 5
    fi
    if ! python3 -c "import json,math,sys; r=json.load(open(sys.argv[1]))['comparisons']['provisional_production_parameters']; ok=math.isclose(float(sys.argv[2]),float(r['reference_dt']),rel_tol=1e-12) and math.isclose(float(sys.argv[3]),float(r['reference_target_stiffness']),rel_tol=1e-12); sys.exit(0 if ok else 1)" "${calibration_report}" "${PRODUCTION_DT}" "${PRODUCTION_TARGET_STIFFNESS}"; then
        echo "The requested common dt/stiffness pair was not certified by calibration." >&2
        echo "Review the pair and set ALLOW_UNCALIBRATED=1 only for an intentional diagnostic run." >&2
        exit 6
    fi
fi

raw_low_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --array=0-2%2 --mem=12gb --time=12:00:00 \
    --output="${STUDY_ROOT}/logs/re128_%A_%a.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_array.slurm")
low_id="${raw_low_id%%;*}"
raw_high_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --array=3-5%2 --mem=24gb --time=36:00:00 \
    --output="${STUDY_ROOT}/logs/re128_%A_%a.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_array.slurm")
high_id="${raw_high_id%%;*}"
raw_analysis_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" --partition=standard \
    --dependency="afterok:${low_id}:${high_id}" \
    --output="${STUDY_ROOT}/logs/re128_analysis_%j.out" \
    --export="${common_export}" \
    "${PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re128_analyze.slurm")
analysis_id="${raw_analysis_id%%;*}"

echo "Production low-resolution array: ${low_id}"
echo "Production high-resolution array: ${high_id}"
echo "Dependent production analysis: ${analysis_id}"
