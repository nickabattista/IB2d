#!/bin/bash
# Create an immutable manifest, submit one Re=64 golden-study stage as a
# two-at-a-time array, and submit a dependent partial-safe analysis job.

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  bash pyIB2d/Examples/Flapping_Plate/docs/hpc/submit_re64_golden_stage.sh PI_ACCOUNT STAGE STUDY_ROOT [MANIFEST_OPTIONS...]

STAGE is one of: stiffness, timestep, box, grid, startup, production

STUDY_ROOT must be an explicit absolute path.  This wrapper keeps each stage
under STUDY_ROOT/STAGE so immutable manifests and results cannot collide.
Extra arguments are passed to the Python driver's `manifest` command.  Common
examples are:
  --target-stiffness 5e8
  --dt 6.25e-7
  --box-chords 16
  --cells-per-chord 32
  --cycles 4
  --stiffness-values 2.5e8,5e8,1e9
  --dt-values 1.25e-6,6.25e-7
  --box-values 12,16,24
  --grid-values 16,24,32

Optional environment overrides:
  PYTHON_MODULE, GOLDEN_PARTITION, GOLDEN_ARRAY_CONCURRENCY,
  GOLDEN_MEMORY, GOLDEN_WALLTIME, GOLDEN_ANALYSIS_MEMORY,
  GOLDEN_ANALYSIS_WALLTIME
EOF
    exit 2
}

if [[ $# -lt 3 ]]; then
    usage
fi

PI_ACCOUNT="$1"
GOLDEN_STAGE="$2"
GOLDEN_STUDY_ROOT="$3"
shift 3
manifest_options=("$@")

case "${GOLDEN_STAGE}" in
    stiffness|timestep|box|grid|startup|production) ;;
    *)
        echo "STAGE must be stiffness, timestep, box, grid, startup, or production." >&2
        usage
        ;;
esac

if [[ -z "${PI_ACCOUNT}" ]]; then
    echo "PI_ACCOUNT cannot be empty." >&2
    exit 2
fi
if [[ "${GOLDEN_STUDY_ROOT}" != /* || "${GOLDEN_STUDY_ROOT}" == "/" ]]; then
    echo "STUDY_ROOT must be an explicit absolute path other than /." >&2
    exit 2
fi
if [[ "${GOLDEN_STUDY_ROOT}" == *','* || "${GOLDEN_STUDY_ROOT}" == *$'\n'* ]]; then
    echo "STUDY_ROOT cannot contain a comma or newline because Slurm exports it." >&2
    exit 2
fi

# These arguments are controlled by this wrapper.  Refusing duplicates avoids
# a manifest whose filename and internal stage disagree.
for argument in "${manifest_options[@]}"; do
    case "${argument}" in
        --stage|--stage=*|--output|--output=*)
            echo "Do not pass ${argument}; STAGE and manifest output are positional." >&2
            exit 2
            ;;
    esac
done

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
IB2D_PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../../../.." && pwd -P)"
RUNNER="${IB2D_PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/run_re64_golden_study.py"
PYTHON_MODULE="${PYTHON_MODULE:-python/3.11/3.11.4}"
GOLDEN_PARTITION="${GOLDEN_PARTITION:-standard}"
GOLDEN_ARRAY_CONCURRENCY="${GOLDEN_ARRAY_CONCURRENCY:-2}"
GOLDEN_ANALYSIS_MEMORY="${GOLDEN_ANALYSIS_MEMORY:-8gb}"
GOLDEN_ANALYSIS_WALLTIME="${GOLDEN_ANALYSIS_WALLTIME:-02:00:00}"

if [[ ! -f "${RUNNER}" ]]; then
    echo "Missing golden-study driver: ${RUNNER}" >&2
    exit 2
fi
if [[ "${IB2D_PROJECT_ROOT}" == *','* || "${IB2D_PROJECT_ROOT}" == *$'\n'* ]]; then
    echo "The project path cannot contain a comma or newline because Slurm exports it." >&2
    exit 2
fi
if ! [[ "${GOLDEN_ARRAY_CONCURRENCY}" =~ ^[1-9][0-9]*$ ]]; then
    echo "GOLDEN_ARRAY_CONCURRENCY must be a positive integer." >&2
    exit 2
fi

# Defaults are intentionally stage-specific.  They can be raised or lowered
# without editing this file by exporting GOLDEN_MEMORY/GOLDEN_WALLTIME.
case "${GOLDEN_STAGE}" in
    stiffness|timestep)
        default_memory="8gb"
        default_walltime="24:00:00"
        ;;
    box|grid)
        default_memory="12gb"
        default_walltime="72:00:00"
        ;;
    startup|production)
        default_memory="8gb"
        default_walltime="72:00:00"
        ;;
esac
GOLDEN_MEMORY="${GOLDEN_MEMORY:-${default_memory}}"
GOLDEN_WALLTIME="${GOLDEN_WALLTIME:-${default_walltime}}"

mkdir -p "${GOLDEN_STUDY_ROOT}"
GOLDEN_STUDY_ROOT="$(cd -- "${GOLDEN_STUDY_ROOT}" && pwd -P)"
GOLDEN_OUTPUT_ROOT="${GOLDEN_STUDY_ROOT}/${GOLDEN_STAGE}"
mkdir -p \
    "${GOLDEN_STUDY_ROOT}/manifests" \
    "${GOLDEN_OUTPUT_ROOT}/logs" \
    "${GOLDEN_OUTPUT_ROOT}/summary"

module purge
if ! module load "${PYTHON_MODULE}"; then
    echo "Python module ${PYTHON_MODULE} was not found." >&2
    module avail python
    exit 3
fi

GOLDEN_MANIFEST="${GOLDEN_STUDY_ROOT}/manifests/${GOLDEN_STAGE}.json"
candidate_dir=$(mktemp -d "${GOLDEN_STUDY_ROOT}/manifests/.${GOLDEN_STAGE}.XXXXXX")
candidate_manifest="${candidate_dir}/manifest.json"
cleanup_candidate() {
    if [[ -f "${candidate_manifest}" ]]; then
        rm -f -- "${candidate_manifest}"
    fi
    if [[ -d "${candidate_dir}" ]]; then
        rmdir -- "${candidate_dir}" 2>/dev/null || true
    fi
}
trap cleanup_candidate EXIT

python3 "${RUNNER}" manifest \
    --stage "${GOLDEN_STAGE}" \
    --output "${candidate_manifest}" \
    "${manifest_options[@]}"

if [[ -f "${GOLDEN_MANIFEST}" ]]; then
    if ! cmp -s -- "${GOLDEN_MANIFEST}" "${candidate_manifest}"; then
        echo "A different immutable manifest already exists:" >&2
        echo "  ${GOLDEN_MANIFEST}" >&2
        echo "Use the same choices, or choose a new STUDY_ROOT for a new study." >&2
        exit 4
    fi
    echo "Reusing matching manifest: ${GOLDEN_MANIFEST}"
else
    mv -- "${candidate_manifest}" "${GOLDEN_MANIFEST}"
    echo "Created manifest: ${GOLDEN_MANIFEST}"
fi

case_count=$(python3 "${RUNNER}" count --manifest "${GOLDEN_MANIFEST}")
if ! [[ "${case_count}" =~ ^[1-9][0-9]*$ ]]; then
    echo "The runner returned an invalid case count: ${case_count}" >&2
    exit 3
fi
last_case=$((case_count - 1))

configuration_file="${GOLDEN_STUDY_ROOT}/manifests/${GOLDEN_STAGE}_configuration.txt"
python3 "${RUNNER}" show --manifest "${GOLDEN_MANIFEST}" > "${configuration_file}"
cat "${configuration_file}"

common_export="ALL,IB2D_PROJECT_ROOT=${IB2D_PROJECT_ROOT},GOLDEN_OUTPUT_ROOT=${GOLDEN_OUTPUT_ROOT},GOLDEN_MANIFEST=${GOLDEN_MANIFEST},GOLDEN_STAGE=${GOLDEN_STAGE},PYTHON_MODULE=${PYTHON_MODULE}"

raw_array_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" \
    --partition="${GOLDEN_PARTITION}" \
    --job-name="re64-${GOLDEN_STAGE}" \
    --array="0-${last_case}%${GOLDEN_ARRAY_CONCURRENCY}" \
    --mem="${GOLDEN_MEMORY}" \
    --time="${GOLDEN_WALLTIME}" \
    --output="${GOLDEN_OUTPUT_ROOT}/logs/${GOLDEN_STAGE}_%A_%a.out" \
    --export="${common_export}" \
    "${IB2D_PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_golden_array.slurm")
array_id="${raw_array_id%%;*}"

# afterany is deliberate: a failed case should still produce a partial report
# identifying what needs to be rerun.  The analyzer must not certify a partial
# stage as complete.
raw_analysis_id=$(sbatch --parsable \
    --account="${PI_ACCOUNT}" \
    --partition="${GOLDEN_PARTITION}" \
    --job-name="re64-${GOLDEN_STAGE}-analysis" \
    --dependency="afterany:${array_id}" \
    --mem="${GOLDEN_ANALYSIS_MEMORY}" \
    --time="${GOLDEN_ANALYSIS_WALLTIME}" \
    --output="${GOLDEN_OUTPUT_ROOT}/logs/${GOLDEN_STAGE}_analysis_%j.out" \
    --export="${common_export}" \
    "${IB2D_PROJECT_ROOT}/pyIB2d/Examples/Flapping_Plate/docs/hpc/puma_re64_golden_analyze.slurm")
analysis_id="${raw_analysis_id%%;*}"

echo
echo "Submitted Re=64 golden-study stage: ${GOLDEN_STAGE}"
echo "Case array (${case_count} cases, at most ${GOLDEN_ARRAY_CONCURRENCY} running): ${array_id}"
echo "Dependent partial-safe analysis: ${analysis_id}"
echo "Study root: ${GOLDEN_STUDY_ROOT}"
echo "Stage output root: ${GOLDEN_OUTPUT_ROOT}"
echo "Manifest: ${GOLDEN_MANIFEST}"
echo "Configuration: ${configuration_file}"
echo "Array logs: ${GOLDEN_OUTPUT_ROOT}/logs/${GOLDEN_STAGE}_${array_id}_<case>.out"
echo "Analysis log: ${GOLDEN_OUTPUT_ROOT}/logs/${GOLDEN_STAGE}_analysis_${analysis_id}.out"
echo
echo "Copy-ready monitoring commands:"
echo "squeue -j ${array_id},${analysis_id} -o \"%.18i %.2t %.10M %.10l %R\""
echo "sacct -j ${array_id},${analysis_id} --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode"
echo "tail -f ${GOLDEN_OUTPUT_ROOT}/logs/${GOLDEN_STAGE}_${array_id}_0.out"
echo
echo "Safe resubmission after this array has stopped: repeat the same command."
echo "Matching manifests are reused and cases already marked COMPLETE are skipped."

