#!/usr/bin/bash -l

######################################################
# Environment setup: Singularity images + conda envs
#
# Downloads all Singularity images required by the
# MAG pipeline and creates the necessary conda
# environments from env/*.yml.
#
# Run from the project root directory before executing
# any pipeline scripts. Requires ~200+ GB free disk space.
#
# Override CONDA/MAMBA paths if your miniforge3 is
# installed elsewhere:
#   CONDA=/path/to/conda MAMBA=/path/to/mamba ./scripts/set_up_environment.sh
######################################################

set -Eeuo pipefail

log() { echo "[$(date +'%F %T')] $*"; }
trap 'rc=$?; log "[ERROR] line ${LINENO}: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    log "SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

module load singularityce

CONDA="${CONDA:-/home/jduque2/miniforge3/bin/conda}"
MAMBA="${MAMBA:-/home/jduque2/miniforge3/bin/mamba}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/../env"

[[ -x "${CONDA}" ]]  || { log "[ERROR] conda not found at ${CONDA}. Set CONDA= to override."; exit 1; }
[[ -x "${MAMBA}" ]]  || { log "[ERROR] mamba not found at ${MAMBA}. Set MAMBA= to override."; exit 1; }
[[ -d "${ENV_DIR}" ]] || { log "[ERROR] env/ directory not found at ${ENV_DIR}"; exit 1; }

# ── Singularity images ──────────────────────────────────────────────────────
# NOTE: verify these Docker sources before pulling on a new cluster.
# bbtools_latest and metabat2 tags may lag upstream; check DockerHub first.

log "Pulling Singularity images..."

declare -A IMAGES=(
    ["bbtools_38.86.sif"]="docker://bryce911/bbtools:38.86"
    ["bbtools_latest.sif"]="docker://bryce911/bbtools:latest"
    ["spades_3.15.2.sif"]="docker://bryce911/spades:3.15.2"
    ["megahit.sif"]="docker://vout/megahit"
    ["metabat2.sif"]="docker://metabat/metabat2"
)

for img in "${!IMAGES[@]}"; do
    if [[ -f "${img}" ]]; then
        log "[SKIP] ${img} already present"
    else
        log "[PULL] ${img}"
        singularity pull "${img}" "${IMAGES[${img}]}"
    fi
done

# ── RQC filter reference data ───────────────────────────────────────────────

RQC_DIR="${PWD}/RQCFilterData"
if [[ -d "${RQC_DIR}" ]]; then
    log "[SKIP] RQCFilterData already present at ${RQC_DIR}"
else
    log "[DOWNLOAD] RQCFilterData (~8 GB) from NERSC portal"
    mkdir -p "${RQC_DIR}"
    wget -O - http://portal.nersc.gov/dna/metagenome/assembly/rqcfilter/RQCFilterData.tar \
        | tar -xf - -C "${RQC_DIR}"
fi

# ── Conda environments ──────────────────────────────────────────────────────

log "Creating conda environments from ${ENV_DIR}/*.yml ..."

for yml in "${ENV_DIR}"/*.yml; do
    env_name="$(grep '^name:' "${yml}" | awk '{print $2}')"
    if "${CONDA}" env list | awk '{print $1}' | grep -qx "${env_name}"; then
        log "[SKIP] conda env '${env_name}' already exists"
    else
        log "[CREATE] conda env '${env_name}' from $(basename "${yml}")"
        "${MAMBA}" env create -f "${yml}"
    fi
done

log "Environment setup complete."
log "Singularity images present in: ${PWD}"
log "Conda environments created from: ${ENV_DIR}"
log "RQCFilterData at: ${RQC_DIR}"
