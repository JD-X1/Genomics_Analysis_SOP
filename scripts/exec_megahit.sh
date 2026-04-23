#!/usr/bin/bash -l

module load singularityce


# MEGAHIT metagenomic assembly script. So arguably you could
# get away with just running exec_assembly.sh, but if you
# want to be thorough and guage the quality of resulting assemblies,
# its best to run exec_megahit.sh and exec_

# Usage: ./exec_megahit.sh <input_dir> <output_dir>
# Input: interleaved filtered fastq files

set -Eeuo pipefail
shopt -s nullglob

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    echo "[INFO] SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

log() {
    echo "[$(date +'%F %T')] $*"
}

trap 'rc=$?; log "[ERROR] line ${LINENO}: command failed: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR

if [[ $# -ne 2 ]]; then
  log "[ERROR] Usage: $0 <input_dir> <output_dir>"
  exit 2
fi

INPUT_DIR="$(realpath -m "$1")"
OUT_DIR="$(realpath -m "$2")"

MEGAHIT_THREADS="${MEGAHIT_THREADS:-40}"
MEGAHIT_MEM="${MEGAHIT_MEM:-858993459200}"
MEGAHIT_IMAGE="${MEGAHIT_IMAGE:-megahit.sif}"

if [[ ! -f "${MEGAHIT_IMAGE}" ]]; then
    log "[ERROR] MEGAHIT singularity image not found at ${MEGAHIT_IMAGE}."
    exit 1
fi

FILTERED_FASTQS=( ${INPUT_DIR}/*_interleaved_filtered.fastq.gz )
mkdir -p "${OUT_DIR}"

if [[ ${#FILTERED_FASTQS[@]} -eq 0 ]]; then
  log "[ERROR] No interleaved filtered files found in ${INPUT_DIR}/"
  exit 1
fi

GLOBAL_TMP_DIR="$(mktemp -d "${PWD}/megahit_tmp.XXXXXXXXXX")"
trap 'rc=$?; rm -rf "${GLOBAL_TMP_DIR}"; log "[ERROR] line ${LINENO}: command failed: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR
trap 'rm -rf "${GLOBAL_TMP_DIR}"' EXIT

for i in "${FILTERED_FASTQS[@]}"; do
    BASE="$(basename "${i}" _interleaved_filtered.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}/megahit_output"
    SAMPLE_OUT_SENS="${SAMPLE_OUT_DIR}/sensitive"
    SAMPLE_OUT_LARGE="${SAMPLE_OUT_DIR}/large"

    RUN_SENS=true
    RUN_LARGE=true

    if [[ -f "${SAMPLE_OUT_SENS}/final.contigs.fa" ]]; then
        log "[SKIP] megahit sensitive already complete for ${BASE}."
        RUN_SENS=false
    elif [[ -d "${SAMPLE_OUT_SENS}" ]]; then
        log "[WARN] Incomplete sensitive run detected for ${BASE}, clearing output dir."
        rm -rf "${SAMPLE_OUT_SENS}"
    fi

    if [[ -f "${SAMPLE_OUT_LARGE}/final.contigs.fa" ]]; then
        log "[SKIP] megahit large already complete for ${BASE}."
        RUN_LARGE=false
    elif [[ -d "${SAMPLE_OUT_LARGE}" ]]; then
        log "[WARN] Incomplete large run detected for ${BASE}, clearing output dir."
        rm -rf "${SAMPLE_OUT_LARGE}"
    fi

    if [[ "${RUN_SENS}" == false && "${RUN_LARGE}" == false ]]; then
        log "[SKIP] Both MEGAHIT runs already complete for ${BASE}, skipping."
    else
        SAMPLE_TMP_DIR="${GLOBAL_TMP_DIR}/${BASE}"
        mkdir -p "${SAMPLE_TMP_DIR}"
        DECOMPRESSED_FASTQ="${SAMPLE_TMP_DIR}/${BASE}_interleaved_filtered.fastq"

        log "[DECOMPRESS] ${BASE}"
        pigz -dc "${i}" > "${DECOMPRESSED_FASTQ}" \
            || gzip -dc "${i}" > "${DECOMPRESSED_FASTQ}"

        CONTAINER_FASTQ="/data/${DECOMPRESSED_FASTQ#${PWD}/}"
        CONTAINER_OUT_SENS="/out/${BASE}/megahit_output/sensitive"
        CONTAINER_OUT_LARGE="/out/${BASE}/megahit_output/large"

        mkdir -p "${SAMPLE_OUT_DIR}"

        if [[ "${RUN_SENS}" == true ]]; then
            log "[RUN] megahit sensitive: ${BASE}"
            singularity exec --cleanenv \
              --bind "${PWD}:/data" \
              --bind "${OUT_DIR}:/out" \
              --pwd /data \
              "${MEGAHIT_IMAGE}" \
              megahit \
              --12 "${CONTAINER_FASTQ}" \
              -t "${MEGAHIT_THREADS}" \
              -m "${MEGAHIT_MEM}" \
              --min-count 2 \
              --k-list 21,29,39,49,59,69,79,89,99,109,129,141 \
              -o "${CONTAINER_OUT_SENS}"
        fi

        if [[ "${RUN_LARGE}" == true ]]; then
            log "[RUN] megahit large: ${BASE}"
            singularity exec --cleanenv \
              --bind "${PWD}:/data" \
              --bind "${OUT_DIR}:/out" \
              --pwd /data \
              "${MEGAHIT_IMAGE}" \
              megahit \
              --12 "${CONTAINER_FASTQ}" \
              -t "${MEGAHIT_THREADS}" \
              -m "${MEGAHIT_MEM}" \
              --k-min 28 \
              --k-max 128 \
              --k-step 10 \
              -o "${CONTAINER_OUT_LARGE}"
        fi

        rm -rf "${SAMPLE_TMP_DIR}"
    fi

    # Symlink each completed assembly into OUT_DIR with the flat names exec_binning.sh expects.
    for _mode in sensitive large; do
        _contigs="${SAMPLE_OUT_DIR}/${_mode}/final.contigs.fa"
        _link="${OUT_DIR}/${BASE}_megahit_${_mode}.fasta"
        if [[ -f "${_contigs}" ]]; then
            ln -sfn "${_contigs}" "${_link}"
            log "[LINK] ${BASE}_megahit_${_mode}.fasta -> megahit_output/${_mode}/final.contigs.fa"
        fi
    done
done

log "MEGAHIT assembly completed."
