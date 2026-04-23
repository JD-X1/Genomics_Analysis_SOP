#!/usr/bin/bash -l

module load singularityce
module load openjdk/17.0.5_8


######################################################
# three phases of execution:
# 1. Bloom Filtering with bbtools
# 2. SPADES
# 3. fungalrelease.sh
######################################################
# Usage: ./exec_assembly.sh <input_dir> <output_dir>
######################################################


# Phase 0: Set Up
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

# Configuration, be sure you know your paths!:
INPUT_DIR="$(realpath -m "$1")"  # expects interleaved filtered files
OUT_DIR="$(realpath -m "$2")"

# java virt max should be set to ~80% of available memory; override with MEM_FOR_BLOOM
MEM_FOR_BLOOM="${MEM_FOR_BLOOM:--Xmx160g}"
# SPAdes memory cap in GB; override with MEM_FOR_SPADES
MEM_FOR_SPADES="${MEM_FOR_SPADES:-800}"
SPADES_THREADS="${SPADES_THREADS:-${SLURM_CPUS_PER_TASK:-48}}"
BBTOOLS_IMAGE="${BBTOOLS_IMAGE:-bbtools_38.86.sif}"
SPADES_IMAGE="${SPADES_IMAGE:-spades_3.15.2.sif}"



FILTERED_FASTQS=( ${INPUT_DIR}/*_interleaved_filtered.fastq.gz )
mkdir -p ${OUT_DIR}


if [[ ${#FILTERED_FASTQS[@]} -eq 0 ]]; then
  log "[ERROR] No interleaved files were produced in ${INPUT_DIR}/"
  exit 1
fi

touch readlen.txt

# Phase 1: Bloom filtering with BBTools

BLOOM_FILTERED=()

# Dependencies: BBTools, Java 17

for i in "${FILTERED_FASTQS[@]}"; do
    log "Running BBCMS & Assembly on ${i} ..."

    BASE="$(basename "${i}" _interleaved_filtered.fastq.gz)"

    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}"
    SAMPLE_COUNTS="${SAMPLE_OUT_DIR}/counts.metadata.json"

    SPADES_IN1="${SAMPLE_OUT_DIR}/bbcms_output.FWD.fastq.gz"
    SPADES_IN2="${SAMPLE_OUT_DIR}/bbcms_output.REV.fastq.gz"

    INTERLEAVED_OUT="${SAMPLE_OUT_DIR}/${BASE}_bbcms_interleaved.fastq.gz"
    if [[ -f "${INTERLEAVED_OUT}" && -f "${SPADES_IN1}" && -f "${SPADES_IN2}" ]]; then
        log "[SKIP] BBCMS outputs already exist for ${BASE}, skipping."
        BLOOM_FILTERED+=( "${INTERLEAVED_OUT}" )
        continue
    elif [[ -f "${SPADES_IN1}" || -f "${SPADES_IN2}" ]] && [[ ! -f "${INTERLEAVED_OUT}" ]]; then
        log "[ERROR] Paired outputs exist but interleaved file ${INTERLEAVED_OUT} is missing for ${BASE}. Possible incomplete run."
        exit 1
    fi
    
    mkdir -p "${SAMPLE_OUT_DIR}"
    log "[RUN] bbcms.sh: $BASE"

    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      readlength.sh -Xmx3g \
      in="${i}" \
      out="${SAMPLE_OUT_DIR}/${BASE}_readlength.txt" \
      overwrite

    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      bbcms.sh ${MEM_FOR_BLOOM} \
        mincount=2 \
        highcountfraction=0.6 \
        in="${i}" \
        out1="${SPADES_IN1}" \
        out2="${SPADES_IN2}" \
        1> "${SAMPLE_OUT_DIR}/bbcms.${BASE}.log" \
        2> "${SAMPLE_OUT_DIR}/bbcms.${BASE}.err"
    
    singularity exec --cleanenv \
      --bind "$PWD:/data" \
      --pwd /data \
      "${BBTOOLS_IMAGE}" \
      reformat.sh -Xmx5g \
        in1="${SPADES_IN1}" \
        in2="${SPADES_IN2}" \
        out="${SAMPLE_OUT_DIR}/${BASE}_bbcms_interleaved.fastq.gz" \
        overwrite=t
    BLOOM_FILTERED+=( "${SAMPLE_OUT_DIR}/${BASE}_bbcms_interleaved.fastq.gz" )
    grep Uniq ${SAMPLE_OUT_DIR}/bbcms.${BASE}.err | awk '{print $NF}' > "${SAMPLE_OUT_DIR}/${BASE}_unique31mer.txt"
    
done

log "Bloom filtering completed. Processed files: ${#BLOOM_FILTERED[@]}"

# PHASE 2: SPADES Assembly

log "Checking for spades singularity image at ${SPADES_IMAGE} ..."

if [[ ! -f "${SPADES_IMAGE}" ]]; then
    log "[ERROR] SPADES singularity image not found at ${SPADES_IMAGE}."
    log "[ERROR] Please ensure the image is built and available."
    log "[ERROR] And check that the singularity image path is correct"
    log "[ERROR] Exiting with code 1."
    exit 1
else
    log "[INFO] SPADES singularity image found at ${SPADES_IMAGE}."
fi

for i in "${BLOOM_FILTERED[@]}"; do
    log "Running SPADES Assembly on ${i} ..."

    BASE="$(basename "${i}" _bbcms_interleaved.fastq.gz)"
    SAMPLE_OUT_DIR="${OUT_DIR}/${BASE}" # should be same as bloom sample dir
    SAMPLE_TMP_DIR="${OUT_DIR}/${BASE}/tmp"
    mkdir -p "${SAMPLE_TMP_DIR}"

    if [[ -f "${SAMPLE_OUT_DIR}/spades_output/scaffolds.fasta" ]]; then
        log "[SKIP] SPAdes output already exists for ${BASE}"
    else
        singularity exec --cleanenv \
          --bind "$PWD:/data" \
          --pwd /data \
          "${SPADES_IMAGE}" \
          spades.py \
          -m ${MEM_FOR_SPADES} \
          -t ${SPADES_THREADS} \
          --tmp-dir ${SAMPLE_TMP_DIR} \
          -o "${SAMPLE_OUT_DIR}/spades_output" \
          --only-assembler \
          -k 33,55,77,99,127 \
          --meta \
          -1 "${SAMPLE_OUT_DIR}/bbcms_output.FWD.fastq.gz" \
          -2 "${SAMPLE_OUT_DIR}/bbcms_output.REV.fastq.gz"
    fi

    # Symlink the assembly into OUT_DIR with the flat name exec_binning.sh expects.
    # Prefers scaffolds.fasta; falls back to contigs.fasta if scaffolding was skipped.
    _spades_link="${OUT_DIR}/${BASE}_spades.fasta"
    if [[ ! -e "${_spades_link}" ]]; then
        if [[ -f "${SAMPLE_OUT_DIR}/spades_output/scaffolds.fasta" ]]; then
            ln -sfn "${SAMPLE_OUT_DIR}/spades_output/scaffolds.fasta" "${_spades_link}"
            log "[LINK] ${BASE}_spades.fasta -> spades_output/scaffolds.fasta"
        elif [[ -f "${SAMPLE_OUT_DIR}/spades_output/contigs.fasta" ]]; then
            ln -sfn "${SAMPLE_OUT_DIR}/spades_output/contigs.fasta" "${_spades_link}"
            log "[LINK] ${BASE}_spades.fasta -> spades_output/contigs.fasta (scaffolding fallback)"
        else
            log "[WARN] No SPAdes FASTA found for ${BASE} — check spades_output/ for errors"
        fi
    fi
done