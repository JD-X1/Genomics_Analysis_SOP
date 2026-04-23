#!/usr/bin/bash -l

# ---------------------       README BEFORE EXEC       ---------------------
# virtual memory setting associated with each of the bbtools commands
# can be adjusted manually in this script. Reccomendations are:
#
# reformat.sh : 10-20g (if interleaving is slow you might have allocated too
#                      much memory)
#
# readlength.sh: 1-4g (generally low memory usage)
#
# rqcfilter2.sh: run entirely in series so highly dependent on dataset 
#                size + avail mem reccomend 70-80% of available memory
#       `        on the node any more and you might run the risk of 
#                out-of-memory kills + significant slowdowns due to swapping.
#                Remember you control memory available on SLURM with --mem flag
#                e.g., --mem=200G for 200 gigabytes of RAM
#
# --------------------------------------------------------------------------------
# Labeling Pattern instructions:
# for paired-end files, provide the common part of 
# the filename that identifies forward and reverse
# reads and all text after it.
# e.g., for files named:
# SampleA_R1_001.fastq.gz
# SampleA_R2_001.fastq.gz
# the labeling patterns would be:
# LABELING_PATTERN_FWD="_R1_001.fastq.gz"
# LABELING_PATTERN_REV="_R2_001.fastq.gz"
# --------------------------------------------------------------------------------
# usage:
# ./exec_rqcFilter.sh <input_dir> \
#   <labeling_pattern_fwd> \
#   <labeling_pattern_rev> \
#   <output_dir>
# --------------------------------------------------------------------------------

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

if [[ $# -ne 4 ]]; then
  log "[ERROR] Usage: $0 <input_dir> <labeling_pattern_fwd> <labeling_pattern_rev> <output_dir>"
  exit 2
fi

INPUT_DIR=$1
LABELING_SUFFIXPATTERN_FWD=$2 # e.g., "_R1_001.fastq.gz"
LABELING_SUFFIXPATTERN_REV=$3 # e.g., "_R2_001.fastq.gz"
OUTPUT_DIR=$4


RQC_INPUT_FILES=()
mkdir -p ${OUTPUT_DIR}

REFORMAT_JOBS=${REFORMAT_JOBS:-${SLURM_CPUS_PER_TASK:-8}}

FWD_FILES=( ${INPUT_DIR}/*${LABELING_SUFFIXPATTERN_FWD} )

if [[ ${#FWD_FILES[@]} -eq 0 ]]; then
    log "[ERROR] No files found in input directory matching pattern: ${LABELING_SUFFIXPATTERN_FWD}"
    exit 1
fi

module load singularityce
module load openjdk/17.0.5_8

log "Interleaving ${#FWD_FILES[@]} read-pairs using ${REFORMAT_JOBS} parallel jobs..."

do_interleave() {
    local FWD_FILE="$1"
    local BASE_NAME
    local REV_FILE
    local OUT_INTERLEAVED

    BASE_NAME=$(basename "$FWD_FILE")
    REV_FILE="${INPUT_DIR}/${BASE_NAME/${LABELING_SUFFIXPATTERN_FWD}/${LABELING_SUFFIXPATTERN_REV}}"
    
    if [[ ! -f "$REV_FILE" ]]; then
        log "[ERROR] Missing reverse pair for: ${FWD_FILE}"
        log "[ERROR] Expected: ${REV_FILE}"
        log "[ERROR] Check your data + labeling patterns. This script only handles paired-end files."
        exit 1
    fi

    OUT_INTERLEAVED="raw_reads/${BASE_NAME/${LABELING_SUFFIXPATTERN_FWD}/_interleaved}.fastq.gz"
    
    if  [[ -f "$OUT_INTERLEAVED" ]]; then
        log "Interleaved file already exists, skipping: ${OUT_INTERLEAVED}"
        return 0
    fi

    log "Interleaving: ${FWD_FILE} + ${REV_FILE} -> ${OUT_INTERLEAVED}"

    singularity exec --cleanenv \
        --bind "$PWD:/data" \
        --pwd /data \
        bbtools_38.86.sif \
        reformat.sh -Xmx20g \
        in="$FWD_FILE" \
        in2="$REV_FILE" \
        out="$OUT_INTERLEAVED" \
        overwrite=t
}

export INPUT_DIR LABELING_SUFFIXPATTERN_FWD LABELING_SUFFIXPATTERN_REV
export -f do_interleave log

printf '%s\0' "${FWD_FILES[@]}" | xargs -0 -n 1 -P "${REFORMAT_JOBS}" bash -c 'do_interleave "$1"' _

RQC_INPUT_FILES=( raw_reads/*_interleaved.fastq.gz )

if [[ ${#RQC_INPUT_FILES[@]} -eq 0 ]]; then
  log "[ERROR] No interleaved files were produced in raw_reads/"
  exit 1
fi

touch readlen.txt


for i in "${RQC_INPUT_FILES[@]}"; do
    log "Running RQC Filter on file: $i"
    BASE_NAME=$(basename $i)

    SAMPLE="${BASE_NAME/.fastq.gz/}" # define samp name
    
    # line below  may need adjustment if the lines defining FINAL_OUT
    # are changed. So double check the rqcfilter2.sh 
    # path flag -> where all final outputs will end up
    # out flag -> from the path flag it will append this suffix
    # to create the final output path
    #
    # the consequence of 
    # path="${OUTPUT_DIR}/filter_${BASE_NAME/.fastq.gz/}" \
    # out="${OUTPUT_DIR}/${BASE_NAME/.fastq.gz/_filtered.fastq.gz}"
    #
    # is that the final output will be:
    # {OUTPUT_DIR}/filter_${BASE_NAME/.fastq.gz/}${OUTPUT_DIR}/${BASE_NAME/.fastq.gz/_filtered.fastq.gz}

    FINAL_OUT="${OUTPUT_DIR}/filter_${SAMPLE}/tempo/${SAMPLE}_filtered.fastq.gz" # line may need adjusttment if 

    if [[ -s "$FINAL_OUT" ]]; then
        log "RQC Filter output already exists, skipping: ${FINAL_OUT}"
        continue
    fi

    singularity exec --cleanenv \
       --bind "$PWD:/data" \
       --pwd /data \
       bbtools_38.86.sif \
       readlength.sh -Xmx3g \
       in=$i \
       out=readlen.txt \
       overwrite


    # this command needs to be adjusted see IMPORTANT NOTE above
    # in this loop
    singularity exec --cleanenv \
        --bind "$PWD:/data" \
        --pwd /data \
        bbtools_38.86.sif \
        rqcfilter2.sh -Xmx96g chastityfilter=f jni=t \
        in=$i \
        rqcfilterdata=/data/RQCFilterData/ \
        path="${OUTPUT_DIR}/filter_${BASE_NAME/.fastq.gz/}" \
        rna=f trimfragadapter=t qtrim=r \
        trimq=0 maxns=3 maq=3 minlen=51 \
        mlf=0.33 phix=t removehuman=t \
        removedog=t removecat=t removemouse=t khist=t \
        removemicrobes=t sketch kapa=t clumpify=t tmpdir="${OUTPUT_DIR}" \
        barcodefilter=f trimpolyg=5 usejni=f \
        out="${OUTPUT_DIR}/${BASE_NAME/.fastq.gz/_filtered.fastq.gz}"
done

 log "That's all folks!"
