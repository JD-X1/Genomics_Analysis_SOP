#!/usr/bin/bash -l

module load singularityce
module load bwa/0.7.19
module load samtools/1.21

######################################################
# Metagenome binning: MaxBin2, CONCOCT, MetaBAT2, QuickBin
# 1. BWA-MEM read mapping + depth profiling
# 2. MaxBin2
# 3. CONCOCT
# 4. MetaBAT2
# 5. QuickBin
######################################################
# Usage: ./exec_binning.sh <assembly_dir> <reads_dir> <output_dir>
######################################################

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

if [[ $# -ne 3 ]]; then
    log "[ERROR] Usage: $0 <assembly_dir> <reads_dir> <output_dir>"
    exit 2
fi

ASSEMBLY_DIR="$(realpath -m "$1")"
READS_DIR="$(realpath -m "$2")"
OUT_DIR="$(realpath -m "$3")"

METABAT2_IMAGE="${METABAT2_IMAGE:-metabat2.sif}"
BBTOOLS_IMAGE="${BBTOOLS_IMAGE:-bbtools_latest.sif}"
CONDA="${CONDA:-/home/jduque2/miniforge3/bin/conda}"
MAXBIN2_ENV="${MAXBIN2_ENV:-maxbin2}"
CONCOCT_ENV="${CONCOCT_ENV:-concoct}"
# Derive MaxBin2 bin directory from the conda installation root + env name.
# Override MAXBIN2_BIN directly if your layout differs.
_CONDA_ROOT="$(dirname "$(dirname "${CONDA}")")"
MAXBIN2_BIN="${MAXBIN2_BIN:-${_CONDA_ROOT}/envs/${MAXBIN2_ENV}/bin}"

for _img in "${METABAT2_IMAGE}" "${BBTOOLS_IMAGE}"; do
    if [[ ! -f "${_img}" ]]; then
        log "[ERROR] Singularity image not found: ${_img}"
        exit 1
    fi
done

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MIN_CONTIG_LEN="${MIN_CONTIG_LEN:-2500}"
CONCOCT_CHUNK="${CONCOCT_CHUNK:-10000}"
MAX_PARALLEL="${MAX_PARALLEL:-1}"
THREADS_PER_SAMPLE=$(( THREADS / MAX_PARALLEL ))

RUN_MAXBIN2="${RUN_MAXBIN2:-true}"
RUN_CONCOCT="${RUN_CONCOCT:-true}"
RUN_METABAT2="${RUN_METABAT2:-true}"
RUN_QUICKBIN="${RUN_QUICKBIN:-true}"

ASSEMBLY_SUFFIXES=( "_megahit_large" "_megahit_sensitive" "_spades" )
ASSEMBLY_SUFFIXES_STR="${ASSEMBLY_SUFFIXES[*]}"

mkdir -p "${OUT_DIR}"

declare -A SAMPLE_CONTIGS=()
declare -A SAMPLE_READS=()

for fasta in "${ASSEMBLY_DIR}"/*.fasta; do
    FASTA_STEM="$(basename "${fasta}" .fasta)"

    SAMPLE_BASE=""
    for suffix in "${ASSEMBLY_SUFFIXES[@]}"; do
        if [[ "${FASTA_STEM}" == *"${suffix}" ]]; then
            SAMPLE_BASE="${FASTA_STEM%${suffix}}"
            break
        fi
    done

    if [[ -z "${SAMPLE_BASE}" ]]; then
        log "[WARN] Could not determine sample name from ${FASTA_STEM}, skipping."
        continue
    fi

    if [[ -n "${SAMPLE_FILTER:-}" && "${FASTA_STEM}" != *"${SAMPLE_FILTER}"* ]]; then
        continue
    fi

    READS="${READS_DIR}/${SAMPLE_BASE}_bbcms_interleaved.fastq.gz"
    if [[ ! -f "${READS}" ]]; then
        log "[WARN] No reads found for ${FASTA_STEM}, skipping."
        continue
    fi

    SAMPLE_CONTIGS["${FASTA_STEM}"]="${fasta}"
    SAMPLE_READS["${FASTA_STEM}"]="${READS}"
    log "[FOUND] ${FASTA_STEM}"
done

if [[ ${#SAMPLE_CONTIGS[@]} -eq 0 ]]; then
    log "[ERROR] No valid (contigs + reads) pairs found. Check assembly_dir and reads_dir."
    exit 1
fi

log "Processing ${#SAMPLE_CONTIGS[@]} sample(s) with ${THREADS_PER_SAMPLE} threads (MAX_PARALLEL=${MAX_PARALLEL})."

GLOBAL_TMP="$(mktemp -d "${PWD}/binning_tmp.XXXXXXXXXX")"
trap 'rm -rf "${GLOBAL_TMP}"' EXIT

cpath() {
    local p="$1"
    p="${p/#${ASSEMBLY_DIR}//assembly}"
    p="${p/#${READS_DIR}//reads}"
    p="${p/#${OUT_DIR}//out}"
    p="${p/#${GLOBAL_TMP}//tmp_run}"
    echo "$p"
}

process_sample() {
    local BASE="$1"
    local CONTIGS="${ASSEMBLY_DIR}/${BASE}.fasta"
    local THREADS=${THREADS_PER_SAMPLE}

    local SAMPLE_BASE=""
    for suffix in ${ASSEMBLY_SUFFIXES_STR}; do
        if [[ "${BASE}" == *"${suffix}" ]]; then
            SAMPLE_BASE="${BASE%${suffix}}"
            break
        fi
    done
    local READS="${READS_DIR}/${SAMPLE_BASE}_bbcms_interleaved.fastq.gz"

    local SAMPLE_OUT="${OUT_DIR}/${BASE}"
    local SAMPLE_TMP="${GLOBAL_TMP}/${BASE}"
    mkdir -p "${SAMPLE_OUT}" "${SAMPLE_TMP}"

    local COMMON_BINDS=(
        --bind "${ASSEMBLY_DIR}:/assembly"
        --bind "${READS_DIR}:/reads"
        --bind "${OUT_DIR}:/out"
        --bind "${GLOBAL_TMP}:/tmp_run"
        --bind "${PWD}:/pwd"
    )

    log "[SAMPLE] ${BASE}"

    # Phase 1: Read mapping + depth profiling

    local MAPPING_DIR="${SAMPLE_OUT}/mapping"
    local SORTED_BAM="${MAPPING_DIR}/${BASE}.sorted.bam"
    local DEPTH_FILE="${MAPPING_DIR}/${BASE}_depth.txt"
    local ABUND_FILE="${MAPPING_DIR}/${BASE}_abund.txt"
    mkdir -p "${MAPPING_DIR}"

    if [[ -f "${SORTED_BAM}" && -f "${SORTED_BAM}.bai" ]]; then
        log "[SKIP] Sorted BAM already exists for ${BASE}"
    else
        log "[RUN] BWA index + mapping: ${BASE}"

        local DECOMPRESSED="${SAMPLE_TMP}/${BASE}.fastq"
        pigz -dc "${READS}" > "${DECOMPRESSED}" \
            || gzip -dc "${READS}" > "${DECOMPRESSED}"

        if [[ ! -f "${CONTIGS}.bwt" ]]; then
            bwa index "${CONTIGS}"
        fi

        bwa mem -p -t "${THREADS}" "${CONTIGS}" "${DECOMPRESSED}" \
            | samtools sort -@ "${THREADS}" -o "${SORTED_BAM}" - \
            && samtools index "${SORTED_BAM}"

        rm -f "${DECOMPRESSED}"
    fi

    if [[ -f "${DEPTH_FILE}" && -f "${ABUND_FILE}" ]]; then
        log "[SKIP] Depth files already exist for ${BASE}"
    else
        log "[RUN] jgi_summarize_bam_contig_depths: ${BASE}"

        singularity exec --cleanenv "${COMMON_BINDS[@]}" --pwd /pwd \
            "${METABAT2_IMAGE}" \
            jgi_summarize_bam_contig_depths \
            --outputDepth "$(cpath "${DEPTH_FILE}")" \
            "$(cpath "${SORTED_BAM}")"

        awk 'NR>1 {print $1"\t"$3}' "${DEPTH_FILE}" > "${ABUND_FILE}"
    fi

    # Phase 2: MaxBin2

    if [[ "${RUN_MAXBIN2}" == true ]]; then
        local MAXBIN2_DIR="${SAMPLE_OUT}/maxbin2"
        local MAXBIN2_SUMMARY="${MAXBIN2_DIR}/bin.summary"
        mkdir -p "${MAXBIN2_DIR}"

        if [[ -f "${MAXBIN2_SUMMARY}" ]]; then
            log "[SKIP] MaxBin2 already complete for ${BASE}"
        else
            log "[RUN] MaxBin2: ${BASE}"

            PATH="${MAXBIN2_BIN}:${PATH}" \
            "${MAXBIN2_BIN}/perl" \
                "${MAXBIN2_BIN}/run_MaxBin.pl" \
                -contig  "${CONTIGS}" \
                -abund   "${ABUND_FILE}" \
                -out     "${MAXBIN2_DIR}/bin" \
                -thread  "${THREADS}" \
                1> "${MAXBIN2_DIR}/maxbin2.${BASE}.log" \
                2> "${MAXBIN2_DIR}/maxbin2.${BASE}.err"
        fi
    fi

    # Phase 3: CONCOCT

    if [[ "${RUN_CONCOCT}" == true ]]; then
        local CONCOCT_DIR="${SAMPLE_OUT}/concoct"
        local CONCOCT_MERGED="${CONCOCT_DIR}/clustering_merged.csv"
        mkdir -p "${CONCOCT_DIR}/fasta_bins"

        if [[ -f "${CONCOCT_MERGED}" ]]; then
            log "[SKIP] CONCOCT already complete for ${BASE}"
        else
            log "[RUN] CONCOCT: ${BASE}"

            local CONCOCT_CUT_FA="${CONCOCT_DIR}/contigs_cut.fa"
            local CONCOCT_CUT_BED="${CONCOCT_DIR}/contigs_cut.bed"
            local CONCOCT_COV="${CONCOCT_DIR}/coverage_table.tsv"

            if [[ ! -s "${CONCOCT_CUT_FA}" ]]; then
                "${CONDA}" run --no-capture-output -n "${CONCOCT_ENV}" \
                    cut_up_fasta.py "${CONTIGS}" \
                    -c "${CONCOCT_CHUNK}" \
                    -o 0 --merge_last \
                    -b "${CONCOCT_CUT_BED}" \
                    > "${CONCOCT_CUT_FA}"
            fi

            if [[ ! -s "${CONCOCT_COV}" ]]; then
                "${CONDA}" run --no-capture-output -n "${CONCOCT_ENV}" \
                    concoct_coverage_table.py \
                    "${CONCOCT_CUT_BED}" \
                    "${SORTED_BAM}" \
                    > "${CONCOCT_COV}"
            fi

            "${CONDA}" run --no-capture-output -n "${CONCOCT_ENV}" \
                concoct \
                --composition_file "${CONCOCT_CUT_FA}" \
                --coverage_file    "${CONCOCT_COV}" \
                --threads          "${THREADS}" \
                --length_threshold 1000 \
                -b                 "${CONCOCT_DIR}/" \
                1> "${CONCOCT_DIR}/concoct.${BASE}.log" \
                2> "${CONCOCT_DIR}/concoct.${BASE}.err"

            local CONCOCT_CLUSTERING=""
            if [[ -f "${CONCOCT_DIR}/clustering_gt1000.csv" ]]; then
                CONCOCT_CLUSTERING="${CONCOCT_DIR}/clustering_gt1000.csv"
            else
                CONCOCT_CLUSTERING="$(ls "${CONCOCT_DIR}"/clustering_gt*.csv 2>/dev/null | head -1)"
            fi
            if [[ -z "${CONCOCT_CLUSTERING}" || ! -f "${CONCOCT_CLUSTERING}" ]]; then
                log "[ERROR] CONCOCT clustering CSV not found in ${CONCOCT_DIR}"
                exit 1
            fi

            "${CONDA}" run --no-capture-output -n "${CONCOCT_ENV}" \
                merge_cutup_clustering.py "${CONCOCT_CLUSTERING}" \
                > "${CONCOCT_MERGED}"

            "${CONDA}" run --no-capture-output -n "${CONCOCT_ENV}" \
                extract_fasta_bins.py \
                "${CONTIGS}" \
                "${CONCOCT_MERGED}" \
                --output_path "${CONCOCT_DIR}/fasta_bins"
        fi
    fi

    # Phase 4: MetaBAT2

    if [[ "${RUN_METABAT2}" == true ]]; then
        local METABAT2_DIR="${SAMPLE_OUT}/metabat2"
        local METABAT2_DONE="${METABAT2_DIR}/.done"
        mkdir -p "${METABAT2_DIR}"

        if [[ -f "${METABAT2_DONE}" ]]; then
            log "[SKIP] MetaBAT2 already complete for ${BASE}"
        else
            log "[RUN] MetaBAT2: ${BASE}"

            singularity exec --cleanenv "${COMMON_BINDS[@]}" --pwd /pwd \
                "${METABAT2_IMAGE}" \
                metabat2 \
                -i "$(cpath "${CONTIGS}")" \
                -a "$(cpath "${DEPTH_FILE}")" \
                -o "$(cpath "${METABAT2_DIR}")/bin" \
                -t "${THREADS}" \
                -m "${MIN_CONTIG_LEN}" \
                --seed 42 \
                1> "${METABAT2_DIR}/metabat2.${BASE}.log" \
                2> "${METABAT2_DIR}/metabat2.${BASE}.err"

            touch "${METABAT2_DONE}"
        fi
    fi

    # Phase 5: QuickBin

    if [[ "${RUN_QUICKBIN}" == true ]]; then
        local QUICKBIN_DIR="${SAMPLE_OUT}/quickbin"
        local QUICKBIN_DONE="${QUICKBIN_DIR}/.done"
        local QUICKBIN_COV="${QUICKBIN_DIR}/coverage.txt"
        mkdir -p "${QUICKBIN_DIR}"

        if [[ -f "${QUICKBIN_DONE}" ]]; then
            log "[SKIP] QuickBin already complete for ${BASE}"
        else
            log "[RUN] QuickBin: ${BASE}"

            local QUICKBIN_COV_ARG
            if [[ -f "${QUICKBIN_COV}" ]]; then
                QUICKBIN_COV_ARG="cov=$(cpath "${QUICKBIN_COV}")"
            else
                QUICKBIN_COV_ARG="covout=$(cpath "${QUICKBIN_COV}")"
            fi

            singularity exec --cleanenv "${COMMON_BINDS[@]}" --pwd /pwd \
                "${BBTOOLS_IMAGE}" \
                quickbin.sh \
                in="$(cpath "${CONTIGS}")" \
                out="$(cpath "${QUICKBIN_DIR}")/bin%.fa" \
                "${QUICKBIN_COV_ARG}" \
                threads="${THREADS}" \
                "$(cpath "${SORTED_BAM}")" \
                1> "${QUICKBIN_DIR}/quickbin.${BASE}.log" \
                2> "${QUICKBIN_DIR}/quickbin.${BASE}.err"

            touch "${QUICKBIN_DONE}"
        fi
    fi

    log "[COMPLETE] ${BASE}"
}

export -f process_sample log cpath
export ASSEMBLY_DIR READS_DIR OUT_DIR GLOBAL_TMP METABAT2_IMAGE BBTOOLS_IMAGE
export CONDA MAXBIN2_ENV MAXBIN2_BIN CONCOCT_ENV
export THREADS_PER_SAMPLE MIN_CONTIG_LEN CONCOCT_CHUNK
export RUN_MAXBIN2 RUN_CONCOCT RUN_METABAT2 RUN_QUICKBIN
export ASSEMBLY_SUFFIXES_STR

printf '%s\0' "${!SAMPLE_CONTIGS[@]}" | \
    xargs -0 -n 1 -P "${MAX_PARALLEL}" bash -c 'process_sample "$1"' _

log "That's all folks!"
