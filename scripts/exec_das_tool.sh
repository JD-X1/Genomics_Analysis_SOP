#!/usr/bin/bash -l

######################################################
# DAS_Tool bin refinement
# 1. Build contig2bin TSVs from each completed binner
# 2. Run DAS_Tool to produce a refined, non-redundant bin set
######################################################
# Usage: ./exec_das_tool.sh <binning_output_dir> <assembly_dir>
######################################################

set -Eeuo pipefail
shopt -s nullglob

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    echo "[INFO] SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

log() { echo "[$(date +'%F %T')] $*"; }
trap 'rc=$?; log "[ERROR] line ${LINENO}: command failed: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR

if [[ $# -ne 2 ]]; then
    log "[ERROR] Usage: $0 <binning_output_dir> <assembly_dir>"
    exit 2
fi

BINNING_DIR="$(realpath -m "$1")"
ASSEMBLY_DIR="$(realpath -m "$2")"

THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MAX_PARALLEL="${MAX_PARALLEL:-1}"
THREADS_PER_SAMPLE=$(( THREADS / MAX_PARALLEL ))
MIN_BINNERS="${MIN_BINNERS:-2}"
SCORE_THRESHOLD="${SCORE_THRESHOLD:-0.5}"
DAS_TOOL_ENV="${DAS_TOOL_ENV:-das-tool}"
CONDA="${CONDA:-/home/jduque2/miniforge3/bin/conda}"

[[ -d "${BINNING_DIR}" ]] || { log "[ERROR] Binning output dir not found: ${BINNING_DIR}"; exit 1; }
[[ -d "${ASSEMBLY_DIR}" ]] || { log "[ERROR] Assembly dir not found: ${ASSEMBLY_DIR}"; exit 1; }

mapfile -t SAMPLE_DIRS < <(
    find "${BINNING_DIR}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
)

if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
    log "[ERROR] No sample directories found in ${BINNING_DIR}"
    exit 1
fi

# Apply optional sample filter
declare -a SAMPLES=()
for s in "${SAMPLE_DIRS[@]}"; do
    if [[ -n "${SAMPLE_FILTER:-}" && "${s}" != *"${SAMPLE_FILTER}"* ]]; then
        continue
    fi
    SAMPLES+=( "${s}" )
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    log "[ERROR] No samples match SAMPLE_FILTER='${SAMPLE_FILTER:-}'"
    exit 1
fi

log "Processing ${#SAMPLES[@]} sample(s) with ${THREADS_PER_SAMPLE} threads (MAX_PARALLEL=${MAX_PARALLEL})."

process_sample() {
    local sample="$1"
    local sample_dir="${BINNING_DIR}/${sample}"
    local das_dir="${sample_dir}/das_tool"
    local das_done="${das_dir}/.done"

    if [[ -f "${das_done}" ]]; then
        log "[SKIP] DAS_Tool already complete for ${sample}"
        return 0
    fi

    # Binner configuration (associative arrays can't be exported; define inline)
    local -a BINNER_ORDER=( maxbin2 concoct metabat2 quickbin )
    declare -A BINNER_SUBDIR=(
        [maxbin2]="maxbin2"
        [concoct]="concoct/fasta_bins"
        [metabat2]="metabat2"
        [quickbin]="quickbin"
    )
    declare -A BINNER_EXT=(
        [maxbin2]="fasta"
        [concoct]="fa"
        [metabat2]="fa"
        [quickbin]="fa"
    )
    # Completion markers written by exec_binning.sh
    declare -A BINNER_MARKER=(
        [maxbin2]="maxbin2/bin.summary"
        [concoct]="concoct/clustering_merged.csv"
        [metabat2]="metabat2/.done"
        [quickbin]="quickbin/.done"
    )

    mkdir -p "${das_dir}"

    local -a c2b_files=()
    local -a active_labels=()

    for binner in "${BINNER_ORDER[@]}"; do
        local marker="${sample_dir}/${BINNER_MARKER[${binner}]}"
        if [[ ! -f "${marker}" ]]; then
            log "[INFO] ${sample}/${binner}: not yet complete, skipping"
            continue
        fi

        local bin_dir="${sample_dir}/${BINNER_SUBDIR[${binner}]}"
        local ext="${BINNER_EXT[${binner}]}"
        local c2b="${das_dir}/${binner}_contig2bin.tsv"

        # Guard: bin_dir must exist and contain bin files
        if [[ ! -d "${bin_dir}" ]]; then
            log "[WARN] ${sample}/${binner}: bin dir missing (${bin_dir}), skipping"
            continue
        fi
        local n_bins
        n_bins=$(find "${bin_dir}" -maxdepth 1 -name "*.${ext}" 2>/dev/null | wc -l)
        if [[ "${n_bins}" -eq 0 ]]; then
            log "[WARN] ${sample}/${binner}: no *.${ext} files in ${bin_dir}, skipping"
            continue
        fi

        local regenerate_c2b=false
        if [[ -s "${c2b}" ]]; then
            local c2b_nf
            c2b_nf="$(awk -F'\t' 'NF {print NF; exit}' "${c2b}")"
            if [[ -z "${c2b_nf}" || "${c2b_nf}" -ne 2 ]]; then
                log "[WARN] ${sample}/${binner}: malformed contig2bin (${c2b_nf:-0} columns), regenerating"
                rm -f "${c2b}"
                regenerate_c2b=true
            fi
        else
            regenerate_c2b=true
        fi

        if [[ "${regenerate_c2b}" == true ]]; then
            log "[RUN] Fasta_to_Contig2Bin: ${sample}/${binner} (${n_bins} bins)"
            # Post-process: strip extra fields from FASTA headers so output is
            # strictly 2-column TSV (contig_id<TAB>bin_id).
            # metabat2 embeds tab-separated depth fields; concoct/quickbin embed
            # space-separated flag/multi/len fields.  In both cases we take the
            # first whitespace/tab-delimited token of col-1 and the last tab-field.
            "${CONDA}" run --no-capture-output -n "${DAS_TOOL_ENV}" \
                Fasta_to_Contig2Bin.sh \
                -i "${bin_dir}" \
                -e "${ext}" \
                | awk -F'\t' '{split($1,a,/[[:space:]]+/); print a[1] "\t" $NF}' \
                > "${c2b}"
        else
            log "[SKIP] contig2bin exists: ${sample}/${binner}"
        fi

        if [[ ! -s "${c2b}" ]]; then
            log "[WARN] Empty contig2bin for ${sample}/${binner}, skipping"
            rm -f "${c2b}"
            continue
        fi

        c2b_files+=( "${c2b}" )
        active_labels+=( "${binner}" )
    done

    local n_active=${#active_labels[@]}
    if [[ "${n_active}" -lt "${MIN_BINNERS}" ]]; then
        log "[SKIP] ${sample}: ${n_active}/${MIN_BINNERS} binner(s) ready (${active_labels[*]:-none})"
        return 0
    fi

    local contigs="${ASSEMBLY_DIR}/${sample}.fasta"
    if [[ ! -f "${contigs}" ]]; then
        log "[ERROR] Contigs not found: ${contigs}"
        return 1
    fi

    local c2b_list label_list
    c2b_list="$(IFS=,; echo "${c2b_files[*]}")"
    label_list="$(IFS=,; echo "${active_labels[*]}")"

    log "[RUN] DAS_Tool: ${sample} (${label_list})"

    local das_rc=0
    "${CONDA}" run --no-capture-output -n "${DAS_TOOL_ENV}" \
        DAS_Tool \
        -i "${c2b_list}" \
        -c "${contigs}" \
        -o "${das_dir}/${sample}" \
        -l "${label_list}" \
        -t "${THREADS_PER_SAMPLE}" \
        --score_threshold "${SCORE_THRESHOLD}" \
        --write_bins \
        --write_bin_evals \
        1> "${das_dir}/das_tool.${sample}.log" \
        2> "${das_dir}/das_tool.${sample}.err" || das_rc=$?

    if [[ "${das_rc}" -ne 0 ]]; then
        log "[ERROR] DAS_Tool failed for ${sample} (exit=${das_rc}); see ${das_dir}/das_tool.${sample}.err"
        return "${das_rc}"
    fi

    touch "${das_done}"
    log "[COMPLETE] ${sample}"
}

export -f process_sample log
export BINNING_DIR ASSEMBLY_DIR THREADS_PER_SAMPLE MIN_BINNERS SCORE_THRESHOLD DAS_TOOL_ENV CONDA

printf '%s\0' "${SAMPLES[@]}" | \
    xargs -0 -n 1 -P "${MAX_PARALLEL}" bash -c 'process_sample "$1"' _

log "That's all folks!"
