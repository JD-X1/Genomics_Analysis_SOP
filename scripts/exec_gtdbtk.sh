#!/usr/bin/bash -l

######################################################
# GTDB-Tk classification workflow with ephemeral DB
# 1. Check /scratch quota headroom
# 2. Download and extract GTDB-Tk reference package
# 3. Run gtdbtk classify_wf
# 4. Remove the downloaded database
######################################################
# Usage:
#   ./exec_gtdbtk.sh --genome_dir <dir> --out_dir <dir> [gtdbtk classify_wf args...]
#   ./exec_gtdbtk.sh --batchfile <tsv> --out_dir <dir> [gtdbtk classify_wf args...]
######################################################

set -Eeuo pipefail

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    echo "[INFO] SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

log() { echo "[$(date +'%F %T')] $*"; }

cleanup() {
    if [[ "${KEEP_DB:-false}" != true && -n "${DB_ROOT:-}" && -d "${DB_ROOT}" ]]; then
        log "[CLEANUP] Removing GTDB-Tk database workspace: ${DB_ROOT}"
        rm -rf "${DB_ROOT}"
    fi
}

trap 'rc=$?; log "[ERROR] line ${LINENO}: ${BASH_COMMAND} (exit=${rc})"; cleanup; exit ${rc}' ERR
trap 'cleanup' EXIT

CONDA="${CONDA:-/home/jduque2/miniforge3/bin/conda}"
GTDBTK_ENV="${GTDBTK_ENV:-gtdbtk-2.6.1}"
REQUIRED_DB_GB="${REQUIRED_DB_GB:-140}"
PRIMARY_DB_URL="${PRIMARY_DB_URL:-https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz}"
MIRROR_DB_URL="${MIRROR_DB_URL:-https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz}"
WORKDIR_REAL="$(pwd -P)"
DB_ROOT="${DB_ROOT:-$(mktemp -d -p "${WORKDIR_REAL}" gtdbtk_db.XXXXXX)}"
DB_ARCHIVE="${DB_ROOT}/gtdbtk_data.tar.gz"
DB_DIR="${DB_ROOT}/db"
RUN_CHECK_INSTALL="${RUN_CHECK_INSTALL:-true}"

[[ $# -gt 0 ]] || { log "[ERROR] No GTDB-Tk arguments provided."; exit 2; }
[[ "${WORKDIR_REAL}" == /scratch/* ]] || { log "[ERROR] Run this script from /scratch so quota checks match the working filesystem. cwd=${PWD} resolved=${WORKDIR_REAL}"; exit 1; }
[[ -x "${CONDA}" ]] || { log "[ERROR] conda executable not found: ${CONDA}"; exit 1; }
command -v quota >/dev/null 2>&1 || { log "[ERROR] quota command not available."; exit 1; }
command -v wget >/dev/null 2>&1 || { log "[ERROR] wget command not available."; exit 1; }
command -v tar >/dev/null 2>&1 || { log "[ERROR] tar command not available."; exit 1; }

declare -a GTDBTK_ARGS=("$@")
GENOME_DIR=""
BATCHFILE=""
OUT_DIR=""
HAS_CPUS=false
HAS_TMPDIR=false
HAS_SCRATCH_DIR=false

for ((i=1; i<=$#; i++)); do
    arg="${!i}"
    next_index=$((i + 1))
    next_value=""
    if [[ ${next_index} -le $# ]]; then
        next_value="${!next_index}"
    fi
    case "${arg}" in
        --genome_dir)
            GENOME_DIR="${next_value}"
            ;;
        --batchfile)
            BATCHFILE="${next_value}"
            ;;
        --out_dir)
            OUT_DIR="${next_value}"
            ;;
        --cpus)
            HAS_CPUS=true
            ;;
        --tmpdir)
            HAS_TMPDIR=true
            ;;
        --scratch_dir)
            HAS_SCRATCH_DIR=true
            ;;
    esac
done

if [[ -n "${GENOME_DIR}" && -n "${BATCHFILE}" ]]; then
    log "[ERROR] Provide only one of --genome_dir or --batchfile."
    exit 2
fi
if [[ -z "${GENOME_DIR}" && -z "${BATCHFILE}" ]]; then
    log "[ERROR] Missing required --genome_dir or --batchfile."
    exit 2
fi
if [[ -z "${OUT_DIR}" ]]; then
    log "[ERROR] Missing required --out_dir."
    exit 2
fi
if [[ -n "${GENOME_DIR}" && ! -d "${GENOME_DIR}" ]]; then
    log "[ERROR] genome_dir not found: ${GENOME_DIR}"
    exit 1
fi
if [[ -n "${BATCHFILE}" && ! -f "${BATCHFILE}" ]]; then
    log "[ERROR] batchfile not found: ${BATCHFILE}"
    exit 1
fi

to_gb() {
    local raw="${1}"
    python3 - "${raw}" <<'PYEOF'
import re
import sys

value = sys.argv[1].strip()
m = re.fullmatch(r'([0-9]*\.?[0-9]+)([KMGTP]?)', value)
if not m:
    raise SystemExit(1)

num = float(m.group(1))
unit = m.group(2)
scale = {
    "": 1 / (1024 ** 3),
    "K": 1 / (1024 ** 2),
    "M": 1 / 1024,
    "G": 1,
    "T": 1024,
    "P": 1024 * 1024,
}
print(num * scale[unit])
PYEOF
}

quota_fields="$(
    quota -vs | awk '
        /\/scratch/ {
            if (NF >= 3) {
                print $2, $3
                exit
            }
            getline
            while ($0 ~ /^[[:space:]]*$/) {
                getline
            }
            print $1, $2
            exit
        }
    '
)"
[[ -n "${quota_fields}" ]] || { log "[ERROR] Could not find /scratch quota line from 'quota -vs'."; exit 1; }

scratch_used_raw="$(awk '{print $1}' <<< "${quota_fields}")"
scratch_quota_raw="$(awk '{print $2}' <<< "${quota_fields}")"
scratch_used_gb="$(to_gb "${scratch_used_raw}")"
scratch_quota_gb="$(to_gb "${scratch_quota_raw}")"
scratch_free_gb="$(python3 - "${scratch_quota_gb}" "${scratch_used_gb}" <<'PYEOF'
import sys
print(float(sys.argv[1]) - float(sys.argv[2]))
PYEOF
)"

log "[INFO] /scratch quota usage: used=${scratch_used_raw}, quota=${scratch_quota_raw}, free=$(printf '%.2f' "${scratch_free_gb}")G"

if ! python3 - "${scratch_free_gb}" "${REQUIRED_DB_GB}" <<'PYEOF'
import sys
free = float(sys.argv[1])
required = float(sys.argv[2])
if free < required:
    raise SystemExit(1)
PYEOF
then
    log "[ERROR] Not enough /scratch space. Need at least ${REQUIRED_DB_GB}G free for the GTDB-Tk database."
    exit 1
fi

mkdir -p "${DB_DIR}"
mkdir -p "${OUT_DIR}"

log "[RUN] Downloading GTDB-Tk database"
download_ok=false
for url in "${PRIMARY_DB_URL}" "${MIRROR_DB_URL}"; do
    log "[INFO] Trying ${url}"
    if wget -O "${DB_ARCHIVE}" "${url}"; then
        download_ok=true
        break
    fi
done
[[ "${download_ok}" == true ]] || { log "[ERROR] Failed to download GTDB-Tk database from both URLs."; exit 1; }

log "[RUN] Extracting GTDB-Tk database into ${DB_DIR}"
tar xvzf "${DB_ARCHIVE}" -C "${DB_DIR}" > /dev/null
rm -f "${DB_ARCHIVE}"

if [[ "${RUN_CHECK_INSTALL}" == true ]]; then
    log "[RUN] gtdbtk check_install"
    GTDBTK_DATA_PATH="${DB_DIR}" \
        "${CONDA}" run --no-capture-output -n "${GTDBTK_ENV}" \
        gtdbtk check_install
fi

if [[ "${HAS_CPUS}" == false ]]; then
    GTDBTK_ARGS+=( --cpus "${SLURM_CPUS_PER_TASK:-16}" )
fi
if [[ "${HAS_TMPDIR}" == false ]]; then
    GTDBTK_ARGS+=( --tmpdir "${DB_ROOT}" )
fi
if [[ "${HAS_SCRATCH_DIR}" == false ]]; then
    GTDBTK_ARGS+=( --scratch_dir "${DB_ROOT}/scratch" )
fi

mkdir -p "${DB_ROOT}/scratch"

log "[RUN] gtdbtk classify_wf"
GTDBTK_DATA_PATH="${DB_DIR}" \
    "${CONDA}" run --no-capture-output -n "${GTDBTK_ENV}" \
    gtdbtk classify_wf "${GTDBTK_ARGS[@]}"

log "[COMPLETE] GTDB-Tk finished successfully"
