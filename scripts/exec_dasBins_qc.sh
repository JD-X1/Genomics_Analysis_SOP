#!/usr/bin/bash -l

module load singularityce

######################################################
# DAS_Tool bin QC
# 1. QuickClade taxonomy for DAS_Tool bins per sample
# 2. BUSCO per DAS_Tool bin
# 3. Merge QC results into das_bins_qc_summary.csv
######################################################
# Usage: ./exec_dasBins_qc.sh <binning_output_dir>
######################################################

set -Eeuo pipefail
shopt -s nullglob

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec 2>&1
    echo "[INFO] SLURM Job ID: ${SLURM_JOB_ID} running on $(hostname)"
fi

log() { echo "[$(date +'%F %T')] $*"; }
trap 'rc=$?; log "[ERROR] line ${LINENO}: ${BASH_COMMAND} (exit=${rc})"; exit ${rc}' ERR

[[ $# -lt 1 ]] && { log "[ERROR] Usage: $0 <binning_output_dir>"; exit 2; }

BINNING_DIR="$(realpath -m "$1")"
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
CONDA="${CONDA:-/home/jduque2/miniforge3/bin/conda}"
BBTOOLS_IMAGE="${BBTOOLS_IMAGE:-bbtools_latest.sif}"
BUSCO_ENV="${BUSCO_ENV:-busco}"
BUSCO_PYTHON="${BUSCO_PYTHON:-/home/jduque2/miniforge3/envs/${BUSCO_ENV}/bin/python}"
BUSCO_MODE="${BUSCO_MODE:-genome}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-}"
BUSCO_BACTERIA_LINEAGE="${BUSCO_BACTERIA_LINEAGE:-bacteria_odb12}"
BUSCO_ARCHAEA_LINEAGE="${BUSCO_ARCHAEA_LINEAGE:-archaea_odb12}"
BUSCO_DOWNLOADS="${BUSCO_DOWNLOADS:-${SCRIPT_DIR}/../busco_downloads}"
BUSCO_ARGS="${BUSCO_ARGS:-}"
RUN_QUICKCLADE="${RUN_QUICKCLADE:-true}"
RUN_BUSCO="${RUN_BUSCO:-true}"
RUN_MERGE="${RUN_MERGE:-true}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MAX_PARALLEL="${MAX_PARALLEL:-1}"
THREADS_PER_SAMPLE=$(( THREADS / MAX_PARALLEL ))

[[ -d "${BINNING_DIR}" ]] || { log "[ERROR] Binning output dir not found: ${BINNING_DIR}"; exit 1; }
[[ -f "${BBTOOLS_IMAGE}" ]] || { log "[ERROR] Image not found: ${BBTOOLS_IMAGE}"; exit 1; }
[[ -x "${BUSCO_PYTHON}" ]] || { log "[ERROR] BUSCO python not executable: ${BUSCO_PYTHON}"; exit 1; }

mapfile -t SAMPLE_DIRS < <(
    find "${BINNING_DIR}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort
)

if [[ ${#SAMPLE_DIRS[@]} -eq 0 ]]; then
    log "[ERROR] No sample directories found in ${BINNING_DIR}"
    exit 1
fi

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
    local bins_dir="${das_dir}/${sample}_DASTool_bins"

    [[ -d "${das_dir}" ]] || { log "[SKIP] ${sample}: missing das_tool dir"; return 0; }
    [[ -f "${das_dir}/.done" ]] || { log "[SKIP] ${sample}: DAS_Tool not complete"; return 0; }
    [[ -d "${bins_dir}" ]] || { log "[SKIP] ${sample}: missing DAS_Tool bins dir"; return 0; }

    mapfile -t bin_files < <(find "${bins_dir}" -maxdepth 1 -name '*.fa' 2>/dev/null | sort)
    [[ ${#bin_files[@]} -gt 0 ]] || { log "[SKIP] ${sample}: no DAS_Tool bin FASTAs"; return 0; }

    local qc_out="${das_dir}/quickclade.tsv"
    local qc_err="${das_dir}/quickclade.err"
    if [[ "${RUN_QUICKCLADE}" == true || "${RUN_BUSCO}" == true ]]; then
        if [[ -s "${qc_out}" ]]; then
            log "[SKIP] QuickClade ${sample}/das_tool"
        else
            local container_in
            container_in="$(
                printf '%s\n' "${bin_files[@]}" | \
                sed "s|^${BINNING_DIR}|/bins|" | paste -sd,
            )"
            log "[RUN] QuickClade: ${sample}/das_tool (${#bin_files[@]} bins)"
            singularity exec --cleanenv \
                --bind "${BINNING_DIR}:/bins" \
                "${BBTOOLS_IMAGE}" \
                quickclade.sh \
                in="${container_in}" \
                oneline \
                out="${qc_out/#${BINNING_DIR}//bins}" \
                2> "${qc_err}"
        fi
    fi

    if [[ "${RUN_BUSCO}" == true ]]; then
        local busco_root="${das_dir}/busco"
        mkdir -p "${busco_root}"

        local -a busco_extra=()
        if [[ -n "${BUSCO_ARGS}" ]]; then
            read -r -a busco_extra <<< "${BUSCO_ARGS}"
        fi

        local qc_map
        qc_map="$(mktemp)"
        if [[ -s "${qc_out}" ]]; then
            python3 - "${qc_out}" > "${qc_map}" <<'PYEOF'
import os
import sys

qc_tsv = sys.argv[1]
with open(qc_tsv) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 16:
            continue
        name = os.path.splitext(os.path.basename(cols[0]))[0]
        lineage = cols[15]
        domain = ""
        if "sk__Archaea" in lineage or "k__Archaea" in lineage:
            domain = "archaea"
        elif "sk__Bacteria" in lineage or "k__Bacteria" in lineage:
            domain = "bacteria"
        print(f"{name}\t{domain}")
PYEOF
        fi

        local fasta bin_name bin_work busco_done busco_log busco_err busco_rc busco_lineage busco_domain
        for fasta in "${bin_files[@]}"; do
            bin_name="$(basename "${fasta}" .fa)"
            bin_work="${busco_root}/${bin_name}"
            busco_done="${bin_work}/.done"
            busco_log="${bin_work}/busco.log"
            busco_err="${bin_work}/busco.err"

            if [[ -f "${busco_done}" ]]; then
                log "[SKIP] BUSCO ${sample}/${bin_name}"
                continue
            fi

            mkdir -p "${bin_work}"
            log "[RUN] BUSCO: ${sample}/${bin_name}"

            busco_domain=""
            if [[ -s "${qc_map}" ]]; then
                busco_domain="$(awk -F'\t' -v name="${bin_name}" '$1==name {print $2; exit}' "${qc_map}")"
            fi

            if [[ -n "${BUSCO_LINEAGE}" ]]; then
                busco_lineage="${BUSCO_LINEAGE}"
            elif [[ "${busco_domain}" == "archaea" ]]; then
                busco_lineage="${BUSCO_ARCHAEA_LINEAGE}"
            elif [[ "${busco_domain}" == "bacteria" ]]; then
                busco_lineage="${BUSCO_BACTERIA_LINEAGE}"
            else
                busco_lineage=""
            fi

            busco_rc=0
            if [[ -n "${busco_lineage}" ]]; then
                "${BUSCO_PYTHON}" -m busco.run_BUSCO \
                    -i "${fasta}" \
                    -o "${bin_name}" \
                    -m "${BUSCO_MODE}" \
                    -l "${busco_lineage}" \
                    --download_path "${BUSCO_DOWNLOADS}" \
                    -c "${THREADS_PER_SAMPLE}" \
                    --out_path "${bin_work}" \
                    "${busco_extra[@]}" \
                    1> "${busco_log}" \
                    2> "${busco_err}" || busco_rc=$?
            else
                "${BUSCO_PYTHON}" -m busco.run_BUSCO \
                    -i "${fasta}" \
                    -o "${bin_name}" \
                    -m "${BUSCO_MODE}" \
                    --auto-lineage-prok \
                    --download_path "${BUSCO_DOWNLOADS}" \
                    -c "${THREADS_PER_SAMPLE}" \
                    --out_path "${bin_work}" \
                    "${busco_extra[@]}" \
                    1> "${busco_log}" \
                    2> "${busco_err}" || busco_rc=$?
            fi

            if [[ "${busco_rc}" -ne 0 ]]; then
                log "[ERROR] BUSCO failed for ${sample}/${bin_name} (exit=${busco_rc}); see ${busco_err}"
                return "${busco_rc}"
            fi

            touch "${busco_done}"
        done
        rm -f "${qc_map}"
    fi

    log "[COMPLETE] ${sample}"
}

export -f process_sample log
export BINNING_DIR CONDA BBTOOLS_IMAGE BUSCO_ENV BUSCO_PYTHON BUSCO_MODE BUSCO_LINEAGE BUSCO_BACTERIA_LINEAGE BUSCO_ARCHAEA_LINEAGE BUSCO_DOWNLOADS BUSCO_ARGS
export RUN_QUICKCLADE RUN_BUSCO THREADS_PER_SAMPLE

printf '%s\0' "${SAMPLES[@]}" | \
    xargs -0 -n 1 -P "${MAX_PARALLEL}" bash -c 'process_sample "$1"' _

if [[ "${RUN_MERGE}" == true ]]; then
    log "[RUN] Merging BUSCO + QuickClade results"
    python3 - "${BINNING_DIR}" <<'PYEOF'
import csv
import glob
import json
import os
import re
import sys

BINNING_DIR = sys.argv[1]
SUFFIX_RE = re.compile(r"_(megahit_large|megahit_sensitive|spades)$")

QC_COLS = ["qc_ref_name", "qc_taxid", "qc_level", "qc_k5dif", "qc_lineage"]
EMPTY_QC = {f: "" for f in QC_COLS}

BUSCO_COLS = [
    "busco_lineage",
    "busco_complete_pct",
    "busco_single_pct",
    "busco_duplicated_pct",
    "busco_fragmented_pct",
    "busco_missing_pct",
    "busco_n_markers",
]
EMPTY_BUSCO = {f: "" for f in BUSCO_COLS}


def parse_quickclade_tsv(tsv_path):
    rows = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 15:
                continue
            bin_name = os.path.splitext(os.path.basename(cols[0]))[0]
            rows[bin_name] = {
                "qc_ref_name": cols[4],
                "qc_taxid": cols[5],
                "qc_level": cols[9],
                "qc_k5dif": cols[14],
                "qc_lineage": cols[15] if len(cols) > 15 else "",
            }
    return rows


def parse_busco_txt(txt_path):
    metrics = EMPTY_BUSCO.copy()
    lineage = ""
    score_re = re.compile(
        r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],"
        r"F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>\d+)"
    )

    with open(txt_path) as fh:
        for raw in fh:
            line = raw.strip()
            m = score_re.search(line)
            if m:
                metrics.update({
                    "busco_complete_pct": m.group("C"),
                    "busco_single_pct": m.group("S"),
                    "busco_duplicated_pct": m.group("D"),
                    "busco_fragmented_pct": m.group("F"),
                    "busco_missing_pct": m.group("M"),
                    "busco_n_markers": m.group("n"),
                })
            if not lineage:
                lm = re.search(r"lineage dataset\s+([A-Za-z0-9_.-]+)", line)
                if lm:
                    lineage = lm.group(1)
    metrics["busco_lineage"] = lineage
    return metrics


def parse_busco_json(json_path):
    metrics = EMPTY_BUSCO.copy()
    with open(json_path) as fh:
        data = json.load(fh)

    lineage = data.get("lineage_dataset", {}) or {}
    if isinstance(lineage, dict):
        metrics["busco_lineage"] = lineage.get("name", "") or lineage.get("creation_date", "")

    results = data.get("results", {}) or {}
    if isinstance(results, dict):
        metrics["busco_complete_pct"] = str(results.get("Complete percentage", "") or results.get("Complete pct", ""))
        metrics["busco_single_pct"] = str(results.get("Single copy percentage", "") or results.get("Single pct", ""))
        metrics["busco_duplicated_pct"] = str(results.get("Multi copy percentage", "") or results.get("Duplicated pct", ""))
        metrics["busco_fragmented_pct"] = str(results.get("Fragmented percentage", "") or results.get("Fragmented pct", ""))
        metrics["busco_missing_pct"] = str(results.get("Missing percentage", "") or results.get("Missing pct", ""))
        metrics["busco_n_markers"] = str(results.get("n_markers", "") or results.get("Total BUSCO groups searched", ""))
    return metrics


def parse_busco_dir(busco_bin_dir):
    json_hits = sorted(glob.glob(os.path.join(busco_bin_dir, "**", "short_summary*.json"), recursive=True))
    if json_hits:
        parsed = parse_busco_json(json_hits[0])
        if any(parsed.values()):
            return parsed

    txt_hits = sorted(glob.glob(os.path.join(busco_bin_dir, "**", "short_summary*.txt"), recursive=True))
    if txt_hits:
        return parse_busco_txt(txt_hits[0])

    return EMPTY_BUSCO.copy()


quickclade = {}
busco = {}

for sample_dir in sorted(glob.glob(os.path.join(BINNING_DIR, "*"))):
    if not os.path.isdir(sample_dir):
        continue
    sample = os.path.basename(sample_dir)
    das_dir = os.path.join(sample_dir, "das_tool")
    qc_tsv = os.path.join(das_dir, "quickclade.tsv")
    if os.path.isfile(qc_tsv):
        quickclade[sample] = parse_quickclade_tsv(qc_tsv)

    busco_root = os.path.join(das_dir, "busco")
    if os.path.isdir(busco_root):
        per_sample = {}
        for bin_dir in sorted(glob.glob(os.path.join(busco_root, "*"))):
            if not os.path.isdir(bin_dir):
                continue
            bin_name = os.path.basename(bin_dir)
            per_sample[bin_name] = parse_busco_dir(bin_dir)
        busco[sample] = per_sample

rows = []
for sample_dir in sorted(glob.glob(os.path.join(BINNING_DIR, "*"))):
    if not os.path.isdir(sample_dir):
        continue
    sample = os.path.basename(sample_dir)
    das_dir = os.path.join(sample_dir, "das_tool")
    bins_dir = os.path.join(das_dir, f"{sample}_DASTool_bins")
    if not os.path.isdir(bins_dir):
        continue

    assembler_match = SUFFIX_RE.search(sample)
    assembler = assembler_match.group(1) if assembler_match else sample
    sample_base = SUFFIX_RE.sub("", sample)

    for fasta in sorted(glob.glob(os.path.join(bins_dir, "*.fa"))):
        bin_name = os.path.splitext(os.path.basename(fasta))[0]
        row = {
            "sample": sample,
            "sample_base": sample_base,
            "assembler": assembler,
            "bin_name": bin_name,
            "fasta_path": fasta,
        }
        row.update(quickclade.get(sample, {}).get(bin_name, EMPTY_QC))
        row.update(busco.get(sample, {}).get(bin_name, EMPTY_BUSCO))
        rows.append(row)

rows.sort(key=lambda r: (r["assembler"], r["sample_base"], r["bin_name"]))

fields = [
    "assembler", "sample_base", "sample", "bin_name", "fasta_path",
    "qc_ref_name", "qc_taxid", "qc_level", "qc_k5dif", "qc_lineage",
    "busco_lineage", "busco_complete_pct", "busco_single_pct",
    "busco_duplicated_pct", "busco_fragmented_pct", "busco_missing_pct",
    "busco_n_markers",
]

out_csv = os.path.join(BINNING_DIR, "das_bins_qc_summary.csv")
with open(out_csv, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

qc_hits = sum(1 for r in rows if r["qc_ref_name"])
busco_hits = sum(1 for r in rows if r["busco_complete_pct"])
print(f"Wrote {len(rows)} rows ({qc_hits} QuickClade hits, {busco_hits} BUSCO hits) -> {out_csv}")
PYEOF
    log "Done. Results: ${BINNING_DIR}/das_bins_qc_summary.csv"
else
    log "Done. QC execution only."
fi
