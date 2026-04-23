#!/usr/bin/bash -l

module load singularityce

######################################################
# Bin analysis: QuickClade taxonomy + assembly stats
# 1. QuickClade per bin
# 2. Assembly statistics via bin_stats.py
# 3. Merge into bin_analysis_results.csv
######################################################
# Usage: ./exec_bin_analysis.sh <binning_output_dir>
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
BBTOOLS_IMAGE="${BBTOOLS_IMAGE:-bbtools_latest.sif}"
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
BINNERS="${BINNERS:-maxbin2,concoct,metabat2,quickbin}"
RUN_STATS="${RUN_STATS:-true}"
RUN_MERGE="${RUN_MERGE:-true}"

[[ -f "${BBTOOLS_IMAGE}" ]] || { log "[ERROR] Image not found: ${BBTOOLS_IMAGE}"; exit 1; }

declare -A BINNER_SUBDIR=(
    [maxbin2]="maxbin2"
    [concoct]="concoct/fasta_bins"
    [metabat2]="metabat2"
    [quickbin]="quickbin"
)
declare -A BINNER_GLOB=(
    [maxbin2]="*.fasta"
    [concoct]="*.fa"
    [metabat2]="bin.*.fa"
    [quickbin]="bin*.fa"
)
IFS=',' read -r -a ACTIVE_BINNERS <<< "${BINNERS}"

# Phase 1: QuickClade per sample per binner

for sample_dir in "${BINNING_DIR}"/*/; do
    [[ -d "${sample_dir}" ]] || continue
    sample="$(basename "${sample_dir}")"

    for binner in "${ACTIVE_BINNERS[@]}"; do
        [[ -n "${BINNER_SUBDIR[${binner}]:-}" ]] || { log "[ERROR] Unknown binner: ${binner}"; exit 1; }
        bin_dir="${sample_dir}${BINNER_SUBDIR[${binner}]}"
        [[ -d "${bin_dir}" ]] || continue

        qc_out="${sample_dir}${binner}/quickclade.tsv"
        if [[ -f "${qc_out}" ]]; then
            log "[SKIP] QuickClade ${sample}/${binner}"
            continue
        fi

        mapfile -t bin_files < <(find "${bin_dir}" -maxdepth 1 -name "${BINNER_GLOB[${binner}]}" 2>/dev/null | sort)
        [[ ${#bin_files[@]} -eq 0 ]] && continue

        container_in="$(printf '%s\n' "${bin_files[@]}" | \
            sed "s|^${BINNING_DIR}|/bins|" | paste -sd,)"

        log "[RUN] QuickClade: ${sample}/${binner} (${#bin_files[@]} bins)"

        singularity exec --cleanenv \
            --bind "${BINNING_DIR}:/bins" \
            "${BBTOOLS_IMAGE}" \
            quickclade.sh \
            in="${container_in}" \
            oneline \
            out="${qc_out/#${BINNING_DIR}//bins}" \
            2> "${sample_dir}${binner}/quickclade.err"
    done
done

# Phase 2: Assembly statistics

if [[ "${RUN_STATS}" == true ]]; then
    log "[RUN] Assembly statistics"
    python3 "${SCRIPT_DIR}/bin_stats.py" "${BINNING_DIR}" "${BINNING_DIR}/bin_stats_summary.csv"
fi

# Phase 3: Merge quickclade + stats into final CSV

if [[ "${RUN_MERGE}" == true ]]; then
log "[RUN] Merging results"
python3 - "${BINNING_DIR}" <<'PYEOF'
import csv, glob, os, re, sys

BINNING_DIR = sys.argv[1]
SUFFIX_RE   = re.compile(r"_(megahit_large|megahit_sensitive|spades)$")

ASSEMBLER_ORDER = {"megahit_large": 0, "megahit_sensitive": 1, "spades": 2}
BINNER_ORDER    = {"maxbin2": 0, "concoct": 1, "metabat2": 2, "quickbin": 3}

QC_COLS = ["qc_ref_name", "qc_taxid", "qc_level", "qc_k5dif", "qc_lineage"]
EMPTY_QC = {f: "" for f in QC_COLS}

# Load quickclade TSVs
# columns: QueryName Q_GC Q_Bases Q_Contigs RefName R_TaxID R_GC R_Bases R_Contigs R_Level GCdif STRdif k3dif k4dif k5dif lineage
qc_data = {}
for sample_dir in sorted(glob.glob(os.path.join(BINNING_DIR, "*"))):
    if not os.path.isdir(sample_dir):
        continue
    sample = os.path.basename(sample_dir)
    for binner in ("maxbin2", "concoct", "metabat2", "quickbin"):
        tsv = os.path.join(sample_dir, binner, "quickclade.tsv")
        if not os.path.isfile(tsv):
            continue
        with open(tsv) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 15:
                    continue
                bin_name = os.path.splitext(os.path.basename(cols[0]))[0]
                qc_data[(sample, binner, bin_name)] = {
                    "qc_ref_name": cols[4],
                    "qc_taxid":    cols[5],
                    "qc_level":    cols[9],
                    "qc_k5dif":    cols[14],
                    "qc_lineage":  cols[15] if len(cols) > 15 else "",
                }

# Load bin stats
stats_csv = os.path.join(BINNING_DIR, "bin_stats_summary.csv")
with open(stats_csv) as fh:
    rows = list(csv.DictReader(fh))

# Merge and sort
out_rows = []
for r in rows:
    m           = SUFFIX_RE.search(r["sample"])
    assembler   = m.group(1) if m else r["sample"]
    sample_base = SUFFIX_RE.sub("", r["sample"])
    qc          = qc_data.get((r["sample"], r["binner"], r["bin_name"]), EMPTY_QC)
    out_rows.append({"assembler": assembler, "sample_base": sample_base, **r, **qc})

out_rows.sort(key=lambda r: (
    ASSEMBLER_ORDER.get(r["assembler"], 99),
    r["sample_base"],
    BINNER_ORDER.get(r["binner"], 99),
    r["bin_name"],
))

FIELDS = [
    "assembler", "sample_base", "sample", "binner", "bin_name",
    "num_contigs", "total_length", "largest_contig",
    "mean_length", "median_length", "n50", "l50", "n90", "l90", "gc_pct",
    "qc_ref_name", "qc_taxid", "qc_level", "qc_k5dif", "qc_lineage",
]

out_csv = os.path.join(BINNING_DIR, "bin_analysis_results.csv")
with open(out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
    w.writeheader()
    w.writerows(out_rows)

matched = sum(1 for r in out_rows if r["qc_ref_name"])
print(f"Wrote {len(out_rows)} rows ({matched} with QuickClade hits) → {out_csv}")
PYEOF

log "Done. Results: ${BINNING_DIR}/bin_analysis_results.csv"
else
    log "Done. QuickClade-only run."
fi
