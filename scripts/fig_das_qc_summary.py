#!/usr/bin/env python3
"""
Quick visual summaries for DAS_Tool QC outputs.

Data source priority:
1. binning_output/das_bins_qc_summary.csv
2. Reconstruct a minimal merged table from per-sample das_tool/quickclade.tsv
   and per-bin BUSCO short summaries

Outputs:
  - das_bins_qc_summary_used.csv
  - busco_completeness_by_assembler.png
  - busco_quality_tiers_by_sample.png
  - quickclade_stacked_bars_phylum.png
  - quickclade_phylum_heatmap.png

Run via:
  ./scripts/make_figures.sh das qc_summary
"""

import csv
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
BINNING_DIR = PROJECT_ROOT / "binning_output"
SUMMARY_CSV = BINNING_DIR / "das_bins_qc_summary.csv"
SUFFIX_RE = re.compile(r"_(megahit_large|megahit_sensitive|spades)$")
TOP_PHYLUM = 12

ASSEMBLERS = ["megahit_large", "megahit_sensitive", "spades"]
ASSEMBLER_LABELS = {
    "megahit_large": "MEGAHIT large",
    "megahit_sensitive": "MEGAHIT sensitive",
    "spades": "SPAdes",
}

QC_COLS = ["qc_ref_name", "qc_taxid", "qc_level", "qc_k5dif", "qc_lineage"]
BUSCO_COLS = [
    "busco_lineage",
    "busco_complete_pct",
    "busco_single_pct",
    "busco_duplicated_pct",
    "busco_fragmented_pct",
    "busco_missing_pct",
    "busco_n_markers",
]
EMPTY_QC = {col: "" for col in QC_COLS}
EMPTY_BUSCO = {col: "" for col in BUSCO_COLS}

output_dir = os.environ.get("FIGURE_OUTPUT_DIR")
if not output_dir:
    print("Error: FIGURE_OUTPUT_DIR not set", file=sys.stderr)
    sys.exit(1)
output_dir = Path(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)


def parse_quickclade_tsv(tsv_path: Path):
    rows = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 16:
                continue
            bin_name = Path(cols[0]).stem
            rows[bin_name] = {
                "qc_ref_name": cols[4],
                "qc_taxid": cols[5],
                "qc_level": cols[9],
                "qc_k5dif": cols[14],
                "qc_lineage": cols[15],
            }
    return rows


def parse_busco_txt(txt_path: Path):
    metrics = EMPTY_BUSCO.copy()
    score_re = re.compile(
        r"C:(?P<C>[0-9.]+)%\[S:(?P<S>[0-9.]+)%,D:(?P<D>[0-9.]+)%\],"
        r"F:(?P<F>[0-9.]+)%,M:(?P<M>[0-9.]+)%,n:(?P<n>\d+)"
    )

    with open(txt_path) as fh:
        for raw in fh:
            line = raw.strip()
            match = score_re.search(line)
            if match:
                metrics.update({
                    "busco_complete_pct": match.group("C"),
                    "busco_single_pct": match.group("S"),
                    "busco_duplicated_pct": match.group("D"),
                    "busco_fragmented_pct": match.group("F"),
                    "busco_missing_pct": match.group("M"),
                    "busco_n_markers": match.group("n"),
                })
            if not metrics["busco_lineage"]:
                lineage_match = re.search(r"lineage dataset is:\s+([A-Za-z0-9_.-]+)", line)
                if lineage_match:
                    metrics["busco_lineage"] = lineage_match.group(1)
    return metrics


def parse_busco_dir(busco_bin_dir: Path):
    json_hits = sorted(busco_bin_dir.glob("**/short_summary*.json"))
    for json_path in json_hits:
        try:
            import json
            with open(json_path) as fh:
                data = json.load(fh)
            results = data.get("results", {}) or {}
            lineage = data.get("lineage_dataset", {}) or {}
            metrics = EMPTY_BUSCO.copy()
            if isinstance(lineage, dict):
                metrics["busco_lineage"] = lineage.get("name", "") or lineage.get("basename", "")
            if isinstance(results, dict):
                metrics["busco_complete_pct"] = str(results.get("Complete percentage", "") or results.get("Complete pct", ""))
                metrics["busco_single_pct"] = str(results.get("Single copy percentage", "") or results.get("Single pct", ""))
                metrics["busco_duplicated_pct"] = str(results.get("Multi copy percentage", "") or results.get("Duplicated pct", ""))
                metrics["busco_fragmented_pct"] = str(results.get("Fragmented percentage", "") or results.get("Fragmented pct", ""))
                metrics["busco_missing_pct"] = str(results.get("Missing percentage", "") or results.get("Missing pct", ""))
                metrics["busco_n_markers"] = str(results.get("n_markers", "") or results.get("Total BUSCO groups searched", ""))
            if any(metrics.values()):
                return metrics
        except Exception:
            continue

    txt_hits = sorted(busco_bin_dir.glob("**/short_summary*.txt"))
    if txt_hits:
        return parse_busco_txt(txt_hits[0])

    return EMPTY_BUSCO.copy()


def build_summary_from_raw():
    rows = []
    sample_dirs = sorted(BINNING_DIR.glob("E17_*"))
    for sample_dir in sample_dirs:
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name
        das_dir = sample_dir / "das_tool"
        bins_dir = das_dir / f"{sample}_DASTool_bins"
        if not bins_dir.is_dir():
            continue

        assembler_match = SUFFIX_RE.search(sample)
        assembler = assembler_match.group(1) if assembler_match else sample
        sample_base = SUFFIX_RE.sub("", sample)

        quickclade = {}
        qc_tsv = das_dir / "quickclade.tsv"
        if qc_tsv.is_file():
            quickclade = parse_quickclade_tsv(qc_tsv)

        busco_root = das_dir / "busco"
        busco = {}
        if busco_root.is_dir():
            for bin_dir in sorted(busco_root.iterdir()):
                if bin_dir.is_dir():
                    busco[bin_dir.name] = parse_busco_dir(bin_dir)

        for fasta_path in sorted(bins_dir.glob("*.fa")):
            bin_name = fasta_path.stem
            row = {
                "assembler": assembler,
                "sample_base": sample_base,
                "sample": sample,
                "bin_name": bin_name,
                "fasta_path": str(fasta_path),
            }
            row.update(EMPTY_QC)
            row.update(EMPTY_BUSCO)
            row.update(quickclade.get(bin_name, {}))
            row.update(busco.get(bin_name, {}))
            rows.append(row)

    return pd.DataFrame(rows)


def load_summary():
    if SUMMARY_CSV.is_file():
        df = pd.read_csv(SUMMARY_CSV)
        source = str(SUMMARY_CSV)
    else:
        df = build_summary_from_raw()
        source = "raw das_tool outputs"
    if df.empty:
        raise SystemExit("No DAS QC rows found in binning_output.")
    return df, source


def extract_taxon(lineage, prefix):
    match = re.search(rf"{prefix}([^;]+)", str(lineage))
    return match.group(1) if match else "Unclassified"


def make_busco_tier(value):
    if pd.isna(value):
        return "No BUSCO result"
    if value >= 90:
        return ">=90%"
    if value >= 70:
        return "70-89%"
    if value >= 50:
        return "50-69%"
    return "<50%"


df, data_source = load_summary()

for col in ["qc_k5dif", "busco_complete_pct", "busco_single_pct", "busco_duplicated_pct", "busco_fragmented_pct", "busco_missing_pct"]:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

df["sample_label"] = (
    df["sample_base"]
    .astype(str)
    .str.replace(r"^E17_", "", regex=True)
    .str.replace(r"_L003$", "", regex=True)
)
df["assembler_label"] = df["assembler"].map(ASSEMBLER_LABELS).fillna(df["assembler"])
df["phylum"] = df["qc_lineage"].apply(lambda x: extract_taxon(x, "p__"))
df["domain"] = df["qc_lineage"].apply(lambda x: extract_taxon(x, "(?:sk__|k__)"))
df["busco_tier"] = df["busco_complete_pct"].apply(make_busco_tier)

used_csv = output_dir / "das_bins_qc_summary_used.csv"
df.to_csv(used_csv, index=False, quoting=csv.QUOTE_MINIMAL)
print(f"Summary used for figures: {used_csv} (source: {data_source})")

samples = sorted(df["sample_label"].dropna().unique().tolist())


def ordered_assemblers(frame):
    present = frame["assembler"].dropna().unique().tolist()
    return [a for a in ASSEMBLERS if a in present] + sorted(a for a in present if a not in ASSEMBLERS)


def build_palette(n):
    tab20 = plt.get_cmap("tab20")
    tab20b = plt.get_cmap("tab20b")
    colors = [tab20(i / 20) for i in range(20)] + [tab20b(i / 20) for i in range(20)]
    return colors[:n]


busco_df = df.dropna(subset=["busco_complete_pct"]).copy()
if not busco_df.empty:
    fig, ax = plt.subplots(figsize=(9, 6))
    order = ordered_assemblers(busco_df)
    palette = plt.get_cmap("Set2")
    positions = np.arange(len(order))
    data = [busco_df.loc[busco_df["assembler"] == assembler, "busco_complete_pct"].dropna().to_numpy() for assembler in order]
    bp = ax.boxplot(data, positions=positions, patch_artist=True, widths=0.55, showfliers=False)
    for patch, idx in zip(bp["boxes"], range(len(order))):
        patch.set_facecolor(palette(idx / max(1, len(order) - 1)))
        patch.set_alpha(0.8)
    rng = np.random.default_rng(7)
    for pos, values in zip(positions, data):
        if len(values) == 0:
            continue
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(np.full(len(values), pos) + jitter, values, s=12, alpha=0.4, color="0.2", linewidths=0)
    ax.set_title("DAS_Tool bin BUSCO completeness by assembler", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("BUSCO complete (%)")
    ax.set_xticks(positions)
    ax.set_xticklabels([ASSEMBLER_LABELS.get(a, a) for a in order], rotation=15, ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8)
    plt.tight_layout()
    out = output_dir / "busco_completeness_by_assembler.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")

    tier_order = [">=90%", "70-89%", "50-69%", "<50%", "No BUSCO result"]
    colors = {
        ">=90%": "#2a9d8f",
        "70-89%": "#8ab17d",
        "50-69%": "#e9c46a",
        "<50%": "#e76f51",
        "No BUSCO result": "#bdbdbd",
    }
    assemblers = ordered_assemblers(df)
    fig, axes = plt.subplots(1, len(assemblers), figsize=(18, 5), sharey=True)
    if len(assemblers) == 1:
        axes = [axes]
    for ax, assembler in zip(axes, assemblers):
        sub = df[df["assembler"] == assembler].copy()
        pivot = (
            sub.groupby(["sample_label", "busco_tier"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=samples, fill_value=0)
        )
        for tier in tier_order:
            if tier not in pivot.columns:
                pivot[tier] = 0
        pivot = pivot[tier_order]
        x = np.arange(len(samples))
        bottom = np.zeros(len(samples))
        for tier in tier_order:
            vals = pivot[tier].to_numpy(dtype=float)
            ax.bar(x, vals, bottom=bottom, color=colors[tier], width=0.75, edgecolor="white", linewidth=0.4)
            bottom += vals
        ax.set_title(ASSEMBLER_LABELS.get(assembler, assembler), fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel("")
    axes[0].set_ylabel("Number of DAS bins")
    handles = [mpatches.Patch(color=colors[t], label=t) for t in tier_order]
    fig.suptitle("BUSCO completeness tiers across DAS_Tool bins", fontsize=13, fontweight="bold", y=1.02)
    fig.legend(handles=handles, loc="lower center", ncol=len(tier_order), frameon=False, bbox_to_anchor=(0.5, -0.05))
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    out = output_dir / "busco_quality_tiers_by_sample.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")


qc_df = df[df["qc_lineage"].fillna("").astype(str) != ""].copy()
if not qc_df.empty:
    top_phyla = qc_df["phylum"].value_counts().head(TOP_PHYLUM).index.tolist()
    all_groups = top_phyla + ["Other"]
    color_map = {group: color for group, color in zip(all_groups, build_palette(len(all_groups)))}

    assemblers = ordered_assemblers(qc_df)
    fig, axes = plt.subplots(1, len(assemblers), figsize=(18, 6), sharey=False)
    if len(assemblers) == 1:
        axes = [axes]
    x = np.arange(len(samples))

    for ax, assembler in zip(axes, assemblers):
        sub = qc_df[qc_df["assembler"] == assembler].copy()
        grouped = sub["phylum"].where(sub["phylum"].isin(top_phyla), other="Other")
        pivot = (
            pd.DataFrame({"sample_label": sub["sample_label"], "phylum_grouped": grouped})
            .groupby(["sample_label", "phylum_grouped"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=samples, fill_value=0)
        )
        col_order = [g for g in all_groups if g in pivot.columns]
        pivot = pivot[col_order]
        props = pivot.div(pivot.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        bottom = np.zeros(len(samples))
        for group in col_order:
            vals = props[group].to_numpy(dtype=float)
            ax.bar(x, vals, bottom=bottom, color=color_map[group], width=0.75, edgecolor="white", linewidth=0.35)
            bottom += vals
        ax.set_title(ASSEMBLER_LABELS.get(assembler, assembler), fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Proportion of DAS bins")
    fig.suptitle("QuickClade phylum composition of DAS_Tool bins", fontsize=13, fontweight="bold", y=1.02)
    legend_handles = [mpatches.Patch(color=color_map[g], label=g) for g in all_groups]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=min(5, len(all_groups)),
        frameon=False,
        bbox_to_anchor=(0.5, -0.12),
        fontsize=8,
        title="Phylum",
        title_fontsize=9,
    )
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    out = output_dir / "quickclade_stacked_bars_phylum.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")

    heat_rows = []
    for assembler in assemblers:
        for sample in samples:
            sub = qc_df[(qc_df["assembler"] == assembler) & (qc_df["sample_label"] == sample)]
            for phylum in top_phyla:
                heat_rows.append({
                    "assembler": ASSEMBLER_LABELS.get(assembler, assembler),
                    "sample": sample,
                    "phylum": phylum,
                    "bin_count": int((sub["phylum"] == phylum).sum()),
                })
    heat_df = pd.DataFrame(heat_rows)
    pivot_heat = heat_df.pivot_table(
        index="phylum",
        columns=["assembler", "sample"],
        values="bin_count",
        aggfunc="sum",
        fill_value=0,
    )
    pivot_heat = pivot_heat.loc[pivot_heat.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(20, 6))
    im = ax.imshow(pivot_heat.to_numpy(), aspect="auto", cmap="YlOrRd")
    ax.set_title("Top QuickClade phyla across samples and assemblers", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks(np.arange(pivot_heat.shape[1]))
    ax.set_xticklabels(
        [f"{assembler}\n{sample}" for assembler, sample in pivot_heat.columns],
        fontsize=7.5,
        rotation=45,
        ha="right",
    )
    ax.set_yticks(np.arange(pivot_heat.shape[0]))
    ax.set_yticklabels(pivot_heat.index.tolist(), fontsize=8)
    n_samples = len(samples)
    for i in range(1, len(assemblers)):
        ax.axvline(x=i * n_samples - 0.5, color="black", linewidth=1.1)
    cbar = fig.colorbar(im, ax=ax, shrink=0.65)
    cbar.set_label("Bin count")
    plt.tight_layout()
    out = output_dir / "quickclade_phylum_heatmap.png"
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")


print("✓ DAS QC figure generation successful")
