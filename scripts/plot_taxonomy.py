#!/usr/bin/env python3
"""
Phylum-level taxonomic breakdown of bins per assembly+binner combination.
Produces a faceted stacked bar chart: rows=binner, cols=assembler, x=sample, stack=phylum.
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# ── config ────────────────────────────────────────────────────────────────────
_binning_dir = os.environ.get("BINNING_DIR", "binning_output")
CSV = Path(_binning_dir) / "bin_analysis_results.csv"
TOP_N = 14          # phyla shown individually; remainder → "Other"
MIN_TOTAL_LEN = 0   # set > 0 to filter low-quality bins by total_length
# ─────────────────────────────────────────────────────────────────────────────

output_dir = os.environ.get("FIGURE_OUTPUT_DIR")
if not output_dir:
    print("Error: FIGURE_OUTPUT_DIR environment variable not set", file=sys.stderr)
    sys.exit(1)
os.makedirs(output_dir, exist_ok=True)

if not CSV.exists():
    print(f"Error: {CSV} not found. Run exec_bin_analysis.sh first.", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(CSV)

if MIN_TOTAL_LEN:
    df = df[df["total_length"] >= MIN_TOTAL_LEN]

# Extract phylum from lineage string
df["phylum"] = df["qc_lineage"].str.extract(r"p__([^;]+)")
df["phylum"] = df["phylum"].fillna("Unclassified")

# Determine top N phyla globally by bin count
top_phyla = (
    df["phylum"].value_counts().head(TOP_N).index.tolist()
)
df["phylum_grouped"] = df["phylum"].where(df["phylum"].isin(top_phyla), other="Other")

# Short sample labels (strip common prefix E17_ and _L003 suffix)
df["sample_label"] = (
    df["sample_base"]
    .str.replace(r"^E17_", "", regex=True)
    .str.replace(r"_L003$", "", regex=True)
)

BINNERS    = ["maxbin2", "concoct", "metabat2", "quickbin"]
ASSEMBLERS = ["megahit_large", "megahit_sensitive", "spades"]
SAMPLES    = sorted(df["sample_label"].unique())

# Consistent colour palette
all_groups = top_phyla + ["Other"]
cmap = plt.get_cmap("tab20", len(all_groups))
color_map = {g: cmap(i) for i, g in enumerate(all_groups)}

# ── figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(
    nrows=len(BINNERS),
    ncols=len(ASSEMBLERS),
    figsize=(22, 18),
    sharey=False,
    sharex=False,
)
fig.suptitle("Phylum-level taxonomic breakdown of bins\nper assembly + binner combination",
             fontsize=14, fontweight="bold", y=0.98)

for row_i, binner in enumerate(BINNERS):
    for col_j, assembler in enumerate(ASSEMBLERS):
        ax = axes[row_i][col_j]

        sub = df[(df["binner"] == binner) & (df["assembler"] == assembler)]

        if sub.empty:
            ax.set_visible(False)
            continue

        # Pivot: rows=sample_label, cols=phylum_grouped, values=bin count
        pivot = (
            sub.groupby(["sample_label", "phylum_grouped"])
            .size()
            .unstack(fill_value=0)
            .reindex(index=SAMPLES, fill_value=0)
        )
        # Order columns: top phyla first in global rank order, then Other
        col_order = [p for p in all_groups if p in pivot.columns]
        pivot = pivot[col_order]

        x = np.arange(len(SAMPLES))
        bottoms = np.zeros(len(SAMPLES))

        for group in col_order:
            vals = pivot[group].values.astype(float)
            ax.bar(x, vals, bottom=bottoms, color=color_map[group],
                   width=0.75, edgecolor="white", linewidth=0.4)
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

        if row_i == 0:
            ax.set_title(assembler.replace("_", " "), fontsize=9, fontweight="bold")
        if col_j == 0:
            ax.set_ylabel(binner, fontsize=9, fontweight="bold")

# Shared legend below the grid
legend_handles = [
    mpatches.Patch(color=color_map[g], label=g) for g in all_groups
]
fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=5,
    fontsize=7.5,
    frameon=False,
    bbox_to_anchor=(0.5, -0.01),
    title="Phylum",
    title_fontsize=8,
)

plt.tight_layout(rect=[0, 0.06, 1, 0.97])
out_path = os.path.join(output_dir, "taxonomy_breakdown.png")
fig.savefig(out_path, bbox_inches="tight", dpi=150)
print(f"Saved: {out_path}")
print("Figure generation successful")
