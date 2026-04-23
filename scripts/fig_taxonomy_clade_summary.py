#!/usr/bin/env python3
"""
Per-sample + assembler taxonomic summary at phylum, class, and family levels.

Outputs:
  - clade_summary_phylum.csv / _class.csv / _family.csv  (bin counts per taxon)
  - stacked_bars_phylum.png / _class.png / _family.png   (relative proportion per sample,
      faceted by assembler — standard metagenomics MAG survey figure)
  - heatmap_phylum.png   (presence/absence heatmap across samples × assemblers)

Aggregates across all four binners; bins are not deduplicated across binners.

Run via: ./scripts/make_figures.sh taxonomy clade_summary
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

# ── config ────────────────────────────────────────────────────────────────────
_binning_dir = os.environ.get("BINNING_DIR", "binning_output")
CSV = Path(_binning_dir) / "bin_analysis_results.csv"
TOP_PHYLUM  = 14
TOP_CLASS   = 14
TOP_FAMILY  = 16
# ─────────────────────────────────────────────────────────────────────────────

output_dir = os.environ.get("FIGURE_OUTPUT_DIR")
if not output_dir:
    print("Error: FIGURE_OUTPUT_DIR not set", file=sys.stderr)
    sys.exit(1)
os.makedirs(output_dir, exist_ok=True)

if not CSV.exists():
    print(f"Error: {CSV} not found. Run exec_bin_analysis.sh first.", file=sys.stderr)
    sys.exit(1)

# ── load & parse ──────────────────────────────────────────────────────────────
df = pd.read_csv(CSV)

def extract(lineage, prefix):
    m = re.search(rf"{prefix}([^;]+)", str(lineage))
    return m.group(1) if m else "Unclassified"

df["phylum"] = df["qc_lineage"].apply(lambda x: extract(x, "p__"))
df["class_"] = df["qc_lineage"].apply(lambda x: extract(x, "c__"))
df["family"] = df["qc_lineage"].apply(lambda x: extract(x, "f__"))

# Short sample labels: strip E17_ prefix and _L003 suffix
df["sample_label"] = (
    df["sample_base"]
    .str.replace(r"^E17_", "", regex=True)
    .str.replace(r"_L003$", "", regex=True)
)

ASSEMBLERS    = ["megahit_large", "megahit_sensitive", "spades"]
ASSEMBLER_LABELS = {"megahit_large": "MEGAHIT large", "megahit_sensitive": "MEGAHIT sensitive", "spades": "SPAdes"}
SAMPLES       = sorted(df["sample_label"].unique())

LEVELS = [
    ("phylum", "phylum",  "Phylum",  TOP_PHYLUM),
    ("class_", "class",   "Class",   TOP_CLASS),
    ("family", "family",  "Family",  TOP_FAMILY),
]

# ── export summary CSVs ───────────────────────────────────────────────────────
for col, tag, label, _ in LEVELS:
    grp = df.groupby(["sample_base", "assembler", col])
    summary = grp.agg(
        bin_count=("qc_k5dif", "count"),
        k5dif_mean=("qc_k5dif", "mean"),
        k5dif_median=("qc_k5dif", "median"),
        k5dif_min=("qc_k5dif", "min"),
        k5dif_max=("qc_k5dif", "max"),
    ).reset_index().rename(columns={col: label.lower()})
    for c in ["k5dif_mean", "k5dif_median", "k5dif_min", "k5dif_max"]:
        summary[c] = summary[c].round(5)
    summary = summary.sort_values(
        ["sample_base", "assembler", "bin_count"], ascending=[True, True, False]
    )
    out_csv = os.path.join(output_dir, f"clade_summary_{tag}.csv")
    summary.to_csv(out_csv, index=False)
    print(f"CSV saved: {out_csv}")

# ── helper: build grouped pivot ───────────────────────────────────────────────
def make_pivot(sub, col, top_n):
    top = df[col].value_counts().head(top_n).index.tolist()
    grouped = sub[col].where(sub[col].isin(top), other="Other")
    pivot = (
        pd.DataFrame({"sample_label": sub["sample_label"], "taxon": grouped})
        .groupby(["sample_label", "taxon"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=SAMPLES, fill_value=0)
    )
    col_order = [t for t in top if t in pivot.columns] + (["Other"] if "Other" in pivot.columns else [])
    return pivot[col_order]

# ── stacked bar charts ────────────────────────────────────────────────────────
# Following Parks et al. 2017 / Tully et al. 2018 style:
# relative proportion per sample, one panel per assembler, one figure per level.

TAB20 = plt.get_cmap("tab20")
TAB20B = plt.get_cmap("tab20b")

def build_palette(n):
    """Combine tab20 + tab20b for up to 40 distinct colours."""
    colors = [TAB20(i / 20) for i in range(20)] + [TAB20B(i / 20) for i in range(20)]
    return colors[:n]

x = np.arange(len(SAMPLES))

for col, tag, label, top_n in LEVELS:
    # global top taxa for a consistent palette across panels
    top_taxa = df[col].value_counts().head(top_n).index.tolist()
    all_groups = top_taxa + ["Other"]
    colors = build_palette(len(all_groups))
    color_map = {g: colors[i] for i, g in enumerate(all_groups)}

    fig, axes = plt.subplots(
        1, len(ASSEMBLERS),
        figsize=(18, 6),
        sharey=False,
        constrained_layout=False,
    )
    fig.suptitle(
        f"{label}-level taxonomic composition of bins per sample\n"
        f"(all binners combined; proportional)",
        fontsize=13, fontweight="bold", y=1.01,
    )

    for j, assembler in enumerate(ASSEMBLERS):
        ax = axes[j]
        sub = df[df["assembler"] == assembler]
        pivot = make_pivot(sub, col, top_n)

        # Convert to proportions row-wise
        row_totals = pivot.sum(axis=1).replace(0, np.nan)
        prop = pivot.div(row_totals, axis=0).fillna(0)

        bottoms = np.zeros(len(SAMPLES))
        for taxon in prop.columns:
            vals = prop[taxon].values
            ax.bar(x, vals, bottom=bottoms,
                   color=color_map.get(taxon, "grey"),
                   width=0.72, edgecolor="white", linewidth=0.3)
            bottoms += vals

        ax.set_title(ASSEMBLER_LABELS[assembler], fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.5, len(SAMPLES) - 0.5)
        ax.spines[["top", "right"]].set_visible(False)
        if j == 0:
            ax.set_ylabel("Proportion of bins", fontsize=10)
        else:
            ax.set_ylabel("")
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))

    legend_handles = [mpatches.Patch(color=color_map[g], label=g) for g in all_groups]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=min(6, len(all_groups)),
        fontsize=7.5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.18),
        title=label,
        title_fontsize=9,
    )

    plt.tight_layout(rect=[0, 0.12, 1, 1])
    out = os.path.join(output_dir, f"stacked_bars_{tag}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")

# ── phylum presence/absence heatmap ──────────────────────────────────────────
# Presence/absence heatmap: rows = phylum, cols = sample × assembler.
# Common in comparative metagenomics (e.g., Nayfach et al. 2021 Earth's microbiome).

top_phyla_global = df["phylum"].value_counts().head(TOP_PHYLUM).index.tolist()

rows = []
for assembler in ASSEMBLERS:
    for sample in SAMPLES:
        sub = df[(df["assembler"] == assembler) & (df["sample_label"] == sample)]
        detected = set(sub["phylum"].dropna().unique())
        for p in top_phyla_global:
            rows.append({
                "phylum": p,
                "sample": sample,
                "assembler": ASSEMBLER_LABELS[assembler],
                "bin_count": sub[sub["phylum"] == p].shape[0],
                "detected": int(p in detected),
            })

heat_df = pd.DataFrame(rows)

# Pivot: rows = phylum, cols = (assembler, sample)
pivot_heat = heat_df.pivot_table(
    index="phylum", columns=["assembler", "sample"],
    values="bin_count", aggfunc="sum", fill_value=0
)

# Order phyla by total bin count descending
phylum_order = pivot_heat.sum(axis=1).sort_values(ascending=False).index
pivot_heat = pivot_heat.loc[phylum_order]

fig_h, ax_h = plt.subplots(figsize=(20, 6))
sns.heatmap(
    pivot_heat,
    ax=ax_h,
    cmap="YlOrRd",
    linewidths=0.3,
    linecolor="white",
    cbar_kws={"label": "Bin count", "shrink": 0.6},
    xticklabels=True,
    yticklabels=True,
)
ax_h.set_title(
    f"Phylum-level bin counts across samples and assemblers\n"
    f"(top {TOP_PHYLUM} phyla by global bin count)",
    fontsize=12, fontweight="bold",
)
ax_h.set_xlabel("")
ax_h.set_ylabel("")
ax_h.tick_params(axis="x", labelsize=7.5, rotation=45)
ax_h.tick_params(axis="y", labelsize=8, rotation=0)

# Add vertical lines to separate assembler groups
n_samples = len(SAMPLES)
for i in range(1, len(ASSEMBLERS)):
    ax_h.axvline(x=i * n_samples, color="black", linewidth=1.5)

# Assembler group labels above the heatmap
for i, assembler in enumerate(ASSEMBLER_LABELS.values()):
    mid = i * n_samples + n_samples / 2 - 0.5
    ax_h.text(mid, -1.2, assembler, ha="center", va="bottom",
              fontsize=9, fontweight="bold", transform=ax_h.transData)

plt.tight_layout()
out_h = os.path.join(output_dir, "heatmap_phylum.png")
fig_h.savefig(out_h, dpi=200, bbox_inches="tight")
plt.close(fig_h)
print(f"Figure saved: {out_h}")

print("Done.")
