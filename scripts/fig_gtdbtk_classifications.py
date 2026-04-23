#!/usr/bin/env python3
"""
GTDB-Tk taxonomic classification summary at phylum and family levels.

Data source: gtdbtk_analysis/ CSVs produced by parse_gtdbtk_results.py

Outputs:
  - stacked_bars_phylum.png   (proportional, one panel per assembler)
  - stacked_bars_family.png
  - heatmap_phylum.png        (bin counts, phylum × sample × assembler)

Run via: ./scripts/make_figures.sh gtdbtk classifications
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
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GTDBTK_DIR = PROJECT_ROOT / "gtdbtk_analysis"
TOP_PHYLUM = 14
TOP_FAMILY = 16

ASSEMBLERS = ["megahit_large", "megahit_sensitive", "spades"]
ASSEMBLER_LABELS = {
    "megahit_large": "MEGAHIT large",
    "megahit_sensitive": "MEGAHIT sensitive",
    "spades": "SPAdes",
}

output_dir = os.environ.get("FIGURE_OUTPUT_DIR")
if not output_dir:
    print("Error: FIGURE_OUTPUT_DIR not set", file=sys.stderr)
    sys.exit(1)
output_dir = Path(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

phylum_csv = GTDBTK_DIR / "gtdbtk_phylum_summary.csv"
family_csv = GTDBTK_DIR / "gtdbtk_family_summary.csv"

if not phylum_csv.exists():
    print(f"Error: {phylum_csv} not found. Run parse_gtdbtk_results.py first.", file=sys.stderr)
    sys.exit(1)
if not family_csv.exists():
    print(f"Error: {family_csv} not found. Run parse_gtdbtk_results.py first.", file=sys.stderr)
    sys.exit(1)

phylum_df = pd.read_csv(phylum_csv)
family_df = pd.read_csv(family_csv)


def clean_sample(s):
    s = re.sub(r"^E17_", "", str(s))
    s = re.sub(r"_L003$", "", s)
    return s


phylum_df["sample_label"] = phylum_df["sample"].apply(clean_sample)
family_df["sample_label"] = family_df["sample"].apply(clean_sample)

SAMPLES = sorted(set(phylum_df["sample_label"].unique()) | set(family_df["sample_label"].unique()))

TAB20 = plt.get_cmap("tab20")
TAB20B = plt.get_cmap("tab20b")


def build_palette(n):
    colors = [TAB20(i / 20) for i in range(20)] + [TAB20B(i / 20) for i in range(20)]
    return colors[:n]


def stacked_bars(df, taxon_col, top_n, level_label, out_stem):
    top_taxa = df[taxon_col].value_counts().head(top_n).index.tolist()
    all_groups = top_taxa + ["Other"]
    colors = build_palette(len(all_groups))
    color_map = {g: colors[i] for i, g in enumerate(all_groups)}

    assemblers_present = [a for a in ASSEMBLERS if a in df["assembly"].values]
    x = np.arange(len(SAMPLES))

    fig, axes = plt.subplots(1, len(assemblers_present), figsize=(18, 6), sharey=False)
    if len(assemblers_present) == 1:
        axes = [axes]
    fig.suptitle(
        f"{level_label}-level GTDB-Tk classification of bins per sample\n(proportional)",
        fontsize=13, fontweight="bold", y=1.01,
    )

    for ax, assembler in zip(axes, assemblers_present):
        sub = df[df["assembly"] == assembler].copy()
        sub["taxon_grouped"] = sub[taxon_col].where(sub[taxon_col].isin(top_taxa), other="Other")
        pivot = (
            sub.groupby(["sample_label", "taxon_grouped"])["count"]
            .sum()
            .unstack(fill_value=0)
            .reindex(index=SAMPLES, fill_value=0)
        )
        col_order = [g for g in all_groups if g in pivot.columns]
        pivot = pivot[col_order]
        row_totals = pivot.sum(axis=1).replace(0, np.nan)
        prop = pivot.div(row_totals, axis=0).fillna(0)

        bottom = np.zeros(len(SAMPLES))
        for taxon in col_order:
            vals = prop[taxon].to_numpy(dtype=float)
            ax.bar(x, vals, bottom=bottom, color=color_map[taxon], width=0.72,
                   edgecolor="white", linewidth=0.3)
            bottom += vals

        ax.set_title(ASSEMBLER_LABELS.get(assembler, assembler), fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.5, len(SAMPLES) - 0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))

    axes[0].set_ylabel("Proportion of bins", fontsize=10)

    legend_handles = [mpatches.Patch(color=color_map[g], label=g) for g in all_groups if g in color_map]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=min(6, len(legend_handles)),
        fontsize=7.5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.18),
        title=level_label,
        title_fontsize=9,
    )
    plt.tight_layout(rect=[0, 0.12, 1, 1])
    out = output_dir / f"stacked_bars_{out_stem}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {out}")


stacked_bars(phylum_df, "phylum", TOP_PHYLUM, "Phylum", "phylum")
stacked_bars(family_df, "family", TOP_FAMILY, "Family", "family")

# Phylum heatmap (bin counts, rows = phylum, cols = assembler × sample)
top_phyla = phylum_df["phylum"].value_counts().head(TOP_PHYLUM).index.tolist()
assemblers_present = [a for a in ASSEMBLERS if a in phylum_df["assembly"].values]

heat_rows = []
for assembler in assemblers_present:
    for sample in SAMPLES:
        sub = phylum_df[(phylum_df["assembly"] == assembler) & (phylum_df["sample_label"] == sample)]
        for p in top_phyla:
            heat_rows.append({
                "phylum": p,
                "sample": sample,
                "assembler": ASSEMBLER_LABELS.get(assembler, assembler),
                "bin_count": int(sub.loc[sub["phylum"] == p, "count"].sum()),
            })

heat_df = pd.DataFrame(heat_rows)
pivot_heat = heat_df.pivot_table(
    index="phylum", columns=["assembler", "sample"],
    values="bin_count", aggfunc="sum", fill_value=0,
)
pivot_heat = pivot_heat.loc[pivot_heat.sum(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(20, 6))
im = ax.imshow(pivot_heat.to_numpy(), aspect="auto", cmap="YlOrRd")
ax.set_title(
    f"GTDB-Tk phylum-level bin counts across samples and assemblers\n(top {TOP_PHYLUM} phyla)",
    fontsize=12, fontweight="bold",
)
ax.set_xticks(np.arange(pivot_heat.shape[1]))
ax.set_xticklabels(
    [f"{assembler}\n{sample}" for assembler, sample in pivot_heat.columns],
    fontsize=7.5, rotation=45, ha="right",
)
ax.set_yticks(np.arange(pivot_heat.shape[0]))
ax.set_yticklabels(pivot_heat.index.tolist(), fontsize=8)
n_samples = len(SAMPLES)
for i in range(1, len(assemblers_present)):
    ax.axvline(x=i * n_samples - 0.5, color="black", linewidth=1.1)
cbar = fig.colorbar(im, ax=ax, shrink=0.65)
cbar.set_label("Bin count")
plt.tight_layout()
out = output_dir / "heatmap_phylum.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"Figure saved: {out}")

print("✓ GTDB-Tk figure generation complete.")
