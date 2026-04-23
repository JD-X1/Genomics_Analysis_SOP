#!/usr/bin/env python3
"""
Template figure: binning method summary comparison.

Convention:
- Script name: scripts/fig_<category>_<name>.py
  → Outputs to: figures/<category>/<name>/
  → Category extracted as: first word after 'fig_'
  → Figure name extracted as: remainder after category

- Figures are always saved as PNG in $FIGURE_OUTPUT_DIR
- Use environment variable: os.environ['FIGURE_OUTPUT_DIR']

Run via: ./scripts/make_figures.sh binning summary
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Get output directory from environment
output_dir = os.environ.get('FIGURE_OUTPUT_DIR')
if not output_dir:
    print("Error: FIGURE_OUTPUT_DIR environment variable not set", file=sys.stderr)
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)

# Load bin statistics
binning_dir = os.environ.get('BINNING_DIR', 'binning_output')
stats_path = os.path.join(binning_dir, 'bin_stats_summary.csv')
if not os.path.exists(stats_path):
    print(f"Error: {stats_path} not found. Run exec_bin_analysis.sh first.", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(stats_path)

# Example: Count bins by method
method_counts = df.groupby('binner').size().reset_index(name='count')

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

sns.barplot(
    data=method_counts,
    x='binner',
    y='count',
    palette='Set2',
    ax=ax
)

ax.set_title('Bin Count by Binning Method', fontsize=14, fontweight='bold')
ax.set_xlabel('Binning Method', fontsize=12)
ax.set_ylabel('Number of Bins', fontsize=12)

plt.tight_layout()

# Save as PNG (not PDF)
output_path = os.path.join(output_dir, 'bin_count_by_method.png')
plt.savefig(output_path, dpi=300, format='png')
print(f"Figure saved: {output_path}")

plt.close()

print("✓ Figure generation successful")
