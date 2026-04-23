#!/usr/bin/env python3
"""
Parse GTDBTk results and aggregate taxonomic classifications by sample + assembly method.
Generates CSV files for phylum and family-level summaries.
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import pandas as pd


def extract_rank(classification_string, rank_prefix):
    """Extract taxonomic rank from GTDB classification string (e.g., 'p__' for phylum)."""
    if not classification_string or pd.isna(classification_string):
        return None

    pattern = f'{rank_prefix}([^;]+)'
    match = re.search(pattern, classification_string)
    if match:
        value = match.group(1).strip()
        return value if value else None
    return None


def parse_gtdbtk_results(base_dir):
    """Parse all GTDBTk summary files and aggregate by sample + assembly method."""
    base_path = Path(base_dir)

    # Collections for aggregation
    phylum_counts = defaultdict(lambda: defaultdict(int))  # {(sample, assembly): {phylum: count}}
    family_counts = defaultdict(lambda: defaultdict(int))   # {(sample, assembly): {family: count}}
    bin_counts = defaultdict(lambda: defaultdict(int))      # {(sample, assembly): count}

    summary_files = list(base_path.glob('*/*/*/gtdbtk.bac120.summary.tsv')) + \
                    list(base_path.glob('*/*/*/gtdbtk.ar53.summary.tsv')) + \
                    list(base_path.glob('*/gtdbtk.bac120.summary.tsv')) + \
                    list(base_path.glob('*/gtdbtk.ar53.summary.tsv'))

    # Remove duplicates (keep root level files, not classify subdirectory)
    summary_files = list(set(summary_files))
    summary_files = [f for f in summary_files if '/classify/' not in str(f)]

    print(f"Found {len(summary_files)} summary files")

    for file_path in sorted(summary_files):
        # Parse the path to extract sample and assembly method
        parts = file_path.parts

        # Find the sample name and assembly method
        # Pattern: gtdbtk_results/{sample}/{assembly}/{library}/gtdbtk.*.summary.tsv
        # or: gtdbtk_results/{sample}_details/gtdbtk.*.summary.tsv

        try:
            gtdbtk_idx = parts.index('gtdbtk_results')
            if gtdbtk_idx + 3 < len(parts):
                sample = parts[gtdbtk_idx + 1]
                assembly = parts[gtdbtk_idx + 2]
            else:
                sample = parts[gtdbtk_idx + 1].replace('_megahit_large', '').replace('_S7_L003', '')
                assembly = 'megahit_large'
        except (ValueError, IndexError):
            print(f"Skipping {file_path}: Could not parse path")
            continue

        key = (sample, assembly)

        try:
            df = pd.read_csv(file_path, sep='\t')
            if 'classification' not in df.columns:
                print(f"Skipping {file_path}: No 'classification' column")
                continue

            for _, row in df.iterrows():
                classification = row['classification']

                phylum = extract_rank(classification, 'p__')
                family = extract_rank(classification, 'f__')

                bin_counts[key][row['user_genome']] = 1

                if phylum and phylum != 'Unclassified':
                    phylum_counts[key][phylum] += 1

                if family and family != 'Unclassified':
                    family_counts[key][family] += 1

        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue

    return phylum_counts, family_counts, bin_counts


def save_aggregated_results(phylum_counts, family_counts, bin_counts, output_dir):
    """Save aggregated results to CSV files."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Phylum summary
    phylum_data = []
    for (sample, assembly), counts in sorted(phylum_counts.items()):
        total = sum(counts.values())
        for phylum, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            phylum_data.append({
                'sample': sample,
                'assembly': assembly,
                'phylum': phylum,
                'count': count,
                'percentage': 100 * count / total if total > 0 else 0
            })

    phylum_df = pd.DataFrame(phylum_data)
    phylum_csv = output_path / 'gtdbtk_phylum_summary.csv'
    phylum_df.to_csv(phylum_csv, index=False)
    print(f"Saved phylum summary: {phylum_csv}")

    # Family summary
    family_data = []
    for (sample, assembly), counts in sorted(family_counts.items()):
        total = sum(counts.values())
        for family, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
            family_data.append({
                'sample': sample,
                'assembly': assembly,
                'family': family,
                'count': count,
                'percentage': 100 * count / total if total > 0 else 0
            })

    family_df = pd.DataFrame(family_data)
    family_csv = output_path / 'gtdbtk_family_summary.csv'
    family_df.to_csv(family_csv, index=False)
    print(f"Saved family summary: {family_csv}")

    # Bin counts summary
    bin_data = []
    for (sample, assembly), bins in sorted(bin_counts.items()):
        bin_data.append({
            'sample': sample,
            'assembly': assembly,
            'total_bins': len(bins)
        })

    bin_df = pd.DataFrame(bin_data)
    bin_csv = output_path / 'gtdbtk_bin_counts.csv'
    bin_df.to_csv(bin_csv, index=False)
    print(f"Saved bin counts: {bin_csv}")


if __name__ == '__main__':
    base_dir = '/scratch/jduque2/oceanc_core_data/gtdbtk_results'
    output_dir = '/scratch/jduque2/oceanc_core_data/gtdbtk_analysis'

    phylum_counts, family_counts, bin_counts = parse_gtdbtk_results(base_dir)
    save_aggregated_results(phylum_counts, family_counts, bin_counts, output_dir)

    print("\nAggregation complete!")
