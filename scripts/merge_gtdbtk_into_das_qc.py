#!/usr/bin/env python3
"""
Merge GTDBTk taxonomic classifications into das_bins_qc_summary.csv.

Adds columns:
  - gtdbtk_classification: full GTDB taxonomic string
  - gtdbtk_source: 'bac120' or 'ar53'
  - gtdbtk_phylum: extracted phylum level
  - gtdbtk_family: extracted family level
  - gtdbtk_ani: ANI to closest reference (if available)
"""

import re
import sys
from pathlib import Path
import pandas as pd


def extract_rank(classification_string, rank_prefix):
    """Extract taxonomic rank from GTDB classification string."""
    if not classification_string or pd.isna(classification_string):
        return None

    pattern = f'{rank_prefix}([^;]+)'
    match = re.search(pattern, str(classification_string))
    if match:
        value = match.group(1).strip()
        return value if value else None
    return None


def load_gtdbtk_results(gtdbtk_dir):
    """Load all GTDBTk summary files and index by (sample, assembly, bin_name)."""
    gtdbtk_dir = Path(gtdbtk_dir)
    gtdbtk_data = {}  # {(sample, assembly, bin_name): {gtdbtk_info}}

    # Find all summary files (avoid duplicates from classify/ subdirs)
    summary_files = []
    summary_files.extend(gtdbtk_dir.glob('*/*/*/gtdbtk.bac120.summary.tsv'))
    summary_files.extend(gtdbtk_dir.glob('*/*/*/gtdbtk.ar53.summary.tsv'))
    summary_files.extend(gtdbtk_dir.glob('*/gtdbtk.bac120.summary.tsv'))
    summary_files.extend(gtdbtk_dir.glob('*/gtdbtk.ar53.summary.tsv'))

    # Remove duplicates and classify subdirectories
    summary_files = list(set(summary_files))
    summary_files = [f for f in summary_files if '/classify/' not in str(f)]

    for file_path in sorted(summary_files):
        parts = file_path.parts

        try:
            gtdbtk_idx = parts.index('gtdbtk_results')
            if gtdbtk_idx + 2 < len(parts):
                sample = parts[gtdbtk_idx + 1]
                assembly = parts[gtdbtk_idx + 2]
            else:
                continue
        except (ValueError, IndexError):
            print(f"Warning: Could not parse {file_path}", file=sys.stderr)
            continue

        source = 'bac120' if 'bac120' in file_path.name else 'ar53'

        try:
            df = pd.read_csv(file_path, sep='\t')
            if 'classification' not in df.columns or 'user_genome' not in df.columns:
                continue

            for _, row in df.iterrows():
                bin_name = str(row['user_genome']).strip()
                classification = row['classification']

                key = (sample, assembly, bin_name)

                # Skip if already have better data (bac120 > ar53)
                if key in gtdbtk_data:
                    if gtdbtk_data[key]['source'] == 'bac120':
                        continue

                ani_val = ''
                try:
                    ani_raw = row.get('closest_genome_ani', '')
                    if pd.notna(ani_raw) and str(ani_raw).strip():
                        ani_val = float(ani_raw)
                except (ValueError, TypeError):
                    ani_val = ''

                gtdbtk_data[key] = {
                    'classification': classification if pd.notna(classification) else '',
                    'source': source,
                    'phylum': extract_rank(classification, 'p__') or '',
                    'family': extract_rank(classification, 'f__') or '',
                    'ani': ani_val,
                }
        except Exception as e:
            print(f"Warning: Error processing {file_path}: {e}", file=sys.stderr)
            continue

    return gtdbtk_data


def merge_gtdbtk_into_qc(qc_csv, gtdbtk_dir, output_csv):
    """Load QC summary and merge in GTDBTk results."""
    if not Path(qc_csv).exists():
        print(f"Error: {qc_csv} not found", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(qc_csv)
    print(f"Loaded {len(df)} bins from {qc_csv}")

    # Load GTDBTk results
    gtdbtk_data = load_gtdbtk_results(gtdbtk_dir)
    print(f"Loaded GTDBTk data for {len(gtdbtk_data)} bin-sample-assembly tuples")

    # Add GTDBTk columns
    df['gtdbtk_classification'] = ''
    df['gtdbtk_source'] = ''
    df['gtdbtk_phylum'] = ''
    df['gtdbtk_family'] = ''
    df['gtdbtk_ani'] = ''

    # Merge GTDBTk data
    matched = 0
    for idx, row in df.iterrows():
        sample_base = row.get('sample_base', '')
        bin_name = row['bin_name']
        assembler = row.get('assembler', '')

        # Try to match using directory structure: gtdbtk_results/{sample}/{assembly}/
        # DAS sample_base is like "E17_1B_S3_L003" and assembler is "megahit_large"
        # GTDBTk sample is like "E17_1B" and assembly is "megahit_large"

        # Extract just the base sample (E17_1B part) from sample_base
        sample_match = re.match(r'^(E\d+_[A-Z0-9]+)', sample_base)
        if sample_match:
            gtdb_sample = sample_match.group(1)
            key = (gtdb_sample, assembler, bin_name)

            if key in gtdbtk_data:
                gtdb = gtdbtk_data[key]
                df.at[idx, 'gtdbtk_classification'] = gtdb['classification']
                df.at[idx, 'gtdbtk_source'] = gtdb['source']
                df.at[idx, 'gtdbtk_phylum'] = gtdb['phylum']
                df.at[idx, 'gtdbtk_family'] = gtdb['family']
                ani_str = str(gtdb['ani']) if gtdb['ani'] else ''
                df.at[idx, 'gtdbtk_ani'] = ani_str
                matched += 1

    print(f"Matched {matched}/{len(df)} bins with GTDBTk results")

    # Save merged result
    df.to_csv(output_csv, index=False)
    print(f"Saved merged table: {output_csv}")

    return df


if __name__ == '__main__':
    project_root = Path(__file__).resolve().parent.parent

    qc_csv = project_root / 'binning_output' / 'das_bins_qc_summary.csv'
    gtdbtk_dir = project_root / 'gtdbtk_results'
    output_csv = project_root / 'binning_output' / 'das_bins_qc_summary_with_gtdbtk.csv'

    if not qc_csv.exists():
        print(f"Error: {qc_csv} not found. Run exec_dasBins_qc.sh first.", file=sys.stderr)
        sys.exit(1)

    if not gtdbtk_dir.exists():
        print(f"Error: {gtdbtk_dir} not found.", file=sys.stderr)
        sys.exit(1)

    df = merge_gtdbtk_into_qc(str(qc_csv), str(gtdbtk_dir), str(output_csv))

    # Print summary statistics
    print("\n" + "="*70)
    print("GTDBTk Integration Summary")
    print("="*70)
    print(f"Total bins: {len(df)}")
    gtdb_classified = (df['gtdbtk_classification'] != '').sum()
    print(f"Bins with GTDBTk classification: {gtdb_classified}")
    print(f"Bins from bac120: {(df['gtdbtk_source'] == 'bac120').sum()}")
    print(f"Bins from ar53: {(df['gtdbtk_source'] == 'ar53').sum()}")

    if gtdb_classified > 0:
        phyla_data = df[df['gtdbtk_phylum'] != '']['gtdbtk_phylum']
        families_data = df[df['gtdbtk_family'] != '']['gtdbtk_family']
        print(f"\nUnique phyla detected: {phyla_data.nunique()}")
        print(f"Unique families detected: {families_data.nunique()}")

        # Show top phyla
        top_phyla = phyla_data.value_counts().head(10)
        print("\nTop 10 phyla:")
        for phylum, count in top_phyla.items():
            print(f"  {phylum}: {count}")
