#!/bin/bash
# Generate pipeline figures from scripts/fig_*.py and scripts/plot_*.py.
# Outputs PNGs into figures/<category>/<figure_name>/.
#
# Set BINNING_DIR to override the default input data directory (default: binning_output).
#
# Usage:
#   ./scripts/generate_pipeline_figures.sh               # Run all figure scripts
#   ./scripts/generate_pipeline_figures.sh binning       # Run only fig_binning_*.py
#   ./scripts/generate_pipeline_figures.sh binning summary  # Run fig_binning_summary.py

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Filter pattern (optional)
PATTERN="${1:-.}"
FIGURE_NAME="${2:-.}"

# Ensure figures directory exists
mkdir -p figures

# Find matching figure scripts
if [[ -n "$FIGURE_NAME" && "$FIGURE_NAME" != "." ]]; then
    # Specific figure: scripts/fig_${PATTERN}_${FIGURE_NAME}.py
    SCRIPTS=("scripts/fig_${PATTERN}_${FIGURE_NAME}.py")
else
    # All figures matching pattern: scripts/fig_${PATTERN}_*.py
    mapfile -t SCRIPTS < <(find scripts -maxdepth 1 -name "fig_${PATTERN}_*.py" -o -name "fig_*.py" 2>/dev/null | grep -E "fig_${PATTERN}" | sort)
fi

if [[ ${#SCRIPTS[@]} -eq 0 ]]; then
    echo "No figure scripts found matching pattern: fig_${PATTERN}_*.py"
    exit 1
fi

echo "Found ${#SCRIPTS[@]} figure script(s) to execute..."

# Execute each figure script
for script in "${SCRIPTS[@]}"; do
    [[ ! -f "$script" ]] && continue

    # Extract category and figure name from script filename
    # Convention: scripts/fig_<category>_<name>.py
    basename=$(basename "$script" .py)
    # Remove 'fig_' prefix
    name_part="${basename#fig_}"

    # Split on first underscore to get category and name
    category="${name_part%%_*}"
    figure_name="${name_part#*_}"

    output_dir="figures/$category/$figure_name"
    mkdir -p "$output_dir"

    echo "▶ Running $script → $output_dir"

    # Run script with output directory as environment variable
    # Script should save figures to $FIGURE_OUTPUT_DIR
    if FIGURE_OUTPUT_DIR="$output_dir" BINNING_DIR="${BINNING_DIR:-binning_output}" python3 "$script"; then
        # Verify PNG(s) were created
        if find "$output_dir" -maxdepth 1 -name "*.png" -type f -exec true \; -print -quit >/dev/null 2>&1; then
            png_count=$(find "$output_dir" -maxdepth 1 -name "*.png" -type f | wc -l)
            echo "  ✓ Created $png_count PNG(s) in $output_dir"
        else
            echo "  ⚠ Warning: No PNG files found in $output_dir after execution"
        fi
    else
        echo "  ✗ Script failed: $script"
        exit 1
    fi
done

echo "✓ Figure generation complete!"
