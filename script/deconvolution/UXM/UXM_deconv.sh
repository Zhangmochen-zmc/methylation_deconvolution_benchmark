#!/bin/bash

# =======================================================
# Script Name: batch_uxm_deconv.sh
# Description:
#   Batch processing script for running `uxm deconv` on multiple
#   .pat.gz files sequentially (no parallel execution).
#
# Usage:
#   bash batch_uxm_deconv.sh <input_dir> <output_dir>
#
# Example:
#   bash batch_uxm_deconv.sh ./pat_files ./results
#
# Requirements:
#   - uxm (installed and available in PATH)
#   - Atlas file (specified below)
# =======================================================

# ================= Configuration =================

# Input directory containing .pat.gz files
INPUT_DIR="$1"

# Output directory (will be created if not exists)
OUTPUT_DIR="$2"

# Absolute path to atlas file
ATLAS="/data/weiyk/work/nature/data/blood/25markers/atlas/25marker_atlas.csv"

# Create output directory if it does not exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# =================================================

echo "Start processing..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Loop through all matching .pat.gz files
for file in "${INPUT_DIR}"/sample_*.pat.gz; do

    # Skip if no files matched
    [ -e "$file" ] || continue

    # Extract filename without extension (e.g., sample_0_1)
    filename=$(basename "$file")
    base="${filename%.pat.gz}"

    echo "Processing: $base"

    # =================================================
    # Core command (sequential execution)
    # =================================================
    uxm deconv \
        -a "$ATLAS" \
        -@ 1 \
        -l 4 \
        "$file" \
        -o "$OUTPUT_DIR/${base}.csv"
    # =================================================

done

echo "All jobs completed!"