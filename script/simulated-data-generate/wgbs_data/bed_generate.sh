#!/bin/bash

# =======================================================
# Script Name: batch_wgbs_pat_to_bed.sh
# Function: Batch process PAT files (.pat.gz) to BED files.
#           1. Index PAT files (generate .beta) if missing
#           2. wgbstools beta2bed + awk (generates merged 1-based BED)
# Usage: bash batch_wgbs_pat_to_bed.sh <input_folder>
# =======================================================

# 1. Argument Check
if [ $# -lt 1 ]; then
    echo "Error: Missing arguments."
    echo "Usage: bash $0 <input_folder_path>"
    echo "Example: bash $0 ./pat_files_folder"
    exit 1
fi

INPUT_DIR=$1
# Note: Threads are rarely needed for beta2bed/awk, but kept variable for consistency
THREADS=1 

# Check if input is a directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input path '$INPUT_DIR' is not a directory!"
    exit 1
fi

# =======================================================
# Dynamic Output Directory Name
# =======================================================
INPUT_FOLDER_NAME=$(basename "${INPUT_DIR}")
OUT_DIR="${INPUT_FOLDER_NAME}_output_pat"

PAT_DIR="${OUT_DIR}/pat"
BED_DIR="${OUT_DIR}/bed"

# Create directories
mkdir -p "$PAT_DIR"
mkdir -p "$BED_DIR"

echo "=== Batch Processing (PAT -> BED) Started ==="
echo "=== Input Directory:  ${INPUT_DIR}"
echo "=== Output Directory: ${OUT_DIR}"

# 3. Iterate through all .pat.gz files in the folder
shopt -s nullglob
PAT_FILES=("${INPUT_DIR}"/*.pat.gz)
total=${#PAT_FILES[@]}
shopt -u nullglob

if [ "$total" -eq 0 ]; then
    echo "Error: No .pat.gz files found in directory '${INPUT_DIR}'."
    exit 1
fi

count=0

for SOURCE_PAT in "${PAT_FILES[@]}"; do
    ((count++))
    # Remove .pat.gz extension to get basename
    FILENAME=$(basename "${SOURCE_PAT}")
    BASENAME="${FILENAME%%.pat.gz}"
    
    echo "------------------------------------------------------"
    echo "Processing file ($count / $total): ${BASENAME}"
    echo "------------------------------------------------------"

    # Define paths in the Output Directory
    # We use symbolic links to keep the output folder self-contained without duplicating large data
    LINKED_PAT="${PAT_DIR}/${BASENAME}.pat.gz"
    BETA_FILE="${PAT_DIR}/${BASENAME}.beta"

    # --- Step 1: Link PAT file and Check/Generate Beta Index ---
    
    # Create a symbolic link to the input PAT file in the output folder
    # This allows wgbstools to write the .beta file in our output folder, keeping input clean
    if [ ! -f "$LINKED_PAT" ]; then
        ln -s "$(realpath "${SOURCE_PAT}")" "$LINKED_PAT"
    fi

    # Check if .beta file exists (either copied or generated previously)
    if [ -f "${INPUT_DIR}/${BASENAME}.beta" ]; then
        # If beta exists in input, link it
        echo "[1/3] Found existing Beta file in input, linking..."
        ln -sf "$(realpath "${INPUT_DIR}/${BASENAME}.beta")" "$BETA_FILE"
    elif [ -f "$BETA_FILE" ]; then
         echo "[1/3] Beta file already exists in output, skipping indexing..."
    else
        # If not, generate it using wgbstools index
        echo "[1/3] Indexing PAT file (generating .beta)..."
        wgbstools index "$LINKED_PAT"
    fi

    if [ ! -f "$BETA_FILE" ]; then
        echo "Error: Beta file generation failed for ${BASENAME}. Skipping."
        continue
    fi

    # --- Step 2: Generate Temporary BEDs ---
    echo "[2/3] Generating temporary BED files..."
    TEMP_MEAN="${BASENAME}_mean.tmp.bed"
    TEMP_READS="${BASENAME}_reads.tmp.bed"

    wgbstools beta2bed --mean -o "$TEMP_MEAN" "$BETA_FILE"
    wgbstools beta2bed -o "$TEMP_READS" "$BETA_FILE"

    # --- Step 3: Merge and Convert (0-based to 1-based) ---
    echo "[3/3] Merging files and converting coordinates..."
    FINAL_BED="${BED_DIR}/${BASENAME}.bed"

    awk -v OFS="\t" '
    NR==FNR {
        key = $1 FS $2 FS $3
        map[key] = $4
        next
    }
    {
        key = $1 FS $2 FS $3
        mean_val = (key in map) ? map[key] : "NA"
        
        # Add 1 to Start coordinate
        $2 = $2 + 1
        
        print $0, mean_val
    }
    ' "$TEMP_MEAN" "$TEMP_READS" > "$FINAL_BED"

    # --- Cleanup Temporary Files ---
    rm "$TEMP_MEAN" "$TEMP_READS"
    
    echo ">> Sample ${BASENAME} finished."
    echo ""

done
