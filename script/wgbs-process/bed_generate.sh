#!/bin/bash

# =======================================================
# Script Name: auto_detect_wgbs.sh
# Function: Automatically detects input file type (BAM or PAT)
#           and processes them to 1-based BED files.
# Usage: bash auto_detect_wgbs.sh <input_folder> [threads]
# =======================================================

# Argument Check
if [ $# -lt 1 ]; then
    echo "Error: Missing arguments."
    echo "Usage: bash $0 <input_folder_path> [threads]"
    exit 1
fi

INPUT_DIR=$1
GENOME="hg38"
THREADS=${2:-1}

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input path '$INPUT_DIR' is not a directory!"
    exit 1
fi

# =======================================================
# Convert Beta to BED (Common Step)
# =======================================================
generate_bed_from_beta() {
    local BASENAME=$1
    local BETA_FILE=$2
    local BED_DIR=$3

    echo "    [BED] Generating temporary files..."
    local TEMP_MEAN="${BASENAME}_mean.tmp.bed"
    local TEMP_READS="${BASENAME}_reads.tmp.bed"

    wgbstools beta2bed --mean -o "$TEMP_MEAN" "$BETA_FILE"
    wgbstools beta2bed -o "$TEMP_READS" "$BETA_FILE"

    # echo "    [BED] Merging and converting to 1-based coordinates..." (已删除)
    local FINAL_BED="${BED_DIR}/${BASENAME}.bed"

    awk -v OFS="\t" '
    NR==FNR {
        key = $1 FS $2 FS $3
        map[key] = $4
        next
    }
    {
        key = $1 FS $2 FS $3
        mean_val = (key in map) ? map[key] : "NA"
        $2 = $2 + 1  # 0-based to 1-based conversion
        print $0, mean_val
    }
    ' "$TEMP_MEAN" "$TEMP_READS" > "$FINAL_BED"

    rm "$TEMP_MEAN" "$TEMP_READS"
    echo "    [BED] Finished: ${BASENAME}.bed"
}

# =======================================================
# Detection Logic
# =======================================================
shopt -s nullglob
BAM_FILES=("${INPUT_DIR}"/*.bam)
PAT_FILES=("${INPUT_DIR}"/*.pat.gz)
shopt -u nullglob

MODE=""
if [ ${#BAM_FILES[@]} -gt 0 ]; then
    MODE="BAM"
    FILES=("${BAM_FILES[@]}")
elif [ ${#PAT_FILES[@]} -gt 0 ]; then
    MODE="PAT"
    FILES=("${PAT_FILES[@]}")
else
    echo "Error: No .bam or .pat.gz files found in '${INPUT_DIR}'."
    exit 1
fi

# =======================================================
# Setup Directories
# =======================================================
INPUT_FOLDER_NAME=$(basename "${INPUT_DIR}")
OUT_DIR="${INPUT_FOLDER_NAME}_output"

PAT_DIR="${OUT_DIR}/pat"
BED_DIR="${OUT_DIR}/bed"
SORTED_DIR="${OUT_DIR}/sorted_bam"

mkdir -p "$PAT_DIR" "$BED_DIR"

if [ "$MODE" == "BAM" ]; then
    mkdir -p "$SORTED_DIR"
fi

echo "=== Auto Detection Result ==="
echo "Input Directory: ${INPUT_DIR}"
echo "Detected Mode:   ${MODE} Processing"
echo "Found Files:     ${#FILES[@]}"
echo "Output Directory: ${OUT_DIR}"
echo "============================="

# =======================================================
# Main Processing Loop
# =======================================================
count=0
total=${#FILES[@]}

for FILE in "${FILES[@]}"; do
    ((count++))
    
    # Get Basename (handle both .bam and .pat.gz extensions)
    FILENAME=$(basename "${FILE}")
    if [ "$MODE" == "BAM" ]; then
        BASENAME="${FILENAME%.bam}"
    else
        BASENAME="${FILENAME%%.pat.gz}"
    fi

    echo "------------------------------------------------------"
    echo "Processing ($count / $total): ${BASENAME}"
    echo "------------------------------------------------------"

    # Define the path where the Beta file will eventually be
    TARGET_BETA="${PAT_DIR}/${BASENAME}.beta"

    # ---------------------------------------------------
    # PATH 1: BAM Input Logic
    # ---------------------------------------------------
    if [ "$MODE" == "BAM" ]; then
        SORTED_BAM="${SORTED_DIR}/${BASENAME}.sorted.bam"
        
        # Sort
        if [ -f "$SORTED_BAM" ]; then
            echo "[1/3] Sorted BAM exists, skipping sort..."
        else
            echo "[1/3] Running samtools sort..."
            samtools sort -@ "$THREADS" -o "$SORTED_BAM" "${FILE}"
        fi

        # Bam2Pat
        echo "[2/3] Running bam2pat..."
        # wgbstools generates output based on input filename
        # input: name.sorted.bam -> output: name.sorted.pat.gz
        wgbstools bam2pat -@ "$THREADS" "${SORTED_BAM}" --genome "$GENOME"

        # Move generated files
        GEN_PAT="${BASENAME}.sorted.pat.gz"
        GEN_CSI="${BASENAME}.sorted.pat.gz.csi"
        GEN_BETA="${BASENAME}.sorted.beta"

        if [ -f "$GEN_PAT" ]; then mv "$GEN_PAT" "${PAT_DIR}/${BASENAME}.pat.gz"; fi
        if [ -f "$GEN_CSI" ]; then mv "$GEN_CSI" "${PAT_DIR}/${BASENAME}.pat.gz.csi"; fi
        if [ -f "$GEN_BETA" ]; then mv "$GEN_BETA" "${TARGET_BETA}"; fi
    
    # ---------------------------------------------------
    # PATH 2: PAT Input Logic
    # ---------------------------------------------------
    else
        echo "[1/2] Processing PAT file..."
        LINKED_PAT="${PAT_DIR}/${BASENAME}.pat.gz"
        
        # Link PAT file
        if [ ! -f "$LINKED_PAT" ]; then
            ln -s "$(realpath "${FILE}")" "$LINKED_PAT"
        fi

        # Check/Generate Index
        if [ -f "${INPUT_DIR}/${BASENAME}.beta" ]; then
            echo "      Found existing beta in input, linking..."
            ln -sf "$(realpath "${INPUT_DIR}/${BASENAME}.beta")" "${TARGET_BETA}"
        elif [ ! -f "${TARGET_BETA}" ]; then
            echo "      Indexing PAT file (generating .beta)..."
            wgbstools index "$LINKED_PAT"
        fi
    fi

    # ---------------------------------------------------
    # Generate BED
    # ---------------------------------------------------
    if [ ! -f "$TARGET_BETA" ]; then
        echo "Error: Beta file missing for ${BASENAME}. Skipping BED generation."
        continue
    fi
    
    # Call the helper function defined at the top
    generate_bed_from_beta "$BASENAME" "$TARGET_BETA" "$BED_DIR"

    echo ">> Sample ${BASENAME} finished."
    echo ""

done

echo "=========================================="
echo "All processing completed!"
echo "Output: ${OUT_DIR}"
echo "=========================================="
