#!/bin/bash

# config
CELFIE_SCRIPT="decon.py"
PERF_WRAPPER="celfie_perf.py"
INPUT_BASE="test_data"
OUTPUT_BASE="celfie_result"
K_VALUE=100

if [ ! -f "$PERF_WRAPPER" ]; then
    echo "Error: Performance monitoring script not found. $PERF_WRAPPER"
    exit 1
fi

# Batch traversal
folders=$(ls -d ${INPUT_BASE}/*/)

echo "Start batch deconvolution processing..."
echo "time: $(date)"

for folder_path in $folders; do
    folder_name=$(basename "$folder_path")
    
    input_file=$(ls ${folder_path}/*.tsv 2>/dev/null | head -n 1)
    
    if [ -z "$input_file" ]; then
        echo "skip: $folder_name (.tsv file not found)"
        continue
    fi

    # Define the output directory for this sample.
    current_out="${OUTPUT_BASE}/${folder_name}"
    mkdir -p "$current_out"

    echo "-------------------------------------------------------"
    echo "Processing task: $folder_name"
    echo "Input path: $input_file"
    echo "Save path: $current_out"

    # Call the performance monitoring wrapper to run Python
    python "$PERF_WRAPPER" "$CELFIE_SCRIPT" "$input_file" "$current_out" "$K_VALUE"

    if [ $? -eq 0 ]; then
        echo "task $folder_name down"
    else
        echo "task $folder_name fail"
    fi
done

echo "-------------------------------------------------------"
echo "All batch tasks completed.！"
echo "Overall completion time: $(date)"
