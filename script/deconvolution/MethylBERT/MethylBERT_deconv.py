#Usage: python script_name.py <input_folder> <output_base_folder>

import os
import subprocess
import sys

# ================= CONFIGURATION AREA =================
# 1. Input folder containing .kmers_fixed.csv files (from previous steps)
INPUT_FOLDER = sys.argv[1] 

# 2. Base directory for results
# The script creates a subfolder for each sample within this directory
OUTPUT_BASE_FOLDER = sys.argv[2]

# 3. Model path
MODEL_PATH = 'output/finetune_b256/'

# 4. Parameters
BATCH_SIZE = 128
THREADS = 1        # Set threads to 1 for single-threaded processing
FILE_SUFFIX = '.kmers_fixed.csv' 
# =====================================================

def run_batch():
    # Check if input directory exists
    if not os.path.exists(INPUT_FOLDER):
        print(f"Error: Input folder '{INPUT_FOLDER}' does not exist.")
        return

    # Filter files by suffix
    files = [f for f in os.listdir(INPUT_FOLDER) if f.endswith(FILE_SUFFIX)]
    
    if not files:
        print(f"No files with suffix '{FILE_SUFFIX}' found in '{INPUT_FOLDER}'.")
        return

    print(f"Detected {len(files)} files. Starting batch processing...\n")

    for filename in files:
        # --- 1. Prepare Paths ---
        input_filepath = os.path.join(INPUT_FOLDER, filename)
        
        # --- 2. Extract sample name and create specific output directory ---
        # Example: 'sample_1.kmers_fixed.csv' -> sample name 'sample_1'
        sample_name = filename.replace(FILE_SUFFIX, '')
        
        # Create directory: e.g., result/real/sample_1/
        specific_output_dir = os.path.join(OUTPUT_BASE_FOLDER, sample_name)
        
        if not os.path.exists(specific_output_dir):
            os.makedirs(specific_output_dir)
            
        print(f">>> Processing sample: {sample_name}")
        print(f"    Output directory: {specific_output_dir}")

        # --- 3. Construct and Execute Command ---
        # Added -t {THREADS} to ensure single-threaded execution if supported by the tool
        # Paths are wrapped in double quotes to handle potential spaces
        cmd = (
            f'methylbert deconvolute '
            f'-i "{input_filepath}" '
            f'-m "{MODEL_PATH}" '
            f'-o "{specific_output_dir}" '
            f'-b {BATCH_SIZE} '
            f'-t {THREADS}'
        )

        try:
            # subprocess.run with check=True will raise an exception if the command fails
            # This allows for immediate error detection during batch runs
            subprocess.run(cmd, shell=True, check=True)
            print("    [Status: Success]\n")
        except subprocess.CalledProcessError as e:
            print(f"    [Status: Error] Failed to process sample {sample_name}")
            print(f"    Error details: {e}\n")

    print("All tasks completed.")

if __name__ == "__main__":
    run_batch()
