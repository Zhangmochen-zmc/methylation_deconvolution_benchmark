#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch Merge TXT Files in Subfolders

This script scans all subfolders in the current directory,
reads all TXT files in each folder, and merges them by 'probeID'.
The merged output for each folder is saved as {folder_name}.txt.

Author: Your Name
Date: 2026-05-08
"""

import os
import glob
import pandas as pd

# --- Configuration ---
# Root directory to scan ('.' means current script directory)
root_directory = 'methylation_deconvolution_benchmark/data/reference_data/450k/original_data'

# --- Script Start ---

# 1. List all subfolders in the root directory
all_items = os.listdir(root_directory)
folders = [f for f in all_items if os.path.isdir(os.path.join(root_directory, f))]

print(f"Detected the following folders to process: {folders}")

# 2. Process each folder
for folder_name in folders:

    folder_path = os.path.join(root_directory, folder_name)
    output_file = f"{folder_name}.txt"  # Output file: folder_name + .txt
    
    # Find all .txt files in the current folder
    file_paths = glob.glob(os.path.join(folder_path, '*.txt'))

    if not file_paths:
        print(f"--- Skipped: No .txt files found in folder '{folder_name}'.")
        continue

    print(f"Processing folder: {folder_name} ({len(file_paths)} files found)...")

    list_of_dfs = []
    for f_path in file_paths:
        try:
            # Extract sample name from file name
            sample_name = os.path.splitext(os.path.basename(f_path))[0]

            # Read TXT file (columns separated by spaces or tabs)
            temp_df = pd.read_csv(f_path, sep='\s+', header=None, names=['probeID', sample_name])

            # Set 'probeID' as index
            temp_df.set_index('probeID', inplace=True)

            list_of_dfs.append(temp_df)
        except Exception as e:
            print(f"  [Warning] Failed to read file {f_path}: {e}")

    # Merge all dataframes if any were successfully read
    if list_of_dfs:
        merged_df = pd.concat(list_of_dfs, axis=1)
        merged_df.to_csv(output_file, sep='\t')
        print(f"--> Done! Merged file saved as: {output_file}")
    else:
        print(f"--- No valid data files were merged in folder '{folder_name}'.")

print("\nAll folders have been processed.")
