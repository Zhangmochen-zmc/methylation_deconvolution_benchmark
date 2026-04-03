# =======================================================
# Script Name: mixdata.py
# Function: Simulate mixed WGBS samples by combining purified cell-type profiles.
#           1. Load cell-type proportions from a CSV file
#           2. Mix PAT files using wgbstools mix_pat
#           3. Generate mixed samples at target coverage
# Usage: python mix_pat.py -i proportions.csv -o output_dir
# =======================================================

import pandas as pd
from multiprocessing import Pool
import subprocess
import os
import argparse

# ================= Argument Parsing =================

parser = argparse.ArgumentParser(description="Generate mixed WGBS samples from cell-type proportions.")
parser.add_argument("-i", "--input", required=True, help="Path to input CSV file (proportions)")
parser.add_argument("-o", "--output", required=True, help="Output directory")
args = parser.parse_args()

csv_file = args.input
output_dir = args.output
os.makedirs(output_dir, exist_ok=True)

# ================= Fixed Configuration =================

pat_files = [
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_bcell.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_cd4.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_cd8.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_monocyte.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_neutrophil.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_nk.sorted.pat.gz"
]

target_cov = 20.0
max_parallel_jobs = 5
threads_per_job = 12

# =====================================================

def process_single_sample(args):
    index, rates = args
    
    output_prefix = os.path.join(output_dir, f"sample_{index}")
    
    cmd = [
        "wgbstools", "mix_pat"
    ] + pat_files + [
        "--rates"
    ] + rates + [
        "-p", output_prefix,
        "-c", str(target_cov),
        "--force",
        "-@", str(threads_per_job)
    ]
    
    print(f"[START] Processing sample {index}...")
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"[DONE] Sample {index} generated successfully.")
        return None
    except subprocess.CalledProcessError as e:
        error_msg = f"[ERROR] Sample {index} failed: {e.stderr}"
        print(error_msg)
        return error_msg

def run_parallel_mixing():
    df = pd.read_csv(csv_file)
    total_samples = len(df)

    print(f"Total samples: {total_samples}")
    print(f"Output directory: {output_dir}")

    tasks = []
    for index, row in df.iterrows():
        rates_str = [str(r) for r in row.values]
        tasks.append((index, rates_str))

    with Pool(processes=max_parallel_jobs) as pool:
        pool.map(process_single_sample, tasks)

if __name__ == "__main__":
    run_parallel_mixing()
