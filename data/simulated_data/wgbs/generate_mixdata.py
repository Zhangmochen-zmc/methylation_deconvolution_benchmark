import pandas as pd
from multiprocessing import Pool
import subprocess
import os

# Path to your CSV file containing mixing proportions
csv_file = '/disk1/weiyk/data/dirichlet/wgbsrandom.csv'

# Source PAT files corresponding to: B cells, CD4+, CD8+, Monocytes, Neutrophils, NK cells
pat_files = [
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_bcell.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_cd4.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_cd8.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_monocyte.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_neutrophil.sorted.pat.gz",
    "/data/weiyk/data/test_data/wgbs_data/blood/merged_data/merged_nk.sorted.pat.gz"
]

# Output directory configuration
output_dir = "../mix_data/simulated_random"
os.makedirs(output_dir, exist_ok=True)

# Target sequencing coverage
target_cov = 20.0

# Parallel resource allocation settings
max_parallel_jobs = 5   # Maximum number of concurrent samples to process
threads_per_job = 12    # Number of threads allocated per individual job

def process_single_sample(args):
    """
    Worker function called by the process pool.
    args is a tuple: (index, rates_list)
    """
    index, rates = args
    
    output_prefix = os.path.join(output_dir, f"sample_{index}")
    
    # Constructing the wgbstools command
    cmd = [
        "wgbstools", "mix_pat"
    ] + pat_files + [
        "--rates"
    ] + rates + [
        "-p", output_prefix,
        "-c", str(target_cov),
        "--force",
        "-@", str(threads_per_job)  # Passing threads per job to the tool
    ]
    
    print(f"[START] Processing Sample {index} in background...")
    
    try:
        # Executing the mixing command
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"[SUCCESS] Sample {index} generated successfully.")
        return None
    except subprocess.CalledProcessError as e:
        error_msg = f"[ERROR] Sample {index} failed: {e.stderr}"
        print(error_msg)
        return error_msg

def run_parallel_mixing():
    # Load mixing proportions from CSV
    df = pd.read_csv(csv_file)
    total_samples = len(df)
    
    print(f"Total samples to process: {total_samples}")
    print(f"Parallel Config: Running {max_parallel_jobs} jobs concurrently, {threads_per_job} threads each.")
    print(f"Total CPU core usage estimate: {max_parallel_jobs * threads_per_job}")

    # Prepare task list for the pool
    tasks = []
    for index, row in df.iterrows():
        rates_str = [str(r) for r in row.values]
        tasks.append((index, rates_str))

    # Initialize process pool and execute tasks
    with Pool(processes=max_parallel_jobs) as pool:
        pool.map(process_single_sample, tasks)

if __name__ == "__main__":
    run_parallel_mixing()
