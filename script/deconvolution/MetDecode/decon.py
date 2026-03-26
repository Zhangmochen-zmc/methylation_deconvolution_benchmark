import time
import gc
import os
import glob
import subprocess
import pandas as pd
from memory_profiler import memory_usage

# config
ATLAS_PATH = "metdecode_ref/atlas.tsv"
DATA_DIR = "test_data"
RESULTS_DIR = "metdecode_result"
LOG_OUTPUT = "metdecode_benchmark.csv"

def run_deconvolution(cmd):
    """Wrapper function for executing subprocess tasks"""
    # shell=False
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

def main():
    # 1. create results dir
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
        print(f"create dir: {RESULTS_DIR}")

    # 2. get all tsv
    all_files = glob.glob(os.path.join(DATA_DIR, "*.tsv"))
    
    performance_records = []

    print(f"find {len(all_files)} files，begin decon...\n")

    for file_path in all_files:
        if os.path.abspath(file_path) == os.path.abspath(ATLAS_PATH):
            continue

        filename = os.path.basename(file_path)
        name = os.path.splitext(filename)[0]
        output_path = os.path.join(RESULTS_DIR, f"{name}.csv")

        # construct cmd
        cmd = ["python3", "MetDecode/run.py", ATLAS_PATH, file_path, output_path]

        print(f"processing: {name} ...")

        # benchmark
        gc.collect()  
        start_time = time.perf_counter()

        try:
            # memory_usage
            # 返回值是采样点内存列表 (单位: MiB)
            mem_history = memory_usage((run_deconvolution, (cmd,)), interval=0.1, timeout=None)
            
            end_time = time.perf_counter()

            MIB_TO_MB = 1.048576
            
            # calculate
            duration = end_time - start_time
            max_mem = (max(mem_history) * MIB_TO_MB) if mem_history else 0
            
            print(f"  > down! time: {duration:.2f}s, peak memory: {max_mem:.2f}MB")

            # save
            performance_records.append({
                "Sample": name,
                "Time_s": round(duration, 4),
                "Memory_MB": round(max_mem, 2),
                "Status": "Success"
            })

        except Exception as e:
            print(f"  > process {name} fail: {e}")
            performance_records.append({
                "Sample": name,
                "Time_s": 0,
                "Memory_MB": 0,
                "Status": f"Failed: {str(e)}"
            })

    # 3. import results
    df_perf = pd.DataFrame(performance_records)
    
    # save results
    df_perf.to_csv(LOG_OUTPUT, index=False)
    
    print("\n" + "="*30)
    print(f"all tasks down！")
    print(f"benchmark has benn saved: {LOG_OUTPUT}")
    print("="*30)
    print(df_perf.to_string(index=False))

if __name__ == "__main__":
    main()
