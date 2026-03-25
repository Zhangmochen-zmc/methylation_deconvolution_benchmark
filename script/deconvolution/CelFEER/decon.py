import time
import gc
import pandas as pd
import subprocess
import os
from memory_profiler import memory_usage

# config
INPUT_PATH = "test_data/test_data.txt"
OUTPUT_DIR = "celfeer_result"
NUM_SAMPLES = "1"
SCRIPT_PATH = "CelFEER/scripts/celfeer.py"

def run_task():
    """
    The core tasks to be performed
    """
    # Build command
    cmd = ["python", SCRIPT_PATH, INPUT_PATH, OUTPUT_DIR, NUM_SAMPLES]
    
    # Running child processes
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    # 1. Clear memory, ready to begin
    gc.collect()
    
    print(f"Start running the task: {SCRIPT_PATH} ...")
    
    # 2. Record start time
    start_perf = time.perf_counter()
    start_time_str = time.strftime("%Y-%m-%d %H:%M:%S")

    # 3. Monitor memory usage
    try:
        mem_history = memory_usage(run_task, interval=0.1, timeout=None)
    except Exception as e:
        print(f"error: {e}")
        mem_history = [0]

    # 4. Record end time
    end_perf = time.perf_counter()
    
    # 5. Calculation indicators

    CONVERSION_FACTOR = 1.048576
    duration = end_perf - start_perf
    max_memory_mb = (max(mem_history) * CONVERSION_FACTOR) if mem_history else 0
    avg_memory_mb = (sum(mem_history) / len(mem_history) * CONVERSION_FACTOR) if mem_history else 0

    # 6. Use Pandas to organize the results
    results = {
        "Task": [os.path.basename(SCRIPT_PATH)],
        "Start_Time": [start_time_str],
        "Duration_Sec": [round(duration, 4)],
        "Max_Memory_MB": [round(max_memory, 2)],
        "Avg_Memory_MB": [round(avg_memory, 2)],
        "Input_File": [INPUT_PATH]
    }
    
    df = pd.DataFrame(results)
    
    # 7. 打印并保存到 CSV
    print("\n--- benchmark ---")
    print(df.to_string(index=False))
    
    log_file = "celfeer_benchmark.csv"
    if not os.path.exists(log_file):
        df.to_csv(log_file, index=False)
    else:
        df.to_csv(log_file, mode='a', header=False, index=False)
        
    print(f"\nPerformance records have been updated to: {log_file}")
