import time
import gc
import os
import sys
import subprocess
import pandas as pd
from memory_profiler import memory_usage

def run_celfie(cmd):
    """run decon"""
    # gc
    gc.collect()
    result = subprocess.run(cmd, shell=True)
    return result.returncode

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python perf_monitor.py <script> <input> <output> <k>")
        sys.exit(1)

    script_py = sys.argv[1]
    input_file = sys.argv[2]
    output_dir = sys.argv[3]
    k_val = sys.argv[4]

    # Construct deconvolution command
    command = f"python {script_py} {input_file} {output_dir} {k_val}"

    start_time = time.time()
    
    # Monitor memory usage (sampled every 0.1 seconds).
    try:
        mem_vals = memory_usage((run_celfie, (command,)), interval=0.1, timeout=None)
        max_mem = max(mem_vals) if mem_vals else 0
    except Exception as e:
        print(f"Profiling Error: {e}")
        max_mem = 0

    end_time = time.time()
    elapsed = end_time - start_time

    # Save performance data to the current task's output directory.
    stats_file = os.path.join(output_dir, "performance_stats.csv")
    df = pd.DataFrame([{
        "folder": os.path.basename(os.path.dirname(input_file)),
        "time_seconds": round(elapsed, 2),
        "peak_memory_mb": round(max_mem, 2)
    }])
    df.to_csv(stats_file, index=False)
    
    print(f"Finished. Time: {elapsed:.2f}s, Max Mem: {max_mem:.2f}MB")
