import os
import glob
import subprocess
import sys
import argparse
import time
import gc
import pandas as pd
from memory_profiler import memory_usage

MAIN_SCRIPT = "deconvolve.py"

def run_benchmark():
    # 1. Define command line arguments
    parser = argparse.ArgumentParser(
        description="Batch run and monitor performance (Time & Memory)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-a', '--atlas', default='methatlas_ref/ref.csv', 
                        help="Reference Atlas CSV")
    parser.add_argument('-i', '--input_dir', default='test_data', 
                        help="Input folder containing the CSV file to be analyzed")
    parser.add_argument('-o', '--out_dir', default='methatlas_result', 
                        help="Output folder")
    
    # The parameters passed to deconvolve.py
    parser.add_argument('--plot', action='store_true', help="plot")
    parser.add_argument('--slim', action='store_true', help="Simplify output")
    parser.add_argument('--residuals', '-r', action='store_true', help="Output residual")

    args = parser.parse_args()

    # 2. check
    if not os.path.exists(MAIN_SCRIPT):
        print(f"[Error] Main program not found: {MAIN_SCRIPT}")
        sys.exit(1)
    if not os.path.exists(args.atlas):
        print(f"[Error] Reference file does not exist: {args.atlas}")
        sys.exit(1)
    if not os.path.exists(args.input_dir):
        print(f"[Error] The input directory does not exist.: {args.input_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # 3. Get file list
    search_pattern = os.path.join(args.input_dir, "*.csv")
    csv_files = glob.glob(search_pattern)

    if not csv_files:
        print(f"[Warning] {args.input_dir} no .csv file。")
        sys.exit(0)

    print(f"==========================================")
    print(f"{len(csv_files)} files detected, batch processing begins....")
    print(f"reference atlas: {args.atlas}")
    print(f"==========================================\n")

    # benchamrk
    benchmark_records = []

    # 4. Looping and monitoring
    for idx, csv_file in enumerate(csv_files):
        file_name = os.path.basename(csv_file)
        print(f"[{idx+1}/{len(csv_files)}] processing: {file_name}")

        # 4.1 Construct command
        cmd = [
            sys.executable, 
            MAIN_SCRIPT,
            "-a", args.atlas,
            csv_file,
            "-o", args.out_dir
        ]
        if args.plot: cmd.append("--plot")
        if args.slim: cmd.append("--slim")
        if args.residuals: cmd.append("--residuals")

        # 4.2 gc
        gc.collect()
        time.sleep(1) 

        # 4.3 Initialize variables
        start_time = time.perf_counter()
        peak_memory_mb = 0
        status = "Success"
        error_msg = "None"

        try:
            # 4.4 memory_usage
            mem_max = memory_usage(
                (subprocess.run, [cmd], {'check': True}), 
                max_usage=True, 
                interval=0.1, 
                include_children=True, 
                retval=False
            )
            
            # return
            if isinstance(mem_max, (list, tuple)):
                peak_memory_mb = max(mem_max)
            else:
                peak_memory_mb = mem_max

        except subprocess.CalledProcessError as e:
            print(f"  >>> [Error] Script execution failed.")
            status = "Failed"
            error_msg = "Subprocess Error"
            peak_memory_mb = 0 
        except Exception as e:
            print(f"  >>> [Error] Monitoring error: {e}")
            status = "Error"
            error_msg = str(e)

        # 4.5 time
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        # 4.6 print
        print(f"  -> Time:   {elapsed_time:.4f} s")
        print(f"  -> Memory: {peak_memory_mb:.2f} MB")
        print("-" * 40)

        # 4.7 save
        benchmark_records.append({
            "File": file_name,
            "Time_Seconds": round(elapsed_time, 4),
            "Peak_Memory_MB": round(peak_memory_mb, 4),
            "Status": status,
            "Error": error_msg
        })

    # 5. save
    log_file = "methatlas_benchmark.csv"
    try:
        df = pd.DataFrame(benchmark_records)
        df.to_csv(log_file, index=False)
        print(f"\nAll done! Performance statistics have been saved to: {log_file}")
    except Exception as e:
        print(f"\nAll done! But saving to CSV failed: {e}")
        print(benchmark_records)

if __name__ == "__main__":
    run_benchmark()
