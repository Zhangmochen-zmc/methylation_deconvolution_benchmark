import os
import glob
import subprocess
import sys
import argparse
import time
import gc
import shutil
import pandas as pd
from memory_profiler import memory_usage

def run_menet(input_file, model_path, output_dir, input_type, menet_path="MEnet", bedtools_path=None):
    """
    Run the MEnet predict command line
    """
    # Build basic commands
    cmd = [
        menet_path, "predict",
        "--input", input_file,
        "--model", model_path,
        "-o", output_dir,
        "--input_type", input_type
    ]

    # If a path to bedtools is provided, add the parameter.
    if bedtools_path:
        cmd.extend(["--bedtools", bedtools_path])

    # Execute commands and capture errors
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result

def main():
    parser = argparse.ArgumentParser(
        description="MEnet",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 1. config
    parser.add_argument('-m', '--model', 
                        default='model/wgbs/CV_best_model.pickle', #select model
                        help="Model file path (.pickle)")

    parser.add_argument('-i', '--input_dir', 
                        default='test_data/wgbs/',
                        help="test_data")

    parser.add_argument('-o', '--output_dir', default='menet_result',
                        help="output_dir")

    # 2. run
    parser.add_argument('-t', '--target_folders', 
                        help="Specify subfolders, separated by commas (e.g., real, random). Leave blank to process all.")

    parser.add_argument('--input_type', default='bismark', help="Input file type (bismark/array)")

    # 3. Critical path repair parameters
    parser.add_argument('--menet_path', default="MEnet")
    
    parser.add_argument('--bedtools', default="bedtools")

    parser.add_argument('--log_name', default='menet_benchmark.csv', help="benchmark")

    args = parser.parse_args()

    # Automatically find MEnet in the current environment
    effective_menet = args.menet_path if args.menet_path else shutil.which("MEnet")
    # Automatically find bedtools in the current environment
    effective_bedtools = args.bedtools if args.bedtools else shutil.which("bedtools")

    # print
    print(f"==================================")
    print(f"Work environment confirmation:")
    print(f"  > MEnet path: {effective_menet}")
    print(f"  > Bedtools path: {effective_bedtools}")
    print(f"  > model file: {args.model}")
    print(f"  > data type: {args.input_type}")
    
    if not effective_menet:
        print("[Error] The MEnet executable cannot be found. Please activate the Conda environment or specify it via --menet_path.")
        sys.exit(1)
    if args.input_type == "bismark" and not effective_bedtools:
        print("[Warning] Processing Bismark format requires bedtools, but it is not found in the system, which may cause the program to report an error.")
    print(f"==================================\n")

    # Determine the subdirectories to be processed.
    if args.target_folders:
        target_list = [s.strip() for s in args.target_folders.split(',')]
        sub_dirs = [os.path.join(args.input_dir, d) for d in target_list]
    else:
        sub_dirs = [os.path.join(args.input_dir, d) for d in os.listdir(args.input_dir) 
                    if os.path.isdir(os.path.join(args.input_dir, d))]

    benchmark_summary = []

    for s_dir in sub_dirs:
        folder_name = os.path.basename(s_dir)
        input_files = glob.glob(os.path.join(s_dir, "*.cov.gz"))
        
        if not input_files:
            continue

        print(f"--- Processing directory: {folder_name} ({len(input_files)} files) ---")

        for f_path in input_files:
            file_full_name = os.path.basename(f_path)
            file_stem = file_full_name.replace(".cov.gz", "")
            
            specific_out_dir = os.path.join(args.output_dir, folder_name, file_stem)
            os.makedirs(specific_out_dir, exist_ok=True)

            print(f"  > Processing: {file_full_name}", end="", flush=True)
            
            gc.collect()
            time.sleep(0.5)

            start_time = time.perf_counter()
            status = "Success"
            error_msg = ""
            peak_memory = 0

            try:
                # Monitor memory and run
                mem_res = memory_usage(
                    (run_menet, (f_path, args.model, specific_out_dir, args.input_type, effective_menet, effective_bedtools)),
                    max_usage=True,
                    interval=0.1,
                    include_children=True
                )
                peak_memory = max(mem_res) if isinstance(mem_res, (list, tuple)) else mem_res

            except subprocess.CalledProcessError as e:
                status = "Failed"
                # Capture specific shell error messages
                error_msg = e.stderr if e.stderr else e.stdout
                print(f"\n    [Error] MEnet failure: {file_full_name}")
                if error_msg:
                    print(f"    reason: {error_msg.strip().splitlines()[-1]}") 
            except Exception as e:
                status = "Failed"
                error_msg = str(e)
                print(f"\n    [Error] System error: {e}")

            end_time = time.perf_counter()
            elapsed = end_time - start_time

            record = {
                "SubFolder": folder_name,
                "Sample": file_stem,
                "Time_Sec": round(elapsed, 2),
                "Peak_Mem_MB": round(peak_memory, 2),
                "Status": status,
                "Error": error_msg.replace('\n', ' ') if error_msg else ""
            }
            benchmark_summary.append(record)
            
            if status == "Success":
                print(f" | Done! {elapsed:.1f}s, {peak_memory:.1f}MB")

    # save bechamrk
    if benchmark_summary:
        summary_df = pd.DataFrame(benchmark_summary)
        summary_df.to_csv(args.log_name, index=False)
        print(f"\nTask completed! The statistics table has been saved.: {args.log_name}")

if __name__ == "__main__":
    main()
