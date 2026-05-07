import os
import glob
import time
import gc 
import pandas as pd 
from memory_profiler import memory_usage
from ARIC import *  

def run_benchmark():
    # 1. path configuration 
    data_dir = "aric_ref"  
    results_dir = "aric_result"
    
    # list used for recording performance metrics
    benchmark_records = []

    # 2. prepare output directory 
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # 3. get list of mix files 
    mix_files = glob.glob(os.path.join(data_dir, "mix_*.csv"))
    
    print(f"==========================================")
    print(f"{len(mix_files)} mix files were detected in '{data_dir}', batch processing begins...")
    print(f"==========================================\n")

    # 4. looping through 
    for i, mix_path in enumerate(mix_files):
        # 4.1 create filename and save path
        file_name = os.path.basename(mix_path) # such as: mix_GSE125105.csv
        base_name = os.path.splitext(file_name)[0] # such as: mix_GSE125105
        
        # 4.1.1 match the corresponding ref file
        if file_name.startswith("mix_"):
            ref_file_name = file_name.replace("mix_", "ref_", 1)
        else:
            ref_file_name = "ref_" + file_name
            
        ref_path = os.path.join(data_dir, ref_file_name)

        # 4.1.2 check if the corresponding ref file exists.
        if not os.path.exists(ref_path):
            print(f"[{i+1}/{len(mix_files)}] skip: {file_name}")
            print(f">>> error: not corresponding ref file: {ref_file_name}")
            benchmark_records.append({
                "File": file_name,
                "Status": "Skipped",
                "Error": f"Ref file {ref_file_name} not found"
            })
            continue

        save_path = os.path.join(results_dir, f"{base_name}_prop.csv")

        print(f"[{i+1}/{len(mix_files)}] pair:")
        print(f"      Mix: {file_name}")
        print(f"      Ref: {ref_file_name}")
        
        # 4.2 gc and time
        gc.collect() 
        time.sleep(1) 

        # 4.3 prepare params
        aric_kwargs = {
            'mix_path': mix_path,
            'ref_path': ref_path, 
            'save_path': save_path,
            'is_methylation': True
        }

        # 4.4 time and memory
        start_time = time.perf_counter()
        
        peak_memory_mb = 0
        error_msg = "None"
        
        try:
            mem_max = memory_usage(
                (ARIC, [], aric_kwargs), 
                max_usage=True, 
                interval=0.1, 
                retval=False
            )
            
            if isinstance(mem_max, (list, tuple)):
                peak_memory_mb = max(mem_max)
            else:
                peak_memory_mb = mem_max

        except Exception as e:
            print(f">>> error: {e}")
            error_msg = str(e)
            peak_memory_mb = 0 

        # ending time
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        # 5. results
        print(f"results: {save_path}")
        print("-" * 40)

        # record
        benchmark_records.append({
            "File": file_name,
            "Ref_File": ref_file_name,
            "Save_Path": save_path,
            "Time_Seconds": round(elapsed_time, 4),
            "Peak_Memory_MB": round(peak_memory_mb, 4),
            "Status": "Success" if error_msg == "None" else "Failed",
            "Error": error_msg
        })

    # 6. benchmark_results
    log_file = "aric_benchmark.csv"
    try:
        df = pd.DataFrame(benchmark_records)
        df.to_csv(log_file, index=False)
        print(f"\n done! bechmark results: {log_file}")
    except Exception as e:
        print(f"\n done! bechmark results fail: {e}")
        print(benchmark_records)

if __name__ == "__main__":
    run_benchmark()
