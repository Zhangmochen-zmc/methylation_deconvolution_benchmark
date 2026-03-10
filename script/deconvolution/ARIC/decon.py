import os
import glob
import time
import gc 
import pandas as pd 
from memory_profiler import memory_usage
from ARIC import *  

def run_benchmark():
    # 1. path configuration 
    data_dir = ""  
    results_dir = "results_random_5"
    
    # 用于记录性能指标的列表
    benchmark_records = []

    # --- 2. 准备输出目录 ---
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # --- 3. 获取 mix 文件列表 ---
    # 只获取 mix_ 开头的文件，然后去寻找对应的 ref_
    mix_files = glob.glob(os.path.join(data_dir, "mix_*.csv"))
    
    print(f"==========================================")
    print(f"在 '{data_dir}' 中检测到 {len(mix_files)} 个 mix 文件，开始批量处理...")
    print(f"==========================================\n")

    # --- 4. 循环遍历 ---
    for i, mix_path in enumerate(mix_files):
        # 4.1 构造文件名和保存路径
        file_name = os.path.basename(mix_path) # 例如: mix_GSE125105.csv
        base_name = os.path.splitext(file_name)[0] # 例如: mix_GSE125105
        
        # 4.1.1 【关键修改】动态匹配对应的 ref 文件
        # 逻辑：将文件名中的 "mix_" 替换为 "ref_"
        if file_name.startswith("mix_"):
            ref_file_name = file_name.replace("mix_", "ref_", 1)
        else:
            # 防止文件名不规范（虽然glob已经过滤了，但为了保险）
            ref_file_name = "ref_" + file_name
            
        ref_path = os.path.join(data_dir, ref_file_name)

        # 4.1.2 检查对应的 Ref 文件是否存在
        if not os.path.exists(ref_path):
            print(f"[{i+1}/{len(mix_files)}] 跳过: {file_name}")
            print(f"      >>> 错误: 找不到对应的 Ref 文件: {ref_file_name}")
            benchmark_records.append({
                "File": file_name,
                "Status": "Skipped",
                "Error": f"Ref file {ref_file_name} not found"
            })
            continue

        # 保持原命名习惯: results/mix_GSE125105_prop.csv
        save_path = os.path.join(results_dir, f"{base_name}_prop.csv")

        print(f"[{i+1}/{len(mix_files)}] 正在处理配对:")
        print(f"      Mix: {file_name}")
        print(f"      Ref: {ref_file_name}")
        
        # 4.2 【关键】强制垃圾回收
        gc.collect() 
        time.sleep(1) 

        # 4.3 准备参数
        # 这里将动态获取的 ref_path 传入
        aric_kwargs = {
            'mix_path': mix_path,
            'ref_path': ref_path, 
            'save_path': save_path,
            'is_methylation': True
        }

        # 4.4 开始监控 (时间和内存)
        start_time = time.perf_counter()
        
        peak_memory_mb = 0
        error_msg = "None"
        
        try:
            # 运行核心函数并监控内存峰值
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
            print(f"      >>> 运行错误: {e}")
            error_msg = str(e)
            peak_memory_mb = 0 

        # 记录结束时间
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        # --- 5. 输出当前文件的结果 ---
        print(f"  -> 耗时 (Time):   {elapsed_time:.4f} s")
        print(f"  -> 峰值内存 (Mem): {peak_memory_mb:.4f} MB")
        print(f"  -> 结果已保存至:   {save_path}")
        print("-" * 40)

        # 记录到列表
        benchmark_records.append({
            "File": file_name,
            "Ref_File": ref_file_name, # 记录一下用了哪个ref
            "Save_Path": save_path,
            "Time_Seconds": round(elapsed_time, 4),
            "Peak_Memory_MB": round(peak_memory_mb, 4),
            "Status": "Success" if error_msg == "None" else "Failed",
            "Error": error_msg
        })

    # --- 6. 保存所有统计结果到 CSV ---
    log_file = "benchmark.csv"
    try:
        df = pd.DataFrame(benchmark_records)
        df.to_csv(log_file, index=False)
        print(f"\n全部完成！性能统计已保存至: {log_file}")
    except Exception as e:
        print(f"\n全部完成！但保存性能统计表失败: {e}")
        print(benchmark_records)

if __name__ == "__main__":
    run_benchmark()
