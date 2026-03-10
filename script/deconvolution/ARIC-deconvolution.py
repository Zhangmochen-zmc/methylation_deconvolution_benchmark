import os
import glob
import time
import gc  
import pandas as pd 
from memory_profiler import memory_usage
from ARIC import *  

def run_benchmark():
    # --- 1. 路径配置 ---
    mix_data_dir = "mix_data"
    ref_file = "ref_raw.csv"
    results_dir = "results"
    
    # 用于记录性能指标的列表
    benchmark_records = []

    # --- 2. 准备输出目录 ---
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    # --- 3. 获取文件列表 ---
    mix_files = glob.glob(os.path.join(mix_data_dir, "*.csv"))
    print(f"==========================================")
    print(f"检测到 {len(mix_files)} 个文件，开始批量处理...")
    print(f"==========================================\n")

    # --- 4. 循环遍历 ---
    for i, mix_path in enumerate(mix_files):
        # 4.1 构造文件名和保存路径
        file_name = os.path.basename(mix_path)
        base_name = os.path.splitext(file_name)[0]
        # 保持原命名习惯: results/文件名_prop.csv
        save_path = os.path.join(results_dir, f"{base_name}_prop.csv")

        print(f"[{i+1}/{len(mix_files)}] 正在处理: {file_name}")
        
        # 4.2 【关键】强制垃圾回收
        # 这步操作是为了清理上一个循环遗留的内存，让本次测量的内存基线更干净
        gc.collect() 
        time.sleep(1) # 稍微暂停一下让系统稳定

        # 4.3 准备参数
        # memory_usage 要求参数以 (函数, 位置参数元组, 关键字参数字典) 传递
        aric_kwargs = {
            'mix_path': mix_path,
            'ref_path': ref_file,
            'save_path': save_path,
            'is_methylation': True
        }

        # 4.4 开始监控 (时间和内存)
        # 记录开始时间 (Wall time)
        start_time = time.perf_counter()
        
        peak_memory_mb = 0
        error_msg = "None"
        
        try:
            # 运行核心函数并监控内存峰值
            # max_usage=True: 返回执行期间达到的最大内存 (MB)
            # interval=0.1: 采样频率 0.1秒
            # include_children=True: 如果ARIC内部调用了子进程，也能尽量捕捉（取决于系统）
            mem_max = memory_usage(
                (ARIC, [], aric_kwargs), 
                max_usage=True, 
                interval=0.1, 
                retval=False
            )
            
            # 处理返回值格式（防止版本差异返回list）
            if isinstance(mem_max, (list, tuple)):
                peak_memory_mb = max(mem_max)
            else:
                peak_memory_mb = mem_max

        except Exception as e:
            print(f"  >>> 错误: {e}")
            error_msg = str(e)
            peak_memory_mb = 0 # 出错则内存数据可能不准确

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
            "Save_Path": save_path,
            "Time_Seconds": round(elapsed_time, 4),
            "Peak_Memory_MB": round(peak_memory_mb, 4),
            "Status": "Success" if error_msg == "None" else "Failed",
            "Error": error_msg
        })

    # --- 6. 保存所有统计结果到 CSV ---
    log_file = "benchmark_summary.csv"
    try:
        df = pd.DataFrame(benchmark_records)
        df.to_csv(log_file, index=False)
        print(f"\n全部完成！性能统计已保存至: {log_file}")
    except Exception as e:
        print(f"\n全部完成！但保存性能统计表失败: {e}")
        # 如果没有pandas，这里只是打印出来供参考
        print(benchmark_records)

if __name__ == "__main__":
    run_benchmark()
