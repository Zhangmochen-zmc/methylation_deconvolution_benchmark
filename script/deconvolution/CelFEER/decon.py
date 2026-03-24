import time
import gc
import pandas as pd
import subprocess
import os
from memory_profiler import memory_usage

# ==================================================
# 配置参数 (与原 Shell 脚本一致)
# ==================================================
INPUT_PATH = "4_data/random_35_fixed.txt"
OUTPUT_DIR = "results/random_35"
NUM_SAMPLES = "100"
SCRIPT_PATH = "../CelFEER/scripts/celfeer.py"

def run_task():
    """
    需要执行的核心任务
    """
    # 构建命令
    cmd = ["python", SCRIPT_PATH, INPUT_PATH, OUTPUT_DIR, NUM_SAMPLES]
    
    # 运行子进程
    # check=True 如果报错会抛出异常
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    # 1. 清理内存，准备开始
    gc.collect()
    
    print(f"开始运行任务: {SCRIPT_PATH} ...")
    
    # 2. 记录起始时间
    start_perf = time.perf_counter()
    start_time_str = time.strftime("%Y-%m-%d %H:%M:%S")

    # 3. 监控内存使用情况
    # memory_usage 会每隔一段时间(interval)采样一次内存
    # retval, mem_history = ... 如果想获取函数返回值可以用这种形式
    try:
        mem_history = memory_usage(run_task, interval=0.1, timeout=None)
    except Exception as e:
        print(f"运行失败: {e}")
        mem_history = [0]

    # 4. 记录结束时间
    end_perf = time.perf_counter()
    
    # 5. 计算指标
    duration = end_perf - start_perf
    max_memory = max(mem_history) if mem_history else 0
    avg_memory = sum(mem_history) / len(mem_history) if mem_history else 0

    # 6. 使用 Pandas 整理结果
    results = {
        "Task": [os.path.basename(SCRIPT_PATH)],
        "Start_Time": [start_time_str],
        "Duration_Sec": [round(duration, 4)],
        "Max_Memory_MiB": [round(max_memory, 2)],
        "Avg_Memory_MiB": [round(avg_memory, 2)],
        "Input_File": [INPUT_PATH]
    }
    
    df = pd.DataFrame(results)
    
    # 7. 打印并保存到 CSV
    print("\n--- 性能报告 ---")
    print(df.to_string(index=False))
    
    # 如果文件不存在则写入，存在则追加（不写表头）
    log_file = "performance_log.csv"
    if not os.path.exists(log_file):
        df.to_csv(log_file, index=False)
    else:
        df.to_csv(log_file, mode='a', header=False, index=False)
        
    print(f"\n性能记录已更新至: {log_file}")
