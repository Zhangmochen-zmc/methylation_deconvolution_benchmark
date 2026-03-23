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
    # 构建基础命令
    cmd = [
        menet_path, "predict",
        "--input", input_file,
        "--model", model_path,
        "-o", output_dir,
        "--input_type", input_type
    ]

    # 如果提供了 bedtools 路径，则加入参数
    if bedtools_path:
        cmd.extend(["--bedtools", bedtools_path])

    # 执行命令并捕获错误
    # stderr=subprocess.STDOUT 将错误信息合并到 stdout，方便 memory_usage 捕获或记录
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result

def main():
    parser = argparse.ArgumentParser(
        description="MEnet WGBS 批量处理工具 (自动兼容 Conda 环境)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 1. 路径配置
    parser.add_argument('-m', '--model', 
                        default='/disk1/yuxy/software/MEnet_8/data/train/wgbs/CV_best_model.pickle',
                        help="模型文件路径 (.pickle)")

    parser.add_argument('-i', '--input_dir', 
                        default='/data/yuxy/data/menet/wgbs/data/',
                        help="输入父目录路径")

    parser.add_argument('-o', '--output_dir', default='results',
                        help="总输出目录")

    # 2. 运行控制
    parser.add_argument('-t', '--target_folders', 
                        help="指定子文件夹，逗号隔开 (如: real,random)。留空则处理所有。")

    parser.add_argument('--input_type', default='bismark', help="输入文件类型 (bismark/array)")

    # 3. 关键路径修复参数
    parser.add_argument('--menet_path', default="/data/yuxy/work/miniforge3/envs/menet/bin/MEnet")
    
    parser.add_argument('--bedtools', default="/data/yuxy/work/miniforge3/envs/menet/bin/bedtools")

    parser.add_argument('--log_name', default='benchmark_summary.csv', help="统计结果保存的文件名")

    args = parser.parse_args()

    # --- 自动路径检测逻辑 ---
    # 自动寻找当前环境下的 MEnet
    effective_menet = args.menet_path if args.menet_path else shutil.which("MEnet")
    # 自动寻找当前环境下的 bedtools
    effective_bedtools = args.bedtools if args.bedtools else shutil.which("bedtools")

    # 打印环境检查
    print(f"==================================")
    print(f"工作环境确认:")
    print(f"  > MEnet 路径: {effective_menet}")
    print(f"  > Bedtools 路径: {effective_bedtools}")
    print(f"  > 模型文件: {args.model}")
    print(f"  > 数据类型: {args.input_type}")
    
    if not effective_menet:
        print("[Error] 找不到 MEnet 可执行程序，请激活 Conda 环境或通过 --menet_path 指定。")
        sys.exit(1)
    if args.input_type == "bismark" and not effective_bedtools:
        print("[Warning] 处理 Bismark 格式需要 bedtools，但系统中未找到，程序可能会报错。")
    print(f"==================================\n")

    # 确定要处理的子目录
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

        print(f"--- 正在处理目录: {folder_name} ({len(input_files)} 个文件) ---")

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
                # 监控内存并运行
                mem_res = memory_usage(
                    (run_menet, (f_path, args.model, specific_out_dir, args.input_type, effective_menet, effective_bedtools)),
                    max_usage=True,
                    interval=0.1,
                    include_children=True
                )
                peak_memory = max(mem_res) if isinstance(mem_res, (list, tuple)) else mem_res

            except subprocess.CalledProcessError as e:
                status = "Failed"
                # 捕捉具体的 shell 报错内容
                error_msg = e.stderr if e.stderr else e.stdout
                print(f"\n    [Error] MEnet 运行失败: {file_full_name}")
                if error_msg:
                    print(f"    原因: {error_msg.strip().splitlines()[-1]}") # 只打印最后一行错误
            except Exception as e:
                status = "Failed"
                error_msg = str(e)
                print(f"\n    [Error] 系统错误: {e}")

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
            
            # 每个文件处理完打个勾，显示耗时
            if status == "Success":
                print(f" | Done! {elapsed:.1f}s, {peak_memory:.1f}MB")

    # 保存最终统计
    if benchmark_summary:
        summary_df = pd.DataFrame(benchmark_summary)
        summary_df.to_csv(args.log_name, index=False)
        print(f"\n任务完成！统计表已保存至: {args.log_name}")

if __name__ == "__main__":
    main()
