#!/bin/bash

# ================= 用户配置 =================
# 原始解卷积脚本
CELFIE_SCRIPT="/data/yuxy/data/wgbs/celfie/celfie/scripts/celfie_new.py"
# 性能监控脚本
PERF_WRAPPER="/data/yuxy/data/wgbs/celfie/celfie/scripts/celfie_perf.py"
# 输入总目录
INPUT_BASE="/data/yuxy/data/wgbs/celfie/pre/2_merge_data/more"
# 输出总目录
OUTPUT_BASE="/data/yuxy/data/wgbs/celfie/data/results"
# 参数 K (Mixtures数量)
K_VALUE=100

# ================= 检查 =================
if [ ! -f "$PERF_WRAPPER" ]; then
    echo "错误: 找不到性能监控脚本 $PERF_WRAPPER"
    exit 1
fi

# ================= 批量遍历 =================
# 查找 INPUT_BASE 下的所有子文件夹
folders=$(ls -d ${INPUT_BASE}/*/)

echo "开始批量解卷积处理..."
echo "时间: $(date)"

for folder_path in $folders; do
    # 提取文件夹名称 (如 real, random_1 等)
    folder_name=$(basename "$folder_path")
    
    # 寻找文件夹下的 .tsv 文件
    input_file=$(ls ${folder_path}/*.tsv 2>/dev/null | head -n 1)
    
    if [ -z "$input_file" ]; then
        echo "跳过: $folder_name (未找到 .tsv 文件)"
        continue
    fi

    # 定义该样本的输出目录
    current_out="${OUTPUT_BASE}/${folder_name}"
    mkdir -p "$current_out"

    echo "-------------------------------------------------------"
    echo "正在处理任务: $folder_name"
    echo "输入路径: $input_file"
    echo "保存路径: $current_out"

    # 调用性能监控包装器运行 Python
    # 参数: 包装器 脚本路径 输入文件 输出目录 K值
    python "$PERF_WRAPPER" "$CELFIE_SCRIPT" "$input_file" "$current_out" "$K_VALUE"

    if [ $? -eq 0 ]; then
        echo "任务 $folder_name 成功完成"
    else
        echo "任务 $folder_name 失败"
    fi
done

echo "-------------------------------------------------------"
echo "所有批量任务完成！"
echo "总体完成时间: $(date)"
