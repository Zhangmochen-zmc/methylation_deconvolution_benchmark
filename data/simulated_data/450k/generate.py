import pandas as pd
import numpy as np
import os
import glob

#比例文件路径
prop_file_path = "../proposion/450klow_nk.csv"

#细胞类型文件路径
test_files = {
    "B cells":              r"/data/weiyk/data/test_data/array_data/blood/GSE88824/bcell.txt",
    "CD4+ T cells":         r"/data/weiyk/data/test_data/array_data/blood/GSE88824/cd4.txt",
    "CD8+ T cells":         r"/data/weiyk/data/test_data/array_data/blood/GSE88824/cd8.txt",
    "Monocytes":            r"/data/weiyk/data/test_data/array_data/blood/GSE88824/monocyte.txt",
    "Neutrophils":          r"/data/weiyk/data/test_data/array_data/blood/GSE88824/neutrophil.txt",
    "Natural Killer cells": r"/data/weiyk/data/test_data/array_data/blood/GSE88824/nk.txt"
}

read_params = {
    "sep": "\t",
    "index_col": 0,      # 第一列作为索引 (CpG ID)
    "header": 0          # 假设第一行是列名 (ID, Beta_Value)，如果没有列名设为 None
}

#加载所有数据
ref_dfs = {}
common_cpgs = None

for cell_type, file_path in test_files.items():
    print(f"  -> 读取 {cell_type}: {file_path}")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"找不到文件: {file_path}")

    # 读取数据
    df = pd.read_csv(file_path, **read_params)

    # 第一次加载时，初始化公共位点集合
    if common_cpgs is None:
        common_cpgs = df.index
    else:
        # 取交集：只保留所有矩阵里都有的 CpG 位点
        common_cpgs = common_cpgs.intersection(df.index)

    ref_dfs[cell_type] = df
    print(f"     样本数: {df.shape[1]}, 当前公共CpG数: {len(common_cpgs)}")

print(f"\n对齐完成！最终使用的公共 CpG 位点数量: {len(common_cpgs)}")

# 统一裁剪所有矩阵，只保留公共行，确保可以直接相加
for cell_type in ref_dfs:
    ref_dfs[cell_type] = ref_dfs[cell_type].loc[common_cpgs]


# 读取比例文件
proportions = pd.read_csv(prop_file_path)
n_samples = len(proportions)
print(f"\n开始生成 {n_samples} 个模拟样本...")

# 创建一个空的 DataFrame 存放结果 (行=CpG, 列=模拟样本名)
simulated_data = pd.DataFrame(
    np.zeros((len(common_cpgs), n_samples)),
    index=common_cpgs,
    columns=[f"Simulated_Sample_{i+1}" for i in range(n_samples)]
)

# 遍历每一个细胞类型进行批量计算 (向量化操作，速度更快)
for cell_type, ref_matrix in ref_dfs.items():
    # 检查该细胞类型在比例文件中是否存在
    if cell_type not in proportions.columns:
        print(f"警告: 比例文件中缺少 {cell_type} 列，跳过。")
        continue
        
    print(f"  正在处理成分: {cell_type} ...")
    
    # 获取该细胞类型的所有可用样本名
    available_samples = ref_matrix.columns
    
    # --- 核心逻辑 ---
    # 为这 100 个模拟样本，每一个都随机挑一个真实样本的列名
    # size=n_samples: 一次性生成 100 个随机选择
    chosen_samples = np.random.choice(available_samples, size=n_samples)
    
    # 提取这些被选中的列，组成一个新的临时矩阵 (CpG x 100)
    # values 属性将其转换为 numpy 数组以便计算
    selected_values = ref_matrix[chosen_samples].values
    
    # 获取对应的比例 (长度为 100 的向量)
    props = proportions[cell_type].values
    
    # 利用广播机制进行加权
    # selected_values 是 (N_CpG, 100)
    # props 是 (100,)
    # Python 会自动把 props 扩展到每一行进行相乘
    weighted_values = selected_values * props
    
    # 累加到最终结果中
    simulated_data.values[:, :] += weighted_values


output_file = "mix_data/simulated_matrix_450klow_nk.csv"
print(f"\n正在保存结果到 {output_file} ...")

simulated_data.to_csv(output_file, float_format='%.4f')

print("全部完成！")
