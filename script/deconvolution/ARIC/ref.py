import pandas as pd
import os

# =========================================
DATA_DIR = "test_data" # pending folder
OUTPUT_DIR = "aric_ref"        
REF_FILE = "ref.csv"  # reference file
# =========================================

# 1.prepare
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Create output directory: {OUTPUT_DIR}")

# 2. read and preprocess Ref 
print(f"[*] Reading the main reference file: {REF_FILE}...")
ref_master = pd.read_csv(REF_FILE, index_col=0)
original_ref_shape = ref_master.shape
ref_master.dropna(inplace=True) # Ref remove NA
print(f"Ref original_ref_shape: {original_ref_shape} -> ref_master.shape: {ref_master.shape}")

# 2. test_data folder
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".csv")]
print(f"[*] 在 '{DATA_DIR}' 下找到 {len(files)} 个 CSV 文件待处理。")

for file_name in files:
    file_path = os.path.join(DATA_DIR, file_name)
    print(f"\n{'='*40}")
    print(f"正在处理文件: {file_name}")
    
    try:
        # 读取 Mix 数据
        mix = pd.read_csv(file_path, index_col=0)
        
        # ----------------------------
        # A. 删除包含 NA 的行
        # ----------------------------
        mix.dropna(inplace=True)
        # Ref 已经在循环外去除了 NA
        
        # ----------------------------
        # B. 找出不在交集中的行并取交集
        # ----------------------------
        # 注意：这里是用当前的 mix 和 主 ref 取交集
        common_cpg = mix.index.intersection(ref_master.index)
        
        # 如果交集为空，跳过
        if len(common_cpg) == 0:
            print(f"  [警告] {file_name} 与参考文件没有公共行(index)，跳过！")
            continue

        # 保留交集
        mix_fixed = mix.loc[common_cpg]
        ref_fixed = ref_master.loc[common_cpg]
        
        # ----------------------------
        # C. 保存文件
        # ----------------------------
        # 为了区分不同 mix 对应的 ref，我们需要重命名
        # 例如: data/sample1.csv -> ref/mix_sample1.csv 和 ref/ref_sample1.csv
        
        mix_out_name = f"mix_{file_name}"
        ref_out_name = f"ref_{file_name}"
        
        mix_out_path = os.path.join(OUTPUT_DIR, mix_out_name)
        ref_out_path = os.path.join(OUTPUT_DIR, ref_out_name)
        
        mix_fixed.to_csv(mix_out_path)
        ref_fixed.to_csv(ref_out_path)
        
        # ----------------------------
        # D. 统计与检查 (针对当前文件)
        # ----------------------------
        print(f"  -> 处理完毕")
        print(f"  -> 最终形状: {mix_fixed.shape}")
        print(f"  -> 保存至: {mix_out_path} 和 {ref_out_path}")
        
        # 检查负值
        neg_mix = (mix_fixed < 0).sum().sum()
        neg_ref = (ref_fixed < 0).sum().sum()
        if neg_mix > 0 or neg_ref > 0:
            print(f"  -> [提示] 负值数量: mix={neg_mix}, ref={neg_ref}")
            
        # 检查行名一致性
        if not mix_fixed.index.equals(ref_fixed.index):
            print("  -> [错误] 行名不一致！")

        # 检查方差 (仅检查 Ref)
        ref_var = ref_fixed.var(axis=1)
        low_var_count = (ref_var < 1e-6).sum()
        if low_var_count > 0:
            print(f"  -> [提示] 方差极低 (<1e-6) 的 CpG 数: {low_var_count}")

    except Exception as e:
        print(f"  [Error] 处理 {file_name} 时发生错误: {e}")

print(f"\n{'='*40}")
print("所有文件处理完成。")
