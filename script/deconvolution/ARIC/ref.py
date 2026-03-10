import pandas as pd
import os

# =========================================
DATA_DIR = "test_data" # pending folder
OUTPUT_DIR = "aric_ref"        
REF_FILE = "ref_raw.csv"  # reference file
# =========================================

# 1. prepare
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Create output directory: {OUTPUT_DIR}")

# 2. read and preprocess Ref 
print(f"[*] Reading the main reference file: {REF_FILE}...")
ref_master = pd.read_csv(REF_FILE, index_col=0)
original_ref_shape = ref_master.shape
ref_master.dropna(inplace=True) # Ref remove NA
print(f"Ref original_ref_shape: {original_ref_shape} -> ref_master.shape: {ref_master.shape}")

# 3. test_data folder
files = [f for f in os.listdir(DATA_DIR) if f.endswith(".csv")]
print(f"[*] '{DATA_DIR}' find {len(files)} CSV files。")

for file_name in files:
    file_path = os.path.join(DATA_DIR, file_name)
    print(f"\n{'='*40}")
    print(f"processing file: {file_name}")
    
    try:
        # read Mix data
        mix = pd.read_csv(file_path, index_col=0)
        
        # ----------------------------
        # A. delete NA
        # ----------------------------
        mix.dropna(inplace=True)
        
        # ----------------------------
        # B. find rows that are not in the intersection and find the intersection.
        # ----------------------------
        common_cpg = mix.index.intersection(ref_master.index)
        
        # if the intersection is empty, skip.
        if len(common_cpg) == 0:
            print(f"  [warning] {file_name} is no common line (index) with the reference file")
            continue

        # Preserve intersection
        mix_fixed = mix.loc[common_cpg]
        ref_fixed = ref_master.loc[common_cpg]
        
        # ----------------------------
        # C. save file
        # ----------------------------
        
        mix_out_name = f"mix_{file_name}"
        ref_out_name = f"ref_{file_name}"
        
        mix_out_path = os.path.join(OUTPUT_DIR, mix_out_name)
        ref_out_path = os.path.join(OUTPUT_DIR, ref_out_name)
        
        mix_fixed.to_csv(mix_out_path)
        ref_fixed.to_csv(ref_out_path)
        
        # ----------------------------
        # D. check
        # ----------------------------
        print(f"  -> done")
        print(f"  -> final shape: {mix_fixed.shape}")
        print(f"  -> save: {mix_out_path} 和 {ref_out_path}")
        
        # check negative values
        neg_mix = (mix_fixed < 0).sum().sum()
        neg_ref = (ref_fixed < 0).sum().sum()
        if neg_mix > 0 or neg_ref > 0:
            print(f"  -> [tip] negative values: mix={neg_mix}, ref={neg_ref}")
            
        # check line name consistency
        if not mix_fixed.index.equals(ref_fixed.index):
            print("  -> [error] Inconsistent line names！")

        # check ref
        ref_var = ref_fixed.var(axis=1)
        low_var_count = (ref_var < 1e-6).sum()
        if low_var_count > 0:
            print(f"  -> [tip] CpG numbers with extremely low variance (<1e-6): {low_var_count}")

    except Exception as e:
        print(f"  [Error] process {file_name} error: {e}")

print(f"\n{'='*40}")
print("Done!")
