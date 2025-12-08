import pandas as pd
import numpy as np
import os
import glob

#path of the proportion file 
prop_file_path = "../proposion/450klow_nk.csv"

#path of input files
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
    "index_col": 0,      
    "header": 0         
}

ref_dfs = {}
common_cpgs = None

for cell_type, file_path in test_files.items():
    print(f"  -> read {cell_type}: {file_path}")

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Can't find the file: {file_path}")

    df = pd.read_csv(file_path, **read_params)

    if common_cpgs is None:
        common_cpgs = df.index
    else:
        common_cpgs = common_cpgs.intersection(df.index)

    ref_dfs[cell_type] = df
    print(f"     Number of samples: {df.shape[1]}, Number of Public CpG: {len(common_cpgs)}")

print(f"\nNumber of Public CpG sites: {len(common_cpgs)}")

# Trim matrices
for cell_type in ref_dfs:
    ref_dfs[cell_type] = ref_dfs[cell_type].loc[common_cpgs]

proportions = pd.read_csv(prop_file_path)
n_samples = len(proportions)
print(f"\nStart generating {n_samples} simulated samples...")

simulated_data = pd.DataFrame(
    np.zeros((len(common_cpgs), n_samples)),
    index=common_cpgs,
    columns=[f"Simulated_Sample_{i+1}" for i in range(n_samples)]
)

for cell_type, ref_matrix in ref_dfs.items():
    if cell_type not in proportions.columns:
        print(f"Warning: The proportion file is missing the {cell_type} column. Skipping.")
        continue
        
    print(f"  Processing: {cell_type} ...")

    available_samples = ref_matrix.columns
    
    chosen_samples = np.random.choice(available_samples, size=n_samples)
    selected_values = ref_matrix[chosen_samples].values
    props = proportions[cell_type].values
    weighted_values = selected_values * props
    
    simulated_data.values[:, :] += weighted_values


output_file = "mix_data/simulated_matrix_450klow_nk.csv"
print(f"\nSaving the result to {output_file} ...")

simulated_data.to_csv(output_file, float_format='%.4f')

print("Completed！")
