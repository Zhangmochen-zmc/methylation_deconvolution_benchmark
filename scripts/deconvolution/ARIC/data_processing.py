import os
import pandas as pd

base_dir = "methylation_deconvolution_benchmark/data/reference_data/450k/original_data"

celltype_avg_list = []
cell_types = []

for cell_type in os.listdir(base_dir):
    cell_path = os.path.join(base_dir, cell_type)
    if os.path.isdir(cell_path):
        print(f"Processing cell type: {cell_type}")
        cell_types.append(cell_type)
        beta_dfs = []
        for txt_file in os.listdir(cell_path):
            if txt_file.endswith(".txt"):
                file_path = os.path.join(cell_path, txt_file)
                try:
                    df = pd.read_csv(file_path, sep="\t", header=None, names=['CpG', 'Beta'])
                    beta_dfs.append(df)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
        if beta_dfs:
            avg_df = pd.concat(beta_dfs).groupby('CpG', as_index=False).mean()
            avg_df.rename(columns={'Beta': cell_type}, inplace=True)
            celltype_avg_list.append(avg_df[['CpG', cell_type]])

# Combined average values ​​of all cell types
if celltype_avg_list:
    ref_df = celltype_avg_list[0]
    for df in celltype_avg_list[1:]:
        ref_df = ref_df.merge(df, on='CpG', how='outer')

    # Set the column names: the first column is ID, and the rest are cell types.
    ref_df.columns = ['ID'] + cell_types

    # save CSV
    ref_df.to_csv("ref_data.csv", index=False)
    print("ref_data.csv generated successfully with column names!")
else:
    print("No data found.")
