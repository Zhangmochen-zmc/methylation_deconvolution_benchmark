import pandas as pd
import argparse
import sys
import numpy as np

def run_pipeline(input_data_path, ref_path, output_path, mode='clean'):
    """
    mode: 
      'full'  - Keeps 'chr' and 'pos' columns at the beginning (for analysis/viz).
      'clean' - Removes all non-numeric columns (for deconvolution computation).
    """
    
    # ---------------------------------------------------------
    # 1. Load Data
    # ---------------------------------------------------------
    print(f"[-] Step 1: Reading methylation data...")
    try:
        # encoding='utf-8-sig' handles potential BOM in CSV files
        df_data = pd.read_csv(input_data_path, index_col=0, encoding='utf-8-sig')
    except Exception as e:
        print(f"[Error] Failed to read data file: {e}")
        return

    print(f"[-] Step 2: Reading reference file (Manifest)...")
    try:
        # Assuming the manifest is a tab-separated file (TSV)
        df_ref = pd.read_csv(ref_path, sep='\t', low_memory=False)
    except Exception as e:
        print(f"[Error] Failed to read reference file: {e}")
        return

    # ---------------------------------------------------------
    # 3. ID Matching and Coordinate Extraction
    # ---------------------------------------------------------
    # Trim whitespace from column names
    df_ref.columns = [c.strip() for c in df_ref.columns]
    target_id_col = 'Probe_ID'
    
    if target_id_col not in df_ref.columns:
        print(f"[Error] Could not find '{target_id_col}' column in reference file.")
        print(f"    Available columns: {list(df_ref.columns)}")
        sys.exit(1)

    df_ref.set_index(target_id_col, inplace=True)
    
    print("[-] Step 3: Extracting coordinates and performing 0-to-1 based conversion...")
    if 'CpG_chrm' not in df_ref.columns or 'CpG_beg' not in df_ref.columns:
        print("[Error] Missing required columns: 'CpG_chrm' or 'CpG_beg'.")
        sys.exit(1)
        
    # Copy coordinate info
    coords = df_ref[['CpG_chrm', 'CpG_beg']].copy()
    # Convert position to numeric, handle errors
    coords['pos'] = pd.to_numeric(coords['CpG_beg'], errors='coerce')
    coords = coords.dropna(subset=['pos'])
    # Convert to 1-based coordinate and ensure integer type
    coords['pos'] = (coords['pos'] + 1).astype(int) 
    coords.rename(columns={'CpG_chrm': 'chr'}, inplace=True)

    # ---------------------------------------------------------
    # 4. Merging and Autosome Filtering (chr1-22)
    # ---------------------------------------------------------
    print("[-] Step 4: Merging data and filtering for autosomes (chr1-22)...")
    df_data.index.name = 'Probe_ID'
    # Inner join ensures only probes with known coordinates are kept
    merged_df = df_data.join(coords[['chr', 'pos']], how='inner')
    
    # List of valid autosomes
    valid_chrs = [f'chr{i}' for i in range(1, 23)] + [str(i) for i in range(1, 23)]
    merged_df['chr'] = merged_df['chr'].astype(str)
    final_df = merged_df[merged_df['chr'].isin(valid_chrs)]
    
    print(f"    - Original probes: {len(df_data)}")
    print(f"    - Remaining probes after filtering: {len(final_df)}")

    # ---------------------------------------------------------
    # 5. Final Formatting and Export
    # ---------------------------------------------------------
    if mode == 'full':
        print("[-] Mode: [FULL] - Preserving coordinate columns...")
        # Move chr and pos to the front
        cols = final_df.columns.tolist()
        cols.remove('chr')
        cols.remove('pos')
        final_df = final_df[['chr', 'pos'] + cols]
    else:
        print("[-] Mode: [CLEAN] - Removing all non-numeric columns (coordinates/text)...")
        # Define keywords for columns to drop
        bad_cols_keywords = ['chr', 'pos', 'chromosome', 'position', 'strand']
        cols_to_drop = [c for c in final_df.columns if c.lower() in bad_cols_keywords]
        final_df.drop(columns=cols_to_drop, inplace=True)
        # Select only numeric data types to ensure a pure matrix
        final_df = final_df.select_dtypes(include=[np.number])

    print(f"[-] Saving processed data to: {output_path}")
    final_df.to_csv(output_path)
    print(f"    Final Dimensions: {final_df.shape}")
    print("[Success] Pipeline completed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Methylation Data Preprocessing Pipeline")
    parser.add_argument("input", default="ref_data.csv",help="Path to the CpG signal value CSV file")
    parser.add_argument("ref", default="HM450.hg38.manifest", help="Path to the Manifest reference TSV file")
    parser.add_argument("-o", "--output", default="marker_ref/ref.csv", help="Output filename")
    parser.add_argument("--mode", choices=['full', 'clean'], default='clean', 
                        help="full: Keep chr/pos columns | clean: Pure numeric matrix for deconvolution")
    
    args = parser.parse_args()
    
    run_pipeline(args.input, args.ref, args.output, args.mode)
