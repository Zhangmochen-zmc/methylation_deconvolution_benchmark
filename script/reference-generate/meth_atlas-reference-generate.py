import pandas as pd
import numpy as np
import argparse
import sys
from scipy.spatial.distance import pdist, squareform

def load_data(file_path):
    """
    读取数据并分离元数据（坐标）和数值数据（甲基化值）。
    适配标准列名: chr, pos
    """
    print(f"[-] 正在读取数据: {file_path} ...")
    try:
        df = pd.read_csv(file_path, index_col=0)
    except Exception as e:
        print(f"[Error] 读取失败: {e}")
        sys.exit(1)
    
    # --- 1. 自动识别坐标列 ---
    # 定义可能的列名映射 (全部转小写对比)
    # 优先级: pos > cpg_beg, chr > cpg_chrm
    coord_map = {
        'chr': 'chr', 'chromosome': 'chr', 'cpg_chrm': 'chr',
        'pos': 'pos', 'position': 'pos', 'start': 'pos', 'cpg_beg': 'pos'
    }
    
    found_cols = {} # 记录找到的列名 { 'chr': '实际列名', 'pos': '实际列名' }
    
    for col in df.columns:
        col_lower = col.lower().strip()
        if col_lower in coord_map:
            std_name = coord_map[col_lower]
            # 如果还没找到，或者找到了更短的(通常更标准)，则更新
            if std_name not in found_cols:
                found_cols[std_name] = col
    
    meta_cols = list(found_cols.values())
    
    # --- 2. 分离数据 ---
    if 'chr' in found_cols and 'pos' in found_cols:
        print(f"[-] 检测到坐标列: {found_cols}")
        meta_df = df[meta_cols].copy()
        # 标准化 meta_df 的列名为 chr, pos，方便后续处理
        meta_df.rename(columns={found_cols['chr']: 'chr', found_cols['pos']: 'pos'}, inplace=True)
        
        # 从主数据中移除坐标列
        data_df = df.drop(columns=meta_cols)
        # 确保剩余全是数值
        data_df = data_df.select_dtypes(include=[np.number])
    else:
        print(f"[Warning] 未能同时检测到 chr 和 pos 列。现有列: {list(df.columns)}")
        print("    -> 将尝试自动分离非数值列作为元数据，但邻居扩展功能可能失效。")
        data_df = df.select_dtypes(include=[np.number])
        meta_df = df.select_dtypes(exclude=[np.number])
        # 如果碰巧有一列叫 pos 但是是数字，它可能还在 data_df 里，这会导致后面出问题
        # 所以这里强制尝试把 pos 提出来
        if 'pos' in df.columns and 'pos' not in meta_df.columns:
            meta_df['pos'] = df['pos']
            data_df = data_df.drop(columns=['pos'], errors='ignore')

    print(f"[-] 数据加载完成: {len(data_df)} 个位点")
    print(f"    - 样本/细胞类型: {len(data_df.columns)} ({list(data_df.columns)})")
    
    return data_df, meta_df

def step1_variance_filtering(df, var_threshold=0.001):
    print("[-] Step 1: 执行缺失值和方差过滤...")
    original_count = len(df)
    df = df.dropna()
    variances = df.var(axis=1)
    df = df[variances >= var_threshold]
    print(f"    - 过滤后剩余: {len(df)} (从 {original_count} 降至 {len(df)})")
    return df

def step2_select_specificity(df, K=100):
    print(f"[-] Step 2: 筛选特异性高/低甲基化位点 (Top {K})...")
    selected_cpgs = set()
    
    # Hyper
    row_sums = df.sum(axis=1) + 1e-9 
    X_norm = df.div(row_sums, axis=0)
    for cell in df.columns:
        selected_cpgs.update(X_norm[cell].nlargest(K).index.tolist())
        
    # Hypo
    X_rev = 1 - df
    row_sums_rev = X_rev.sum(axis=1) + 1e-9
    X_rev_norm = X_rev.div(row_sums_rev, axis=0)
    for cell in df.columns:
        selected_cpgs.update(X_rev_norm[cell].nlargest(K).index.tolist())
        
    print(f"    - 初步筛选位点数量: {len(selected_cpgs)}")
    return selected_cpgs

def step3_add_neighbors(selected_cpgs, all_cpgs_index, meta_df, dist_threshold=50):
    print("[-] Step 3: 检查并添加邻居位点 (50bp)...")
    
    if meta_df is None or 'chr' not in meta_df.columns or 'pos' not in meta_df.columns:
        print("    [Warning] 缺少标准化的 'chr' 或 'pos' 数据，跳过邻居扩展。")
        return selected_cpgs
    
    current_set = set(selected_cpgs)
    
    # 只处理当前过滤后剩下的位点
    # 注意：all_cpgs_index 是经过方差过滤后的 data_df 的索引
    # 我们需要用这个索引去 meta_df 里取坐标
    valid_meta = meta_df.loc[all_cpgs_index].copy()
    
    # 确保 pos 是整数
    try:
        valid_meta['pos'] = valid_meta['pos'].astype(int)
    except:
        print("    [Warning] pos 列包含非数字字符，跳过邻居扩展。")
        return selected_cpgs

    added_count = 0
    # 排序加速查找
    candidates = valid_meta.sort_values(['chr', 'pos'])
    
    # 获取已选位点的坐标
    selected_info = candidates.loc[list(current_set & set(candidates.index))]
    
    for cpg in selected_info.index:
        chrom = selected_info.loc[cpg, 'chr']
        pos = selected_info.loc[cpg, 'pos']
        
        # 查找范围内的邻居
        neighbors = candidates[
            (candidates['chr'] == chrom) & 
            (candidates['pos'] >= pos - dist_threshold) & 
            (candidates['pos'] <= pos + dist_threshold)
        ].index.tolist()
        
        for n in neighbors:
            if n not in current_set:
                current_set.add(n)
                added_count += 1
                
    print(f"    - 添加了 {added_count} 个邻居位点。当前总数: {len(current_set)}")
    return current_set

def step4_pairwise_iterative(df, current_cpgs, max_iterations=50):
    print(f"[-] Step 4: 迭代优化配对差异 (Max Iter: {max_iterations})...")
    
    s_list = list(current_cpgs)
    all_cpgs = df.index.tolist()
    s_set = set(s_list)
    
    for i in range(max_iterations):
        current_atlas = df.loc[s_list]
        dists = pdist(current_atlas.T, metric='euclidean')
        dist_matrix = squareform(dists)
        np.fill_diagonal(dist_matrix, np.inf)
        
        min_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
        cell_a, cell_b = df.columns[min_idx[0]], df.columns[min_idx[1]]
        min_dist_val = dist_matrix[min_idx]
        
        remaining_cpgs = list(set(all_cpgs) - s_set)
        if not remaining_cpgs:
            break
            
        diffs = (df.loc[remaining_cpgs, cell_a] - df.loc[remaining_cpgs, cell_b]).abs()
        best_cpg = diffs.idxmax()
        
        s_list.append(best_cpg)
        s_set.add(best_cpg)
        
        print(f"    Iter {i+1}: ({cell_a} vs {cell_b}), Dist={min_dist_val:.3f}. Add {best_cpg} (Diff={diffs.max():.3f})")
    
    return s_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("-o", "--output", default="final_reference_atlas.csv")
    parser.add_argument("-k", type=int, default=100)
    parser.add_argument("--max_iter", type=int, default=50)
    args = parser.parse_args()
    
    # 1. Load (Fixed logic)
    data_df, meta_df = load_data(args.input_file)
    
    # 2. Filter
    filtered_df = step1_variance_filtering(data_df)
    
    # 3. Specificity
    selected_set = step2_select_specificity(filtered_df, K=args.k)
    
    # 4. Neighbors (Now uses standardized meta_df)
    if meta_df is not None:
        # 对齐索引
        aligned_meta = meta_df.loc[filtered_df.index]
        selected_set = step3_add_neighbors(selected_set, filtered_df.index, aligned_meta)
    
    # 5. Iterative
    final_cpg_list = step4_pairwise_iterative(filtered_df, selected_set, max_iterations=args.max_iter)
    
    # 6. Save
    print(f"[-] 保存结果... 最终位点: {len(final_cpg_list)}")
    final_atlas = data_df.loc[final_cpg_list]
    
    if meta_df is not None:
        final_output = pd.concat([meta_df.loc[final_cpg_list], final_atlas], axis=1)
    else:
        final_output = final_atlas
        
    final_output.to_csv(args.output)
    print(f"[Success] Done: {args.output}")

if __name__ == "__main__":
    main()
