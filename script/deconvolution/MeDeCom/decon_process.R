library(MeDeCom)

# Configure parameters (modify here based on the PDF results from script 1).
best_K <- 6
best_lambda <- 0.001 # <--- 请根据图表调整此值

# 3. 加载脚本1的结果和原始数据
medecom_result <- readRDS("rds/simulated_matrix_450krandom_missing_0.95.rds")
bulk_data <- read.csv("/data/yuxy/data/data/450k/miss/simulated_matrix_450krandom_missing_0.95.csv", row.names = 1, check.names = FALSE)
# 4. 提取 LMCs 和 Proportions
cat("Extracting LMCs and proportions...\n")
lmcs <- getLMCs(medecom_result, K=best_K, lambda=best_lambda)
rownames(lmcs) <- rownames(bulk_data)

props <- getProportions(medecom_result, K=best_K, lambda=best_lambda)

# 5. 加载参考数据进行细胞类型识别
cat("Processing reference data for cell type matching...\n")
meth_data <- read.table("/data/weiyk/R/EDec/450k/ref_merge_data/refdata.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
meth_data <- as.matrix(meth_data)
meta_data <- read.csv("/data/weiyk/R/EDec/450k/ref_merge_data/refmeta.csv")

# 计算参考数据每种细胞类型的平均值 (Tref)
unique_cell_types <- unique(meta_data$Tissue)
Tref <- sapply(unique_cell_types, function(ct) {
  samples_for_ct <- meta_data$FileID[meta_data$Tissue == ct]
  subset_data <- meth_data[, samples_for_ct, drop = FALSE]
  rowMeans(subset_data, na.rm = TRUE)
})
Tref <- as.matrix(Tref)

# 6. 相关性分析与标签匹配 (使用匈牙利算法确保 1-to-1 匹配)
cat("Matching components to cell types using optimal assignment...\n")
subset_markers <- intersect(rownames(Tref), rownames(lmcs))
cor_matrix <- cor(lmcs[subset_markers, ], Tref[subset_markers, ])

# 检查维度
if(nrow(cor_matrix) != ncol(cor_matrix)) {
  cat("Warning: Number of components (K) and reference cell types are not equal!\n")
}

# 匈牙利算法默认是最小化损失，所以我们用 (1 - 相关性) 作为代价矩阵
# 安装并加载 clue 包
if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
library(clue)

# 计算最优分配 (使得总相关性最大)
# solve_LSAP 要求矩阵是非负的代价矩阵，所以用 1 - cor
cost_matrix <- 1 - cor_matrix
# 处理可能存在的负相关（如果有的话，强制转为 0 附近）
cost_matrix[cost_matrix < 0] <- 0 

assignment <- solve_LSAP(as.matrix(cost_matrix), maximum = FALSE)

# assignment 存储了每个 LMC 对应的 Tref 列索引
assigned_labels <- colnames(cor_matrix)[as.numeric(assignment)]

# 打印匹配结果查看
matching_summary <- data.frame(
  Component = paste0("LMC_", 1:best_K),
  Assigned_CellType = assigned_labels,
  Correlation = sapply(1:best_K, function(i) cor_matrix[i, assignment[i]])
)
print(matching_summary)

# 7. 整理并保存结果
# 按照匹配到的标签重命名比例矩阵
props_final <- t(props)
colnames(props_final) <- assigned_labels

# 按照参考细胞类型的标准顺序重新排序列（可选，方便后续对比）
# 比如你想让列的顺序始终是固定的：
target_order <- unique_cell_types 
props_final <- props_final[, target_order, drop = FALSE]

# 保存
write.csv(props_final, "results/simulated_matrix_450krandom_missing_0.95.csv")
cat("Success! 1-to-1 matching complete. Final proportions saved.\n")
