##### 步骤1：从每种细胞类型中提取所有用于构建reference matrix的数据，对其进行合并

# 加载必要的包
library(limma)
library(dplyr)
library(tools)
library(tibble)
library(EpiDISH)

# 封装数据处理为函数
process_methylation_data <- function(base_path) {
  all_data <- list()  # 用来存放数据
  
  # 获取所有 .txt 文件
  files <- list.files(base_path, pattern = "\\.txt$", full.names = TRUE)  # 只列出 .txt 文件
  
  cat("正在处理文件夹：", base_path, "\n")  # 打印正在处理的文件夹
  
  # 循环读取每个文件并整合数据
  for (file in files) {
    cat("正在读取文件：", file, "\n")  # 打印正在读取的文件名
    
    # 读取数据
    data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # 获取文件名（去掉路径和扩展名）
    file_name <- file_path_sans_ext(basename(file))
    
    # 为列名赋予“文件夹名-文件名”的格式
    # 根据 cell type 动态赋予列名
    col_name <- paste(basename(base_path), file_name, sep = "-")
    
    # 设置列名：cg_probe 和 文件夹-文件名
    colnames(data) <- c("cg_probe", col_name)
    
    # 合并数据：按cg_probe合并
    if (length(all_data) == 0) {
      all_data <- list(data)
      cat("初始化 all_data\n")  # 打印初始化提示
    } else {
      cat("正在合并数据...\n")  # 打印合并提示
      all_data <- mapply(function(x, y) merge(x, y, by = "cg_probe", all = TRUE),
                         all_data, list(data), SIMPLIFY = FALSE)
    }
  }
  
  # 将所有数据合并为一个data.frame
  cat("开始合并数据...\n")  # 打印合并提示
  final_data <- Reduce(function(x, y) merge(x, y, by = "cg_probe", all = TRUE), all_data)
  
  cat("数据合并完成，行数：", nrow(final_data), " 列数：", ncol(final_data), "\n")  # 打印合并后的数据维度
  
  # 输出合并后的结果
  output_path <- paste0(base_path, "_merged_methylation_data.txt")
  write.table(final_data, file = output_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("数据整合完成，已保存为", output_path, "\n")  # 打印文件保存路径
  
  cat("数据处理完成！\n")
  
  return(final_data)  # 返回合并后的数据框
}

# 调用函数进行数据处理
cell_types <- c("bcell", "cd4", "cd8", "nk", "monocyte", "neutrophil")
for (cell_type in cell_types) {
  base_path <- paste0("/data/yuxy/data/data/wgbs_850k/ref_data/", cell_type)
  processed_data <- process_methylation_data(base_path)
  
  # 将生成的所有结果文件移动到指定路径
  output_dir <- "ref"
  file.rename(paste0(base_path, "_merged_methylation_data.txt"), 
              file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")))
  cat("文件已移动至：", file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")), "\n")
}

##### 步骤2：整合数据，生成merged_data数据框和merged_data_matrix矩阵

# 封装数据处理为函数
process_methylation_matrix <- function(base_path, output_file) {
  # 读取各个细胞类型的合并数据
  bcell <- read.table(file.path(base_path, "bcell_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd4 <- read.table(file.path(base_path, "cd4_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd8 <- read.table(file.path(base_path, "cd8_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  neutrophil <- read.table(file.path(base_path, "neutrophil_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  monocyte <- read.table(file.path(base_path, "monocyte_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  nk <- read.table(file.path(base_path, "nk_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 修改列名以反映细胞类型
  colnames(cd4) <- gsub("bcell", "cd4", colnames(cd4))
  colnames(cd8) <- gsub("bcell", "cd8", colnames(cd8))
  colnames(nk) <- gsub("bcell", "nk", colnames(nk))
  colnames(neutrophil) <- gsub("bcell", "neutrophil", colnames(neutrophil))
  colnames(monocyte) <- gsub("bcell", "monocyte", colnames(monocyte))
  
  # 合并数据框
  merged_data <- merge(bcell, cd4, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, cd8, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, neutrophil, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, monocyte, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, nk, by = "cg_probe", all = TRUE)
  
  # 打印合并后的数据维度
  cat("合并后的数据维度：", nrow(merged_data), "行，", ncol(merged_data), "列\n")
  
  # 删除属性值为NA的行
  merged_data <- merged_data %>% filter(!apply(merged_data, 1, function(row) any(row == "NA")))
  
  # 打印删除NA行后的数据维度
  cat("删除包含NA的条目后的数据维度：", nrow(merged_data), "行，", ncol(merged_data), "列\n")
  
  # 输出合并后的数据
  write.table(merged_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("数据已保存为：", output_file, "\n")
  
  # 将 merged_data 转换为 matrix 形式，并设置行名为 cg_probe
  if ("cg_probe" %in% colnames(merged_data)) {
    merged_data_matrix <- as.matrix(merged_data[, -which(colnames(merged_data) == "cg_probe")])  # 删除 cg_probe 列
    rownames(merged_data_matrix) <- merged_data$cg_probe  # 设置行名为 cg_probe
    
    # 返回转换后的矩阵
    return(list(merged_data = merged_data, merged_data_matrix = merged_data_matrix))
  } else {
    cat("警告：数据框中不存在 'cg_probe' 列！\n")
    return(NULL)
  }
}

# 调用函数进行数据处理
base_path <- "ref"
output_file <- "ref/merged_data.txt"
result <- process_methylation_matrix(base_path, output_file)

merged_data <- result$merged_data
merged_data_matrix <- result$merged_data_matrix

# 如果需要查看矩阵的前几行
if (!is.null(merged_data_matrix)) {
  cat("合并后的数据矩阵的前几行：\n")
  print(head(merged_data_matrix))
}


####### 步骤3：针对每个细胞类型生成对应的样本数量，然后每个都跑一遍差异分析，筛选DMC

# 封装数据处理为函数
process_cell_type_data_dynamic <- function(merged_data_matrix, cell_types, output_dir) {
  
  # 存储所有结果
  results <- list()
  
  # 遍历每个细胞类型
  for (cell_type in cell_types) {
    
    # 统计当前细胞类型的样本数
    case_count <- length(grep(cell_type, colnames(merged_data_matrix)))  # case 组数量是该细胞类型的样本数
    control_count <- sum(sapply(cell_types, function(ct) length(grep(ct, colnames(merged_data_matrix))))) - case_count  # control 组是其他细胞类型的样本总数
    
    # 获取包含当前细胞类型的列名
    cell_columns <- grep(cell_type, colnames(merged_data_matrix), value = TRUE)
    
    # 获取不包含当前细胞类型的列名
    other_columns <- setdiff(colnames(merged_data_matrix), cell_columns)
    
    # 将当前细胞类型的列放到前面，其他列放到后面
    cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]
    
    # 创建设计矩阵，表示组别
    group <- factor(c(rep("case", case_count), rep("control", control_count)))
    group <- factor(group, levels = c("case", "control"), ordered = FALSE)
    design <- model.matrix(~ group)
    colnames(design) <- levels(group)
    
    # 使用 lmFit 拟合线性模型
    fit <- lmFit(cell_matrix, design)
    
    # 进行对比，并应用贝叶斯调整
    contrast.matrix <- makeContrasts(case - control, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # 提取所有显著基因
    top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")
    
    # 保存结果到指定路径
    output_path <- file.path(output_dir, paste0(cell_type, "_DMC_results.txt"))
    write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
    cat("数据已保存为：", output_path, "\n")
    
    # 筛选 p 值小于 0.05 的行
    DMC <- top_genes_all[top_genes_all[, 5] < 0.05, ]  # 筛选第六列小于0.05的行
    DMC <- na.omit(DMC)  # 删除任何包含 NA 的行
    
    # 将结果存储到结果列表中
    results[[cell_type]] <- DMC
  }
  
  # 返回所有细胞类型的差异基因结果
  return(results)
}

# 使用示例
output_dir <- "ref"

# 定义细胞类型列表
cell_types <- c("bcell", "cd4", "cd8", "neutrophil", "monocyte", "nk")

# 调用函数进行数据处理
cell_type_results <- process_cell_type_data_dynamic(merged_data_matrix, cell_types, output_dir)

# 查看每个细胞类型的 DMC 结果的维度
lapply(cell_type_results, dim)


####### 步骤4：和位于DHS中的cg探针编号进行合并，筛选每个细胞类型中所有备选的DHS-DMC
# 封装数据处理为函数
filter_and_save_DMC <- function(DHS_cpg_path, DMC_list, output_dir, cell_types) {
  # 遍历 cell_types 中的每个细胞类型
  for (cell_type in cell_types) {
    # 构建DHS文件的路径（根据细胞类型命名）
    # 下面的是所有的all.bed
    DHS_cpg_file <- DHS_cpg_path
    
    # 下面是分开cell_type做DHS的合并
    # DHS_cpg_file <- paste(DHS_cpg_path, "/", cell_type, "-850k.bed", sep = "")
    
    # 读取 DHS_cpg 文件
    if (file.exists(DHS_cpg_file)) {
      DHS_cpg <- read.table(DHS_cpg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      
      # 提取 DHS_cpg 第四列作为筛选条件
      DHS_cpg_ids <- DHS_cpg[, 4]
      
      # 获取当前细胞类型的 DMC 数据框
      DMC_data <- DMC_list[[cell_type]]
      
      # 筛选 DMC 数据框中行名为 DHS_cpg_ids 的条目
      filtered_DMC <- DMC_data[rownames(DMC_data) %in% DHS_cpg_ids, ]
      
      # 将筛选后的结果保存回 DMC_list 中
      DMC_list[[cell_type]] <- filtered_DMC  # 更新 DMC_list
      
      # 设置输出文件路径
      output_path <- file.path(output_dir, paste0(cell_type, "_DHS_DMC.txt"))
      
      # 保存筛选后的结果
      write.table(filtered_DMC, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
      
      # 打印保存结果的提示
      cat(paste(cell_type, "筛选结果已保存。\n"))
    } else {
      cat(paste(cell_type, "未查找到相关DHS文件。\n"))
    }
  }
  return(DMC_list)
}

# 设定 DHS_cpg 文件路径和输出文件夹
# 下面的是DHS合并起来的结果
DHS_cpg_path <- "/data/zhangmch/deconvolution_benchmark/requirement/DHS/epigenomic-roadmap/filtered_cg_info_850k.bed"
# 下面的是DHS单独的结果
# DHS_cpg_path <- "/data/zhangmch/deconvolution_benchmark/requirement/DHS/epigenomic-roadmap/850k_DHS_cell_type"
output_dir <- "ref"

# 调用函数进行筛选并保存结果
cell_type_results_updated <- filter_and_save_DMC(DHS_cpg_path, cell_type_results, output_dir, cell_types)


###### 步骤5：计算不同细胞类型中case组和control组delta beta-value值，并按照DNA甲基化差值进行排序
# 封装数据处理为函数
process_multiple_DHS_DMC_optimized <- function(DHS_DMC_list, merged_data, cell_types) {
  
  # 存储每个细胞类型的处理结果
  results_list <- list()
  
  # 遍历每个细胞类型进行处理
  for (cell_type in cell_types) {
    
    # 提取对应的 DHS_DMC 数据
    DHS_DMC_data <- DHS_DMC_list[[cell_type]]
    
    # 调用处理每个细胞类型的函数
    result <- process_DHS_DMC_optimized(DHS_DMC_data, merged_data, cell_type)
    
    # 存储结果
    results_list[[cell_type]] <- result
  }
  
  # 返回所有细胞类型的结果列表
  return(results_list)
}

# 原来的单细胞类型处理函数
process_DHS_DMC_optimized <- function(DHS_DMC_data, merged_data, cell_type) {
  
  # 将 DHS_DMC_data 的行名转换为一个新的列 cg_probe
  DHS_DMC_data$cg_probe <- rownames(DHS_DMC_data)
  
  # 移除含有 NA 的 cg_probe 行
  DHS_DMC_data <- DHS_DMC_data[!is.na(rownames(DHS_DMC_data)), ]
  merged_data <- merged_data[!is.na(merged_data$cg_probe), ]
  
  # 获取包含 cell_type 的列
  cell_columns <- grep(cell_type, colnames(merged_data), value = TRUE)
  
  # 获取不包含 cell_type 的列，排除 cg_probe 列
  other_columns <- setdiff(colnames(merged_data), c(cell_columns, "cg_probe"))
  
  # 使用 merge 将 DHS_DMC_data 和 merged_data 合并（按 cg_probe）
  merged_data_with_DMC <- merge(DHS_DMC_data, merged_data, by = "cg_probe", all.x = TRUE)
  
  # 检查合并后的数据
  if (any(is.na(merged_data_with_DMC))) {
    warning("合并后数据中存在 NA 值，可能是 cg_probe 匹配失败。")
  }
  
  # 强制转换为数值型数据（以防止字符型数据导致错误）
  merged_data_with_DMC[cell_columns] <- sapply(merged_data_with_DMC[cell_columns], as.numeric)
  merged_data_with_DMC[other_columns] <- sapply(merged_data_with_DMC[other_columns], as.numeric)
  
  # 检查转换后的数据
  if (any(is.na(merged_data_with_DMC[cell_columns]))) {
    warning("cell_columns 中包含 NA 值，检查数据源。")
  }
  if (any(is.na(merged_data_with_DMC[other_columns]))) {
    warning("other_columns 中包含 NA 值，检查数据源。")
  }
  
  # 计算 cell_average 和 other_average
  merged_data_with_DMC$cell_average <- rowMeans(merged_data_with_DMC[, cell_columns], na.rm = TRUE)
  merged_data_with_DMC$other_average <- rowMeans(merged_data_with_DMC[, other_columns], na.rm = TRUE)
  
  # 检查计算的平均值是否为 NA
  if (any(is.na(merged_data_with_DMC$cell_average)) || any(is.na(merged_data_with_DMC$other_average))) {
    warning("计算 cell_average 或 other_average 时出现 NA 值。")
  }
  
  # 计算 cell_diff 和 cell_diff+
  merged_data_with_DMC$cell_diff <- merged_data_with_DMC$cell_average - merged_data_with_DMC$other_average
  merged_data_with_DMC$cell_diff_plus <- abs(merged_data_with_DMC$cell_diff)
  
  # 提取需要的结果，并按 cell_diff+ 排序
  result <- merged_data_with_DMC[, c("cg_probe", "cell_average", "other_average", "cell_diff", "cell_diff_plus")]
  colnames(result) <- c("cg_probe", "cell_average", "other_average", "cell_diff", "cell_diff_plus")
  
  # 排序：按 cell_diff+ 值排序
  result_sorted <- result[order(result$cell_diff_plus, decreasing = TRUE), ]
  
  # 设置行名为 cg_probe
  rownames(result_sorted) <- result_sorted$cg_probe
  
  # 移除 cg_probe 列
  result_sorted$cg_probe <- NULL
  
  # 返回排序后的结果
  return(result_sorted)
}

# 假设 merged_data 已经是处理好的数据框

# 调用处理多个细胞类型的函数
results <- process_multiple_DHS_DMC_optimized(cell_type_results_updated, merged_data, cell_types)

# 查看 bcell 结果
head(results$bcell)


##### 步骤6：针对每一个细胞类型筛选前50个差异性最大的cg探针编号，然后计算和EpiDISH给出的reference中cg探针的重叠程度
##### 步骤6续：根据非重复探针列表生成对应的reference matrix
# 封装数据处理为函数
process_and_calculate_methylation_centroids <- function(merged_data, bcell_results_optimized, nk_results_optimized, 
                                                        cd8_results_optimized, cd4_results_optimized, 
                                                        monocyte_results_optimized, neutrophil_results_optimized) {
  # 提取前50个 DMC 的行名
  bcell_names <- rownames(bcell_results_optimized)[1:50]
  nk_names <- rownames(nk_results_optimized)[1:50]
  cd8_names <- rownames(cd8_results_optimized)[1:50]
  cd4_names <- rownames(cd4_results_optimized)[1:50]
  monocyte_names <- rownames(monocyte_results_optimized)[1:50]
  neutrophil_names <- rownames(neutrophil_results_optimized)[1:50]
  
  # 合并所有行名并提取唯一值
  EpiDISH_DHS_DMC <- unique(c(bcell_names, nk_names, cd8_names, cd4_names, monocyte_names, neutrophil_names))
  
  # 提取 merged_data 中对应 cg_probe 的行
  subset_data <- merged_data[merged_data$cg_probe %in% EpiDISH_DHS_DMC, ]
  
  # 提取包含细胞类型的列
  bcell_columns <- grep("bcell", colnames(subset_data), value = TRUE)
  cd4_columns <- grep("cd4", colnames(subset_data), value = TRUE)
  cd8_columns <- grep("cd8", colnames(subset_data), value = TRUE)
  nk_columns <- grep("nk", colnames(subset_data), value = TRUE)
  monocyte_columns <- grep("monocyte", colnames(subset_data), value = TRUE)
  neutrophil_columns <- grep("neutrophil", colnames(subset_data), value = TRUE)
  
  # 计算每种细胞类型的平均值
  bcell_avg <- rowMeans(subset_data[, bcell_columns], na.rm = TRUE)
  cd4_avg <- rowMeans(subset_data[, cd4_columns], na.rm = TRUE)
  cd8_avg <- rowMeans(subset_data[, cd8_columns], na.rm = TRUE)
  nk_avg <- rowMeans(subset_data[, nk_columns], na.rm = TRUE)
  monocyte_avg <- rowMeans(subset_data[, monocyte_columns], na.rm = TRUE)
  neutrophil_avg <- rowMeans(subset_data[, neutrophil_columns], na.rm = TRUE)
  
  # 创建新的数据框 EpiDISH_450k_reference
  EpiDISH_850k_reference <- data.frame(
    b = bcell_avg,
    cd4 = cd4_avg,
    cd8 = cd8_avg,
    nk = nk_avg,
    monocyte = monocyte_avg,
    neutrophil = neutrophil_avg
  )
  
  # 设置行名为 EpiDISH_DHS_DMC
  rownames(EpiDISH_850k_reference) <- EpiDISH_DHS_DMC
  
  # 返回结果
  return(EpiDISH_850k_reference)
}

# 示例调用
# 假设 bcell_results_optimized, nk_results_optimized 等数据框已经存在
EpiDISH_850k_reference <- process_and_calculate_methylation_centroids(
  merged_data,results$bcell, results$nk, results$cd8, results$cd4,
  results$monocyte, results$neutrophil)

# 查看结果
print(EpiDISH_850k_reference)

# 设置保存路径
output_result_file <- "ref/EpiDISH_850k_reference_result.csv"

# 保存数据框为 CSV 文件
write.csv(EpiDISH_850k_reference, file = output_result_file, row.names = TRUE)

# 打印提示信息
cat("数据已保存为：", output_result_file, "\n")
