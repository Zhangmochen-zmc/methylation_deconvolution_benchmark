#################################### 一、读取450k数据的reference data情况

#################################### 步骤1：从每种细胞类型中提取所有用于构建reference matrix的数据，对其进行合并
library(limma)
library(dplyr)
library(tools)
library(tibble)

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

################################################# 步骤2：整合数据，生成merged_data数据框和merged_data_matrix矩阵

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


################################################# 步骤3：整合数据，将每个细胞类型生成一个平均值
# 定义细胞类型关键词
cell_types <- c('bcell', 'cd4', 'cd8', 'neutrophil', 'monocyte', 'nk')

# 创建一个空的数据框来存储每个细胞类型的平均值
averaged_df <- data.frame(cg_id = rownames(merged_data_matrix))

# 遍历每个细胞类型并计算每个细胞类型的平均值
for (cell_type in cell_types) {
  # 提取包含当前细胞类型关键词的列
  selected_columns <- grep(cell_type, colnames(merged_data_matrix), value = TRUE)
  
  # 提取相应的列数据
  selected_data <- merged_data_matrix[, selected_columns]
  
  # 计算该细胞类型每行的平均值
  averaged_data <- rowMeans(selected_data, na.rm = TRUE)
  
  # 将计算结果添加到新的数据框中
  averaged_df[[cell_type]] <- averaged_data
}

# 显示结果
print(averaged_df)

# 如果你需要对列名进行重命名，已经在 averaged_df 中包含了 cg_id 和细胞类型的平均值列
rownames(averaged_df) <- averaged_df$cg_id
averaged_df <- averaged_df[,-1]

# 定义保存路径，将生成的reference matrix进行保存
output_path <- "ref/wgbs_850k_reference_output_PRMeth.csv"

# 将 averaged_df 保存为 CSV 文件
write.csv(averaged_df, file = output_path, row.names = TRUE, col.names = TRUE)

# 提示保存成功
cat("File saved to:", output_path)
