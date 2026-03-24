library(dplyr)
library(RefFreeEWAS)
library(FactoMineR)
library(factoextra)
library(brgedata)
library(minfi)
library(medepir)
library(peakRAM) # 新增：用于监控时间和内存

# ================= 配置区域 (请根据实际情况修改路径) =================

# Reference Matrix 文件路径
ref_file_path <- "ref/wgbs_850k_reference_output_RefFreeEWAS.csv"

# 待分析样本的文件夹路径 (存放待处理CSV文件的文件夹)
# 假设这里存放的是已经整理好的单样本CSV文件
input_folder <- "/data/yuxy/data/data/wgbs_850k/random_1" 

# 结果输出文件夹路径
output_folder <- "results"

# 性能日志保存路径
benchmark_log_file <- file.path(output_folder, "benchmark.csv")

# ====================================================================

# 创建输出目录
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# --------------------------------------------------------------------
# 步骤1：导入并预处理 Reference Matrix (只需加载一次)
# --------------------------------------------------------------------
cat("正在加载 Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# 处理 Reference 行名 (假设第一列是X或者ID)
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_850k <- as.matrix(df_ref)
cat("Reference Matrix 加载完成，维度:", dim(ref_850k), "\n")


# --------------------------------------------------------------------
# 步骤2：批量读取样本并运行解卷积
# --------------------------------------------------------------------

# 获取所有 .csv 文件
sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(sample_files) == 0) {
  stop("在输入文件夹中未找到 .csv 文件！")
}

# 初始化日志数据框
benchmark_results <- data.frame()

cat("开始批量处理，共找到", length(sample_files), "个文件...\n\n")

for (file_path in sample_files) {
  
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("正在处理样本:", file_name, "\n")
  
  # 强制垃圾回收
  gc()
  
  # 使用 peakRAM 监控
  monitor_metrics <- peakRAM({
    
    # --- A. 读取数据 ---
    # 假设 CSV 第一列是 Probe ID
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # --- B. 数据格式转换 (保留你要求的核心逻辑) ---
    # 1. 将矩阵转换为数值型
    # apply(..., 2, ...) 表示按列处理，as.numeric 将字符转为数字
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)
    
    # 2. 关键一步：apply 转换后会丢失行名，需要把原来的行名（cg...）找回来
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    
    first_data_matrix <- first_data_matrix_num
    
    # --- C. 选择和reference能对应上的位点值 ---
    # 获取行名交集
    common_rows <- intersect(rownames(first_data_matrix), rownames(ref_850k))
    
    # 如果交集为空，跳过该样本
    if (length(common_rows) == 0) {
      warning(paste("样本", file_name, "与 Reference 没有重合的 CpG 位点，跳过。"))
    } else {
      
      # 从 first_data_matrix 和 ref 中提取这些行
      first_data_matrix_updated <- first_data_matrix[common_rows, , drop=FALSE]
      ref_850k_updated <- ref_850k[common_rows, , drop=FALSE]
      
      # --- D. 运行 RefFreeEWAS ---
      # 这里的逻辑是 Reference-Based (传入了 mu0)
      cat("  Running RefFreeCellMix...\n")
      cell_mix <- RefFreeEWAS::RefFreeCellMix(
        Y = first_data_matrix_updated, 
        mu0 = ref_850k_updated, 
        iters = 10, 
        verbose = FALSE # 批量运行时建议关闭 verbose 以减少刷屏，或者设为 FALSE
      )
      
      # --- E. 保存结果 ---
      output_txt_path <- file.path(output_folder, paste0(file_name, "_RefFreeEWAS_result.txt"))
      
      # 保存 Omega (细胞混合比例)
      write.table(cell_mix$Omega, file = output_txt_path, sep = "\t", col.names = NA, quote = FALSE)
    }
    
  }) # peakRAM 结束
  
  # 提取性能指标
  elapsed_time <- monitor_metrics$Elapsed_Time_sec
  peak_mem <- monitor_metrics$Peak_RAM_Used_MiB
  
  cat("处理完成: ", file_name, "\n")
  cat("耗时: ", elapsed_time, " 秒 | 内存峰值: ", peak_mem, " MiB\n")
  
  # 记录到日志表
  benchmark_results <- rbind(benchmark_results, data.frame(
    Sample_Name = file_name,
    Elapsed_Time_Sec = elapsed_time,
    Peak_RAM_MiB = peak_mem,
    Common_CpGs = length(common_rows) # 记录一下使用了多少个位点
  ))
  
  # 实时保存日志
  write.csv(benchmark_results, file = benchmark_log_file, row.names = FALSE)
}

cat("\n所有样本处理完毕！\n")
cat("Benchmark 日志已保存至: ", benchmark_log_file, "\n")
