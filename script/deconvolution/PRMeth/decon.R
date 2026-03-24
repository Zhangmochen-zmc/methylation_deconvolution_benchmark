library(matrixStats)
library(quadprog)
library(peakRAM) # 用于监控时间和内存

# ================= 配置区域 (请根据实际情况修改路径) =================

# PRMeth R函数所在的目录
rdir <- "/data/zhangmch/deconvolution_benchmark/script_WGBS/PRMeth/PRMeth-main/PRMeth/R"

# Reference Matrix 文件路径
ref_file_path <- "ref/wgbs_850k_reference_output_PRMeth.csv"

# 待分析样本的文件夹路径 (存放待处理CSV文件的文件夹)
input_folder <- "/data/yuxy/data/data/wgbs_850k/random_10" 

# 结果输出文件夹路径
output_folder <- "results"

# 性能日志保存路径
benchmark_log_file <- file.path(output_folder, "more_benchmark.csv")

# ====================================================================

# 创建输出目录（如果不存在）
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# 2. 导入 PRMeth 运行所需的所有 R 函数
cat("正在加载 PRMeth 函数...\n")
rfiles <- list.files(rdir, pattern="\\.R$", full.names=TRUE)
# 跳过出错文件
rfiles <- rfiles[basename(rfiles) != "test.R"]
invisible(lapply(rfiles, source))

# 3. 导入并处理 Reference Matrix (只需加载一次)
cat("正在加载 Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# 处理 Reference 行名 (假设第一列是X或者ID)
# 如果第一列名为 "X"，则设为行名并删除该列
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  # 否则假设第一列即为ID
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_850k <- as.matrix(df_ref)
cat("Reference Matrix 加载完成，维度:", dim(ref_850k), "\n")


# 4. 批量处理循环

# 获取所有 .csv 文件
sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(sample_files) == 0) {
  stop("在输入文件夹中未找到 .csv 文件！")
}

# 初始化日志数据框
benchmark_results <- data.frame()

cat("开始批量处理，共找到", length(sample_files), "个文件...\n\n")

for (file_path in sample_files) {
  
  # 获取文件名（不带后缀），用于命名输出文件
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("正在处理样本:", file_name, "\n")
  
  # 强制进行垃圾回收，释放内存
  gc()
  
  # 使用 peakRAM 监控代码块的执行
  monitor_metrics <- peakRAM({
    
    # ---------------------------------------------------------
    # (A) 读取并预处理样本数据
    # ---------------------------------------------------------
    # 假设输入CSV的第一列是 CpG Probe 名称
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # 确保矩阵为数值型 (处理可能存在的字符型转换问题)
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    first_data_matrix <- first_data_matrix_num
    
    # ---------------------------------------------------------
    # (B) 情况1：Partial Reference-Based
    # ---------------------------------------------------------
    # selecting CpG sites by the coefficient of variation (cv)
    feat <- select_feature(first_data_matrix, 1, 1000) # 前1000个变异最大的位点
    
    # 提取交集
    feat_in_ref <- feat[feat %in% rownames(ref_850k)]
    
    first_data_matrix_1000 <- first_data_matrix[feat_in_ref, , drop=FALSE]
    ref_850k_1000 <- ref_850k[feat_in_ref, , drop=FALSE]
    
    # Determining the total number of cell types by λ_BIC
    # 注意：getCellTypeNumber 比较耗时，这里一并计入时间
    optimalK <- getCellTypeNumber(first_data_matrix_1000, ref_850k_1000, 10)
    
    out_partial <- prmeth(Y = first_data_matrix_1000, W1 = ref_850k_1000, 
                          K = optimalK$optimal_K, iters = 1000, rssDiffStop = 1e-10)
    
    # 保存结果 1
    out_path_partial <- file.path(output_folder, paste0(file_name, "_PRMeth_partial.txt"))
    write.table(out_partial$H, file = out_path_partial, sep = "\t", col.names = NA, quote = FALSE)
    
    # ---------------------------------------------------------
    # (C) 情况2：Fully Reference-Based
    # ---------------------------------------------------------
    # 取交集（CpG行名）
    common_cpg <- intersect(rownames(first_data_matrix), rownames(ref_850k))
    
    first_data_matrix_common  <- first_data_matrix[common_cpg, , drop = FALSE]
    ref_850k_common <- ref_850k[common_cpg, , drop = FALSE]
    
    out_ref_based <- prmeth(Y = first_data_matrix_common, W1 = ref_850k_common, 
                            K = ncol(ref_850k_common), iters = 1000, rssDiffStop = 1e-10)
    
    # 保存结果 2
    out_path_ref <- file.path(output_folder, paste0(file_name, "_PRMeth_ref_based.txt"))
    write.table(out_ref_based$H, file = out_path_ref, sep = "\t", col.names = NA, quote = FALSE)
    
    # ---------------------------------------------------------
    # (D) 情况3：Reference-Free
    # ---------------------------------------------------------
    # 使用情况1中确定的 optimalK
    current_K <- optimalK$optimal_K
    if (ncol(first_data_matrix_1000) < current_K) {
      current_K <- 6
    }
    
    out_ref_free <- prmeth(Y = first_data_matrix_1000, W1 = NULL, 
                           K = current_K, iters = 1000, rssDiffStop = 1e-10)
    
    # 保存结果 3
    out_path_free <- file.path(output_folder, paste0(file_name, "_PRMeth_ref_free.txt"))
    write.table(out_ref_free$H, file = out_path_free, sep = "\t", col.names = NA, quote = FALSE)
    
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
    Optimal_K = optimalK$optimal_K, # 记录一下算法算出的K值
    Lambda = optimalK$lambda
  ))
  
  # 实时保存日志（防止程序中途崩溃数据丢失）
  write.csv(benchmark_results, file = benchmark_log_file, row.names = FALSE)
}

cat("\n所有样本处理完毕！\n")
cat("Benchmark 日志已保存至: ", benchmark_log_file, "\n")
