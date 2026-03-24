library(TOAST)
library(dplyr)
library(stringr)
library(peakRAM) # 用于监控

# ================= 配置区域 =================

# Reference Matrix 文件路径
ref_file_path <- "ref/wgbs_850k_reference_output_Tsisal.csv"

# 待分析样本文件夹 (存放 CSV 文件)
input_folder <- "/data/yuxy/data/data/wgbs_850k/random_1" 

# 结果输出文件夹
output_folder <- "results"

# 性能日志路径
benchmark_log_file <- file.path(output_folder, "benchmark.csv")

# ============================================

# 创建输出目录
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# --------------------------------------------------------------------
# 步骤1：导入并预处理 Reference Matrix
# --------------------------------------------------------------------
cat("正在加载 Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# 处理行名
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_850k <- as.matrix(df_ref)

# [优化]：Reference 的基础清洗可以先做一次
# 对应原逻辑：检查参考集 (knowRef) 是否也有异常值
ref_clean_master <- as.matrix(ref_850k)
# 剔除 Reference 中含有 NA, NaN, Inf 的行
ref_clean_master <- ref_clean_master[rowSums(!is.finite(ref_clean_master)) == 0, ]
cat("Reference Matrix 加载及基础清洗完成。维度:", dim(ref_clean_master), "\n")


# --------------------------------------------------------------------
# 步骤2：批量处理样本
# --------------------------------------------------------------------

sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
if (length(sample_files) == 0) stop("未找到 .csv 文件")

benchmark_results <- data.frame()

cat("开始批量处理", length(sample_files), "个文件...\n\n")

for (file_path in sample_files) {
  
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("正在处理样本:", file_name, "\n")
  
  gc() # 内存清理
  
  # 监控开始
  monitor_metrics <- peakRAM({
    
    # --- A. 读取数据 ---
    # 假设 CSV 第一列是 Probe ID
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # --- B. 核心数据格式转换 (保留你的原逻辑) ---
    
    # 1. 将矩阵转换为数值型 (apply 按列处理)
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)
    
    # [补丁]：如果样本只有1列，apply可能会返回向量，导致丢失维度，这里强制转回矩阵
    if (is.vector(first_data_matrix_num)) {
      first_data_matrix_num <- as.matrix(first_data_matrix_num)
      colnames(first_data_matrix_num) <- colnames(first_data_matrix)
    }
    
    # 2. 关键一步：找回行名
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    first_data_matrix <- first_data_matrix_num
    
    # 3. 删除含有NA的位点 (na.omit)
    first_data_matrix <- na.omit(first_data_matrix)
    
    # --- C. 严格的数据清洗 (保留 Tsisal 报错处理逻辑) ---
    
    # 1. 再次确保为矩阵和数值型
    data_clean <- as.matrix(first_data_matrix)
    class(data_clean) <- "numeric"
    
    # 2. 剔除含有 NA, NaN, Inf 的行
    # rowSums(!is.finite(...)) == 0 表示这一行所有值都是正常的
    if (nrow(data_clean) > 0) {
      data_clean <- data_clean[rowSums(!is.finite(data_clean)) == 0, , drop=FALSE]
    }
    
    # 3. Reference 已经预清洗过 (ref_clean_master)
    
    # 4. 取两者的交集位点，确保匹配
    common_probes <- intersect(rownames(data_clean), rownames(ref_clean_master))
    
    if (length(common_probes) < 100) {
      warning("  [Warning] 有效位点过少 (<100)，跳过此样本")
      out_estProp <- NULL # 标记为失败
    } else {
      
      data_input <- data_clean[common_probes, , drop=FALSE]
      ref_input <- ref_clean_master[common_probes, , drop=FALSE]
      
      # 5. 计算行方差，过滤掉变异度极小的位点
      # 注意：如果只有一个样本，方差计算结果为NA或0，Tsisal可能无法运行
      if (ncol(data_input) > 1) {
        row_vars <- apply(data_input, 1, var)
        # 只保留方差大于 1e-8 的位点
        keep_idx <- which(row_vars > 1e-8)
        data_input <- data_input[keep_idx, , drop=FALSE]
        ref_input <- ref_input[rownames(data_input), , drop=FALSE]
      } else {
        # 如果是单样本，无法计算方差(var返回NA)，则跳过方差过滤或仅保留非零行
        # 这里为了保持脚本稳健，不做额外方差过滤，直接使用
        cat("  [Info] 单样本无法计算行方差，跳过方差过滤步骤。\n")
      }
      
      # 检查过滤后是否还有数据
      if (nrow(data_input) == 0) {
        stop("经过清洗和方差过滤后，没有剩余的有效位点。")
      }
      
      # --- D. 运行 Tsisal ---
      cat("  Running Tsisal (K=6)...\n")
      
      # 运行核心函数
      out <- Tsisal(data_input, K = 6, knowRef = ref_input)
      out_estProp <- out$estProp
      
      # --- E. 保存结果 ---
      output_txt_path <- file.path(output_folder, paste0(file_name, "_Tsisal_result.txt"))
      write.table(out_estProp, file = output_txt_path, sep = "\t", col.names = NA, quote = FALSE)
    }
    
  }) # peakRAM End
  
  # 提取指标
  elapsed_time <- monitor_metrics$Elapsed_Time_sec
  peak_mem <- monitor_metrics$Peak_RAM_Used_MiB
  
  cat("处理完成: ", file_name, "\n")
  cat("耗时: ", elapsed_time, "s | 内存: ", peak_mem, "MB\n")
  
  benchmark_results <- rbind(benchmark_results, data.frame(
    Sample = file_name,
    Time_Sec = elapsed_time,
    Peak_Mem_MB = peak_mem,
    Result = ifelse(exists("out_estProp") && !is.null(out_estProp), "Success", "Failed")
  ))
  
  # 实时保存日志
  write.csv(benchmark_results, benchmark_log_file, row.names = FALSE)
}

cat("\n所有任务结束。日志已保存。\n")
