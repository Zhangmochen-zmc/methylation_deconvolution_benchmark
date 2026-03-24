#############################################
# Tsisal 批量解卷积运行脚本 (修正版)
#############################################

library(TOAST) 
# 1. 加载必要的包
#library(Tsisal)  <--- 删除这行，因为报错说没有这个包
library(peakRAM)   # 确保这个包已经安装 (install.packages("peakRAM"))

# ================= 重要：加载 Tsisal 函数 =================
# 请确认 Tsisal.R 文件的真实路径！
# 如果 Tsisal 是一个文件夹里的多个脚本，可能需要 source 那个主文件
#tsisal_script_path <- "/data/yuxy/data/tsisal/Tsisal.R" # <--- 请修改这里！！！

#if (file.exists(tsisal_script_path)) {
#  source(tsisal_script_path)
#  cat("成功加载 Tsisal 脚本文件。\n")
#} else {
#  stop(paste("找不到 Tsisal 脚本文件，请检查路径:", tsisal_script_path))
#}

# ================= 2. 配置路径 =================
ref_path <- "/data/yuxy/data/tsisal/850k/ref/850k_reference_output_Tsisal.csv"
input_folder <- "/data/yuxy/data/data/850k/miss" 
output_folder <- "/data/yuxy/data/tsisal/850k/miss_results"
log_file_path <- file.path(output_folder, "miss_benchmark.csv")

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ================= 3. 加载参考矩阵 =================
cat(">>> 正在加载参考矩阵...\n")
ref_df <- read.csv(ref_path)
rownames(ref_df) <- ref_df[, 1] # 假设第一列是行名
ref_850k <- as.matrix(ref_df[, -1])
cat("    参考矩阵加载完成。\n\n")

# ================= 4. 批量循环处理 =================
benchmark_log <- data.frame()
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("未找到 .csv 文件，请检查路径。")

cat(paste0(">>> 找到 ", length(csv_files), " 个文件。开始批量处理...\n\n"))

for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] 处理文件: %s\n", i, length(csv_files), file_name))
  gc(verbose = FALSE) # 清理内存
  
  # === 监控开始 ===
  monitor_res <- peakRAM({
    
    # 读取数据
    raw_data <- read.csv(file_path, header = TRUE, check.names = FALSE)
    
    # 转换为矩阵 (假设第1列是探针名)
    rownames(raw_data) <- raw_data[, 1]
    data_matrix <- as.matrix(raw_data[, -1])
    
    # 运行 Tsisal
    # 注意：确保 loaded 的 Tsisal 函数名确实叫 "Tsisal"
    out <- Tsisal(data_matrix, K = 6, knowRef = ref_850k)
    
    # 保存结果 (使用 write.table 避免 write.csv 的警告)
    output_file <- file.path(output_folder, paste0(file_name, "_result_Tsisal.txt"))
    write.table(out$estProp, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
    
  })
  # === 监控结束 ===
  
  # 记录日志
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB
  
  cat("    -> 耗时:", elapsed_time, "s | 内存:", peak_mem, "MB\n")
  
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem
  ))
}

# 保存日志
write.csv(benchmark_log, log_file_path, row.names = FALSE)
cat("完成。\n")
