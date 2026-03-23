library(EpiDISH)
library(peakRAM) # 引入用于记录时间和内存的包

# config
input_folder <- "test_data"
# 输出文件夹：存放解卷积结果
output_folder <- "results"
# 性能日志文件保存路径
log_file_path <- "benchmark.csv" 

# 如果输出目录不存在，则创建
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ================= 加载参考矩阵 (只运行一次) =================
cat("正在加载参考矩阵...\n")

# 1. 加载自己的 Reference Matrix (850k)
reference_EpiDISH_850k_path <- "ref/EpiDISH_850k_reference_result.csv"
own_ref_data <- read.csv(reference_EpiDISH_850k_path, row.names = 1)
own_ref_matrix <- as.matrix(own_ref_data)

# 2. 加载 EpiDISH 自带的参考矩阵 (cent12CT.m)
data("cent12CT.m") 

# ================= 初始化日志 =================
benchmark_log <- data.frame()

# ================= 循环处理 CSV 文件 =================

# 获取所有 .csv 文件
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("在输入文件夹中未找到CSV文件，请检查路径: ", input_folder)
}

cat(paste0(">>> 找到 ", length(csv_files), " 个文件。开始批量处理...\n\n"))

for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  
  # 获取文件名（不带后缀），用于命名输出文件
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] 正在处理文件: %s\n", i, length(csv_files), file_name))
  
  # [性能优化] 强制清理内存，确保本次测量的内存不受上一个循环影响
  gc(verbose = FALSE)
  
  # ================= 核心运行与监控 (peakRAM) =================
  # peakRAM 会执行大括号内的代码，并返回时间和内存消耗
  monitor_res <- peakRAM({
    
    # --- 1. 读取数据 ---
    # 假设CSV的第一列是探针ID (cg_probe)，将其作为行名
    beta_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
    beta_matrix <- as.matrix(beta_data)
    
    cat("    数据维度: ", dim(beta_matrix)[1], "行, ", dim(beta_matrix)[2], "列\n")
    
    # -----------------------------------------------------------------
    # 第一部分：使用【自己的】参考矩阵 (Own Reference - 850k)
    # -----------------------------------------------------------------
    
    # 1. RPC 方法
    # cat("    运行 Own Reference (850k) - RPC...\n")
    res_own_rpc <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "RPC")$estF
    write.table(res_own_rpc, 
                file = file.path(output_folder, paste0(file_name, "_result_own_RPC.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 2. CBS 方法
    # cat("    运行 Own Reference (850k) - CBS...\n")
    res_own_cbs <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CBS")$estF
    write.table(res_own_cbs, 
                file = file.path(output_folder, paste0(file_name, "_result_own_CBS.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 3. CP 方法
    # cat("    运行 Own Reference (850k) - CP...\n")
    res_own_cp <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CP")$estF
    write.table(res_own_cp, 
                file = file.path(output_folder, paste0(file_name, "_result_own_CP.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    
    # -----------------------------------------------------------------
    # 第二部分：使用【EpiDISH自带】参考矩阵 (cent12CT.m)
    # -----------------------------------------------------------------
    
    # 1. RPC 方法
    # cat("    运行 EpiDISH Reference (cent12CT) - RPC...\n")
    res_epi_rpc <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "RPC")$estF
    write.table(res_epi_rpc, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_RPC.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 2. CBS 方法
    # cat("    运行 EpiDISH Reference (cent12CT) - CBS...\n")
    res_epi_cbs <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "CBS")$estF
    write.table(res_epi_cbs, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_CBS.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 3. CP 方法
    # cat("    运行 EpiDISH Reference (cent12CT) - CP...\n")
    res_epi_cp <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "CP")$estF
    write.table(res_epi_cp, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_CP.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
  })
  
  # ================= 记录性能指标 =================
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB
  
  cat("    -> 处理耗时 :", elapsed_time, "秒\n")
  cat("    -> 内存峰值 :", peak_mem, "MB\n")
  cat("------------------------------------------------------\n")
  
  # 添加到日志 DataFrame
  # 这里的 Sample_Count 记录的是列数(样本数)
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Sample_Count = ncol(read.csv(file_path, row.names = 1, check.names = FALSE, nrows = 5)),
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem,
    Status = "Success"
  ))
}

# ================= 保存性能日志 =================
write.csv(benchmark_log, log_file_path, row.names = FALSE)

cat("\n所有文件处理完毕！\n")
cat("性能监控日志已保存至:", log_file_path, "\n")
