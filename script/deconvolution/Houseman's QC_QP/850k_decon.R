library(minfi)
library(FlowSorted.Blood.EPIC) 
library(matrixStats)
library(quadprog)
library(peakRAM) 


# 1. config
ref_dir <- "ref"

# 2. 待分析样本的文件夹路径 (存放已生成好的 CSV 文件)
input_folder <- "/data/yuxy/data/data/wgbs_850k/random_1"

# 3. 结果输出文件夹路径
output_folder <- "results"

# ==============================================================================

# 创建输出目录
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# ------------------------------------------------------------------------------
# A. 定义核心函数 (保持原逻辑不变)
# ------------------------------------------------------------------------------

# 1. Signature Probe 选择函数
select_signature <- function(beta_ref, cell_ref, cellTypes,
                             probeSelect = c("any","both"),
                             p_cutoff = 1e-8, n_any = 100, n_each = 50) {
  probeSelect <- match.arg(probeSelect)
  
  N <- length(cell_ref)
  K <- length(cellTypes)
  overall <- rowMeans(beta_ref, na.rm = TRUE)
  
  means <- sapply(cellTypes, function(ct) rowMeans(beta_ref[, cell_ref == ct, drop = FALSE], na.rm = TRUE))
  ns    <- sapply(cellTypes, function(ct) sum(cell_ref == ct))
  
  vars <- sapply(cellTypes, function(ct) rowVars(beta_ref[, cell_ref == ct, drop = FALSE], na.rm = TRUE))
  SSW  <- rowSums(t(t(vars) * (ns - 1)), na.rm = TRUE)
  SSB  <- rowSums(t(t((means - overall)^2) * ns), na.rm = TRUE)
  
  df1 <- K - 1
  df2 <- N - K
  Fstat <- (SSB/df1) / (SSW/df2)
  pval  <- pf(Fstat, df1, df2, lower.tail = FALSE)
  pval[is.na(pval) | !is.finite(pval)] <- 1
  
  rng <- apply(means, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  rng[!is.finite(rng)] <- -Inf
  
  ok <- which(pval < p_cutoff & is.finite(rng) & rng > 0)
  
  if (length(ok) == 0) return(NULL) # 如果没有筛选到，返回 NULL
  
  if (probeSelect == "any") {
    ix <- ok[order(rng[ok], decreasing = TRUE)]
    return(rownames(beta_ref)[ix[seq_len(min(n_any, length(ix)))]])
  }
  
  # both
  maxMean <- apply(means, 1, max, na.rm = TRUE)
  minMean <- apply(means, 1, min, na.rm = TRUE)
  
  ix_hyper <- ok[order(maxMean[ok], decreasing = TRUE)]
  ix_hypo  <- ok[order(minMean[ok], decreasing = FALSE)]
  
  pick_hyper <- ix_hyper[seq_len(min(n_each, length(ix_hyper)))]
  pick_hypo  <- ix_hypo [seq_len(min(n_each, length(ix_hypo )))]
  
  feat <- unique(c(pick_hyper, pick_hypo))
  rownames(beta_ref)[feat]
}

# 2. QP 求解函数
houseman_qp <- function(Y, W, sum_to_one = TRUE) {
  # 确保行名一致
  common <- intersect(rownames(Y), rownames(W))
  Y <- Y[common, , drop=FALSE]
  W <- W[common, , drop=FALSE]
  
  K <- ncol(W)
  nS <- ncol(Y)
  ct <- colnames(W)
  
  P <- matrix(NA_real_, nrow = nS, ncol = K, dimnames = list(colnames(Y), ct))
  
  # 约束矩阵准备
  Amat <- diag(K)
  bvec <- rep(0, K)
  meq  <- 0
  
  if (sum_to_one) {
    Amat <- cbind(rep(1, K), Amat)
    bvec <- c(1, bvec)
    meq  <- 1
  } else {
    Amat <- cbind(-rep(1, K), Amat)
    bvec <- c(-1, bvec)
    meq  <- 0
  }
  
  for (j in seq_len(nS)) {
    y <- Y[, j]
    ok <- is.finite(y)
    
    if(sum(ok) < K) { P[j,] <- NA; next }
    
    Wj <- W[ok, , drop = FALSE]
    y  <- y[ok]
    
    D <- crossprod(Wj) + diag(1e-8, K)
    d <- crossprod(Wj, y)
    
    try({
      sol <- solve.QP(Dmat = D, dvec = as.vector(d), Amat = Amat, bvec = bvec, meq = meq)$solution
      p <- pmax(sol, 0)
      if (sum_to_one && sum(p) > 0) p <- p / sum(p)
      P[j, ] <- p
    }, silent = TRUE)
  }
  
  P
}

# ------------------------------------------------------------------------------
# B. 加载参考数据 (只做一次)
# ------------------------------------------------------------------------------
cat(">>> [Step 1] Loading Reference Data (850k/EPIC)...\n")

# 加载 850k 参考数据
cell_ref_global <- readRDS(file.path(ref_dir, "cell_ref_own_850k.rds"))
beta_ref_global <- readRDS(file.path(ref_dir, "beta_ref_own_850k.rds"))
cellTypes <- unique(cell_ref_global)

cat("    Reference Loaded. Cell Types:", paste(cellTypes, collapse=", "), "\n\n")


# ------------------------------------------------------------------------------
# C. 批量处理流程
# ------------------------------------------------------------------------------

# 获取输入文件夹下的所有 CSV
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("未找到CSV文件，请检查 input_folder 路径！")

benchmark_log <- data.frame()
cat(paste0(">>> [Step 2] Start Batch Processing: ", length(csv_files), " files...\n\n"))

for (i in seq_along(csv_files)) {
  
  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(csv_files), file_name))
  gc(verbose = FALSE) # 内存清理
  
  # === 监控开始 ===
  monitor_res <- peakRAM({
    
    # 1. 读取 CSV 数据
    # ------------------
    sample_df <- read.csv(file_path, header = TRUE, check.names = FALSE)
    # 假设第一列是 CpG Probe ID
    rownames(sample_df) <- sample_df[, 1]
    beta_matrix <- as.matrix(sample_df[, -1])
    
    stopifnot(is.matrix(beta_matrix) || is.data.frame(beta_matrix))
    
    # 2. 取 CpG 交集
    # ------------------
    common_cpg <- intersect(rownames(beta_matrix), rownames(beta_ref_global))
    
    if (length(common_cpg) < 200) {
      # 错误处理：交集过少
      write.table("Error: common_cpg < 200", file.path(output_folder, paste0(file_name, "_ERROR.txt")))
    } else {
      
      # 3. 数据准备 (切片)
      # ------------------
      beta_use <- beta_matrix[common_cpg, , drop = FALSE]
      beta_ref_current <- beta_ref_global[common_cpg, , drop = FALSE]
      
      # 4. 计算 Reference Centroids
      # ------------------
      ref_centroids <- sapply(cellTypes, function(ct) {
        rowMeans(beta_ref_current[, cell_ref_global == ct, drop = FALSE], na.rm = TRUE)
      })
      colnames(ref_centroids) <- cellTypes
      
      # 5. Select Signature
      # ------------------
      probeSelect <- "any"
      feat <- select_signature(beta_ref_current, cell_ref_global, cellTypes, probeSelect = probeSelect)
      
      if (is.null(feat)) {
        write.table("Error: No signature probes found", file.path(output_folder, paste0(file_name, "_NO_SIG.txt")))
      } else {
        
        # 6. 构建矩阵 W 和 Y
        W <- ref_centroids[feat, , drop = FALSE]
        Y <- beta_use[feat, , drop = FALSE]
        
        # 7. 运行 QP 解卷积
        # ------------------
        sum_to_one <- TRUE
        P <- houseman_qp(Y, W, sum_to_one = sum_to_one)
        
        # 8. 保存结果
        # ------------------
        output_file_txt <- file.path(output_folder, paste0(file_name, "_Houseman_QP_result.txt"))
        write.table(P, file = output_file_txt, sep = "\t", col.names = NA, quote = FALSE)
      }
    }
  })
  # === 监控结束 ===
  
  # 记录日志
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB
  cat(sprintf("    -> Time: %.2f s | Peak RAM: %.2f MB\n", elapsed_time, peak_mem))
  
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem
  ))
}

# 保存性能日志
write.csv(benchmark_log, file.path(output_folder, "benchmark.csv"), row.names = FALSE)

cat("==============================================\n")
cat("All tasks finished.\n")
