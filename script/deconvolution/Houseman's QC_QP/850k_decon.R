library(minfi)
library(FlowSorted.Blood.EPIC) 
library(matrixStats)
library(quadprog)
library(peakRAM) 


# config
ref_dir <- "houseman_ref"
input_folder <- "test_data"
output_folder <- "hosueman_result"

# Create output directory
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# A. Define the core function

# 1. Signature Probe
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
  
  if (length(ok) == 0) return(NULL) 
  
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

# 2. QP 
houseman_qp <- function(Y, W, sum_to_one = TRUE) {
  common <- intersect(rownames(Y), rownames(W))
  Y <- Y[common, , drop=FALSE]
  W <- W[common, , drop=FALSE]
  
  K <- ncol(W)
  nS <- ncol(Y)
  ct <- colnames(W)
  
  P <- matrix(NA_real_, nrow = nS, ncol = K, dimnames = list(colnames(Y), ct))
  
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

# B. Load reference data (do only once)
cat(">>> [Step 1] Loading Reference Data (850k/EPIC)...\n")

# load 850k ref
cell_ref_global <- readRDS(file.path(ref_dir, "cell_ref_own_850k.rds"))
beta_ref_global <- readRDS(file.path(ref_dir, "beta_ref_own_850k.rds"))
cellTypes <- unique(cell_ref_global)

cat("    Reference Loaded. Cell Types:", paste(cellTypes, collapse=", "), "\n\n")

# C. Batch processing workflow

# get CSV
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("CSV file not found. Please check the input_folder path!")

benchmark_log <- data.frame()
cat(paste0(">>> [Step 2] Start Batch Processing: ", length(csv_files), " files...\n\n"))

for (i in seq_along(csv_files)) {
  
  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(csv_files), file_name))
  gc(verbose = FALSE) 
  
  # peakRAM
  monitor_res <- peakRAM({
    
    # 1. read CSV
    # ------------------
    sample_df <- read.csv(file_path, header = TRUE, check.names = FALSE)
    rownames(sample_df) <- sample_df[, 1]
    beta_matrix <- as.matrix(sample_df[, -1])
    
    stopifnot(is.matrix(beta_matrix) || is.data.frame(beta_matrix))
    
    # 2. Find the intersection of CpG
    common_cpg <- intersect(rownames(beta_matrix), rownames(beta_ref_global))
    
    if (length(common_cpg) < 200) {
      # Error handling: Insufficient intersection
      write.table("Error: common_cpg < 200", file.path(output_folder, paste0(file_name, "_ERROR.txt")))
    } else {
      
      # 3. data preparation
      # ------------------
      beta_use <- beta_matrix[common_cpg, , drop = FALSE]
      beta_ref_current <- beta_ref_global[common_cpg, , drop = FALSE]
      
      # 4. calculate Reference Centroids
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
        
        # 6. Construct matrices W and Y
        W <- ref_centroids[feat, , drop = FALSE]
        Y <- beta_use[feat, , drop = FALSE]
        
        # 7. Run QP deconvolution
        # ------------------
        sum_to_one <- TRUE
        P <- houseman_qp(Y, W, sum_to_one = sum_to_one)
        
        # 8. save results
        # ------------------
        output_file_txt <- file.path(output_folder, paste0(file_name, "_Houseman_QP_result.txt"))
        write.table(P, file = output_file_txt, sep = "\t", col.names = NA, quote = FALSE)
      }
    }
  })
  # Monitoring Ended
  
  # log
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB* 1.048576
  cat(sprintf("    -> Time: %.2f s | Peak RAM: %.2f MB\n", elapsed_time, peak_mem))
  
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem
  ))
}

# save benchamrk
write.csv(benchmark_log, file.path(output_folder, "houseman_benchmark.csv"), row.names = FALSE)

cat("==============================================\n")
cat("All tasks finished.\n")
