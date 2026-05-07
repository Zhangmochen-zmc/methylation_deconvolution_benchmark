library(matrixStats)
library(quadprog)
library(peakRAM) 

# config
rdir <- "PRMeth-main/PRMeth/R"
ref_file_path <- "prmeth_ref/reference_output_PRMeth.csv"
input_folder <- "test_data" 
output_folder <- "prmeth_result"
benchmark_log_file <- file.path(output_folder, "prmeth_benchmark.csv")

# Create the output directory.
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Import all the R functions required for PRMeth to run.
cat("Loading the PRMeth function...\n")
rfiles <- list.files(rdir, pattern="\\.R$", full.names=TRUE)
# Skip error files
rfiles <- rfiles[basename(rfiles) != "test.R"]
invisible(lapply(rfiles, source))

# Import and process the Reference Matrix (loads only once)
cat("Loading Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# Process the Reference row name (assuming the first column is X or ID).
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  # Otherwise, assume the first column is the ID.
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_prmeth <- as.matrix(df_ref)
cat("Reference Matrix loading complete, dimensions:", dim(ref_prmeth), "\n")

# Batch processing loop
# Get all .csv files
sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(sample_files) == 0) {
  stop("No .csv file found in the input folder.！")
}

# Initialize log data frame
benchmark_results <- data.frame()

cat("Batch processing started, find", length(sample_files), "files...\n\n")

for (file_path in sample_files) {
  
  # Get the filename (without the extension) to name the output file.
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("Sample being processed:", file_name, "\n")
  
  # gc
  gc()
  
  # peakRAM
  monitor_metrics <- peakRAM({
    
    # (A) Read and preprocess sample data
    # Assume the first column of the input CSV is the CpG Probe name.
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # Ensure the matrix is ​​a numeric type (handle potential character type conversion issues).
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    first_data_matrix <- first_data_matrix_num
    
    # (B) Case 1: Partial Reference-Based
    # selecting CpG sites by the coefficient of variation (cv)
    feat <- select_feature(first_data_matrix, 1, 1000)
    
    # Extracting the intersection
    feat_in_ref <- feat[feat %in% rownames(ref_prmeth)]
    
    first_data_matrix_1000 <- first_data_matrix[feat_in_ref, , drop=FALSE]
    ref_prmeth_1000 <- ref_prmeth[feat_in_ref, , drop=FALSE]
    
    # Determining the total number of cell types by λ_BIC
    optimalK <- getCellTypeNumber(first_data_matrix_1000, ref_prmeth_1000, 10)
    
    out_partial <- prmeth(Y = first_data_matrix_1000, W1 = ref_prmeth_1000, 
                          K = optimalK$optimal_K, iters = 1000, rssDiffStop = 1e-10)
    
    # Save Result 1
    out_path_partial <- file.path(output_folder, paste0(file_name, "_PRMeth_partial.txt"))
    write.table(out_partial$H, file = out_path_partial, sep = "\t", col.names = NA, quote = FALSE)
    
    # (C) Scenario 2: Fully Reference-Based
    # Find the intersection (CpG line name)
    common_cpg <- intersect(rownames(first_data_matrix), rownames(ref_prmeth))
    
    first_data_matrix_common  <- first_data_matrix[common_cpg, , drop = FALSE]
    ref_prmeth_common <- ref_prmeth[common_cpg, , drop = FALSE]
    
    out_ref_based <- prmeth(Y = first_data_matrix_common, W1 = ref_prmeth_common, 
                            K = ncol(ref_prmeth_common), iters = 1000, rssDiffStop = 1e-10)
    
    # Save result 2
    out_path_ref <- file.path(output_folder, paste0(file_name, "_PRMeth_ref_based.txt"))
    write.table(out_ref_based$H, file = out_path_ref, sep = "\t", col.names = NA, quote = FALSE)
    
    # (D) Case 3: Reference-Free
    # The optimalK determined in use case 1
    current_K <- optimalK$optimal_K
    if (ncol(first_data_matrix_1000) < current_K) {
      current_K <- 6
    }
    
    out_ref_free <- prmeth(Y = first_data_matrix_1000, W1 = NULL, 
                           K = current_K, iters = 1000, rssDiffStop = 1e-10)
    
    # Save results 3
    out_path_free <- file.path(output_folder, paste0(file_name, "_PRMeth_ref_free.txt"))
    write.table(out_ref_free$H, file = out_path_free, sep = "\t", col.names = NA, quote = FALSE)
    
  }) 
  
  # Extract performance indicators
  elapsed_time <- monitor_metrics$Elapsed_Time_sec
  peak_mem <- monitor_metrics$Peak_RAM_Used_MiB * 1.048576
  
  cat("down: ", file_name, "\n")
  cat("time: ", elapsed_time, " s | peak memory: ", peak_mem, " MB\n")
  
  # Recorded in the log table
  benchmark_results <- rbind(benchmark_results, data.frame(
    Sample_Name = file_name,
    Elapsed_Time_Sec = elapsed_time,
    Peak_RAM_MB = peak_mem,
    Optimal_K = optimalK$optimal_K, 
    Lambda = optimalK$lambda
  ))
  
  # Real-time log saving (to prevent data loss in case of program crash).
  write.csv(benchmark_results, file = benchmark_log_file, row.names = FALSE)
}

cat("\nAll samples have been processed.！\n")
cat("Benchmark logs have been saved to: ", benchmark_log_file, "\n")
