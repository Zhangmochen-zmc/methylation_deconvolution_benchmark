library(TOAST)
library(dplyr)
library(stringr)
library(peakRAM) 

# config
ref_file_path <- "tsisal_ref/reference_output_Tsisal.csv"
input_folder <- "test_data" 
output_folder <- "tsisal_result"
benchmark_log_file <- file.path(output_folder, "tsisal_benchmark.csv")

# Create output directory
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Step 1: Import and preprocess the Reference Matrix
cat("loading Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# Process line name
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_tsisal <- as.matrix(df_ref)

# Check if the reference set also contains outliers.
ref_clean_master <- as.matrix(ref_tsisal)
# Remove rows containing NA, NaN, or Inf from the Reference list.
ref_clean_master <- ref_clean_master[rowSums(!is.finite(ref_clean_master)) == 0, ]
cat("Reference Matrix loading and basic cleaning complete. Dimensions:", dim(ref_clean_master), "\n")

# Step 2: Batch processing of samples

sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
if (length(sample_files) == 0) stop("csv file not found")

benchmark_results <- data.frame()

cat("Start batch processing", length(sample_files), "files...\n\n")

for (file_path in sample_files) {
  
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("Sample being processed:", file_name, "\n")
  
  gc() 
  
  # peakRAM
  monitor_metrics <- peakRAM({
    
    # A. Read data
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # B. Core data format conversion
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)

    if (is.vector(first_data_matrix_num)) {
      first_data_matrix_num <- as.matrix(first_data_matrix_num)
      colnames(first_data_matrix_num) <- colnames(first_data_matrix)
    }
    
    # Find the line name
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    first_data_matrix <- first_data_matrix_num
    
    # Delete sites containing NA
    first_data_matrix <- na.omit(first_data_matrix)
    
    # C. Strict data cleaning
    
    # 1. Please ensure again that it is a matrix and a numerical type.
    data_clean <- as.matrix(first_data_matrix)
    class(data_clean) <- "numeric"
    
    # 2. Remove rows containing NA, NaN, or Inf.
    if (nrow(data_clean) > 0) {
      data_clean <- data_clean[rowSums(!is.finite(data_clean)) == 0, , drop=FALSE]
    }
    
    # 3. Reference has been pre-cleaned
    
    # 4. Find the intersection point of the two to ensure a match.
    common_probes <- intersect(rownames(data_clean), rownames(ref_clean_master))
    
    if (length(common_probes) < 100) {
      warning("  [Warning] Insufficient valid loci (<100), skip this sample.")
      out_estProp <- NULL 
    } else {
      
      data_input <- data_clean[common_probes, , drop=FALSE]
      ref_input <- ref_clean_master[common_probes, , drop=FALSE]
      
      # 5. Calculate row variance and filter out sites with extremely low variability.
      if (ncol(data_input) > 1) {
        row_vars <- apply(data_input, 1, var)
        # Only retain sites with variance greater than 1e-8.
        keep_idx <- which(row_vars > 1e-8)
        data_input <- data_input[keep_idx, , drop=FALSE]
        ref_input <- ref_input[rownames(data_input), , drop=FALSE]
      } else {
        cat("  [Info] Row variance cannot be calculated for a single sample; therefore, the variance filtering step is skipped.\n")
      }
      
      # Check if there is still data after filtering.
      if (nrow(data_input) == 0) {
        stop("After cleaning and variance filtering, no valid sites remain.")
      }
      
      # --- D. run Tsisal ---
      cat("  Running Tsisal (K=6)...\n")
      
      # run core function
      out <- Tsisal(data_input, K = 6, knowRef = ref_input)
      out_estProp <- out$estProp
      
      # E. save
      output_txt_path <- file.path(output_folder, paste0(file_name, "_Tsisal_result.txt"))
      write.table(out_estProp, file = output_txt_path, sep = "\t", col.names = NA, quote = FALSE)
    }
    
  }) 
  
  # Extraction indicators
  elapsed_time <- monitor_metrics$Elapsed_Time_sec
  peak_mem <- monitor_metrics$Peak_RAM_Used_MiB * 1.048576
  
  cat("down: ", file_name, "\n")
  cat("time: ", elapsed_time, "s | peak memory: ", peak_mem, "MB\n")
  
  benchmark_results <- rbind(benchmark_results, data.frame(
    Sample = file_name,
    Time_Sec = elapsed_time,
    Peak_Mem_MB = peak_mem,
    Result = ifelse(exists("out_estProp") && !is.null(out_estProp), "Success", "Failed")
  ))
  
  # Real-time log saving
  write.csv(benchmark_results, benchmark_log_file, row.names = FALSE)
}

cat("\nAll tasks completed. Log saved.\n")
