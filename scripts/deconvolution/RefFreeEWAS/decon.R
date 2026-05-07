library(dplyr)
library(RefFreeEWAS)
library(FactoMineR)
library(factoextra)
library(brgedata)
library(minfi)
library(medepir)
library(peakRAM) 

# config
ref_file_path <- "reffreeewas_ref/reference_output_RefFreeEWAS.csv"
input_folder <- "test_data" 
output_folder <- "reffreeewas_result"
benchmark_log_file <- file.path(output_folder, "reffreeewas_benchmark.csv")

# Create output directory
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Step 1: Import and preprocess the Reference Matrix (load only once)
cat("loadins Reference Matrix...\n")
df_ref <- read.csv(ref_file_path)

# Process the Reference row name (assuming the first column is X or ID).
if (colnames(df_ref)[1] == "X") {
  rownames(df_ref) <- df_ref$X
  df_ref <- df_ref[,-1]
} else {
  rownames(df_ref) <- df_ref[,1]
  df_ref <- df_ref[,-1]
}

ref_reffreeewas <- as.matrix(df_ref)
cat("Reference Matrix loading complete, dimensions...:", dim(ref_reffreeewas), "\n")

# Step 2: Read samples in batches and run deconvolution.
sample_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(sample_files) == 0) {
  stop("No .csv file found in the input folder!")
}

# Initialize log data frame
benchmark_results <- data.frame()

cat("Batch processing started, find", length(sample_files), "files...\n\n")

for (file_path in sample_files) {
  
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cat("======================================================\n")
  cat("Sample being processed:", file_name, "\n")
  
  # gc
  gc()
  
  # peakRAM 
  monitor_metrics <- peakRAM({
    
    # read data
    raw_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
    first_data_matrix <- as.matrix(raw_data)
    
    # Convert data 
    # Convert the matrix to a numerical value.
    first_data_matrix_num <- apply(first_data_matrix, 2, as.numeric)
    
    # The `apply` command will lose line names after conversion; you need to retrieve the original line names (cg...).
    rownames(first_data_matrix_num) <- rownames(first_data_matrix)
    
    first_data_matrix <- first_data_matrix_num
    
    # Select the site value that corresponds to the reference
    common_rows <- intersect(rownames(first_data_matrix), rownames(ref_reffreeewas))
    
    # If the intersection is empty, skip the sample.
    if (length(common_rows) == 0) {
      warning(paste("sample", file_name, "Skip any CpG sites that do not overlap with the Reference."))
    } else {
      
      # Extract these lines from first_data_matrix and ref.
      first_data_matrix_updated <- first_data_matrix[common_rows, , drop=FALSE]
      ref_reffreeewas_updated <- ref_reffreeewas[common_rows, , drop=FALSE]
      
      # run RefFreeEWAS 
      cat("  Running RefFreeCellMix...\n")
      cell_mix <- RefFreeEWAS::RefFreeCellMix(
        Y = first_data_matrix_updated, 
        mu0 = ref_reffreeewas_updated, 
        iters = 10, 
        verbose = FALSE 
      )
      
      # save results
      output_txt_path <- file.path(output_folder, paste0(file_name, "_RefFreeEWAS_result.txt"))
      
      # save Omega 
      write.table(cell_mix$Omega, file = output_txt_path, sep = "\t", col.names = NA, quote = FALSE)
    }
    
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
    Common_CpGs = length(common_rows) 
  ))
  
  # Real-time log saving
  write.csv(benchmark_results, file = benchmark_log_file, row.names = FALSE)
}

cat("\nAll samples have been processed!\n")
cat("Benchmark logs have been saved to: ", benchmark_log_file, "\n")
