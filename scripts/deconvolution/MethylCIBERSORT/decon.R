library(devtools)
devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
library(peakRAM) 

mix_data_dir  <- "test_data"         
results_dir   <- "methylcibersort_result"           
ref_file_path <- "test_ref_Signature.txt" 
log_file_path <- "methylcibersort_benchmark.csv" 

# 1. check results_dir
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  cat("Created results directory:", results_dir, "\n")
}

# 2. prepare Signature Matrix
cat(">>> Loading Reference Data...\n")
sig_matrix <- read.table(
  ref_file_path,
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)
sig_no_na <- na.omit(sig_matrix)
cat("    Original Ref rows:", nrow(sig_matrix), "\n")
cat("    Cleaned  Ref rows:", nrow(sig_no_na), "\n")
temp_sig_path <- "temp_sig_processed.txt"
write.table(sig_no_na, temp_sig_path, sep = "\t", quote = FALSE, col.names = NA)

# 3. get file list
mix_files <- list.files(path = mix_data_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(mix_files) == 0) {
  stop("Error: No .csv files found in 'mix_data' directory!")
}
cat(paste0(">>> Found ", length(mix_files), " files. Starting batch processing...\n\n"))
benchmark_log <- data.frame()

# 4. process
for (i in seq_along(mix_files)) {
  mix_path <- mix_files[i]
  
  # 4.1 file name and path
  file_name <- basename(mix_path)
  base_name <- tools::file_path_sans_ext(file_name)
  save_path <- file.path(results_dir, paste0(base_name, "_res.txt"))
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(mix_files), file_name))
  
  # gc
  gc(verbose = FALSE)
  
  # 4.3 peakRAM
  monitor_res <- peakRAM({
    
    # A. read Mix data
    exp <- read.csv(mix_path, header = TRUE, row.names = 1, check.names = FALSE)
    
    # B. delete NA
    exp_clean <- na.omit(exp)
    
    if (nrow(exp_clean) == 0) {
      warning(paste("Skipping", file_name, "- All rows removed due to NA."))
      next # skip
    }
    
    # C. CIBERSORT
    res <- cibersort(temp_sig_path, exp_clean, perm = 1000, QN = FALSE)
    
    # D. save results
    write.table(res, save_path, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    
  })
  
  # --- 4.4 benchmark
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB* 1.048576
  
  cat("    -> Save Path    :", save_path, "\n")
  cat("    -> Running Time :", elapsed_time, "s\n")
  cat("    -> Memory Usage :", peak_mem, "MB (Peak)\n")
  cat("--------------------------------------------------\n")
  
  # log
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Save_Path = save_path,
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem,
    Status = "Success"
  ))
}

# 5. save benchmark
write.csv(benchmark_log, log_file_path, row.names = FALSE)

cat("\n>>> Batch processing completed!\n")
cat(">>> Benchmark summary saved to:", log_file_path, "\n")

# delete temp_reference 
if (file.exists(temp_sig_path)) {
  file.remove(temp_sig_path)
}
