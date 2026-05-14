library(EpiDISH)
library(peakRAM)

# ================= CONFIG =================
input_folder <- "test_data"
output_folder <- "epidish_result"
log_file_path <- "epidish_benchmark.csv" 

# create output folder if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# ================= LOAD REFERENCE MATRICES =================
cat("Loading reference matrices...\n")

# 1. Load custom 450k reference matrix
reference_EpiDISH_450k_path <- "marker_ref/450k/EpiDISH_450k_reference_result.csv"
own_ref_data <- read.csv(reference_EpiDISH_450k_path, row.names = 1)
own_ref_matrix <- as.matrix(own_ref_data)

# 2. Load built-in EpiDISH reference matrix for 450k
data("centDHSbloodDMC.m") 

# ================= INITIALIZE BENCHMARK LOG =================
benchmark_log <- data.frame()

# ================= LOOP THROUGH CSV FILES =================
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No CSV files found in input folder. Please check the path: ", input_folder)
}

cat(paste0(">>> Found ", length(csv_files), " files. Starting batch processing...\n\n"))

for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] Processing: %s\n", i, length(csv_files), file_name))
  
  # Clear memory before processing each file
  gc(verbose = FALSE)
  
  # ================= RUN DECONVOLUTION AND MEASURE PERFORMANCE =================
  monitor_res <- peakRAM({
    
    # read CSV
    beta_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
    beta_matrix <- as.matrix(beta_data)
    
    cat("Data dimensions: ", dim(beta_matrix)[1], "row, ", dim(beta_matrix)[2], "col\n")

    # --- Custom 450k Reference ---
    res_own_rpc <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "RPC")$estF
    write.table(res_own_rpc, file.path(output_folder, paste0(file_name, "_result_own_RPC.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    res_own_cbs <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CBS")$estF
    write.table(res_own_cbs, file.path(output_folder, paste0(file_name, "_result_own_CBS.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    res_own_cp <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CP")$estF
    write.table(res_own_cp, file.path(output_folder, paste0(file_name, "_result_own_CP.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    # --- Built-in EpiDISH Reference (centDHSbloodDMC.m) ---
    res_epi_rpc <- epidish(beta.m = beta_matrix, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
    write.table(res_epi_rpc, file.path(output_folder, paste0(file_name, "_result_EpiDISH_RPC.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    res_epi_cbs <- epidish(beta.m = beta_matrix, ref.m = centDHSbloodDMC.m, method = "CBS")$estF
    write.table(res_epi_cbs, file.path(output_folder, paste0(file_name, "_result_EpiDISH_CBS.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    res_epi_cp <- epidish(beta.m = beta_matrix, ref.m = centDHSbloodDMC.m, method = "CP")$estF
    write.table(res_epi_cp, file.path(output_folder, paste0(file_name, "_result_EpiDISH_CP.txt")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
  })
  
  # ================= RECORD PERFORMANCE =================
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB
  
  cat("    -> Elapsed time :", elapsed_time, "s\n")
  cat("    -> Peak memory :", peak_mem, "MB\n")
  cat("------------------------------------------------------\n")
  
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Sample_Count = ncol(read.csv(file_path, row.names = 1, check.names = FALSE, nrows = 5)),
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem,
    Status = "Success"
  ))
}

# ================= SAVE PERFORMANCE LOG =================
write.csv(benchmark_log, log_file_path, row.names = FALSE)

cat("\nAll files processed!\n")
cat("Performance log saved to:", log_file_path, "\n")
