library(EpiDISH)
library(peakRAM)

# config
input_folder <- "test_data"
output_folder <- "epidish_result"
log_file_path <- "epidish_benchmark.csv" 

# If the output directory does not exist, create it.
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Load the reference matrix (run only once)
cat("Loading reference matrix...\n")

# 1. load Reference Matrix (850k)
reference_EpiDISH_850k_path <- "epidish_ref/EpiDISH_850k_reference_result.csv"
own_ref_data <- read.csv(reference_EpiDISH_850k_path, row.names = 1)
own_ref_matrix <- as.matrix(own_ref_data)

# 2. load EpiDISH Reference Matrix (cent12CT.m)
data("cent12CT.m") 

# Initialization log
benchmark_log <- data.frame()

# Loop processing of CSV files

# get csv file
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("The CSV file was not found in the input folder. Please check the path.: ", input_folder)
}

cat(paste0(">>> find ", length(csv_files), " file。Start batch processing...\n\n"))

for (i in seq_along(csv_files)) {
  file_path <- csv_files[i]
  
  # Get the filename (without the extension) to name the output file.
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] processing: %s\n", i, length(csv_files), file_name))
  
  # Force memory cleanup to ensure that the memory measured in this cycle is not affected by the previous cycle.
  gc(verbose = FALSE)
  
  # peakRAM
  monitor_res <- peakRAM({
    
    # read data
    beta_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
    beta_matrix <- as.matrix(beta_data)
    
    cat("Data Dimensions: ", dim(beta_matrix)[1], "row, ", dim(beta_matrix)[2], "col\n")

    # Part 1: Using a custom reference matrix (Own Reference - 850k)

    # 1. RPC
    # cat("    run Own Reference (850k) - RPC...\n")
    res_own_rpc <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "RPC")$estF
    write.table(res_own_rpc, 
                file = file.path(output_folder, paste0(file_name, "_result_own_RPC.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 2. CBS 
    # cat("    run Own Reference (850k) - CBS...\n")
    res_own_cbs <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CBS")$estF
    write.table(res_own_cbs, 
                file = file.path(output_folder, paste0(file_name, "_result_own_CBS.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 3. CP 
    # cat("    run Own Reference (850k) - CP...\n")
    res_own_cp <- epidish(beta.m = beta_matrix, ref.m = own_ref_matrix, method = "CP")$estF
    write.table(res_own_cp, 
                file = file.path(output_folder, paste0(file_name, "_result_own_CP.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # Part Two: Using the built-in EpiDISH reference matrix (cent12CT.m)
    
    # 1. RPC
    # cat("    run EpiDISH Reference (cent12CT) - RPC...\n")
    res_epi_rpc <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "RPC")$estF
    write.table(res_epi_rpc, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_RPC.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 2. CBS 
    # cat("    run EpiDISH Reference (cent12CT) - CBS...\n")
    res_epi_cbs <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "CBS")$estF
    write.table(res_epi_cbs, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_CBS.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # 3. CP 
    # cat("    run EpiDISH Reference (cent12CT) - CP...\n")
    res_epi_cp <- epidish(beta.m = beta_matrix, ref.m = cent12CT.m, method = "CP")$estF
    write.table(res_epi_cp, 
                file = file.path(output_folder, paste0(file_name, "_result_EpiDISH_CP.txt")), 
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
  })
  
  # Record performance metrics
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB
  
  cat("    -> time :", elapsed_time, "s\n")
  cat("    -> peak memory :", peak_mem, "MB\n")
  cat("------------------------------------------------------\n")
  
  # DataFrame
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = file_name,
    Sample_Count = ncol(read.csv(file_path, row.names = 1, check.names = FALSE, nrows = 5)),
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem,
    Status = "Success"
  ))
}

# save benchmark
write.csv(benchmark_log, log_file_path, row.names = FALSE)

cat("\nDown!\n")
cat("save benchmark:", log_file_path, "\n")
