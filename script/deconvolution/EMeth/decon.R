library(tools)
library(peakRAM) 
library(quadprog)
library(Matrix)

# 0. prepare

initial_wd <- getwd()
log_file_path <- "emeth_benchmark.csv" 

# get file
csv_files_list <- list.files(path = "test_data", pattern = "\\.csv$", full.names = TRUE)

# initialize log data frame
benchmark_log <- data.frame()

cat("the following files have been detected and require processing.：\n")
print(csv_files_list)
cat("--------------------------------------------------\n")

# decon
for (current_csv_file in csv_files_list) {
  
  # gc
  gc(verbose = FALSE)
  
  cat(paste0("processing file : ", current_csv_file, "\n"))
  
  # peakRAM
  monitor_res <- peakRAM({

    ##########################################step 1: read the merged CSV and convert its format.
    
    # 1. file path
    input_csv <- current_csv_file
    
    # generate output_file path
    file_base_name <- file_path_sans_ext(basename(input_csv)) 
    output_file <- file.path("txt", paste0(file_base_name, ".txt"))
    
    # check dir
    if(!dir.exists(dirname(output_file))) {
      dir.create(dirname(output_file), recursive = TRUE)
    }
    
    # 2. read CSV file
    first_data <- read.csv(input_csv, header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)
    
    # 3. normalized column names
    colnames(first_data)[1] <- "cg_probe"
    
    # 4. Output converted to TXT format
    write.table(first_data, file = output_file, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    # 5. generating matrix
    rownames(first_data) <- first_data$cg_probe
    first_data <- first_data[, -1] 
    first_data_matrix <- as.matrix(first_data)
    
    ##########################################step 2: fill missing values
    
    first_data_imputed_row_median <- first_data
    first_data_imputed_row_median[] <- t(apply(first_data_imputed_row_median, 1, function(x) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
      return(x)
    }))
    
    first_data_imputed_row_median_matrix <- as.matrix(first_data_imputed_row_median)
    colnames(first_data_imputed_row_median_matrix) <- NULL
    
    ##########################################step 3: set multiple environment variables
    
    # change working directory
    setwd("EMeth")
    
    # reload the library
    library(quadprog)
    library(Matrix)
    
    purity <- rep(0.0000000001, ncol(first_data))
    
    load("emeth_ref/avg_data_matrix.RData")
    mu <- avg_data_matrix 
    
    first_data_nona_matrix_EMeth_reference <- first_data_imputed_row_median_matrix[rownames(first_data_imputed_row_median_matrix) %in% rownames(mu), ]
    
    na_rows <- apply(is.na(first_data_nona_matrix_EMeth_reference), 1, all)
    first_data_nona_matrix_EMeth_reference[na_rows, ] <- 0
    
    check_missing_rows_and_remove <- function(mu, reference_data) {
      mu_rows <- rownames(mu)
      reference_rows <- rownames(reference_data)
      missing_rows <- mu_rows[!mu_rows %in% reference_rows]
      missing_mu_data <- mu[missing_rows, ]
      mu <- mu[!rownames(mu) %in% missing_rows, ]
      return(list(missing_mu_data = missing_mu_data, updated_mu = mu))
    }
    
    result <- check_missing_rows_and_remove(mu, first_data_nona_matrix_EMeth_reference)
    mu <- result$updated_mu
    
    penalty= (dim(mu)[1])*(10^seq(-2,1,0.5)) 
    source('./EMeth-master/R/utils.R')
    source('./EMeth-master/R/cv.emeth.R')
    source('./EMeth-master/R/emeth.R')
    cellTypes = colnames(mu)
    
    ##########################################step 4: run EMeth（LaplaceEM）
    
    print('LaplaceEM')
    hundrediter_laplace = cv.emeth(first_data_nona_matrix_EMeth_reference,
                                   purity,
                                   mu,
                                   aber = FALSE, V='c', init = 'default',
                                   family = 'laplace', nu = penalty, folds = 5, 
                                   maxiter = 50, verbose = TRUE)
    rho.laplace = hundrediter_laplace[[1]]$rho
    colnames(rho.laplace) <- colnames(mu)
    rownames(rho.laplace) <- colnames(first_data)
    
    ##########################################step 5: run EMeth（NormalEM）
    
    print('NormalEM')
    hundrediter = cv.emeth(first_data_nona_matrix_EMeth_reference,
                           purity,
                           mu,
                           aber = FALSE, V='c', init = 'default',
                           family = 'normal', nu = penalty, folds = 5, 
                           maxiter = 50, verbose = TRUE)
    rho.normal = hundrediter[[1]]$rho
    colnames(rho.normal) <- colnames(mu)
    rownames(rho.normal) <- colnames(first_data)
    
    ############################################# step6: save results
    
    # Here we choose to switch back to the initial directory before saving to ensure the path is correct.
    setwd(initial_wd) 
    
    core_name <- file_path_sans_ext(basename(output_file)) 
    save_dir <- "emeth_results" 
    if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    file_path_laplace <- file.path(save_dir, paste0(core_name, "_EMeth_rho_laplace.txt"))
    file_path_normal  <- file.path(save_dir, paste0(core_name, "_EMeth_rho_normal.txt"))
    
    write.table(rho.laplace, file = file_path_laplace, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    write.table(rho.normal, file = file_path_normal, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
  }) 
  
  # benchmark
  elapsed_time <- monitor_res$Elapsed_Time_sec
  peak_mem     <- monitor_res$Peak_RAM_Used_MiB * 1.048576
  
  cat("    -> Running Time :", elapsed_time, "s\n")
  cat("    -> Memory Usage :", peak_mem, "MB (Peak)\n")
  cat("--------------------------------------------------\n")
  
  # log
  benchmark_log <- rbind(benchmark_log, data.frame(
    File_Name = basename(current_csv_file),
    Time_Seconds = elapsed_time,
    Peak_Memory_MB = peak_mem,
    Status = "Success"
  ))

  setwd(initial_wd)
}

# save benchmark
write.csv(benchmark_log, log_file_path, row.names = FALSE)
cat(">>> Done！\n")
cat(">>> benchmark :", log_file_path, "\n")
