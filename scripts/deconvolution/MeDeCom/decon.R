suppressPackageStartupMessages({
  library(MeDeCom)
  library(peakRAM)
})

# config
input_dir    <- "test_data"  
output_rds   <- "medecom_rds"                             
output_pdf   <- "medecom_pdf"                            
output_bench <- "medecom_benchmark"                       

# create dir
dir.create(output_rds,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_pdf,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_bench, showWarnings = FALSE, recursive = TRUE)

# Parameter settings
k_values     <- 6
lambda_grid  <- c(0, 10^(-5:-1)) # 0, 0.00001, ..., 0.1

# Get the list of CSV files to be processed
file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# print
cat(sprintf("Checking directories... Done.\n"))
cat(sprintf("Found %d CSV files to process.\n", length(file_list)))
cat("========================================================\n")

# 2. Loop processing
for (i in seq_along(file_list)) {
  
  file_path <- file_list[i]
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  cat(sprintf("[%d/%d] Starting: %s\n", i, length(file_list), file_name))
  
  # check
  target_rds <- file.path(output_rds, paste0(file_name, ".rds"))
  if (file.exists(target_rds)) {
    cat("  Notice: RDS file already exists. Skipping to next...\n")
    next
  }

  # 3. load data
  cat("  -> Loading CSV data...")
  bulk_data <- tryCatch({
    as.matrix(read.csv(file_path, row.names = 1, check.names = FALSE))
  }, error = function(e) {
    cat(sprintf("\n  ERROR: Could not read %s: %s\n", file_name, e$message))
    return(NULL)
  })
  
  if (is.null(bulk_data)) next
  cat(" Done.\n")

  # 4. run MeDeCom
  cat("  -> Running MeDeCom (this may take a long time)... ")
  
  # time and memory
  ram_log <- peakRAM({
    medecom_result <- tryCatch({
      runMeDeCom(
        bulk_data, 
        Ks = k_values, 
        lambdas = lambda_grid, 
        NINIT = 10,       
        NFOLDS = 10,      
        ITERMAX = 300,    
        NCORES = 1     
      )
    }, error = function(e) {
      cat(sprintf("\n  ERROR: runMeDeCom failed on %s: %s\n", file_name, e$message))
      return(NULL)
    })
  })

  # If the calculation fails, proceed to the next step.
  if (is.null(medecom_result)) {
    rm(bulk_data); gc()
    next
  }
  cat(" Done.\n")

  # 5. save 
  cat("  -> Saving RDS and benchmark log...")
  write.csv(ram_log, file.path(output_bench, paste0(file_name, "_bench.csv")), row.names = FALSE)
  saveRDS(medecom_result, target_rds)
  cat(" Done.\n")

  # 6. PDF
  cat("  -> Generating diagnostic plots...")
  pdf_path <- file.path(output_pdf, paste0(file_name, ".pdf"))
  
  tryCatch({
    pdf(pdf_path, width = 10, height = 7)
      # Plot the CVE/LMC error curves for all parameters.
      plotParameters(medecom_result)
      # Lambda refinement curve for K=6
      plotParameters(medecom_result, K=k_values, lambdaScale="log")
    dev.off()
  }, error = function(e) {
    if (!is.null(dev.list())) dev.off()
    cat(sprintf("\n  Warning: Plotting failed for %s", file_name))
  })
  cat(" Done.\n")

  # 7. Memory cleanup
  cat(sprintf("  -> Cleaning memory for %s...\n", file_name))
  rm(bulk_data, medecom_result, ram_log)
  gc() 
  
  cat(sprintf(">>> Finished %s <<<\n\n", file_name))
}

cat("========================================================\n")
cat("All tasks completed!\n")
