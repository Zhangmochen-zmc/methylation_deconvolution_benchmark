library(MASS)
library(peakRAM) 

# 1. load data
if(!file.exists("episcore_ref/episcore.RData"))
load("episcore_ref/episcore.RData")

cat(paste(">>> Data loading complete. Number of projects pending processing:", length(project_data_list), "\n"))

# 2. Create output directory
output_dir <- "episcore_result"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
  cat(paste(">>> Results folder has been created:", output_dir, "\n"))
}

# 3. Initialization performance log table
metrics_df <- data.frame()

# 4. Define the deconvolution function (Matrix Input)
wRPC_Matrix <- function(data_mat, ref_mat, maxit=100){
  # data_mat: Rows=Genes, Cols=Samples
  # ref_mat:  Rows=Genes, Cols=CellTypes + Weight
  
  # Forcefully remove the last column (weight) of Ref.
  actual_ref <- ref_mat[, -ncol(ref_mat), drop=FALSE]
  
  # Initialization result matrix
  est_mat <- matrix(nrow=ncol(data_mat), ncol=ncol(actual_ref))
  
  # Key step: Setting row and column names
  colnames(est_mat) <- colnames(actual_ref) 
  rownames(est_mat) <- colnames(data_mat)   
  
  # Loop through all samples under this Project
  for(s in 1:ncol(data_mat)){
    sample_data_vector <- data_mat[, s]
    rlm.o <- rlm(sample_data_vector ~ actual_ref, maxit=maxit)
    coef.v <- summary(rlm.o)$coef[2:(ncol(actual_ref)+1), 1]
    coef.v[which(coef.v < 0)] <- 0
    total <- sum(coef.v)
    if(total == 0) final_fracs <- rep(0, length(coef.v)) else final_fracs <- coef.v/total
    est_mat[s,] <- final_fracs
  }
  
  return(est_mat)
}

# 5. Process each Project in a loop
cat(">>> Start analysis by Project...\n")

project_names <- names(project_data_list)

for(i in 1:length(project_names)){
  
  p_name <- project_names[i]
  p_data <- project_data_list[[p_name]] 
  
  # A. Start analysis by Project
  common_genes <- intersect(rownames(p_data), rownames(ref_gene_full))
  
  if(length(common_genes) < 100) {
    warning(paste("Project", p_name, "Insufficient effective genes, skip!"))
    next
  }
  
  # Cut data
  curr_proj_data <- p_data[common_genes, , drop=FALSE]
  curr_ref       <- ref_gene_full[common_genes, , drop=FALSE]
  
  # Construct Ref 
  ref_for_run <- cbind(curr_ref, Weight = rep(1, nrow(curr_ref)))
  
  # B. Core deconvolution + peakRAM
  
  final_result_mat <- NULL
  
  perf_log <- peakRAM({
     # wRPC
     final_result_mat <<- wRPC_Matrix(curr_proj_data, ref_for_run, maxit=100)
  })
  
  # --- C. save results
  
  # Filename construction: Result_original_filename.csv
  clean_name <- sub("\\.csv$", "", p_name)
  out_file_name <- paste0("Result_", clean_name, ".csv")
  out_file_path <- file.path(output_dir, out_file_name)
  
  # save CSV
  write.csv(final_result_mat, file = out_file_path, row.names = TRUE)
  
  # D. save benchmark
  
  perf_log$Project_Name <- clean_name
  perf_log$Sample_Count <- ncol(p_data)  
  perf_log$Gene_Count <- length(common_genes)
  perf_log$Function_Call <- NULL 
  
  metrics_df <- rbind(metrics_df, perf_log)
  
  # print
  cat(paste0("[", i, "/", length(project_names), "] ", clean_name, 
             " (Samples: ", ncol(p_data), ")",
             " | Time: ", round(perf_log$Elapsed_Time_sec, 3), "s\n"))
}

# 6. save all benchmarks
write.csv(metrics_df, "low_benchmark.csv", row.names = FALSE)

cat("\n>>> Down！\n")
cat(paste("The results of each project (Samples x CellTypes) are saved in:", output_dir, "\n"))
cat("Performance statistics tables are saved in: episcore_benchmark.csv\n")
