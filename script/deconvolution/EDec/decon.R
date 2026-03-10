library(devtools)
devtools::install_github("BRL-BCM/EDec")
library(EDec)
library(gplots)
library(RColorBrewer)
library(peakRAM)
library(clue)

input_folder <- "test_data"  
output_folder <- "edec_result" 
k_num <- 6
# load markers 
markers <- readRDS("edec_stage0_markers.rds") 

# 1. preprocess ref
cat(">>> prepare Tref...\n")
meth_data <- read.table("ref_data/refdata.txt", 
                        row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
meth_data <- as.matrix(meth_data)
meta_data <- read.csv("ref_data/refmeta.csv")

unique_cell_types <- unique(meta_data$Tissue)
Tref <- sapply(unique_cell_types, function(ct) {
  samples_for_ct <- meta_data$FileID[meta_data$Tissue == ct]
  subset_data <- meth_data[, samples_for_ct, drop = FALSE]
  rowMeans(subset_data, na.rm = TRUE)
})
Tref <- as.matrix(Tref)

# 2. prepare cycle
if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

file_list <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
metrics_df <- data.frame()

cat(paste(">>> find", length(file_list), "file...\n"))

# 3. decon
for (i in 1:length(file_list)) {
  current_file_path <- file_list[i]
  file_name_raw <- basename(current_file_path)
  file_basename <- sub("\\.csv$", "", file_name_raw)
  
  cat(paste0("\n[", i, "/", length(file_list), "] processing: ", file_name_raw, "\n"))
  
  # peakRAM
  perf_log <- peakRAM({
    
    # (1) read data
    bulk_data <- read.csv(current_file_path, row.names = 1, check.names = FALSE)
    
    # (2) run EDec Stage 1
    edec_result <- run_edec_stage_1(
      meth_bulk_samples = bulk_data,
      informative_loci = markers,
      num_cell_types = k_num,
      max_its = 2000,
      rss_diff_stop = 1e-10
    )
    
    # (3) cell type
    meth_profiles <- edec_result$methylation
    subset_markers <- intersect(markers, rownames(meth_profiles))
    subset_markers <- intersect(subset_markers, rownames(Tref)) 
    
    # Rows: EDec components, Cols: Reference types
    cor_matrix <- cor(meth_profiles[subset_markers, ], Tref[subset_markers, ])
    
    # cost_matrix
    cost_matrix <- 1 - cor_matrix
    cost_matrix[cost_matrix < 0] <- 0 
    
    # solve_LSAP
    assignment <- solve_LSAP(as.matrix(cost_matrix), maximum = FALSE)
    
    # labels
    assigned_labels <- colnames(cor_matrix)[as.numeric(assignment)]
    
    cat("    Matching results:\n")
    for(j in 1:length(assigned_labels)){
       cat(paste0("    - Component ", j, " -> ", assigned_labels[j], 
                  " (r=", round(cor_matrix[j, assignment[j]], 3), ")\n"))
    }

    # edec results
    colnames(edec_result$proportions) <- assigned_labels
    colnames(edec_result$methylation) <- assigned_labels
    
    # (4) save results
    write.csv(edec_result$proportions, 
              file.path(output_folder, paste0(file_basename, "_Proportions_Labeled.csv")))
    write.csv(edec_result$methylation, 
              file.path(output_folder, paste0(file_basename, "_Methylation_Labeled.csv")))
  })

  gc()

  # benchamrk
  perf_log$File_Name <- file_name_raw
  metrics_df <- rbind(metrics_df, perf_log)
  
  cat(paste(">>> Done！ time:", round(perf_log$Elapsed_Time_sec, 2), "s | memory:", 
            round(perf_log$Total_RAM_Used_MiB* 1.048576, 2), "MB\n"))
}

# 4. save benchmarks
write.csv(metrics_df, file.path(output_folder, "low_benchmark.csv"), row.names = FALSE)
cat("\n>>> Done! result:", output_folder, "\n")
