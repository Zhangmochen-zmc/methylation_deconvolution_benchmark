##### Step 1: Extract all data used to construct the reference matrix from each cell type and merge them.
library(limma)
library(dplyr)
library(tools)
library(tibble)
library(EpiDISH)

# Encapsulate data processing as functions
process_methylation_data <- function(base_path) {
  all_data <- list() 
  
  # get .txt file
  files <- list.files(base_path, pattern = "\\.txt$", full.names = TRUE) 
  
  cat("Folder being processed：", base_path, "\n") 
  
  # Loop through each file and integrate the data.
  for (file in files) {
    cat("Reading file：", file, "\n") 
    
    # read data
    data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # get file name
    file_name <- file_path_sans_ext(basename(file))
    
    # Assign column names the format "folder name-file name".
    # Dynamically assign column names based on cell type
    col_name <- paste(basename(base_path), file_name, sep = "-")
    
    # Set column names: cg_probe and folder-filename.
    colnames(data) <- c("cg_probe", col_name)
    
    # Merge data: Merge by cg_probe
    if (length(all_data) == 0) {
      all_data <- list(data)
      cat("initialization all_data\n")  
    } else {
      cat("Data being merged...\n")  
      all_data <- mapply(function(x, y) merge(x, y, by = "cg_probe", all = TRUE),
                         all_data, list(data), SIMPLIFY = FALSE)
    }
  }
  
  # Merge all data into a single data.frame
  cat("Start merging data...\n")  
  final_data <- Reduce(function(x, y) merge(x, y, by = "cg_probe", all = TRUE), all_data)
  
  cat("Data merging complete，rows ：", nrow(final_data), " cols：", ncol(final_data), "\n")  
  
  # output results
  output_path <- paste0(base_path, "_merged_methylation_data.txt")
  write.table(final_data, file = output_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("Data integration is complete and has been saved as ", output_path, "\n")  
  
  cat("Data processing completed！\n")
  
  return(final_data) 
}

# Calling functions to process data
cell_types <- c("bcell", "cd4", "cd8", "nk", "monocyte", "neutrophil")
for (cell_type in cell_types) {
  base_path <- paste0("ref_data/", cell_type)
  processed_data <- process_methylation_data(base_path)
  
  # Move all generated result files to the specified path.
  output_dir <- "epidish_ref"
  file.rename(paste0(base_path, "_merged_methylation_data.txt"), 
              file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")))
  cat("The file has been moved to：", file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")), "\n")
}

##### Step 2: Integrate the data to generate a merged_data data frame and a merged_data_matrix matrix.

# Encapsulate data processing as functions
process_methylation_matrix <- function(base_path, output_file) {
  # Read merged data from various cell types
  bcell <- read.table(file.path(base_path, "bcell_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd4 <- read.table(file.path(base_path, "cd4_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd8 <- read.table(file.path(base_path, "cd8_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  neutrophil <- read.table(file.path(base_path, "neutrophil_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  monocyte <- read.table(file.path(base_path, "monocyte_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  nk <- read.table(file.path(base_path, "nk_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Modify column names to reflect cell type
  colnames(cd4) <- gsub("bcell", "cd4", colnames(cd4))
  colnames(cd8) <- gsub("bcell", "cd8", colnames(cd8))
  colnames(nk) <- gsub("bcell", "nk", colnames(nk))
  colnames(neutrophil) <- gsub("bcell", "neutrophil", colnames(neutrophil))
  colnames(monocyte) <- gsub("bcell", "monocyte", colnames(monocyte))
  
  # Merge Data Frames
  merged_data <- merge(bcell, cd4, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, cd8, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, neutrophil, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, monocyte, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, nk, by = "cg_probe", all = TRUE)
  
  # Print merged data dimensions
  cat("Merged data dimensions：", nrow(merged_data), "rows，", ncol(merged_data), "cols\n")
  
  # delete NA
  merged_data <- merged_data %>% filter(!apply(merged_data, 1, function(row) any(row == "NA")))
  
  # Print the data dimensions after deleting NA rows.
  cat("Data dimensions after deleting entries containing NA：", nrow(merged_data), "rows，", ncol(merged_data), "cols\n")
  
  # Output merged data
  write.table(merged_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("The data has been saved as：", output_file, "\n")
  
  # Convert merged_data to matrix form and set the row name to cg_probe.
  if ("cg_probe" %in% colnames(merged_data)) {
    merged_data_matrix <- as.matrix(merged_data[, -which(colnames(merged_data) == "cg_probe")]) 
    rownames(merged_data_matrix) <- merged_data$cg_probe 
    
    # Return the transformed matrix
    return(list(merged_data = merged_data, merged_data_matrix = merged_data_matrix))
  } else {
    cat("Warning: The column 'cg_probe' does not exist in the data frame!\n")
    return(NULL)
  }
}

# Calling functions to process data
base_path <- "epidish_ref"
output_file <- "epidish_ref/merged_data.txt"
result <- process_methylation_matrix(base_path, output_file)

merged_data <- result$merged_data
merged_data_matrix <- result$merged_data_matrix

# If you need to view the first few rows of the matrix
if (!is.null(merged_data_matrix)) {
  cat("The first few rows of the merged data matrix：\n")
  print(head(merged_data_matrix))
}


####### Step 3: Generate the corresponding number of samples for each cell type, and then run differential analysis on each one to screen DMCs.

# Encapsulate data processing as functions
process_cell_type_data_dynamic <- function(merged_data_matrix, cell_types, output_dir) {
  
  # Store all results
  results <- list()
  
  # Traverse each cell type
  for (cell_type in cell_types) {
    
    # Count the number of samples for the current cell type
    case_count <- length(grep(cell_type, colnames(merged_data_matrix)))  # case : current cell type
    control_count <- sum(sapply(cell_types, function(ct) length(grep(ct, colnames(merged_data_matrix))))) - case_count  # control : other cell type
    
    # Retrieve the column names containing the current cell type
    cell_columns <- grep(cell_type, colnames(merged_data_matrix), value = TRUE)
    
    # Get the column names not containing the current cell type
    other_columns <- setdiff(colnames(merged_data_matrix), cell_columns)
    
    # Place the column containing the current cell type first, and place other columns last
    cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]
    
    # Create a design matrix to represent the groups
    group <- factor(c(rep("case", case_count), rep("control", control_count)))
    group <- factor(group, levels = c("case", "control"), ordered = FALSE)
    design <- model.matrix(~ group)
    colnames(design) <- levels(group)
    
    # lmFit
    fit <- lmFit(cell_matrix, design)
    
    # Perform comparisons and apply Bayesian adjustment
    contrast.matrix <- makeContrasts(case - control, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # Extract all significant genes
    top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")
    
    # save results
    output_path <- file.path(output_dir, paste0(cell_type, "_DMC_results.txt"))
    write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
    cat("The data has been saved as...：", output_path, "\n")
    
    # Filter rows with p-values ​​less than 0.05
    DMC <- top_genes_all[top_genes_all[, 5] < 0.05, ] 
    DMC <- na.omit(DMC)  
    
    # save results
    results[[cell_type]] <- DMC
  }
  
  # Return the differential gene results for all cell types.
  return(results)
}

# usage
output_dir <- "epidish_ref"

# define cell type
cell_types <- c("bcell", "cd4", "cd8", "neutrophil", "monocyte", "nk")

# Calling functions for data processing
cell_type_results <- process_cell_type_data_dynamic(merged_data_matrix, cell_types, output_dir)

# View the dimensions of the DMC results for each cell type
lapply(cell_type_results, dim)


####### Step 4: Merge with the cg probe numbers located in DHS, and filter all candidate DHS-DMCs for each cell type
# Encapsulate data processing into functions
filter_and_save_DMC <- function(DHS_cpg_path, DMC_list, output_dir, cell_types) {
  # Loop through each cell type in cell_types
  for (cell_type in cell_types) {
    # Construct the path to the DHS file (named according to cell type)
    DHS_cpg_file <- DHS_cpg_path
    
    # Merging of DHSs separately for each cell_type
    # DHS_cpg_file <- paste(DHS_cpg_path, "/", cell_type, "-450k.bed", sep = "")
    
    #read DHS_cpg file
    if (file.exists(DHS_cpg_file)) {
      DHS_cpg <- read.table(DHS_cpg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      
      # Extract the fourth column of DHS_cpg as the filter criteria.
      DHS_cpg_ids <- DHS_cpg[, 4]
      
      # Get the DMC data frame for the current cell type
      DMC_data <- DMC_list[[cell_type]]
      
      # Filter entries in the DMC data frame whose row name is DHS_cpg_ids
      filtered_DMC <- DMC_data[rownames(DMC_data) %in% DHS_cpg_ids, ]
      
      # Save the filtered results back to DMC_list
      DMC_list[[cell_type]] <- filtered_DMC 
      
      # Set the output file path
      output_path <- file.path(output_dir, paste0(cell_type, "_DHS_DMC.txt"))
      
      # save results
      write.table(filtered_DMC, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
      
      # print
      cat(paste(cell_type, "The screening results have been saved.\n"))
    } else {
      cat(paste(cell_type, "No relevant DHS file was found.\n"))
    }
  }
  return(DMC_list)
}

# Set the DHS_cpg file path and output folder
# DHS merger
DHS_cpg_path <- "DHS/epigenomic-roadmap/filtered_cg_info_450k.bed"
# DHS alone
# DHS_cpg_path <- "DHS/epigenomic-roadmap/450k_DHS_cell_type"
output_dir <- "epidish_ref"

# Call the function to filter and save the results
cell_type_results_updated <- filter_and_save_DMC(DHS_cpg_path, cell_type_results, output_dir, cell_types)


###### Step 5: Calculate the delta-beta-values ​​of the case and control groups in different cell types, and sort them according to the DNA methylation difference.
# Encapsulate data processing as functions
process_multiple_DHS_DMC_optimized <- function(DHS_DMC_list, merged_data, cell_types) {
  
  # Store the processing results for each cell type
  results_list <- list()
  
  # Process each cell type
  for (cell_type in cell_types) {
    
    # Extract the corresponding DHS_DMC data
    DHS_DMC_data <- DHS_DMC_list[[cell_type]]
    
    # Calling the function that processes each cell type
    result <- process_DHS_DMC_optimized(DHS_DMC_data, merged_data, cell_type)
    
    # Save results
    results_list[[cell_type]] <- result
  }
  
  # Returns a list of results for all cell types.
  return(results_list)
}

# Original single-cell type processing function
process_DHS_DMC_optimized <- function(DHS_DMC_data, merged_data, cell_type) {
  
  # Convert the row names of DHS_DMC_data to a new column cg_probe
  DHS_DMC_data$cg_probe <- rownames(DHS_DMC_data)
  
  # Remove cg_probe lines containing NA
  DHS_DMC_data <- DHS_DMC_data[!is.na(rownames(DHS_DMC_data)), ]
  merged_data <- merged_data[!is.na(merged_data$cg_probe), ]
  
  # Get the column containing cell_type
  cell_columns <- grep(cell_type, colnames(merged_data), value = TRUE)
  
  # Retrieve columns that do not contain cell_type, excluding the cg_probe column.
  other_columns <- setdiff(colnames(merged_data), c(cell_columns, "cg_probe"))
  
  # Use merge to combine DHS_DMC_data and merged_data (by cg_probe)
  merged_data_with_DMC <- merge(DHS_DMC_data, merged_data, by = "cg_probe", all.x = TRUE)
  
  # Check the merged data
  if (any(is.na(merged_data_with_DMC))) {
    warning("The presence of NA values ​​in the merged data may indicate that cg_probe failed to match.")
  }
  
  # Force conversion to numeric data (to prevent errors caused by character data)
  merged_data_with_DMC[cell_columns] <- sapply(merged_data_with_DMC[cell_columns], as.numeric)
  merged_data_with_DMC[other_columns] <- sapply(merged_data_with_DMC[other_columns], as.numeric)
  
  # Check the converted data
  if (any(is.na(merged_data_with_DMC[cell_columns]))) {
    warning("The cell_columns contains NA values; check the data source.")
  }
  if (any(is.na(merged_data_with_DMC[other_columns]))) {
    warning("The other_columns contains NA values; check the data source.")
  }
  
  # Calculate cell_average and other_average
  merged_data_with_DMC$cell_average <- rowMeans(merged_data_with_DMC[, cell_columns], na.rm = TRUE)
  merged_data_with_DMC$other_average <- rowMeans(merged_data_with_DMC[, other_columns], na.rm = TRUE)
  
  # Check if the calculated average is NA.
  if (any(is.na(merged_data_with_DMC$cell_average)) || any(is.na(merged_data_with_DMC$other_average))) {
    warning("NA values ​​appear when calculating cell_average or other_average.")
  }
  
  # Calculate cell_diff and cell_diff+
  merged_data_with_DMC$cell_diff <- merged_data_with_DMC$cell_average - merged_data_with_DMC$other_average
  merged_data_with_DMC$cell_diff_plus <- abs(merged_data_with_DMC$cell_diff)
  
  # Extract the desired results and sort them by cell_diff+
  result <- merged_data_with_DMC[, c("cg_probe", "cell_average", "other_average", "cell_diff", "cell_diff_plus")]
  colnames(result) <- c("cg_probe", "cell_average", "other_average", "cell_diff", "cell_diff_plus")
  
  # Sort by cell_diff+ value
  result_sorted <- result[order(result$cell_diff_plus, decreasing = TRUE), ]
  
  # Set the line name to cg_probe
  rownames(result_sorted) <- result_sorted$cg_probe
  
  # Remove the cg_probe column
  result_sorted$cg_probe <- NULL
  
  # Return the sorted results
  return(result_sorted)
}

# Calling functions that handle multiple cell types
results <- process_multiple_DHS_DMC_optimized(cell_type_results_updated, merged_data, cell_types)

# View bcell results
head(results$bcell)


##### Step 6: For each cell type, select the top 50 most differentially expressed CG probes and then calculate their overlap with the CG probes in the EpiDISH reference.
# Encapsulate data processing as functions
process_and_calculate_methylation_centroids <- function(merged_data, bcell_results_optimized, nk_results_optimized, 
                                                        cd8_results_optimized, cd4_results_optimized, 
                                                        monocyte_results_optimized, neutrophil_results_optimized) {
  # Extract the row names of the first 50 DMCs
  bcell_names <- rownames(bcell_results_optimized)[1:50]
  nk_names <- rownames(nk_results_optimized)[1:50]
  cd8_names <- rownames(cd8_results_optimized)[1:50]
  cd4_names <- rownames(cd4_results_optimized)[1:50]
  monocyte_names <- rownames(monocyte_results_optimized)[1:50]
  neutrophil_names <- rownames(neutrophil_results_optimized)[1:50]
  
  # Merge all row names and extract unique values
  EpiDISH_DHS_DMC <- unique(c(bcell_names, nk_names, cd8_names, cd4_names, monocyte_names, neutrophil_names))
  
  # Extract the row corresponding to cg_probe from merged_data.
  subset_data <- merged_data[merged_data$cg_probe %in% EpiDISH_DHS_DMC, ]
  
  # Extract columns containing cell types
  bcell_columns <- grep("bcell", colnames(subset_data), value = TRUE)
  cd4_columns <- grep("cd4", colnames(subset_data), value = TRUE)
  cd8_columns <- grep("cd8", colnames(subset_data), value = TRUE)
  nk_columns <- grep("nk", colnames(subset_data), value = TRUE)
  monocyte_columns <- grep("monocyte", colnames(subset_data), value = TRUE)
  neutrophil_columns <- grep("neutrophil", colnames(subset_data), value = TRUE)
  
  # Calculate the average value for each cell type
  bcell_avg <- rowMeans(subset_data[, bcell_columns], na.rm = TRUE)
  cd4_avg <- rowMeans(subset_data[, cd4_columns], na.rm = TRUE)
  cd8_avg <- rowMeans(subset_data[, cd8_columns], na.rm = TRUE)
  nk_avg <- rowMeans(subset_data[, nk_columns], na.rm = TRUE)
  monocyte_avg <- rowMeans(subset_data[, monocyte_columns], na.rm = TRUE)
  neutrophil_avg <- rowMeans(subset_data[, neutrophil_columns], na.rm = TRUE)
  
  # Create a new data frame EpiDISH_450k_reference
  EpiDISH_450k_reference <- data.frame(
    b = bcell_avg,
    cd4 = cd4_avg,
    cd8 = cd8_avg,
    nk = nk_avg,
    monocyte = monocyte_avg,
    neutrophil = neutrophil_avg
  )
  
  # Set the line name to EpiDISH_DHS_DMC
  rownames(EpiDISH_450k_reference) <- EpiDISH_DHS_DMC
  
  # return results
  return(EpiDISH_450k_reference)
}

# usage
EpiDISH_450k_reference <- process_and_calculate_methylation_centroids(
  merged_data,results$bcell, results$nk, results$cd8, results$cd4,
  results$monocyte, results$neutrophil)

# pring
print(EpiDISH_450k_reference)

# Set save path
output_result_file <- "epidish_ref/EpiDISH_450k_reference_result.csv"

# save
write.csv(EpiDISH_450k_reference, file = output_result_file, row.names = TRUE)

# print
cat("The data has been saved as：", output_result_file, "\n")

