library(limma)
library(dplyr)
library(tools)
library(tibble)

# Encapsulate data processing as functions
process_methylation_data <- function(base_path) {
  all_data <- list() 
  
  # get .txt file
  files <- list.files(base_path, pattern = "\\.txt$", full.names = TRUE) 
  
  cat("processing：", base_path, "\n")  
  
  # Loop through each file and integrate the data
  for (file in files) {
    cat("processing：", file, "\n")  
    
    # read data
    data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # get filename
    file_name <- file_path_sans_ext(basename(file))
    
    # Assign column names the format "folder name-file name".
    col_name <- paste(basename(base_path), file_name, sep = "-")
    
    # Set column names: cg_probe and folder-filename.
    colnames(data) <- c("cg_probe", col_name)
    
    # Merge data
    if (length(all_data) == 0) {
      all_data <- list(data)
      cat("initialization\n") 
    } else {
      cat("Data being merged...\n")  
      all_data <- mapply(function(x, y) merge(x, y, by = "cg_probe", all = TRUE),
                         all_data, list(data), SIMPLIFY = FALSE)
    }
  }
  
  # Merge all data into a single data.frame
  cat("Start merging data...\n") 
  final_data <- Reduce(function(x, y) merge(x, y, by = "cg_probe", all = TRUE), all_data)
  
  cat("Data merging complete, number of rows：", nrow(final_data), " cols：", ncol(final_data), "\n") 
  
  # Output the merged result
  output_path <- paste0(base_path, "_merged_methylation_data.txt")
  write.table(final_data, file = output_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("Data integration is complete and has been saved as", output_path, "\n") 
  
  cat("Data processing completed！\n")
  
  return(final_data)  
}

# Calling functions to process data
cell_types <- c("bcell", "cd4", "cd8", "nk", "monocyte", "neutrophil")
for (cell_type in cell_types) {
  base_path <- paste0("ref_data/", cell_type)
  processed_data <- process_methylation_data(base_path)
  
  # Move all generated result files to the specified path.
  output_dir <- "houseman_ref"
  file.rename(paste0(base_path, "_merged_methylation_data.txt"), 
              file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")))
  cat("The file has been moved to：", file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")), "\n")
}

# Integrate the data to generate a merged_data data frame and a merged_data_matrix.

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
  cat("Merged data dimensions：", nrow(merged_data), "row，", ncol(merged_data), "col\n")
  
  # Delete rows with the attribute value NA
  merged_data <- merged_data %>% filter(!apply(merged_data, 1, function(row) any(row == "NA")))
  
  # Print the data dimensions after deleting NA rows.
  cat("Data dimensions after deleting entries containing NA：", nrow(merged_data), "row，", ncol(merged_data), "col\n")
  
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
base_path <- "housemman_ref"
output_file <- "houseman_ref/merged_data.txt"
result <- process_methylation_matrix(base_path, output_file)

merged_data <- result$merged_data
merged_data_matrix <- result$merged_data_matrix

# If you need to view the first few rows of the matrix
if (!is.null(merged_data_matrix)) {
  cat("The first few rows of the merged data matrix：\n")
  print(head(merged_data_matrix))
}

beta_ref <- merged_data_matrix
# Retrieve the relevant cell types from merged_data_matrix.
cell_ref <- sub("\\..*$", "", colnames(merged_data_matrix))

# Save cell_ref and beta_ref
out_dir <- "houseman_ref"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1) Save as RDS
saveRDS(cell_ref, file = file.path(out_dir, "cell_ref_own.rds"))
saveRDS(beta_ref, file = file.path(out_dir, "beta_ref_own.rds"))

cat("Saved to:", out_dir, "\n")

# Define a save path to save the generated reference matrix.
output_path <- "hosueman_ref/reference_output_houseman.csv"

# Save averaged_df as a CSV file
write.csv(merged_data_matrix, file = output_path, row.names = TRUE, col.names = TRUE)

# print
cat("File saved to:", output_path)
