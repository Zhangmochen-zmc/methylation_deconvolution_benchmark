##### step 1: Extract all data used to construct the reference matrix from each cell type and merge them.

library(limma)
library(dplyr)
library(tools)
library(tibble)

# Encapsulate data processing as functions
process_methylation_data <- function(base_path) {
  all_data <- list() 
  
  # get txt file
  files <- list.files(base_path, pattern = "\\.txt$", full.names = TRUE) 
  
  cat("processing file：", base_path, "\n")  
  
  # loop through each file and integrate the data.
  for (file in files) {
    cat("processing file：", file, "\n")
    
    # read data
    data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # get file_name
    file_name <- file_path_sans_ext(basename(file))
    
    # assign column names the format "folder name-file name".
    # assign column names based on cell type
    col_name <- paste(basename(base_path), file_name, sep = "-")
    
    # set column names: cg_probe and folder-filename.
    colnames(data) <- c("cg_probe", col_name)
    
    # merge data: merge by cg_probe
    if (length(all_data) == 0) {
      all_data <- list(data)
      cat("initialization all_data\n")
    } else {
      cat("data merging in progress...\n") 
      all_data <- mapply(function(x, y) merge(x, y, by = "cg_probe", all = TRUE),
                         all_data, list(data), SIMPLIFY = FALSE)
    }
  }
  
  # merge all data into a single data.frame
  cat("data merging in progress...\n") 
  final_data <- Reduce(function(x, y) merge(x, y, by = "cg_probe", all = TRUE), all_data)
  
  cat("data merging complete, number of rows：", nrow(final_data), " rows ：", ncol(final_data), "\n") 
  
  # output the merged result
  output_path <- paste0(base_path, "_merged_methylation_data.txt")
  write.table(final_data, file = output_path, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("data integration is complete and has been saved as", output_path, "\n")  
  
  cat("data processing done！\n")
  
  return(final_data) 
}

# calling functions to process data
cell_types <- c("bcell", "cd4", "cd8", "nk", "monocyte", "neutrophil")
for (cell_type in cell_types) {
  base_path <- paste0("ref_data/", cell_type)
  processed_data <- process_methylation_data(base_path)
  
  # move all generated result files to the specified path.
  output_dir <- "emeth_ref"
  file.rename(paste0(base_path, "_merged_methylation_data.txt"), 
              file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")))
  cat("The file has been moved to ：", file.path(output_dir, paste0(cell_type, "_merged_methylation_data.txt")), "\n")
}



##### Step 2: Integrate the data to generate a merged_data data frame and a merged_data_matrix matrix.
                       
# Encapsulate data processing as functions
process_methylation_matrix <- function(base_path, output_file) {
  # read merged data from various cell types
  bcell <- read.table(file.path(base_path, "bcell_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd4 <- read.table(file.path(base_path, "cd4_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cd8 <- read.table(file.path(base_path, "cd8_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  neutrophil <- read.table(file.path(base_path, "neutrophil_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  monocyte <- read.table(file.path(base_path, "monocyte_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  nk <- read.table(file.path(base_path, "nk_merged_methylation_data.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # modify column names to reflect cell type
  colnames(cd4) <- gsub("bcell", "cd4", colnames(cd4))
  colnames(cd8) <- gsub("bcell", "cd8", colnames(cd8))
  colnames(nk) <- gsub("bcell", "nk", colnames(nk))
  colnames(neutrophil) <- gsub("bcell", "neutrophil", colnames(neutrophil))
  colnames(monocyte) <- gsub("bcell", "monocyte", colnames(monocyte))
  
  # merge_data
  merged_data <- merge(bcell, cd4, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, cd8, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, neutrophil, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, monocyte, by = "cg_probe", all = TRUE)
  merged_data <- merge(merged_data, nk, by = "cg_probe", all = TRUE)
  
  # print
  cat("merged data dimensions：", nrow(merged_data), "row，", ncol(merged_data), "col\n")
  
  # delete NA
  merged_data <- merged_data %>% filter(!apply(merged_data, 1, function(row) any(row == "NA")))
  
  # print the data dimensions after deleting NA rows.
  cat("data dimensions after deleting entries containing NA：", nrow(merged_data), "row，", ncol(merged_data), "col\n")
  
  # output merged data
  write.table(merged_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  cat("the data has been saved as：", output_file, "\n")
  
  # convert merged_data to matrix form and set the row name to cg_probe.
  if ("cg_probe" %in% colnames(merged_data)) {
    merged_data_matrix <- as.matrix(merged_data[, -which(colnames(merged_data) == "cg_probe")])  
    rownames(merged_data_matrix) <- merged_data$cg_probe  
    
    # return the transformed matrix
    return(list(merged_data = merged_data, merged_data_matrix = merged_data_matrix))
  } else {
    cat("warning: The column 'cg_probe' does not exist in the data frame!\n")
    return(NULL)
  }
}

# calling functions to process data
base_path <- "emeth_ref"
output_file <- "emeth_ref/merged_data.txt"
result <- process_methylation_matrix(base_path, output_file)

merged_data <- result$merged_data
merged_data_matrix <- result$merged_data_matrix

# if you need to view the first few rows of the matrix
if (!is.null(merged_data_matrix)) {
  cat("the first few rows of the merged data matrix：\n")
  print(head(merged_data_matrix))
}



####### step 3: Generate the corresponding number of samples for each cell type, and then run differential analysis on each one to screen DMCs.

# Encapsulate data processing as functions
process_cell_type_data_dynamic <- function(merged_data_matrix, cell_types, output_dir) {
  
  # save results
  results <- list()
  
  # traverse each cell type
  for (cell_type in cell_types) {
    
    # count the number of samples for the current cell type
    case_count <- length(grep(cell_type, colnames(merged_data_matrix)))  # case: current cell type
    control_count <- sum(sapply(cell_types, function(ct) length(grep(ct, colnames(merged_data_matrix))))) - case_count  # control: other cell type
    
    # get the column name containing the current cell type
    cell_columns <- grep(cell_type, colnames(merged_data_matrix), value = TRUE)
    
    # get column names that do not contain the current cell type
    other_columns <- setdiff(colnames(merged_data_matrix), cell_columns)
    
    # move the column for the current cell type to the front, and all other columns to the back.
    cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]
    
    # create a design matrix to represent groups.
    group <- factor(c(rep("case", case_count), rep("control", control_count)))
    group <- factor(group, levels = c("case", "control"), ordered = FALSE)
    design <- model.matrix(~ group)
    colnames(design) <- levels(group)
    
    # lmFit
    fit <- lmFit(cell_matrix, design)
    
    # compare and apply bayesian adjustment.
    contrast.matrix <- makeContrasts(case - control, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # extract all significant genes
    top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")
    
    # save results
    output_path <- file.path(output_dir, paste0(cell_type, "_DMC_results.txt"))
    write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
    cat("The data has been saved as ：", output_path, "\n")
    
    # filter rows where the p-value is less than 0.01.
    DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ] 
    DMC <- na.omit(DMC)  
    
    # save
    results[[cell_type]] <- DMC
  }
  
  # return differential gene results for all cell types
  return(results)
}

# usage
output_dir <- "emeth_ref"

# define cell type list
cell_types <- c("bcell", "cd4", "cd8", "neutrophil", "monocyte", "nk")

# calling functions to process data
cell_type_results <- process_cell_type_data_dynamic(merged_data_matrix, cell_types, output_dir)

# view the dimensions of DMC results for each cell type
lapply(cell_type_results, dim)


##### screening for differential methylation sites for the remaining cell types
# 1、B cell、CD4+ Tcell、CD8+ Tcell to Monocyte、neutrophil、NK cell
# 2、B cell to NK 
# 3、B cell、NK cell to CD4+ Tcell、CD8+ Tcell
# 4、CD4+ T cell to CD8+ Tcell
# 5、neutrophil VS Monocyte



########################### B cell and NK cell DMC

cell_types_b_nk <- c("bcell", "nk")

# count the number of samples for the current cell type
case_count <- length(grep(cell_types_b_nk[1], colnames(merged_data_matrix)))  # case : current cell type
control_count <- length(grep(cell_types_b_nk[2], colnames(merged_data_matrix)))  # control : other cell type

# get the column name containing the current cell type
cell_columns <- grep(cell_types_b_nk[1], colnames(merged_data_matrix), value = TRUE)

# get column names that do not contain the current cell type
other_columns <- grep(cell_types_b_nk[2], colnames(merged_data_matrix), value = TRUE)

# move the column for the current cell type to the front, and all other columns to the back.
cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]

# create a design matrix to represent groups.
group <- factor(c(rep("case", case_count), rep("control", control_count)))
group <- factor(group, levels = c("case", "control"), ordered = FALSE)
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# lmFit
fit <- lmFit(cell_matrix, design)

# compare and apply Bayesian adjustment.
contrast.matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract all significant genes
top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")

# save results
output_path <- file.path(output_dir, "bcell_nkcell_DMC_results.txt")
write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
cat("the data has been saved as：", output_path, "\n")

# filter rows where the p-value is less than 0.01.
DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ]  
DMC <- na.omit(DMC)  

# save results
cell_type_results[["bcell_nkcell"]] <- DMC


########################### CD4+ T cell and CD8+ T cell DMC

cell_types_cd4_cd8 <- c("cd4", "cd8")

# count the number of samples for the current cell type
case_count <- length(grep(cell_types_cd4_cd8[1], colnames(merged_data_matrix)))  # case : current cell type
control_count <- length(grep(cell_types_cd4_cd8[2], colnames(merged_data_matrix)))  # control : other cell type

# get the column name containing the current cell type
cell_columns <- grep(cell_types_cd4_cd8[1], colnames(merged_data_matrix), value = TRUE)

# get column names that do not contain the current cell type
other_columns <- grep(cell_types_cd4_cd8[2], colnames(merged_data_matrix), value = TRUE)

# move the column for the current cell type to the front, and all other columns to the back.
cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]

# create a design matrix to represent groups.
group <- factor(c(rep("case", case_count), rep("control", control_count)))
group <- factor(group, levels = c("case", "control"), ordered = FALSE)
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# lmFit
fit <- lmFit(cell_matrix, design)

# compare and apply Bayesian adjustment.
contrast.matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract all significant genes
top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")

# save results
output_path <- file.path(output_dir, "cd4tcell_cd8tcell_DMC_results.txt")
write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
cat("the data has been saved as：", output_path, "\n")

# filter rows where the p-value is less than 0.01.
DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ]  
DMC <- na.omit(DMC)  

# save results
cell_type_results[["cd4tcell_cd8tcell"]] <- DMC




########################### neutrophil cell and monocyte cell DMC

cell_types_neutrophil_monocyte <- c("neutrophil", "monocyte")

# count the number of samples for the current cell type
case_count <- length(grep(cell_types_neutrophil_monocyte[1], colnames(merged_data_matrix)))  # case : current cell type
control_count <- length(grep(cell_types_neutrophil_monocyte[2], colnames(merged_data_matrix)))  # control : other cell type

# get the column name containing the current cell type
cell_columns <- grep(cell_types_neutrophil_monocyte[1], colnames(merged_data_matrix), value = TRUE)

# get column names that do not contain the current cell type
other_columns <- grep(cell_types_neutrophil_monocyte[2], colnames(merged_data_matrix), value = TRUE)

# move the column for the current cell type to the front, and all other columns to the back.
cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]

# create a design matrix to represent groups.
group <- factor(c(rep("case", case_count), rep("control", control_count)))
group <- factor(group, levels = c("case", "control"), ordered = FALSE)
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# lmFit
fit <- lmFit(cell_matrix, design)

# compare and apply Bayesian adjustment.
contrast.matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract all significant genes
top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")

# save results
output_path <- file.path(output_dir, "neutrophilcell_monocytecell_DMC_results.txt")
write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
cat("The data has been saved as ：", output_path, "\n")

# filter rows where the p-value is less than 0.01.
DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ]  
DMC <- na.omit(DMC)  

# save results
cell_type_results[["neutrophilcell_monocytecell"]] <- DMC




########################### b cell and nk cell and cd4+ t cell and cd8+ t cell DMC

cell_types_bcell_nk <- c("bcell", "nk")
cell_types_cd4_cd8 <- c("cd4", "cd8")

# count the number of samples for the current cell type
case_count <- length(grep(cell_types_bcell_nk[1], colnames(merged_data_matrix))) + length(grep(cell_types_bcell_nk[2], colnames(merged_data_matrix)))  # case : current cell type
control_count <- length(grep(cell_types_cd4_cd8[1], colnames(merged_data_matrix))) + length(grep(cell_types_cd4_cd8[2], colnames(merged_data_matrix))) # control : other cell type

# get the column names containing the current cell type
cell_columns_1 <- grep(cell_types_bcell_nk[1], colnames(merged_data_matrix), value = TRUE)
cell_columns_2 <- grep(cell_types_bcell_nk[2], colnames(merged_data_matrix), value = TRUE)
cell_columns <- c(cell_columns_1, cell_columns_2)

# get the column names not containing the current cell type
other_columns_1 <- grep(cell_types_cd4_cd8[1], colnames(merged_data_matrix), value = TRUE)
other_columns_2 <- grep(cell_types_cd4_cd8[2], colnames(merged_data_matrix), value = TRUE)
other_columns <- c(other_columns_1, other_columns_2)

# place the column for the current cell type first, and the other columns last.
cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]

# create a design matrix to represent the groups
group <- factor(c(rep("case", case_count), rep("control", control_count)))
group <- factor(group, levels = c("case", "control"), ordered = FALSE)
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# lmFit
fit <- lmFit(cell_matrix, design)

# perform comparisons and apply Bayesian adjustment
contrast.matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract all significant genes
top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")

# save results
output_path <- file.path(output_dir, "bnkcell_cd4cd8cell_DMC_results.txt")
write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
cat("Data has been saved as ：", output_path, "\n")

# filter rows with p-values ​​less than 0.01
DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ]  
DMC <- na.omit(DMC)  

# save results
cell_type_results[["bnkcell_cd4cd8cell"]] <- DMC


########################### b cell, cd4+ t cell, and cd8+ t cell and monocyte, neutrophil, nk DMC

cell_types_bcell_cd4_cd8 <- c("bcell", "cd4", "cd8")
cell_types_nk_monocyte_neutrophil <- c("nk", "monocyte", "neutrophil")

# count the number of samples for the current cell type
case_count <- length(grep(cell_types_bcell_cd4_cd8[1], colnames(merged_data_matrix))) + length(grep(cell_types_bcell_cd4_cd8[2], colnames(merged_data_matrix))) + + length(grep(cell_types_bcell_cd4_cd8[3], colnames(merged_data_matrix))) # case 组数量是该细胞类型的样本数
control_count <- length(grep(cell_types_nk_monocyte_neutrophil[1], colnames(merged_data_matrix))) + length(grep(cell_types_nk_monocyte_neutrophil[2], colnames(merged_data_matrix))) + length(grep(cell_types_nk_monocyte_neutrophil[3], colnames(merged_data_matrix)))# case 组数量是该细胞类型的样本数

# get the column names containing the current cell type
cell_columns_1 <- grep(cell_types_bcell_cd4_cd8[1], colnames(merged_data_matrix), value = TRUE)
cell_columns_2 <- grep(cell_types_bcell_cd4_cd8[2], colnames(merged_data_matrix), value = TRUE)
cell_columns_3 <- grep(cell_types_bcell_cd4_cd8[3], colnames(merged_data_matrix), value = TRUE)
cell_columns <- c(cell_columns_1, cell_columns_2, cell_columns_3)

# get the column names not containing the current cell type
other_columns_1 <- grep(cell_types_nk_monocyte_neutrophil[1], colnames(merged_data_matrix), value = TRUE)
other_columns_2 <- grep(cell_types_nk_monocyte_neutrophil[2], colnames(merged_data_matrix), value = TRUE)
other_columns_3 <- grep(cell_types_nk_monocyte_neutrophil[3], colnames(merged_data_matrix), value = TRUE)
other_columns <- c(other_columns_1, other_columns_2, other_columns_3)

# place the column for the current cell type first, and the other columns last
cell_matrix <- merged_data_matrix[, c(cell_columns, other_columns)]

# create a design matrix to represent the groups
group <- factor(c(rep("case", case_count), rep("control", control_count)))
group <- factor(group, levels = c("case", "control"), ordered = FALSE)
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# lmFit
fit <- lmFit(cell_matrix, design)

# perform comparisons and apply Bayesian adjustment
contrast.matrix <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# extract all significant genes
top_genes_all <- topTable(fit2, number = Inf, adjust = "BH")

# save results
output_path <- file.path(output_dir, "bcd4cd8cell_monocyteneutrophilnk_DMC_results.txt")
write.table(top_genes_all, file = output_path, sep = "\t", row.names = TRUE, col.names = TRUE)
cat("Data has been saved as：", output_path, "\n")

# filter rows with p-values ​​less than 0.01
DMC <- top_genes_all[top_genes_all[, 5] < 0.01 & abs(top_genes_all[, "logFC"]) > 1, ]  
DMC <- na.omit(DMC)  

# save results
cell_type_results[["bcd4cd8cell_monocyteneutrophilnk"]] <- DMC

#################################### from these, select the top 100 significant comparisons for each comparison.

# sort each data frame in cell_type_results by the logFC column
cell_type_results_sorted <- lapply(cell_type_results, function(df) {
  df_sorted <- df[order(df$logFC, decreasing = TRUE), ]  
  return(df_sorted)
})

all_comparisons_DMC <- do.call(rbind, lapply(cell_type_results_sorted, function(df) {
  head(df, 100)  
}))

# filter for the names of non-repeating CG sites (rownames)
cg_names <- rownames(all_comparisons_DMC)


# split the string and extract the part after "cg"
cg_only <- sapply(cg_names, function(name) {
  parts <- strsplit(name, "\\.")[[1]]
  cg_part <- parts[length(parts)] 
  return(cg_part)
})

# filter for unique CG sites
unique_cgs <- unique(cg_only)

# count the number of unique CG sites
length(unique_cgs)

#################################### filter for the DNA methylation corresponding to these specific sites
# extract rows from merged_data whose row names are in unique_cgs
filtered_data <- merged_data[merged_data$cg_probe %in% unique_cgs, ]

# extract the required column names
bcell_cols <- grep("bcell", colnames(filtered_data), value = TRUE)
cd4_cols <- grep("cd4", colnames(filtered_data), value = TRUE)
cd8_cols <- grep("cd8", colnames(filtered_data), value = TRUE)
neutrophil_cols <- grep("neutrophil", colnames(filtered_data), value = TRUE)
monocyte_cols <- grep("monocyte", colnames(filtered_data), value = TRUE)
nk_cols <- grep("nk", colnames(filtered_data), value = TRUE)

# calculate the mean for each cell type
bcell_avg <- rowMeans(filtered_data[, bcell_cols], na.rm = TRUE)
cd4_avg <- rowMeans(filtered_data[, cd4_cols], na.rm = TRUE)
cd8_avg <- rowMeans(filtered_data[, cd8_cols], na.rm = TRUE)
neutrophil_avg <- rowMeans(filtered_data[, neutrophil_cols], na.rm = TRUE)
monocyte_avg <- rowMeans(filtered_data[, monocyte_cols], na.rm = TRUE)
nk_avg <- rowMeans(filtered_data[, nk_cols], na.rm = TRUE)

# store the mean as a new data frame
avg_data <- data.frame(
  cg_probe = filtered_data$cg_probe,  # 保留 cg_probe 列
  bcell_avg = bcell_avg,
  cd4_avg = cd4_avg,
  cd8_avg = cd8_avg,
  neutrophil_avg = neutrophil_avg,
  monocyte_avg = monocyte_avg,
  nk_avg = nk_avg
)

# cg_probe
rownames(avg_data) <- avg_data$cg_probe

# delete cg_probe
avg_data <- avg_data[, -which(names(avg_data) == "cg_probe")]

# delete _avg
colnames(avg_data) <- gsub("_avg", "", colnames(avg_data))

# show results
# head(avg_data)

# save avg_data as an EMeth reference format that can be used as input
avg_data_matrix <- as.matrix(avg_data)

# save results
save(avg_data_matrix, file = "emeth_ref/avg_data_matrix.RData")
