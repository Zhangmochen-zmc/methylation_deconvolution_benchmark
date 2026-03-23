library(MeDeCom)

# Configure parameters (modify here based on the PDF results from script decon.R).
best_K <- 6
best_lambda <- 0.001 # Please adjust this value according to the chart.

# Load the results and raw data
medecom_result <- readRDS("medecomrds/test_data.rds")
bulk_data <- read.csv("test_data/test_data.csv", row.names = 1, check.names = FALSE)
# LMCs 和 Proportions
cat("Extracting LMCs and proportions...\n")
lmcs <- getLMCs(medecom_result, K=best_K, lambda=best_lambda)
rownames(lmcs) <- rownames(bulk_data)

props <- getProportions(medecom_result, K=best_K, lambda=best_lambda)

# Load reference data for cell type identification.
cat("Processing reference data for cell type matching...\n")
meth_data <- read.table("ref_data/refdata.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
meth_data <- as.matrix(meth_data)
meta_data <- read.csv("ref_data/refmeta.csv")

# Calculate the mean (Tref) for each cell type in the reference data.
unique_cell_types <- unique(meta_data$Tissue)
Tref <- sapply(unique_cell_types, function(ct) {
  samples_for_ct <- meta_data$FileID[meta_data$Tissue == ct]
  subset_data <- meth_data[, samples_for_ct, drop = FALSE]
  rowMeans(subset_data, na.rm = TRUE)
})
Tref <- as.matrix(Tref)

# Correlation analysis and label matching (using the Hungarian algorithm to ensure 1-to-1 matching)
cat("Matching components to cell types using optimal assignment...\n")
subset_markers <- intersect(rownames(Tref), rownames(lmcs))
cor_matrix <- cor(lmcs[subset_markers, ], Tref[subset_markers, ])

# check
if(nrow(cor_matrix) != ncol(cor_matrix)) {
  cat("Warning: Number of components (K) and reference cell types are not equal!\n")
}

# The Hungarian algorithm defaults to minimizing the loss, so we use (1 - correlation) as the cost matrix.
if (!requireNamespace("clue", quietly = TRUE)) install.packages("clue")
library(clue)

# Calculate the optimal allocation (which maximizes the overall correlation).
cost_matrix <- 1 - cor_matrix
cost_matrix[cost_matrix < 0] <- 0 

assignment <- solve_LSAP(as.matrix(cost_matrix), maximum = FALSE)

# assignment 
assigned_labels <- colnames(cor_matrix)[as.numeric(assignment)]

# print
matching_summary <- data.frame(
  Component = paste0("LMC_", 1:best_K),
  Assigned_CellType = assigned_labels,
  Correlation = sapply(1:best_K, function(i) cor_matrix[i, assignment[i]])
)
print(matching_summary)

# Organize and save the results
props_final <- t(props)
colnames(props_final) <- assigned_labels

# The sequences were rearranged according to the standard order of the reference cell type
target_order <- unique_cell_types 
props_final <- props_final[, target_order, drop = FALSE]

# save
write.csv(props_final, "medecom_result/test_data.csv") # Modify by yourself
cat("Success! 1-to-1 matching complete. Final proportions saved.\n")
