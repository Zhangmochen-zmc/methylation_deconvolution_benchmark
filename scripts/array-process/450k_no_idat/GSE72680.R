setwd("/data/zhangmch/ewas_array/data/GSE72680")
library(GEOquery)
library(gmqn) 
library(minfi)
library(gplots)
library('limma')
library('WGCNA')
library(impute)
library(IlluminaHumanMethylationEPICmanifest)
library(wateRmelon)
library(ChAMP)
library(R.utils)

################################################ Read data provided by the authors
file_path <- "./GSE72680_signals.txt"

GSE72680 <- read.table(
  file = file_path,
  header = TRUE,          # First row contains column names
  sep = "\t",             # Tab-separated
  stringsAsFactors = FALSE,
  check.names = FALSE     # Keep original column names (e.g., 81136_UnMethylated_signal)
)

# Check first few rows
head(GSE72680)
str(GSE72680)

########################### Construct mapping relationship
# Read mapping file
mapping_data <- read.table("/data/zhangmch/ewas_array/data/GSE72680/mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract GSM_ID and Column_Prefix vectors
GSM_ID <- mapping_data$V1
Column_Prefix <- mapping_data$V2

# Create gsm_mapping data frame
gsm_mapping <- data.frame(GSM_ID, Column_Prefix)

# Print mapping table
print(gsm_mapping)

# Check column names of GSE72680 and remove prefix 'X'
colnames(GSE72680) <- gsub("^X", "", colnames(GSE72680))

# Check if column names are correct
print(colnames(GSE72680))

# Ensure matches exist during replacement
for (i in 1:nrow(gsm_mapping)) {
  # Check for matching prefixes
  matching_cols <- grep(paste0("^", gsm_mapping$Column_Prefix[i]), colnames(GSE72680), value = TRUE)
  if (length(matching_cols) > 0) {
    print(paste("Matched columns:", matching_cols))
  }
  
  # Replace Column_Prefix with corresponding GSM_ID
  colnames(GSE72680) <- gsub(paste0("^", gsm_mapping$Column_Prefix[i]), 
                              gsm_mapping$GSM_ID[i], 
                              colnames(GSE72680))
}

# Print replaced column names
print("Replaced column names:")
print(colnames(GSE72680))

GSE72680 <- as.data.frame(GSE72680)
rownames(GSE72680) <- GSE72680[,1]

################################### Split into M, U, and P matrices
# Remove ID_REF column
GSE72680[,1] <- NULL

# Column names for UnMethylated / Methylated / Pval
u_cols <- grep("_Unmethylated signal$", colnames(GSE72680), value = TRUE)
m_cols <- grep("_Methylated signal$",  colnames(GSE72680), value = TRUE)
p_cols <- grep("_Detection PVal$",     colnames(GSE72680), value = TRUE)

# Extract sample IDs (e.g., 81136, 81084)
u_ids <- sub("_Unmethylated signal$", "", u_cols)
m_ids <- sub("_Methylated signal$",   "", m_cols)
p_ids <- sub("_Detection PVal$",      "", p_cols)

# Ensure consistent order (avoid column mismatch)
sample_ids <- sort(unique(u_ids))

u_cols <- u_cols[match(sample_ids, u_ids)]
m_cols <- m_cols[match(sample_ids, m_ids)]
p_cols <- p_cols[match(sample_ids, p_ids)]

# Combine into matrices
U <- as.matrix(GSE72680[, u_cols, drop = FALSE])
M <- as.matrix(GSE72680[, m_cols, drop = FALSE])
P <- as.matrix(GSE72680[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE72680)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids


################################### Calculate Beta values and perform QC using P-values
offset <- 100

beta <- M / (M + U + offset)

# QC: set values with P > 0.01 to NA
beta[P > 0.01] <- NA


beta_na = data.frame(beta)

############################# BMIQ normalization for Type I and Type II probes
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="450K",cores=24)
beta_BMIQ <- myNorm

colnames(beta_BMIQ) = gsub("_.*?$","", colnames(beta_BMIQ)) 

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 


############################### Beta QC: remove low-quality samples and probes, then KNN imputation
beta = beta_BMIQ
NA_r = rowSums(is.na(beta)) 
NA_c = colSums(is.na(beta))
beta_remain = beta[, which(NA_c <= 485512*0.15) ]

beta_remain = beta_remain[which(NA_r <=  dim(beta_remain)[2]/10),] 
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)

########################## Remove SNP probes
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

######################## Remove sex chromosome probes
load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")
gene_annotation <- b
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrX"),]
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrY"),]
beta = beta[ intersect(row.names(beta), gene_annotation$X1),]

###################### Save beta values

# Function to save data
save_columns_as_files <- function(data, output_dir = getwd()) {
  # Get column names (sample names)
  column_names <- colnames(data)
  
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each column
  for (col_name in column_names) {
    # Create a data frame with cg probe IDs and values
    result <- data.frame(cg_id = rownames(data), attribute_value = data[[col_name]], stringsAsFactors = FALSE)
    
    # Define file name
    file_name <- paste0(col_name, ".txt")
    
    # Create full file path
    file_path <- file.path(output_dir, file_name)
    
    # Write data using tab separation
    write.table(result, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Print saved file name (optional)
    cat("File saved as:", file_path, "\n")
  }
}

# Call function to save results
output_dir <- "/data/zhangmch/ewas_array/result/GSE72680"
save_columns_as_files(beta, output_dir)
