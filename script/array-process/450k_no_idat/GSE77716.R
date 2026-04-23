setwd("/data/zhangmch/ewas_array/data/GSE77716")
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
file_path <- "./GSE77716_Matrix_signal_intensities_updated_2018-Jan08.tsv"

GSE77716 <- read.table(
  file = file_path,
  header = TRUE,          # First row contains column names
  sep = "\t",             # Tab-separated
  stringsAsFactors = FALSE,
  check.names = FALSE     # Keep original column names (e.g., 81136_UnMethylated_signal)
)

# Check first few rows
head(GSE77716)
str(GSE77716)

########################### Construct mapping relationship

mapping_data <- read.table("/data/zhangmch/ewas_array/data/GSE77716/mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract GSM_ID and Column_Prefix vectors
GSM_ID <- mapping_data$V1
Column_Prefix <- mapping_data$V2

gsm_mapping <- data.frame(GSM_ID, Column_Prefix)

# Print mapping table
print(gsm_mapping)

######################### Replace column names with GSM
# Create a function to replace column prefixes with corresponding GSM IDs
replace_column_names <- function(colnames, gsm_mapping) {
  new_colnames <- colnames
  # Iterate through each row in gsm_mapping
  for (i in 1:nrow(gsm_mapping)) {
    # Modify regex in gsub to ensure exact matching of prefix
    new_colnames <- gsub(paste0("\\b", gsm_mapping$Column_Prefix[i], "\\b"), gsm_mapping$GSM_ID[i], new_colnames)
  }
  return(new_colnames)
}

# Replace column names
colnames(GSE77716) <- replace_column_names(colnames(GSE77716), gsm_mapping)

# Print new column names
print(colnames(GSE77716))

GSE77716 <- as.data.frame(GSE77716)
rownames(GSE77716) <- GSE77716$ID_REF


################################### Split into M, U, and P matrices
# Remove ID_REF column
GSE77716$ID_REF <- NULL

# Extract column names for UnMethylated, Methylated, and Pval
u_cols <- grep("Unmethylated Signal$", colnames(GSE77716), value = TRUE)
m_cols <- grep("Methylated Signal$",  colnames(GSE77716), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE77716), value = TRUE)

# Extract sample IDs (e.g., SAMPLE 1) and remove suffix
u_ids <- sub("^(GSM\\d+).*", "\\1", u_cols)
m_ids <- sub("^(GSM\\d+).*", "\\1", m_cols)
p_ids <- sub("^(GSM\\d+).*", "\\1", p_cols)

# Ensure consistent order (avoid column mismatch)
sample_ids <- sort(unique(u_ids))

# Match columns by sample IDs
u_cols <- u_cols[match(sample_ids, u_ids)]
m_cols <- m_cols[match(sample_ids, m_ids)]
p_cols <- p_cols[match(sample_ids, p_ids)]

# Combine into matrices
U <- as.matrix(GSE77716[, u_cols, drop = FALSE])
M <- as.matrix(GSE77716[, m_cols, drop = FALSE])
P <- as.matrix(GSE77716[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE77716)
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

# Define function to save data
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

# Call function to save
output_dir <- "/data/zhangmch/ewas_array/result/GSE77716"
save_columns_as_files(beta, output_dir)
