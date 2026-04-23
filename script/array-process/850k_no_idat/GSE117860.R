setwd("/data/zhangmch/ewas_array/data/GSE117860")
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
file_path <- "./GSE117860_MethylatedSignal.txt"

GSE117860 <- read.table(
  file = file_path,
  header = TRUE,          # First row contains column names
  sep = "\t",             # Tab-separated
  stringsAsFactors = FALSE,
  check.names = FALSE     # Keep original column names (e.g., 81136_UnMethylated_signal)
)

colnames(GSE117860)[1] <- "ID_REF"

# Check first few rows
head(GSE117860)
str(GSE117860)

########################### Construct mapping relationship
# Read mapping file
mapping_data <- read.table("/data/zhangmch/ewas_array/data/GSE117860/mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract GSM_ID and Column_Prefix vectors
GSM_ID <- mapping_data$V2
Column_Prefix <- mapping_data$V1

# Create gsm_mapping data frame
gsm_mapping <- data.frame(GSM_ID, Column_Prefix)

# Print mapping table
print(gsm_mapping)

# Check column names of GSE166844 and remove prefix 'X'
colnames(GSE117860) <- gsub("^X", "", colnames(GSE117860))

# Check if column names are correct
print(colnames(GSE117860))

# Ensure matches exist during replacement
#for (i in 1:nrow(gsm_mapping)) {
#  # Check for matching prefixes
#  matching_cols <- grep(paste0("^", gsm_mapping$Column_Prefix[i]), colnames(GSE117860), value = TRUE)
#  if (length(matching_cols) > 0) {
#    print(paste("Matched columns:", matching_cols))
#  }
#  
#  # Replace Column_Prefix with corresponding GSM_ID
#  colnames(GSE117860) <- gsub(paste0("^", gsm_mapping$Column_Prefix[i]), 
#                              gsm_mapping$GSM_ID[i], 
#                              colnames(GSE117860))
#}

# Avoid matching "sample 11" with "sample 1"
for (i in 1:nrow(gsm_mapping)) {
  prefix <- gsm_mapping$Column_Prefix[i]
  gsmi   <- gsm_mapping$GSM_ID[i]
  
  # Prefix must be followed by whitespace or end of string
  pattern <- paste0("^", prefix, "(\\s|$)")
  
  # Find matching columns
  matching_cols <- grep(pattern, colnames(GSE117860), value = TRUE, perl = TRUE)
  if (length(matching_cols) > 0) {
    print(paste("Matched columns:", paste(matching_cols, collapse = ", ")))
  }
  
  # Replace only the prefix while keeping suffix
  colnames(GSE117860) <- sub(
    pattern,
    paste0(gsmi, "\\1"),   # \\1 keeps (\\s|$)
    colnames(GSE117860),
    perl = TRUE
  )
}

# Print replaced column names
print("Replaced column names:")
print(colnames(GSE117860))


GSE117860 <- as.data.frame(GSE117860)
rownames(GSE117860) <- GSE117860$ID_REF


################################### Split into M, U, and P matrices
# Remove ID_REF column
GSE117860$ID_REF <- NULL

# Extract column names for UnMethylated, Methylated, and Pval
u_cols <- grep("Unmethylated signal$", colnames(GSE117860), value = TRUE)
m_cols <- grep("Methylated signal$",  colnames(GSE117860), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE117860), value = TRUE)

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
U <- as.matrix(GSE117860[, u_cols, drop = FALSE])
M <- as.matrix(GSE117860[, m_cols, drop = FALSE])
P <- as.matrix(GSE117860[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE117860)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids


################################### Calculate Beta values and perform QC using P-values
offset <- 100

beta <- M / (M + U + offset)

# QC: set values with P > 0.01 to NA
beta[P > 0.01] <- NA


beta_na = data.frame(beta)

############################# BMIQ normalization for Type I and Type II probes
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="EPIC",cores=24)
beta_BMIQ <- myNorm

colnames(beta_BMIQ) = gsub("_.*?$","", colnames(beta_BMIQ)) 

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 


############################### Beta QC: remove low-quality samples and probes, then KNN imputation
beta = beta_BMIQ
NA_r = rowSums(is.na(beta)) 
NA_c = colSums(is.na(beta))
beta_remain = beta[, which(NA_c <= 865859*0.15) ]

beta_remain = beta_remain[which(NA_r <=  dim(beta_remain)[2]/10),] 
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)

########################## Remove SNP probes
load("/data/zhangmch/ewas_array/script/850k/snp_cg_850K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

######################## Remove sex chromosome probes
load("/data/zhangmch/ewas_array/script/850k/850K_cg_annotation.RData")
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
output_dir <- "/data/zhangmch/ewas_array/result/GSE117860"
save_columns_as_files(beta, output_dir)
