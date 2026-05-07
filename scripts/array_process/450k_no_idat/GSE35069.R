setwd("/data/zhangmch/ewas_array/data/GSE35069")
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
file_path <- "./GSE35069_Matrix_signal_intensities.txt"

GSE35069 <- read.table(
  file = file_path,
  header = TRUE,          # First row contains column names
  sep = "\t",             # Tab-separated values
  stringsAsFactors = FALSE,
  check.names = FALSE     # Keep original column names (e.g., 81136_UnMethylated_signal)
)

# View first few rows for confirmation
head(GSE35069)
str(GSE35069)

########################### Construct mapping relationship
gsm_mapping <- data.frame(
  GSM_ID = c("GSM861635", "GSM861636", "GSM861637", "GSM861638", "GSM861639", 
             "GSM861640", "GSM861641", "GSM861642", "GSM861643", "GSM861644", 
             "GSM861645", "GSM861646", "GSM861647", "GSM861648", "GSM861649", 
             "GSM861650", "GSM861651", "GSM861652", "GSM861653", "GSM861654", 
             "GSM861655", "GSM861656", "GSM861657", "GSM861658", "GSM861659", 
             "GSM861660", "GSM861661", "GSM861662", "GSM861663", "GSM861664", 
             "GSM861665", "GSM861666", "GSM861667", "GSM861668", "GSM861669", 
             "GSM861670", "GSM861671", "GSM861672", "GSM861673", "GSM861674", 
             "GSM861675", "GSM861676", "GSM861677", "GSM861678", "GSM861679", 
             "GSM861680", "GSM861681", "GSM861682", "GSM861683", "GSM861684", 
             "GSM861685", "GSM861686", "GSM861687", "GSM861688", "GSM861689", 
             "GSM861690", "GSM861691", "GSM861692", "GSM861693", "GSM861694"),
  Column_Prefix = c("SAMPLE 1", "SAMPLE 2", "SAMPLE 3", "SAMPLE 4", "SAMPLE 5", 
                    "SAMPLE 6", "SAMPLE 7", "SAMPLE 8", "SAMPLE 9", "SAMPLE 10", 
                    "SAMPLE 11", "SAMPLE 12", "SAMPLE 13", "SAMPLE 14", "SAMPLE 15", 
                    "SAMPLE 16", "SAMPLE 17", "SAMPLE 18", "SAMPLE 19", "SAMPLE 20", 
                    "SAMPLE 21", "SAMPLE 22", "SAMPLE 23", "SAMPLE 24", "SAMPLE 25", 
                    "SAMPLE 26", "SAMPLE 27", "SAMPLE 28", "SAMPLE 29", "SAMPLE 30", 
                    "SAMPLE 31", "SAMPLE 32", "SAMPLE 33", "SAMPLE 34", "SAMPLE 35", 
                    "SAMPLE 36", "SAMPLE 37", "SAMPLE 38", "SAMPLE 39", "SAMPLE 40", 
                    "SAMPLE 41", "SAMPLE 42", "SAMPLE 43", "SAMPLE 44", "SAMPLE 45", 
                    "SAMPLE 46", "SAMPLE 47", "SAMPLE 48", "SAMPLE 49", "SAMPLE 50", 
                    "SAMPLE 51", "SAMPLE 52", "SAMPLE 53", "SAMPLE 54", "SAMPLE 55", 
                    "SAMPLE 56", "SAMPLE 57", "SAMPLE 58", "SAMPLE 59", "SAMPLE 60")
)

# Print mapping table
print(gsm_mapping)

######################### Replace column names with GSM IDs
# Function to replace column prefixes with corresponding GSM IDs
replace_column_names <- function(colnames, gsm_mapping) {
  new_colnames <- colnames
  # Iterate over each row in gsm_mapping
  for (i in 1:nrow(gsm_mapping)) {
    # Use regex to ensure exact matching of SAMPLE numbers
    new_colnames <- gsub(paste0("\\b", gsm_mapping$Column_Prefix[i], "\\b"), gsm_mapping$GSM_ID[i], new_colnames)
  }
  return(new_colnames)
}

# Replace column names
colnames(GSE35069) <- replace_column_names(colnames(GSE35069), gsm_mapping)

# Print new column names
print(colnames(GSE35069))

GSE35069 <- as.data.frame(GSE35069)
rownames(GSE35069) <- GSE35069$ID_REF


################################### Split into M, U, and P matrices
# Remove ID_REF column
GSE35069$ID_REF <- NULL

# Extract column names for UnMethylated, Methylated, and P values
u_cols <- grep("Unmethylated Signal$", colnames(GSE35069), value = TRUE)
m_cols <- grep("Methylated signal$",  colnames(GSE35069), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE35069), value = TRUE)

# Extract sample IDs (e.g., GSMxxxx)
u_ids <- sub("^(GSM\\d+).*", "\\1", u_cols)
m_ids <- sub("^(GSM\\d+).*", "\\1", m_cols)
p_ids <- sub("^(GSM\\d+).*", "\\1", p_cols)

# Ensure consistent ordering
sample_ids <- sort(unique(u_ids))

# Match columns by sample IDs
u_cols <- u_cols[match(sample_ids, u_ids)]
m_cols <- m_cols[match(sample_ids, m_ids)]
p_cols <- p_cols[match(sample_ids, p_ids)]

# Combine into matrices
U <- as.matrix(GSE35069[, u_cols, drop = FALSE])
M <- as.matrix(GSE35069[, m_cols, drop = FALSE])
P <- as.matrix(GSE35069[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE35069)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids


################################### Calculate Beta values and apply QC using P-values
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

########################## Remove SNP-related probes
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

######################## Remove sex chromosome probes
load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")
gene_annotation <- b
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrX"),]
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrY"),]
beta = beta[ intersect(row.names(beta), gene_annotation$X1),]

###################### Save beta values

# Function to save each column as a file
save_columns_as_files <- function(data, output_dir = getwd()) {
  # Get column names (sample names)
  column_names <- colnames(data)
  
  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each column
  for (col_name in column_names) {
    # Create a data frame with cg probe ID and values
    result <- data.frame(cg_id = rownames(data), attribute_value = data[[col_name]], stringsAsFactors = FALSE)
    
    # Define file name
    file_name <- paste0(col_name, ".txt")
    
    # Create full file path
    file_path <- file.path(output_dir, file_name)
    
    # Write data to file using tab separation
    write.table(result, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Print saved file name (optional)
    cat("File saved as:", file_path, "\n")
  }
}

# Call function to save files
output_dir <- "/data/zhangmch/ewas_array/result/GSE35069"
save_columns_as_files(beta, output_dir)
