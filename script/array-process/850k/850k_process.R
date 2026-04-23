setwd("/data/zhangmch/ewas_array/data")
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

################################################  Convert idat files to beta and preprocess

# Decompress .idat.gz files
gz_files <- list.files(path = "./GSE184269", pattern = "\\.idat\\.gz$", full.names = TRUE)

# Decompress all .idat.gz files
for (file in gz_files) {
  gunzip(file, remove = FALSE)  # Decompress and keep original .gz files
}

rgset =  read.metharray.exp("./GSE184269")
mset = preprocessIllumina(rgset)
beta_na = getBeta(mset, offset = 100)  # Extract methylation proportion, formula adds +100
beta_na = round(beta_na,3)

### Remove low-quality data
pvalue = detectionP(rgset)              
for (i in 1:dim(pvalue)[2]){                      
  beta_na[which(pvalue[,i] > 0.01),i] = NA
}

beta_na = data.frame(beta_na)

############################# Perform BMIQ normalization for Type I and Type II probes
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="EPIC",cores=12)
beta_BMIQ <- myNorm

colnames(beta_BMIQ) = gsub("_.*?$","", colnames(beta_BMIQ)) 

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 


############################### Beta quality control: remove low-quality samples and probes, then KNN imputation
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
output_dir <- "/data/zhangmch/ewas_array/result/GSE184269"
save_columns_as_files(beta, output_dir)
