# DNA Methylation Data Processing and Quality Control

## Project Overview
This script is used for processing and analyzing Illumina 450K DNA methylation data, covering the following major steps:

- Unzipping `.idat.gz` files and reading the data
- Preprocessing, including background correction and normalization
- Removing low-quality probes and samples
- BMIQ normalization for probe type bias adjustment
- Imputing missing data using KNN
- Removing SNP and sex chromosome probes
- Saving the final methylation beta values into text files

## Install Dependencies

```bash
install.packages(c("GEOquery", "gmqn", "minfi", "gplots", "limma", "WGCNA", "impute", "IlluminaHumanMethylationEPICmanifest", "wateRmelon", "ChAMP"))

Usage
1-Set the working directory and load the required libraries.
2-Unzip .idat.gz files and read the methylation data.
3-Perform data preprocessing, including removing low-quality probes and samples.
4-Apply BMIQ normalization to adjust for probe type differences.
5-Impute missing methylation data using the KNN method.
6-Remove SNP and sex chromosome probes.
7-Save the processed data to individual text files.

## Output

Each sample's methylation beta values will be saved as a separate .txt file.

## Code Explanation: DNA Methylation Array Data (450k) Processing and Quality Control

1. Set Working Directory and Load Required Libraries

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

1) setwd: Sets the working directory to the folder containing the data files.
2) library: Loads the required R packages for methylation data processing, quality control, statistical analysis, and visualization.

2. Unzip .idat.gz Files and Read Data

gz_files <- list.files(path = "./GSE125105", pattern = "\\.idat\\.gz$", full.names = TRUE)
for (file in gz_files) {
  gunzip(file, remove = FALSE)  # Unzip and keep the original .gz files
}
rgset =  read.metharray.exp("./GSE125105")
mset = preprocessIllumina(rgset)
beta_na = getBeta(mset, offset = 100)  # Extract methylation beta values, adding offset of 100
beta_na = round(beta_na, 3)

1) list.files: Lists all .idat.gz files in the specified directory.
2) gunzip: Unzips the .idat.gz files and retains the original compressed files.
3) read.metharray.exp: Reads the Illumina methylation array data from the .idat files.
4) preprocessIllumina: Preprocesses the data, including background correction and normalization.
5) getBeta: Extracts methylation beta values, with an offset of 100 to avoid division by zero errors.

3. Remove Low-Quality Probes

pvalue = detectionP(rgset)              
for (i in 1:dim(pvalue)[2]){                      
  beta_na[which(pvalue[,i] > 0.01),i] = NA
}
beta_na = data.frame(beta_na)

1) detectionP: Calculates detection P-values for each probe, indicating whether a probe was successfully detected in the sample.
2) for loop: Iterates over each sample, and for probes with a P-value greater than 0.01 (indicating poor detection), their values are set to NA.

4. BMIQ Normalization

beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean, arraytype="450K", cores=24)
beta_BMIQ <- myNorm

1) na.omit: Removes rows containing NA values.
2) champ.norm: Normalizes the data using the BMIQ method to adjust for type 1 and type 2 probe biases. The method is specifically for 450K arrays.
3)cores=24: Specifies using 24 CPU cores for parallel computation, optimizing processing time.

5. Process Beta Values and Impute Missing Data

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 
beta = beta_BMIQ
NA_r = rowSums(is.na(beta)) 
NA_c = colSums(is.na(beta))
beta_remain = beta[, which(NA_c <= 485512*0.15) ]
beta_remain = beta_remain[which(NA_r <=  dim(beta_remain)[2]/10),] 
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)


1) gsub: Removes the extra parts of the sample names, keeping only the main part.
2)rowSums and colSums: Count the number of missing values in each row and column.
3) beta_remain: Filters out samples with more than 15% missing data and probes with more than 10% missing data.
4) impute.knn: Uses the k-Nearest Neighbors (KNN) method to impute missing beta values.

6. Remove SNP Probes

load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

1) load: Loads a file containing the SNP probes.
2) setdiff: Removes SNP probes from the beta data, as these probes may interfere with methylation analysis.

7. Remove Sex Chromosome Probes

load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")
gene_annotation <- b
gene_annotation = gene_annotation[which(gene_annotation$X3 != "chrX"),]
gene_annotation = gene_annotation[which(gene_annotation$X3 != "chrY"),]
beta = beta[intersect(row.names(beta), gene_annotation$X1),]

1) load: Loads the file containing the gene annotations.
2) which: Filters out probes located on the X and Y chromosomes.
3) intersect: Keeps only the probes that are not on sex chromosomes (X and Y).

8. Save Beta Values to Files

save_columns_as_files <- function(data, output_dir = getwd()) {
  column_names <- colnames(data)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  for (col_name in column_names) {
    result <- data.frame(cg_id = rownames(data), attribute_value = data[[col_name]], stringsAsFactors = FALSE)
    file_name <- paste0(col_name, ".txt")
    file_path <- file.path(output_dir, file_name)
    write.table(result, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("File saved as:", file_path, "\n")
  }
}
output_dir <- "/data/zhangmch/ewas_array/result/GSE125105"
save_columns_as_files(beta, output_dir)

1) save_columns_as_files: A custom function to save each sample's beta values into individual .txt files.
2) column_names: Extracts the column names (i.e., sample names).
3) write.table: Writes each sample's methylation data (CG probe IDs and attribute values) into a .txt file, using tab as the separator.
4) output_dir: Specifies the directory to save the output files.

