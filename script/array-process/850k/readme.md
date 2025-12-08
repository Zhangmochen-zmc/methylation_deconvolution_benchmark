# DNA Methylation Data Processing and Quality Control Pipeline for 850K (EPIC) Array Data (GSE123914, GSE112618, GSE110554, GSE184269) with IDAT Files

This project provides a complete **R language pipeline** for processing raw Illumina 450K DNA methylation array data (`.idat` files). The script covers the entire workflow from reading raw data to generating high-quality $\beta$ values, adhering strictly to best practices for bioinformatics data preprocessing and quality control (QC).

## Project Overview

This pipeline is designed to transform raw Illumina IDAT files into a normalized and filtered $\beta$ value matrix, suitable for downstream Epigenome-Wide Association Studies (EWAS) or deconvolution analyses.

| Step | Purpose | Key Technique/Library |
| :--- | :--- | :--- |
| **Data Reading** | Reads and unzips raw `.idat` files. | `R.utils`, `minfi` |
| **Basic Preprocessing** | Performs background correction and initial normalization. | `minfi::preprocessIllumina` |
| **Quality Control (QC)** | Removes low-quality probes (based on P-value) and filters samples/probes with excessive missing data. | `minfi::detectionP`, Missing value filtering |
| **Probe Bias Adjustment** | Applies the **BMIQ** (Beta Mixture Quantile) normalization method to correct for systemic differences between Type I and Type II probes. | `ChAMP::champ.norm` |
| **Missing Data Imputation** | Fills remaining missing $\beta$ values using the **k-Nearest Neighbors (KNN)** algorithm. | `impute::impute.knn` |
| **Probe Filtering** | Removes **SNP**-related probes and probes located on sex chromosomes (`chrX/chrY`). | `setdiff`, `intersect` |
| **Output Generation** | Saves the final processed $\beta$ values per sample into individual text files. | Custom `R` function |

## Installation Dependencies

To run this script, the following R packages are required. It is highly recommended to install them within a **Conda** or **renv** environment to avoid dependency conflicts.

```R
# Run the following commands in the R console:
install.packages(c("GEOquery", "gmqn", "minfi", "gplots", "limma", "WGCNA", "impute", "wateRmelon", "ChAMP", "R.utils"))

# Note: Some dependency packages (like minfi or manifest files) must be installed from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("minfi", "impute", "IlluminaHumanMethylation450kmanifest"))
```

## Usage Guide

### Running the Pipeline:

1. **Prepare Data**: Place all raw IDAT files (often initially in `.idat.gz` format) inside the subdirectory specified in the script (e.g., `GSE125105`).
2. **Set Paths**: Modify the input (`setwd`) and output (`output_dir`) paths in Sections 1 and 8 of the code to match your environment.
3. **Execute**: Run the complete R script line by line or as a batch job.

### Critical Paths in the Script (Please Modify to Your Local Environment):

```r
# 1. Set Working Directory
setwd("/data/zhangmch/ewas_array/data")

# 6. Path to SNP Probe List file
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")

# 7. Path to Probe Annotation file
load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")

# 8. Output Directory
output_dir <- "/data/zhangmch/ewas_array/result/GSE125105"
```

## Code Explanation

### 1. Set Working Directory and Load Libraries
```r
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
```

### 2. Read and Initial Preprocessing
```r
gz_files <- list.files(path = "./GSE125105", pattern = "\\.idat\\.gz$", full.names = TRUE)
for (file in gz_files) {
  gunzip(file, remove = FALSE) # Unzip and retain original files
}
rgset = read.metharray.exp("./GSE125105")
mset = preprocessIllumina(rgset)
beta_na = getBeta(mset, offset = 100) # Extract beta values with offset
beta_na = round(beta_na, 3)
```

### 3. Remove Low-Quality Probes (Detection P-Value QC)
```r
pvalue = detectionP(rgset) 
for (i in 1:dim(pvalue)[2]){
  # Set beta values to NA for probes with P-value > 0.01
  beta_na[which(pvalue[,i] > 0.01),i] = NA
}
beta_na = data.frame(beta_na)
```

### 4. BMIQ Normalization (Bias Adjustment)
```r
beta_na_clean <- na.omit(beta_na) # Temporarily remove NA rows for champ.norm input
myNorm <- champ.norm(beta=beta_na_clean, arraytype="450K", cores=24)
beta_BMIQ <- myNorm
```

### 5. Missing Data Filtering and KNN Imputation
```r
beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) # Clean up sample names
beta = beta_BMIQ
NA_r = rowSums(is.na(beta))
NA_c = colSums(is.na(beta))

# Filter samples with > 15% missing data
beta_remain = beta[, which(NA_c <= 485512*0.15)] 
# Filter probes with > 10% missing data
beta_remain = beta_remain[which(NA_r <= dim(beta_remain)[2]/10),] 

# KNN imputation for remaining missing values
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)
```

### 6. Remove SNP Probes
```
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]
```

### 7. Remove Sex Chromosome Probes (chrX/chrY)
```
load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")
gene_annotation <- b
# Filter out X chromosome probes
gene_annotation = gene_annotation[which(gene_annotation$X3 != "chrX"),]
# Filter out Y chromosome probes
gene_annotation = gene_annotation[which(gene_annotation$X3 != "chrY"),]
# Keep the intersection of existing probes and filtered annotation
beta = beta[intersect(row.names(beta), gene_annotation$X1),]
```

### 8. Save $\beta$ Values to Files
```r
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
```

## Output Files

A separate `.txt` file will be generated for each sample in the specified `output_dir` directory.

### File Naming Format:
`[Sample_ID].txt` (e.g., `GSM123456.txt`)

### File Content:
Each file contains two tab-separated columns (with no headers): the CpG ID and the sample's processed $\beta$ value.

| cg_id      | attribute_value |
| ---------- | --------------- |
| cg00000029 | 0.845           |
| cg00000109 | 0.122           |
| ...        | ...             |
