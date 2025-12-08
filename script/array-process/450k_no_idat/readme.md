# DNA Methylation Data Processing and Quality Control Pipeline for 450K Array Data (GSE125105, GSE127824, GSE88824) without IDAT Files

This project provides a complete **R language pipeline** for processing raw Illumina 450K DNA methylation array data (`.idat` files). The script covers the entire workflow from reading raw data to generating high-quality $\beta$ values, adhering strictly to best practices for bioinformatics data preprocessing and quality control (QC).

## Project Overview

This pipeline is designed to process raw Illumina DNA methylation data in tab-delimited text format, transforming it into a normalized and filtered $\beta$ value matrix suitable for downstream Epigenome-Wide Association Studies (EWAS) or deconvolution analyses. The pipeline performs a series of preprocessing steps, including data reading, format conversion, quality control, BMIQ normalization, missing data imputation, and probe filtering. The final processed data is saved into individual files for each sample.

One of the key steps includes replacing the sample column prefixes with corresponding **GSM IDs** to enhance the clarity and consistency of the dataset. This is achieved using the custom `replace_column_names()` function, which ensures that the column names accurately reflect the sample identifiers.

| Step                      | Purpose                                                                 | Key Technique/Library                  |
| ------------------------- | ----------------------------------------------------------------------- | -------------------------------------- |
| **Data Reading & Format Changing** | Reads raw tab-delimited data files and converts the format for further processing. | `read.table`                          |
| **Column Name Replacement** | Replaces sample column prefixes (e.g., `SAMPLE 1`, `SAMPLE 2`) with corresponding **GSM IDs**. | Custom `replace_column_names()` function |
| **Basic Preprocessing**    | Performs background correction and initial normalization.              | `minfi::preprocessIllumina`           |
| **Quality Control (QC)**   | Removes low-quality probes (based on P-value) and filters samples/probes with excessive missing data. | `minfi::detectionP`, Missing value filtering |
| **Probe Bias Adjustment**  | Applies the **BMIQ** (Beta Mixture Quantile) normalization method to correct for systemic differences between Type I and Type II probes. | `ChAMP::champ.norm`                   |
| **Missing Data Imputation**| Fills remaining missing $\beta$ values using the **k-Nearest Neighbors (KNN)** algorithm. | `impute::impute.knn`                  |
| **Probe Filtering**        | Removes **SNP**-related probes and probes located on sex chromosomes (`chrX/chrY`). | `setdiff`, `intersect`                |
| **Output Generation**      | Saves the final processed $\beta$ values per sample into individual text files. | Custom `R` function                   |

---

### Explanation of the Key Steps:

1. **Data Reading & Format Changing**: 
   - The pipeline reads raw tab-delimited data files (e.g., `GSE35069_Matrix_signal_intensities.txt`) and converts it into a usable format for further processing.

2. **Column Name Replacement**: 
   - The custom `replace_column_names()` function is crucial in this pipeline. It replaces the generic sample column prefixes (e.g., `SAMPLE 1`, `SAMPLE 2`, ...) with **GSM IDs** from a pre-defined mapping (`gsm_mapping`). This ensures that the dataset has clear, informative column names reflecting the actual sample identifiers.

   ```r
   replace_column_names <- function(colnames, gsm_mapping) {
     new_colnames <- colnames
     # Iterate over each GSM mapping row
     for (i in 1:nrow(gsm_mapping)) {
       # Replace the sample prefix with the corresponding GSM ID
       new_colnames <- gsub(paste0("\\b", gsm_mapping$Column_Prefix[i], "\\b"), gsm_mapping$GSM_ID[i], new_colnames)
     }
     return(new_colnames)
   }

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

This pipeline processes raw Illumina DNA methylation data in tab-delimited text format and generates a normalized and filtered $\beta$ value matrix for downstream Epigenome-Wide Association Studies (EWAS) or deconvolution analyses. The following steps outline the process:

### 1. **Set Working Directory**
First, set the working directory to the location of your dataset:

```r
setwd("/data/zhangmch/ewas_array/data/GSE35069")
```

### 2. **Load Required Libraries**
The following R packages are required to run the pipeline:
```r
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

### 3. **Read Raw Data**
The raw data in tab-delimited format (e.g., GSE35069_Matrix_signal_intensities.txt) is read into R:
```r
file_path <- "./GSE35069_Matrix_signal_intensities.txt"
GSE35069 <- read.table(
  file = file_path,
  header = TRUE,          
  sep = "\t",             
  stringsAsFactors = FALSE,
  check.names = FALSE    
)
```

### 4. **Map GSM IDs to Column Prefixes**
The GSM IDs are mapped to the sample column prefixes (e.g., SAMPLE 1, SAMPLE 2, etc.) using a custom mapping:
```r
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
```

### 5. **Replace Column Names with GSM IDs**
The function replace_column_names() replaces the sample column prefixes (e.g., SAMPLE 1) with the corresponding GSM IDs from the mapping:
```r
replace_column_names <- function(colnames, gsm_mapping) {
  new_colnames <- colnames
  for (i in 1:nrow(gsm_mapping)) {
    new_colnames <- gsub(paste0("\\b", gsm_mapping$Column_Prefix[i], "\\b"), gsm_mapping$GSM_ID[i], new_colnames)
  }
  return(new_colnames)
}
# Replace column names in the dataset
colnames(GSE35069) <- replace_column_names(colnames(GSE35069), gsm_mapping)
```

### 6. **Preprocessing: Extract and Normalize Data**
Data is split into Unmethylated, Methylated, and P-value columns. Quality control (QC) steps are applied to remove low-quality probes and perform BMIQ normalization:
```r
# Split the data into U (Unmethylated), M (Methylated), and P (P-value) matrices
u_cols <- grep("Unmethylated Signal$", colnames(GSE35069), value = TRUE)
m_cols <- grep("Methylated signal$",  colnames(GSE35069), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE35069), value = TRUE)
# Process the data for Beta value calculation
beta <- M / (M + U + offset)
beta[P > 0.01] <- NA
```

### 7. **BMIQ Normalization**
The BMIQ normalization method is applied to correct for probe bias between Type I and Type II probes:
```r
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="450K",cores=24)
beta_BMIQ <- myNorm
```

### 8. **Missing Data Imputation**
The pipeline uses k-Nearest Neighbors (KNN) to impute missing $\beta$ values:
```r
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)
```

### 9. **Probe Filtering**
SNP-related probes and probes on sex chromosomes are removed from the dataset:
```r
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]
```

### 10. **Saving the Processed Data**
The final processed data is saved as .txt files for each sample:
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
    cat("文件保存为:", file_path, "\n")
  }
}
output_dir <- "/data/zhangmch/ewas_array/result/GSE35069"
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
