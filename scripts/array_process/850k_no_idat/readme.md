# DNA Methylation Data Processing and Quality Control Pipeline for EPIC 850K Array Data (GSE117860, GSE166844) without IDAT Files

This project provides a complete **R language pipeline** for processing raw Illumina EPIC 850K DNA methylation array data without IDAT files. The scripts cover the entire workflow from reading raw signal intensity data to generating high-quality $\beta$ values, adhering strictly to best practices for bioinformatics data preprocessing and quality control (QC).

## Project Overview

This pipeline is designed to process raw Illumina EPIC 850K DNA methylation data in tab-delimited or comma-separated text format, transforming it into a normalized and filtered $\beta$ value matrix suitable for downstream Epigenome-Wide Association Studies (EWAS) or deconvolution analyses. The pipeline performs a series of preprocessing steps, including data reading, format conversion, GSM ID mapping, quality control, BMIQ normalization, missing data imputation, and probe filtering. The final processed data is saved into individual files for each sample.

A critical feature of this pipeline is the use of an **external mapping file** (`mapping.txt`) to associate raw column prefixes with **GSM IDs**, rather than hardcoding sample mappings within the script. For GSE117860, an enhanced prefix replacement strategy employing **boundary detection** is used to prevent partial prefix matching (e.g., preventing `sample 1` from matching `sample 11`).

| Step                      | Purpose                                                                 | Key Technique/Library                  |
| ------------------------- | ----------------------------------------------------------------------- | -------------------------------------- |
| **Data Reading & Format Changing** | Reads raw signal intensity files (tab-delimited `.txt` or `.csv`) and standardizes the format for further processing. | `read.table`                          |
| **Column Name Mapping**   | Loads GSM ID–to–column prefix mappings from an external `mapping.txt` file and replaces raw column prefixes with GSM IDs. | External `mapping.txt`, regex replacement |
| **Boundary-Safe Replacement** (GSE117860) | Replaces column prefixes using word-boundary pattern matching to avoid partial prefix collisions (e.g., `sample 1` vs. `sample 11`). | `sub()` with `perl = TRUE`, `(\\s|$)` anchor |
| **Matrix Construction**   | Splits the combined signal table into separate Unmethylated (U), Methylated (M), and Detection P-value (P) matrices. | `grep`, `as.matrix`                   |
| **Basic Preprocessing**   | Calculates raw $\beta$ values using `M / (M + U + 100)` and applies P-value–based QC to flag unreliable measurements as `NA`. | Beta formula, `detectionP` threshold  |
| **Probe Bias Adjustment** | Applies the **BMIQ** (Beta Mixture Quantile) normalization method to correct for systematic differences between Type I and Type II probes, using EPIC array type. | `ChAMP::champ.norm(arraytype="EPIC")` |
| **Missing Data Imputation** | Removes samples and probes with excessive missing values (>15% and >10%, respectively), then imputes remaining `NA` values using the **k-Nearest Neighbors (KNN)** algorithm. | `impute::impute.knn`                  |
| **Probe Filtering**       | Removes **SNP**-related probes (using 850K-specific SNP list) and probes located on sex chromosomes (`chrX/chrY`) using EPIC annotation. | `snp_cg_850K.RData`, `850K_cg_annotation.RData` |
| **Output Generation**     | Saves the final processed $\beta$ values per sample into individual text files. | Custom `R` function                   |

---

### Explanation of the Key Steps:

1. **Data Reading & Format Changing**:
   - GSE117860 reads a tab-delimited signal intensity file (`GSE117860_MethylatedSignal.txt`), and GSE166844 reads a comma-separated file (`GSE166844_Variance_raw_Signal.csv`). Both are read with `check.names = FALSE` to preserve original column names.

   ```r
   # GSE117860 (tab-delimited)
   GSE117860 <- read.table(file = "./GSE117860_MethylatedSignal.txt",
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE)

   # GSE166844 (comma-separated)
   GSE166844 <- read.table(file = "./GSE166844_Variance_raw_Signal.csv",
                           header = TRUE, sep = ",",
                           stringsAsFactors = FALSE, check.names = FALSE)
   ```

2. **External Mapping File**:
   - GSM ID mappings are loaded from an external `mapping.txt` file rather than hardcoded in the script. This makes the pipeline flexible and reusable across datasets.

   ```r
   mapping_data <- read.table("./mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
   GSM_ID         <- mapping_data$V2   # Column 2: GSM ID
   Column_Prefix  <- mapping_data$V1   # Column 1: raw column prefix
   gsm_mapping    <- data.frame(GSM_ID, Column_Prefix)
   ```

3. **Boundary-Safe Column Replacement** (GSE117860):
   - To avoid partial prefix matching (e.g., `sample 1` matching `sample 10` or `sample 11`), GSE117860 uses a regular expression with a trailing `(\\s|$)` anchor, ensuring that only exact prefix matches followed by whitespace or end-of-string are replaced.

   ```r
   for (i in 1:nrow(gsm_mapping)) {
     prefix  <- gsm_mapping$Column_Prefix[i]
     gsmi    <- gsm_mapping$GSM_ID[i]
     pattern <- paste0("^", prefix, "(\\s|$)")   # boundary-safe pattern
     colnames(GSE117860) <- sub(pattern,
                                paste0(gsmi, "\\1"),
                                colnames(GSE117860),
                                perl = TRUE)
   }
   ```

4. **Matrix Construction**:
   - Column names are parsed to extract Unmethylated, Methylated, and Detection P-value columns. Sample ordering is synchronized across all three matrices before combining.

   ```r
   u_cols     <- grep("Unmethylated signal$", colnames(GSE117860), value = TRUE)
   m_cols     <- grep("Methylated signal$",   colnames(GSE117860), value = TRUE)
   p_cols     <- grep("Detection Pval$",      colnames(GSE117860), value = TRUE)
   sample_ids <- sort(unique(sub("^(GSM\\d+).*", "\\1", u_cols)))
   U <- as.matrix(GSE117860[, u_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", u_cols))]])
   M <- as.matrix(GSE117860[, m_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", m_cols))]])
   P <- as.matrix(GSE117860[, p_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", p_cols))]])
   colnames(U) <- colnames(M) <- colnames(P) <- sample_ids
   ```

## Installation Dependencies

To run these scripts, the following R packages are required. It is highly recommended to install them within a **Conda** or **renv** environment to avoid dependency conflicts.

```R
# Run the following commands in the R console:
install.packages(c("GEOquery", "gmqn", "minfi", "gplots", "limma", "WGCNA", "impute", "wateRmelon", "ChAMP", "R.utils"))

# Note: Some packages must be installed from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("minfi", "impute", "IlluminaHumanMethylationEPICmanifest"))
```

## Usage Guide

This pipeline processes raw Illumina EPIC 850K DNA methylation signal data and generates a normalized and filtered $\beta$ value matrix for downstream EWAS or deconvolution analyses. The following steps outline the process.

### 1. **Set Working Directory**

```r
setwd("/data/zhangmch/ewas_array/data/GSE117860")
# or
setwd("/data/zhangmch/ewas_array/data/GSE166844")
```

### 2. **Load Required Libraries**

```r
library(GEOquery)
library(gmqn)
library(minfi)
library(gplots)
library(limma)
library(WGCNA)
library(impute)
library(IlluminaHumanMethylationEPICmanifest)
library(wateRmelon)
library(ChAMP)
library(R.utils)
```

### 3. **Prepare the Mapping File**

Place a tab-separated `mapping.txt` file (no header) in the dataset directory. Each row maps a raw column prefix to its corresponding GSM ID. Note that the column order differs between the two datasets:

```
# GSE117860: Column_Prefix <TAB> GSM_ID
# GSE166844: GSM_ID <TAB> Column_Prefix
```

### 4. **Read Raw Signal Data**

```r
# GSE117860 (tab-delimited)
GSE117860 <- read.table(file = "./GSE117860_MethylatedSignal.txt",
                        header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, check.names = FALSE)
colnames(GSE117860)[1] <- "ID_REF"

# GSE166844 (comma-separated)
GSE166844 <- read.table(file = "./GSE166844_Variance_raw_Signal.csv",
                        header = TRUE, sep = ",",
                        stringsAsFactors = FALSE, check.names = FALSE)
colnames(GSE166844)[1] <- "ID_REF"
```

### 5. **Load Mapping and Replace Column Names with GSM IDs**

```r
mapping_data <- read.table("./mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gsm_mapping  <- data.frame(GSM_ID = mapping_data$V2, Column_Prefix = mapping_data$V1)

# Remove leading 'X' that R may prepend to numeric column names
colnames(GSE117860) <- gsub("^X", "", colnames(GSE117860))

# Boundary-safe replacement (GSE117860)
for (i in 1:nrow(gsm_mapping)) {
  pattern <- paste0("^", gsm_mapping$Column_Prefix[i], "(\\s|$)")
  colnames(GSE117860) <- sub(pattern, paste0(gsm_mapping$GSM_ID[i], "\\1"),
                             colnames(GSE117860), perl = TRUE)
}
```

### 6. **Construct U, M, P Matrices**

```r
u_cols <- grep("Unmethylated signal$", colnames(GSE117860), value = TRUE)
m_cols <- grep("Methylated signal$",   colnames(GSE117860), value = TRUE)
p_cols <- grep("Detection Pval$",      colnames(GSE117860), value = TRUE)

sample_ids <- sort(unique(sub("^(GSM\\d+).*", "\\1", u_cols)))
u_cols <- u_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", u_cols))]
m_cols <- m_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", m_cols))]
p_cols <- p_cols[match(sample_ids, sub("^(GSM\\d+).*", "\\1", p_cols))]

U <- as.matrix(GSE117860[, u_cols])
M <- as.matrix(GSE117860[, m_cols])
P <- as.matrix(GSE117860[, p_cols])
rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE117860)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids
```

### 7. **Calculate Beta Values and Apply P-value QC**

```r
offset <- 100
beta         <- M / (M + U + offset)
beta[P > 0.01] <- NA
beta_na      <- data.frame(beta)
```

### 8. **BMIQ Normalization (EPIC)**

```r
beta_na_clean <- na.omit(beta_na)
myNorm        <- champ.norm(beta = beta_na_clean, arraytype = "EPIC", cores = 24)
beta_BMIQ     <- myNorm
colnames(beta_BMIQ) <- gsub("_.*?$", "", colnames(beta_BMIQ))
```

### 9. **Missing Data Filtering and KNN Imputation**

```r
beta        <- beta_BMIQ
beta_remain <- beta[, which(colSums(is.na(beta)) <= 865859 * 0.15)]
beta_remain <- beta_remain[which(rowSums(is.na(beta_remain)) <= dim(beta_remain)[2] / 10), ]
beta_knn    <- impute.knn(as.matrix(beta_remain))
beta        <- data.frame(beta_knn$data)
```

### 10. **Probe Filtering**

```r
# Remove SNP probes (850K-specific list)
load("/data/zhangmch/ewas_array/script/850k/snp_cg_850K.RData")
beta <- beta[setdiff(row.names(beta), snp_cg$cg_snp), ]

# Remove sex chromosome probes (EPIC 850K annotation)
load("/data/zhangmch/ewas_array/script/850k/850K_cg_annotation.RData")
gene_annotation <- b
gene_annotation <- gene_annotation[gene_annotation$X3 != "chrX", ]
gene_annotation <- gene_annotation[gene_annotation$X3 != "chrY", ]
beta <- beta[intersect(row.names(beta), gene_annotation$X1), ]
```

### 11. **Saving the Processed Data**

```r
save_columns_as_files <- function(data, output_dir = getwd()) {
  column_names <- colnames(data)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  for (col_name in column_names) {
    result    <- data.frame(cg_id = rownames(data), attribute_value = data[[col_name]], stringsAsFactors = FALSE)
    file_path <- file.path(output_dir, paste0(col_name, ".txt"))
    write.table(result, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("File saved as:", file_path, "\n")
  }
}

output_dir <- "/data/zhangmch/ewas_array/result/GSE117860"
save_columns_as_files(beta, output_dir)
```

## Output Files

A separate `.txt` file will be generated for each sample in the specified `output_dir` directory.

### File Naming Format:
`[Sample_ID].txt` (e.g., `GSM123456.txt`)

### File Content:
Each file contains two tab-separated columns (with no headers): the CpG probe ID and the sample's processed $\beta$ value.

| cg_id      | attribute_value |
| ---------- | --------------- |
| cg00000029 | 0.845           |
| cg00000109 | 0.122           |
| ...        | ...             |
