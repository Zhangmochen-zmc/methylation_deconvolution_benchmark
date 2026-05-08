if (!requireNamespace("this.path", quietly = TRUE)) {
  install.packages("this.path")
}
library(this.path)
library(data.table)
library(tools) 

script_dir <- this.dir()                
parent_dir <- dirname(script_dir)    

# Path configuration
root_dir <- "ref_data"
output_file <- file.path(parent_dir, "ref_data.RData")

# Get subfolder names for all cell types
cell_types <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
print(paste("Detected cell types:", paste(cell_types, collapse = ", ")))

# -------------------------------------------------------
# 2. Batch data reading
# -------------------------------------------------------
data_list <- list()   # To store Beta value vectors for each sample
sample_names <- c()   # To store sample names (e.g., GSM IDs)
pheno_labels <- c()   # To store group labels (Cell Types)

# To record probe IDs common to all samples (for intersection)
common_probes <- NULL 

cat("Start reading files...\n")

for (ct in cell_types) {
  # Construct path for the current cell type
  ct_path <- file.path(root_dir, ct)
  
  # Get all .txt files in the current folder
  files <- list.files(ct_path, pattern = "\\.txt$", full.names = TRUE)
  
  cat(paste("Processing", ct, ": found", length(files), "files.\n"))
  
  for (f in files) {
    # Extract file name as sample ID (remove .txt extension)
    s_name <- file_path_sans_ext(basename(f))
    
    # Rapidly read file with fread (assuming 1st column is probe ID, 2nd is Beta value)
    # Note: Adjust header=TRUE/FALSE based on your specific file structure
    tmp <- fread(f, header = TRUE, sep = "\t")
    
    # Assuming inconsistent column names, force 1st column as ID and 2nd as Value
    # If the txt file contains only values without IDs, modify this logic
    probes <- as.character(tmp[[1]])
    betas <- as.numeric(tmp[[2]])
    
    # Handle common probes (initialize on first loop, then intersect)
    if (is.null(common_probes)) {
      common_probes <- probes
    } else {
      common_probes <- intersect(common_probes, probes)
    }
    
    # Temporarily store in list (avoid merging here to save memory)
    # Stored as a named vector for easier subsequent alignment
    names(betas) <- probes
    
    data_list[[s_name]] <- betas
    sample_names <- c(sample_names, s_name)
    pheno_labels <- c(pheno_labels, ct)
  }
}

cat("\nFiles read finished. Merging data...\n")
cat(paste("Number of common probes:", length(common_probes), "\n"))

# -------------------------------------------------------
# 3. Data alignment and merging
# -------------------------------------------------------
# Create an empty matrix
ref <- matrix(NA, nrow = length(common_probes), ncol = length(sample_names))
rownames(ref) <- common_probes
colnames(ref) <- sample_names

# Loop to fill the matrix, ensuring strict row alignment
for (i in 1:length(sample_names)) {
  s_name <- sample_names[i]
  # Extract data for the sample, ordered by common_probes
  ref[, i] <- data_list[[s_name]][common_probes]
  
  if(i %% 10 == 0) cat(".") # Progress indicator
}
cat("\nMatrix construction done.\n")

# -------------------------------------------------------
# 4. Construct Phenotype and Save
# -------------------------------------------------------
# Construct ref.pheno (must be a vector)
ref.pheno <- pheno_labels
names(ref.pheno) <- sample_names

# Check dimensions (columns should match the total number of txt files)
print(dim(ref))
print(table(ref.pheno))

# Save as RData
# Note: Keeping variable names as 'ref' and 'ref.pheno' for downstream compatibility
save(ref, ref.pheno, file = output_file)

cat(paste("Data saved to:", output_file, "\n"))
