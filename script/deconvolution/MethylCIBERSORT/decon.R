library(tools)           
library("MethylCIBERSORT")

beta_file <- "/data/yuxy/data/data/wgbs_850k/data/real.csv" 
ref_file  <- "/data/yuxy/data/data/wgbs_850k/wgbs_850k_ref_raw.RData"
csv_name <- file_path_sans_ext(basename(beta_file))
output_dir <- "ref"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
project <- file.path(output_dir, csv_name)

print(paste("Processing Project:", csv_name))
print(paste("Output will be saved to prefix:", project))


Mat <- read.table(file=beta_file, sep = ',', header = TRUE, row.names=1, check.names=FALSE)
load(ref_file)  
RefData <- ref
RefPheno <- ref.pheno

RefData <- RefData[sort(rownames(RefData)), ]
Mat <- Mat[sort(rownames(Mat)), ]
Int <- intersect(rownames(Mat), rownames(RefData))
Mat <- Mat[match(Int, rownames(Mat)), ]
RefData <- RefData[match(Int, rownames(RefData)), ]

Signature <- FeatureSelect.V4(
  CellLines.matrix = NULL,
  Heatmap = TRUE,
  export = TRUE,
  sigName = paste0(csv_name, "_ref"), 
  Stroma.matrix = RefData,
  deltaBeta = 0.2,
  FDR = 0.01,
  MaxDMRs = 1000,
  Phenotype.stroma = RefPheno
)

save(Mat, RefData, RefPheno, Signature, file = paste0(project, "_processed.RData"))
