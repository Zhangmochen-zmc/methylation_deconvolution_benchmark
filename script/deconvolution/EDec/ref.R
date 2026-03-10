library(EDec)

# screening sites
meth_data <- read.table("ref/refdata.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
meth_data <- as.matrix(meth_data)
meta_data <- read.csv("ref/refmeta.csv")

common_samples <- intersect(colnames(meth_data), meta_data$FileID)
meth_data <- meth_data[, common_samples] 

# reorder meta_data to match the column order of meth_data
meta_data <- meta_data[match(common_samples, meta_data$FileID), ]

# extracting cell type vectors
my_classes <- meta_data$Tissue

if(ncol(meth_data) != length(my_classes)) {
  stop("error：sample numbers not match！")
}

# markers
markers <- run_edec_stage_0(
  reference_meth = as.matrix(meth_data), 
  reference_classes = my_classes,
  max_p_value = 1e-5,      # 0.00001
  num_markers = 500,       
  version = "one.vs.rest"  
)

# save
# print(head(markers))

saveRDS(markers, file = "edec_stage0_markers.rds")


