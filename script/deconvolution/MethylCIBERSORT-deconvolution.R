library(devtools)
#devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(pheatmap)
library(tibble)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggthemes)

#load refdata 
sig_matrix <- read.table(
  "450k_ref_Signature.txt",
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

#refdata-na
sig_no_na <- na.omit(sig_matrix)
cat("Signature:", nrow(sig_matrix), "\n")
cat("delete:", nrow(sig_no_na), "\n")
cat("invaild", nrow(sig_matrix) - nrow(sig_no_na), "NA。\n")

#temporary filr
write.table(sig_no_na, "temp_sig.txt", sep = "\t", quote = FALSE, col.names = NA)

#load_real
exp <- read.table("beta_values.tsv",
                  sep = "\t",
                  header = TRUE,
                  row.names = 1,
                  fill = TRUE,
                  check.names = FALSE)
#real-0 
exp_clean <- exp
exp_clean[is.na(exp_clean)] <- 0

res <- cibersort("temp_sig.txt", exp_clean, perm = 1000, QN = F)
write.table(res, "res.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

#QN array:T，sequence:F
