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

################################################  idat文件转换为beta，预处理

# 解压缩 .idat.gz 文件
gz_files <- list.files(path = "./GSE184269", pattern = "\\.idat\\.gz$", full.names = TRUE)

# 解压所有的 .idat.gz 文件
for (file in gz_files) {
  gunzip(file, remove = FALSE)  # 解压并保留原始 .gz 文件
}

rgset =  read.metharray.exp("./GSE184269")
mset = preprocessIllumina(rgset)
beta_na = getBeta(mset, offset = 100)  #提取甲基化比例，公式下面+100
beta_na = round(beta_na,3)
### 去低质量
pvalue = detectionP(rgset)              
for (i in 1:dim(pvalue)[2]){                      
  beta_na[which(pvalue[,i] > 0.01),i] = NA
}

beta_na = data.frame(beta_na)

############################# BMIQ进行一类二类探针校正
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="EPIC",cores=12)
beta_BMIQ <- myNorm

colnames(beta_BMIQ) = gsub("_.*?$","", colnames(beta_BMIQ)) 

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 


############################### beta质控 去低质量的样本与位点，并KNN填补
beta = beta_BMIQ
NA_r = rowSums(is.na(beta)) 
NA_c = colSums(is.na(beta))
beta_remain = beta[, which(NA_c <= 865859*0.15) ]
beta_remain = beta_remain[which(NA_r <=  dim(beta_remain)[2]/10),] 
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)

########################## 去除SNP位点
load("/data/zhangmch/ewas_array/script/850k/snp_cg_850K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

######################## 去除性染色体位点
load("/data/zhangmch/ewas_array/script/850k/850K_cg_annotation.RData")
gene_annotation <- b
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrX"),]
gene_annotation= gene_annotation[which( gene_annotation$X3  != "chrY"),]
beta = beta[ intersect(row.names(beta), gene_annotation$X1),]

###################### 保存beta值

# 定义保存数据的函数
save_columns_as_files <- function(data, output_dir = getwd()) {
  # 获取列名（即样本名）
  column_names <- colnames(data)
  
  # 如果指定的输出目录不存在，则创建它
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 循环遍历每一列
  for (col_name in column_names) {
    # 创建一个新的数据框，将该列的cg探针ID和属性值保存到一起
    result <- data.frame(cg_id = rownames(data), attribute_value = data[[col_name]], stringsAsFactors = FALSE)
    
    # 定义文件名，保存为 col_name.txt
    file_name <- paste0(col_name, ".txt")
    
    # 创建完整的文件路径
    file_path <- file.path(output_dir, file_name)
    
    # 将数据写入文件，使用 tab 分隔
    write.table(result, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # 打印保存的文件名（可选）
    cat("文件保存为:", file_path, "\n")
  }
}

# 调用函数进行保存
output_dir <- "/data/zhangmch/ewas_array/result/GSE184269"
save_columns_as_files(beta, output_dir)
