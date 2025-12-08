setwd("/data/zhangmch/ewas_array/data/GSE77716")
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

################################################读取作者提供的数据
file_path <- "./GSE77716_Matrix_signal_intensities_updated_2018-Jan08.tsv"

GSE77716 <- read.table(
  file = file_path,
  header = TRUE,          # 第一行是列名
  sep = "\t",             # 制表符分隔
  stringsAsFactors = FALSE,
  check.names = FALSE     # 不改列名，保留 81136_UnMethylated_signal 这种
)

# 看前几行确认一下
head(GSE77716)
str(GSE77716)

########################### 构建映射关系

mapping_data <- read.table("/data/zhangmch/ewas_array/data/GSE77716/mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 提取 GSM_ID 和 Column_Prefix 向量
GSM_ID <- mapping_data$V1
Column_Prefix <- mapping_data$V2


gsm_mapping <- data.frame(GSM_ID, Column_Prefix)

# 打印映射表
print(gsm_mapping)

######################### 替换列名为GSM
# 创建一个函数，用于替换列名前缀为对应的 GSM ID
replace_column_names <- function(colnames, gsm_mapping) {
  new_colnames <- colnames
  # 遍历 gsm_mapping 中的每一行
  for (i in 1:nrow(gsm_mapping)) {
    # 修改 gsub 函数的正则表达式，确保只替换准确的 SAMPLE 数字组合
    new_colnames <- gsub(paste0("\\b", gsm_mapping$Column_Prefix[i], "\\b"), gsm_mapping$GSM_ID[i], new_colnames)
  }
  return(new_colnames)
}

# 替换列名
colnames(GSE77716) <- replace_column_names(colnames(GSE77716), gsm_mapping)


# 打印新的列名
print(colnames(GSE77716))

GSE77716 <- as.data.frame(GSE77716)
rownames(GSE77716) <- GSE77716$ID_REF


################################### 拆分M、U、P三列
# 除去 ID_REF 这一列
GSE77716$ID_REF <- NULL

# 提取 UnMethylated、Methylated 和 Pval 的列名
u_cols <- grep("Unmethylated Signal$", colnames(GSE77716), value = TRUE)
m_cols <- grep("Methylated Signal$",  colnames(GSE77716), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE77716), value = TRUE)

# 提取样本 ID（例如 SAMPLE 1），并去掉后面的部分
u_ids <- sub("^(GSM\\d+).*", "\\1", u_cols)
m_ids <- sub("^(GSM\\d+).*", "\\1", m_cols)
p_ids <- sub("^(GSM\\d+).*", "\\1", p_cols)

# 确保顺序一致（防止列顺序乱）
sample_ids <- sort(unique(u_ids))

# 按照样本 ID 匹配列
u_cols <- u_cols[match(sample_ids, u_ids)]
m_cols <- m_cols[match(sample_ids, m_ids)]
p_cols <- p_cols[match(sample_ids, p_ids)]

# 组合成矩阵
U <- as.matrix(GSE77716[, u_cols, drop = FALSE])
M <- as.matrix(GSE77716[, m_cols, drop = FALSE])
P <- as.matrix(GSE77716[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE77716)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids


################################### 计算Beta值，并用P值做QC
offset <- 100

beta <- M / (M + U + offset)

# 质控：P 值 > 0.01 的设为 NA
beta[P > 0.01] <- NA


beta_na = data.frame(beta)

############################# BMIQ进行一类二类探针校正
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="450K",cores=24)
beta_BMIQ <- myNorm

colnames(beta_BMIQ) = gsub("_.*?$","", colnames(beta_BMIQ)) 

beta_na = data.frame(beta_na)
colnames(beta_na) = gsub("_.*?$","", colnames(beta_na)) 


############################### beta质控 去低质量的样本与位点，并KNN填补
beta = beta_BMIQ
NA_r = rowSums(is.na(beta)) 
NA_c = colSums(is.na(beta))
beta_remain = beta[, which(NA_c <= 485512*0.15) ]

beta_remain = beta_remain[which(NA_r <=  dim(beta_remain)[2]/10),] 
beta_knn = impute.knn(as.matrix(beta_remain))
beta = data.frame(beta_knn$data)

########################## 去除SNP位点
load("/data/zhangmch/ewas_array/script/450k/snp_cg_450K.RData")
beta = beta[ setdiff( row.names(beta), snp_cg$cg_snp),]

######################## 去除性染色体位点
load("/data/zhangmch/ewas_array/script/450k/450K_cg_annotation.RData")
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
output_dir <- "/data/zhangmch/ewas_array/result/GSE77716"
save_columns_as_files(beta, output_dir)
