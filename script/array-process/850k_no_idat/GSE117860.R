setwd("/data/zhangmch/ewas_array/data/GSE117860")
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
file_path <- "./GSE117860_MethylatedSignal.txt"

GSE117860 <- read.table(
  file = file_path,
  header = TRUE,          # 第一行是列名
  sep = "\t",             # 制表符分隔
  stringsAsFactors = FALSE,
  check.names = FALSE     # 不改列名，保留 81136_UnMethylated_signal 这种
)

colnames(GSE117860)[1] <- "ID_REF"

# 看前几行确认一下
head(GSE117860)
str(GSE117860)

########################### 构建映射关系
# 读取映射文件
mapping_data <- read.table("/data/zhangmch/ewas_array/data/GSE117860/mapping.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# 提取 GSM_ID 和 Column_Prefix 向量
GSM_ID <- mapping_data$V2
Column_Prefix <- mapping_data$V1

# 创建 gsm_mapping 数据框
gsm_mapping <- data.frame(GSM_ID, Column_Prefix)

# 打印映射表
print(gsm_mapping)

# 检查 GSE166844 的列名，去掉前缀 'X'
colnames(GSE117860) <- gsub("^X", "", colnames(GSE117860))

# 检查列名是否正确
print(colnames(GSE117860))

# 确保替换时有匹配
#for (i in 1:nrow(gsm_mapping)) {
#  # 检查是否有匹配的前缀
#  matching_cols <- grep(paste0("^", gsm_mapping$Column_Prefix[i]), colnames(GSE117860), value = TRUE)
#  if (length(matching_cols) > 0) {
#    print(paste("匹配列：", matching_cols))
#  }
#  
#  # 替换 Column_Prefix 为对应的 GSM_ID
#  colnames(GSE117860) <- gsub(paste0("^", gsm_mapping$Column_Prefix[i]), 
#                              gsm_mapping$GSM_ID[i], 
#                              colnames(GSE117860))
#}

# 避免出现sample 11匹配上了sample 1的情况
for (i in 1:nrow(gsm_mapping)) {
  prefix <- gsm_mapping$Column_Prefix[i]
  gsmi   <- gsm_mapping$GSM_ID[i]
  
  # 前缀后面必须是空白字符或字符串结尾：Sample 1 + 空格 / 结束
  pattern <- paste0("^", prefix, "(\\s|$)")
  
  # 找匹配的列
  matching_cols <- grep(pattern, colnames(GSE117860), value = TRUE, perl = TRUE)
  if (length(matching_cols) > 0) {
    print(paste("匹配列：", paste(matching_cols, collapse = ", ")))
  }
  
  # 只替换前缀部分，保留后面的空格或其他内容
  colnames(GSE117860) <- sub(
    pattern,
    paste0(gsmi, "\\1"),   # \\1 把 (\\s|$) 捕获回去
    colnames(GSE117860),
    perl = TRUE
  )
}

# 打印替换后的列名
print("替换后的列名:")
print(colnames(GSE117860))


GSE117860 <- as.data.frame(GSE117860)
rownames(GSE117860) <- GSE117860$ID_REF


################################### 拆分M、U、P三列
# 除去 ID_REF 这一列
GSE117860$ID_REF <- NULL

# 提取 UnMethylated、Methylated 和 Pval 的列名
u_cols <- grep("Unmethylated signal$", colnames(GSE117860), value = TRUE)
m_cols <- grep("Methylated signal$",  colnames(GSE117860), value = TRUE)
p_cols <- grep("Detection Pval$",     colnames(GSE117860), value = TRUE)

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
U <- as.matrix(GSE117860[, u_cols, drop = FALSE])
M <- as.matrix(GSE117860[, m_cols, drop = FALSE])
P <- as.matrix(GSE117860[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(GSE117860)
colnames(U) <- colnames(M) <- colnames(P) <- sample_ids


################################### 计算Beta值，并用P值做QC
offset <- 100

beta <- M / (M + U + offset)

# 质控：P 值 > 0.01 的设为 NA
beta[P > 0.01] <- NA


beta_na = data.frame(beta)

############################# BMIQ进行一类二类探针校正
beta_na_clean <- na.omit(beta_na)
myNorm <- champ.norm(beta=beta_na_clean,arraytype="EPIC",cores=24)
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
output_dir <- "/data/zhangmch/ewas_array/result/GSE117860"
save_columns_as_files(beta, output_dir)
