setwd("/data/zhangmch/ewas_array/data/GSE58888")
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
file_path <- "./GSE58888_methylated_unmethylated_signal_intensities.txt"

gse58888 <- read.table(
  file = file_path,
  header = TRUE,          # 第一行是列名
  sep = "\t",             # 制表符分隔
  stringsAsFactors = FALSE,
  check.names = FALSE     # 不改列名，保留 81136_UnMethylated_signal 这种
)

# 看前几行确认一下
head(gse58888)
str(gse58888)

########################### 构建映射关系
gsm_mapping <- data.frame(
  GSM_ID = c(
    "GSM1421611", "GSM1421612", "GSM1421613", "GSM1421614", "GSM1421615",
    "GSM1421616", "GSM1421617", "GSM1421618", "GSM1421619", "GSM1421620",
    "GSM1421621", "GSM1421622", "GSM1421623", "GSM1421624", "GSM1421625",
    "GSM1421626", "GSM1421627", "GSM1421628", "GSM1421629", "GSM1421630",
    "GSM1421631", "GSM1421632", "GSM1421633", "GSM1421634", "GSM1421635",
    "GSM1421636", "GSM1421637", "GSM1421638", "GSM1421639", "GSM1421640",
    "GSM1421641", "GSM1421642", "GSM1421643", "GSM1421644", "GSM1421645",
    "GSM1421646", "GSM1421647", "GSM1421648", "GSM1421649", "GSM1421650",
    "GSM1421651", "GSM1421652", "GSM1421653", "GSM1421654", "GSM1421655",
    "GSM1421656", "GSM1421657", "GSM1421658", "GSM1421659", "GSM1421660",
    "GSM1421661", "GSM1421662", "GSM1421663", "GSM1421664", "GSM1421665",
    "GSM1421666", "GSM1421667", "GSM1421668", "GSM1421669", "GSM1421670",
    "GSM1421671", "GSM1421672", "GSM1421673", "GSM1421674", "GSM1421675",
    "GSM1421676", "GSM1421677", "GSM1421678", "GSM1421679", "GSM1421680",
    "GSM1421681", "GSM1421682", "GSM1421683", "GSM1421684", "GSM1421685",
    "GSM1421686", "GSM1421687", "GSM1421688", "GSM1421689", "GSM1421690",
    "GSM1421691", "GSM1421692", "GSM1421693", "GSM1421694", "GSM1421695",
    "GSM1421696", "GSM1421697", "GSM1421698", "GSM1421699", "GSM1421700",
    "GSM1421701", "GSM1421702", "GSM1421703", "GSM1421704", "GSM1421705",
    "GSM1421706", "GSM1421707", "GSM1421708", "GSM1421709", "GSM1421710",
    "GSM1421711", "GSM1421712", "GSM1421713", "GSM1421714", "GSM1421715",
    "GSM1421716", "GSM1421717", "GSM1421718", "GSM1421719", "GSM1421720",
    "GSM1421721", "GSM1421722", "GSM1421723", "GSM1421724", "GSM1421725",
    "GSM1421726", "GSM1421727", "GSM1421728", "GSM1421729", "GSM1421730",
    "GSM1421731", "GSM1421732", "GSM1421733", "GSM1421734", "GSM1421735",
    "GSM1421736", "GSM1421737", "GSM1421738", "GSM1421739", "GSM1421740",
    "GSM1421741", "GSM1421742", "GSM1421743", "GSM1421744", "GSM1421745",
    "GSM1421746", "GSM1421747", "GSM1421748", "GSM1421749", "GSM1421750",
    "GSM1421751", "GSM1421752", "GSM1421753"
  ),
  Column_Prefix = c(
    "28028", "28048", "28165", "28172", "28191", 
    "28200", "28242", "28258", "28272", "28277", 
    "28284", "28301", "28306", "28381", "28396", 
    "28399", "28487", "28489", "28502", "28542", 
    "28544", "28559", "28581", "28582", "28583", 
    "28592", "28611", "28663", "28675", "28704", 
    "28752", "28778", "28806", "28015", "28020", 
    "28032", "28035", "28042", "28058", "28073", 
    "28119", "28123", "28131", "28140", "28148", 
    "28175", "28179", "28183", "28188", "28189", 
    "28198", "28202", "28204", "28232", "28247", 
    "28269", "28286", "28288", "28289", "28290", 
    "28291", "28298", "28303", "28319", "28330", 
    "28341", "28346", "28358", "28376", "28393", 
    "28395", "28401", "28406", "28438", "28448", 
    "28453", "28456", "28472", "28498", "28501", 
    "28575", "28588", "28591", "28600", "28603", 
    "28612", "28614", "28615", "28623", "28629", 
    "28630", "28647", "28649", "28655", "28657", 
    "28664", "28667", "28698", "28706", "28716", 
    "28720", "28726", "28735", "28736", "28760", 
    "28782", "28798", "28803", "28807", "28812", 
    "28816", "28861", "28895", "28897", "81012", 
    "81018", "81019", "81059", "81064", "81074", 
    "81084", "81136", "42008", "42017", "42020", 
    "42023", "42029", "42034", "42038", "42005", 
    "42007", "42009", "42016", "42018", "42019", 
    "42021", "42022", "42025", "42028", "42030", 
    "42035", "42036", "42037"
  )
)

# 打印映射表
print(gsm_mapping)

######################### 替换列名为GSM
# 创建一个函数，用于替换列名前缀为对应的 GSM ID
replace_column_names <- function(colnames, gsm_mapping) {
  new_colnames <- colnames
  for (i in 1:nrow(gsm_mapping)) {
    new_colnames <- gsub(gsm_mapping$Column_Prefix[i], gsm_mapping$GSM_ID[i], new_colnames)
  }
  return(new_colnames)
}

# 替换列名
colnames(gse58888) <- replace_column_names(colnames(gse58888), gsm_mapping)

# 打印新的列名
print(colnames(gse58888))

gse58888 <- as.data.frame(gse58888)
rownames(gse58888) <- gse58888$ID_REF


################################### 拆分M、U、P三列
# 除去 ID_REF 这一列
gse58888$ID_REF <- NULL

# UnMethylated / Methylated / Pval 列名
u_cols <- grep("_UnMethylated_signal$", colnames(gse58888), value = TRUE)
m_cols <- grep("_Methylated_Signal$",  colnames(gse58888), value = TRUE)
p_cols <- grep("_Detection_Pval$",     colnames(gse58888), value = TRUE)

# 提取样本 ID（81136, 81084 这种）
u_ids <- sub("_UnMethylated_signal$", "", u_cols)
m_ids <- sub("_Methylated_Signal$",   "", m_cols)
p_ids <- sub("_Detection_Pval$",      "", p_cols)

# 确保顺序一致（防止列顺序乱）
sample_ids <- sort(unique(u_ids))

u_cols <- u_cols[match(sample_ids, u_ids)]
m_cols <- m_cols[match(sample_ids, m_ids)]
p_cols <- p_cols[match(sample_ids, p_ids)]

# 组合成矩阵
U <- as.matrix(gse58888[, u_cols, drop = FALSE])
M <- as.matrix(gse58888[, m_cols, drop = FALSE])
P <- as.matrix(gse58888[, p_cols, drop = FALSE])

rownames(U) <- rownames(M) <- rownames(P) <- rownames(gse58888)
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
output_dir <- "/data/zhangmch/ewas_array/result/GSE58888"
save_columns_as_files(beta, output_dir)
