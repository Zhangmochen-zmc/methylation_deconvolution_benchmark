# 1. config
ref_file_path <- "episcore_ref/wgbs_850k_ref_raw.csv"          # 你的 Ref 文件
sample_folder_path <- "/data/yuxy/data/data/wgbs_850k/episcore" # 存放 Project CSV 的文件夹
chip_type <- "850k"                      # "450k" 或"850k"

# --------------------------
# 2. 加载环境与函数
# --------------------------
library(EpiSCORE)
if(chip_type == "450k") { data("probeInfo450k") } else { data("probeInfo850k") }

# 转换函数 (保持不变，功能是 CpG -> Gene)
constAvBetaTSS <- function(beta.m, type=c("450k","850k","EPICv2")){
  if(type=="450k"){ if(exists("probeInfo450k.lv")) probeInfoALL.lv <- probeInfo450k.lv } 
  else if (type=="850k"){ probeInfoALL.lv <- probeInfo850k.lv } 
  else if (type=="EPICv2"){ data("probeInfoEPICv2"); probeInfoALL.lv <- probeInfoEPICv2.lv }
  
  map.idx <- match(rownames(beta.m), probeInfoALL.lv$probeID)
  probeInfo.lv <- lapply(probeInfoALL.lv, function(tmp.v,ext.idx){return(tmp.v[ext.idx]);}, map.idx);
  beta.lm <- list();
  
  for (g in 1:6) {
    group.idx <- which(probeInfo.lv$GeneGroup == g)
    if(length(group.idx) == 0) next;
    tmp.m <- beta.m[group.idx, , drop=FALSE]
    rownames(tmp.m) <- probeInfo.lv$EID[group.idx];
    sel.idx <- which(!is.na(rownames(tmp.m)));
    tmp.m <- tmp.m[sel.idx, , drop=FALSE];
    if(nrow(tmp.m) == 0) next;
    nL <- length(factor(rownames(tmp.m)));
    nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
    beta.lm[[g]] <- rowsum(tmp.m, group=rownames(tmp.m))/nspg.v;
  }
  genes2 <- if(!is.null(beta.lm[[2]])) rownames(beta.lm[[2]]) else c()
  genes4 <- if(!is.null(beta.lm[[4]])) rownames(beta.lm[[4]]) else c()
  unqEID.v <- unique(c(genes2, genes4));
  avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(beta.m))
  colnames(avbeta.m) <- colnames(beta.m)
  rownames(avbeta.m) <- unqEID.v
  if(!is.null(beta.lm[[4]])) avbeta.m[match(rownames(beta.lm[[4]]), rownames(avbeta.m)), ] <- beta.lm[[4]]
  if(!is.null(beta.lm[[2]])) avbeta.m[match(rownames(beta.lm[[2]]), rownames(avbeta.m)), ] <- beta.lm[[2]]
  return(avbeta.m);
}

# --------------------------
# 3. 处理 Reference
# --------------------------
cat(">>> 处理 Reference 矩阵...\n")
if(!file.exists(ref_file_path)) stop("找不到 Ref 文件！")
raw_ref <- read.csv(ref_file_path, row.names=1, check.names=FALSE)
# 转换后：行=Gene, 列=CellType
ref_gene_full <- constAvBetaTSS(as.matrix(raw_ref), type=chip_type) 
cat(paste("Ref 处理完毕，基因数:", nrow(ref_gene_full), "\n"))

# --------------------------
# 4. 批量处理 Project CSV
# --------------------------
project_files <- list.files(path = sample_folder_path, pattern = "\\.csv$", full.names = TRUE)
if(length(project_files) == 0) stop("文件夹为空！")

cat(paste(">>> 发现", length(project_files), "个项目文件，开始处理...\n"))

# 列表用于存储处理后的矩阵
# 每个元素的结构：Row=Genes, Col=Samples (Patient1, Patient2...)
project_data_list <- list()

for(i in 1:length(project_files)){
  f_path <- project_files[i]
  f_name <- basename(f_path)
  
  cat(paste0(">>> [", i, "/", length(project_files), "] 处理 Project: ", f_name, " ... "))
  
  # 1. 读取原始矩阵 (Row=CpG, Col=Samples)
  raw_proj <- read.csv(f_path, row.names=1, check.names=FALSE)
  mat_proj <- as.matrix(raw_proj)
  
  # 2. 转换为 Gene-level (Row=Gene, Col=Samples)
  # 注意：列名保持不变，即原本是 Patient1, 转换后还是 Patient1
  gene_proj <- constAvBetaTSS(mat_proj, type=chip_type)
  
  # 3. 存入列表
  project_data_list[[f_name]] <- gene_proj
  
  cat(paste("完成 (包含样本数:", ncol(gene_proj), ")\n"))
}

# --------------------------
# 5. 保存中间结果
# --------------------------
save(ref_gene_full, project_data_list, file = "ref/episcore.RData")
cat(">>> 数据已保存至 'Project_Level_Data.RData'\n")
