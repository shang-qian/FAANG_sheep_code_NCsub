
library(data.table)
####Enh_location
file1="~/FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt"
P1=fread(file1,header=T)
Pair12= paste(P1$promoter,P1$enhancer,sep="_")

split_result <- strsplit(Pair12, "_")
split_df <- do.call(rbind, split_result)
colnames(split_df) <- c("GeneName", "Promoter_Chrom", "Promoter_Start", "Promoter_End", "Enhancer_Chrom", "Enhancer_Start", "Enhancer_End")

gene_promoter=unique(split_df[,1:4])
gene_enhancer=unique(split_df[,c(1,5:7)])


###methylation level
meth_files <- list.files(path ="~/FAANG/analysis/11_2Methy/11_2Methy/01Methy_region/01Methy_loc", pattern = "*methy.bed", full.names = TRUE, recursive = TRUE)
print(meth_files)

read_bed <- function(file) {
gr=fread(file, header = F)
return(gr)
}

ME_list <- lapply(meth_files, read_bed)
tissue_name_ME <- gsub(".*/([^/]+)\\.methy\\.bed$", "\\1", meth_files)


args=c("~/FAANG/analysis/09Promoter_TPM/03pro_TPM")
input_files <- list.files(path =args[1], pattern = "*_01final_promoter_TPM.txt", full.names = TRUE, recursive = TRUE)
print(input_files)

read_bed1 <- function(file) {
gr=fread(file, header = T)
return(gr)
}


GE_list <- lapply(input_files, read_bed1)
tissue_names_GE <- gsub(".*/03pro_TPM/([^/]+)/.*", "\\1", input_files)


colnames(gene_promoter)=c("geneID","Pchr","Pstart","Pend")
gene_promoter <- as.data.frame(gene_promoter)

for (i in 1:length(tissue_name_ME))
{
tis=tissue_name_ME[i]

Tis_ME=ME_list[i]
Tis_GE=GE_list[which(tissue_names_GE==tis)]

Tis_GE_df <- as.data.frame(Tis_GE[[1]])

library(data.table)

GE_prom_df=Tis_GE_df[,c(1,2,9,10,17)]   ##gene, chr, start, end, TPM

# 确保数据是 data.table 格式
merged_data_dt <- as.data.table(GE_prom_df)
colnames(merged_data_dt)=c("geneID","Pchr","Pstart","Pend","TPM")

Tis_ME_dt <- as.data.table(Tis_ME[[1]])
colnames(Tis_ME_dt) <- c("Pchr", "Position", "End", "Value")

# 确保 Pstart 和 Pend 是整数类型
merged_data_dt[, Pstart := as.integer(Pstart)]
merged_data_dt[, Pend := as.integer(Pend)]

# 添加两个新列并初始化
merged_data_dt[, `:=`(Num_Sites = 0, Mean_Value = NA)]

# 设置键以加速区间联接
setkey(Tis_ME_dt, Pchr, Position, End)

# 区间联接并计算统计信息
results <- Tis_ME_dt[merged_data_dt, 
                     on = .(Pchr = Pchr, Position >= Pstart, Position <= Pend), 
                     .(Num_Sites = .N, Mean_Value = mean(Value, na.rm = TRUE)), 
                     by = .EACHI]

# 将结果合并回 merged_data_dt
merged_data_dt[, Num_Sites := results$Num_Sites]
merged_data_dt[, Mean_Value := results$Mean_Value]


merged_data_clean <- na.omit(merged_data_dt)


# 保存绘图到 PDF 文件
pdf(paste0(tis,"_Promoter_all.TPM_ML.pdf"))

plot(
  merged_data_clean$Mean_Value, 
  log10(merged_data_clean$TPM+0.0001),
  xlab = "Mean methylation level in promoter region", 
  ylab = "TPM",  # 对数变换的 Y 轴标签
  main = paste0(tis,"_Promoter_ML_TPM"),
  pch = 19,       # 实心圆点
  col = "orangered4",    # 蓝色点
  axes=F
)

axis(2, at = c(-4,-2,0,2,4), labels = c("0","0.01","1","100","10000"))
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1"))

# 
lm_model <- lm(log10(TPM+0.0001) ~ Mean_Value, data = merged_data_clean)

#
abline(lm_model, col = "red", lwd = 2)  # 红色线，线宽 2

# 
dev.off()

}



