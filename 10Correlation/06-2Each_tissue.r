# parameters
args <- commandArgs(trailingOnly = TRUE)

args=c("~/FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt",
"~/FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt")

library(data.table)

#library(parallel)
#library(GenomicRanges)


P1=fread(args[1], header = TRUE)
P2=fread(args[2],header=TRUE)

P13=P1[,c("promoter","enhancer","enh_tis#")]
colnames(P13)=c("genename","enh_tmp","enh_tis")
P23=P2[,c("genename","enh_tmp","enh_tis")]
Pall=rbind(P13,P23)

Pall$promoter=sub(".*?(CM.*)", "\\1", Pall$genename)
Pall$pair=paste(Pall$genename,Pall$enh_tmp,sep="_")



#########each tissue
tisname=colnames(P2)[3:26]
Pairs_num=matrix(,24,3)
rownames(Pairs_num)=tisname
colnames(Pairs_num)=c("1total_no","2total_enh","3total_pro")


labels_combined <- c(1:20, ">20")
Prom_hub=matrix(,24,21)
rownames(Prom_hub)=tisname
colnames(Prom_hub)=labels_combined

Enh_hub=matrix(,24,6)
rownames(Enh_hub)=tisname
colnames(Enh_hub)=c(1:5,">5")


for (i in 1:length(tisname) )
{
TisN=tisname[i]
print(TisN)
Tis_p=Pall[grepl(TisN,Pall$enh_tis),]
print(nrow(Tis_p))

dir.create(TisN)
write.table(as.data.frame(Tis_p), file=paste0(TisN,"/",TisN,"_01Enh_prom_total.txt"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


####Total pairs from Coexp CTCF, CAGE
total_no= length(unique(Tis_p$pair))
total_enh= length(unique(Tis_p$enh_tmp))
total_pro= length(unique(Tis_p$promoter))
Pairs_num[i,]=c(total_no,total_enh,total_pro)

# Promoter hub frequency
frequency <- table(table(unique(Tis_p)$promoter))

greater_than_30 <- sum(frequency[21:length(frequency)])
frequency_combined <- c(frequency[1:20], greater_than_30)
frequency_combined[is.na(frequency_combined)]=0
Prom_hub[i,]=frequency_combined

#Enhancer
frequency <- table(table(unique(Tis_p)$enh_tmp))
frequency_combined <- c(frequency[1:5], sum(frequency[6:length(frequency)]))
frequency_combined[is.na(frequency_combined)]=0
Enh_hub[i,]=frequency_combined

}

AllP=cbind(Pairs_num,Prom_hub,Enh_hub)
write.table(as.data.frame(AllP), file=paste0("01Enh_prom_hub_24_each_tissue.txt"), sep="\t", row.names=TRUE, quote=FALSE, col.names=TRUE)


All2=cbind(Prom_hub,Enh_hub)
#write.table(as.data.frame(All2), file=paste0("05Enh_prom_hub_2_each_tissue.txt"), sep="\t", row.names=TRUE, quote=FALSE, col.names=TRUE)
Ename=paste0("E",1:21)
Pname=paste0("P",1:6)
colnames(All2)=c(Ename,Pname)


library(pheatmap)
#
data <- as.data.frame(All2)

# 按身体部位顺序对组织进行排序
body_order <- c("CerebralCortex", "Cerebellum", "SpinalCord", "Lung", "HeartRightAtrium", "HeartRightVentricle", 
                "Gallbladder", "Duodenum", "Bladder", "AdrenalCortex", "AdrenalMedulla", 
                "Ileum", "IleumPeyersPatch", "DescendingColon", "SpiralColon", "RumenAtrium", "Reticulum", 
                "Ovary", "Oviduct", "Uterus", "MuscleSM", "Tongue", "Tonsil", "LymphNodeMesenteric")

# 按身体顺序进行排序
data_sorted <- data[match(body_order, rownames(data)),]
# 将 0 替换为 0.0001
data_sorted[data_sorted == 0] <- 0.1
# 将数据转换为矩阵
data_matrix <- as.matrix(data_sorted) 

# 对数据进行 log10 转换
data_matrix_log10 <- log2(data_matrix)

# 绘制热图，不进行列聚类
pdf(paste0("02Heatmap_hub.pdf")) # 打开PDF设备

pheatmap(data_matrix_log10, 
         cluster_cols = FALSE,   # 不对列进行聚类
         cluster_rows = FALSE,   # 不对行进行聚类
         color = colorRampPalette(c("white", "orange"))(1000),  # 颜色渐变
         border_color = NA,      # 去掉格线
         main = "Heatmap of Tissue Distributions (log2)",
         show_rownames = TRUE,   # 显示行名（组织名称在左侧）
         show_colnames = TRUE)  # 隐藏列名（如果你不需要列名）
dev.off()


# 创建柱状图
pdf(paste0("03Promoter_hub.pdf")) # 打开PDF设备

means <- apply(Prom_hub,2,mean)
std_dev <- apply(Prom_hub,2,sd)
labels <- c(1:20, "20+")

plot(1:21, means, type = "o", pch = 16, xaxt = "n", ylim = c(0, max(means + std_dev) * 1.2),col="orange",
     main = "Line Plot of Means with Standard Deviation", xlab = "Enhancer # in promoter hub", ylab = "Promoter hub #", las = 2)

# 绘制柱状图
#barplot_heights <- barplot(means, names.arg = labels, ylim = c(0, max(means + std_dev) * 1.2),
#                           col = "orange", main = "Bar Plot of Means with Standard Deviation",
#                           xlab = "Enhancer # in promoter hub", ylab = "Promoter hub #", las = 2)
# 添加标准差（误差条）
arrows(1:21, means - std_dev, 1:21, means + std_dev, angle = 90, code = 3, length = 0.1, col="orange")

# 添加自定义的 X 轴标签
axis(1, at = 1:21, labels = labels)

# 在每个点上显示整数数值
text(1:21, means + 10, labels = round(means, 1), cex = 0.8, pos = 3)
dev.off()

#####Enhancer
# 创建数据
means <- apply(Enh_hub,2,mean)
std_dev <- apply(Enh_hub,2,sd)
labels <- c("1", "2", "3", "4", "5", "5+")
# 绘制柱状图
pdf(paste0("04Enhancer_hub.pdf")) # 打开PDF设备

plot(1:6, means, type = "o", pch = 16, xaxt = "n", ylim = c(0, max(means + std_dev) * 1.2),col="red",
     main = "Line Plot of Means with Standard Deviation", xlab = "Enhancer # in promoter hub", ylab = "Promoter hub #", las = 2)

# 绘制柱状图
#barplot_heights <- barplot(means, names.arg = labels, ylim = c(0, max(means + std_dev) * 1.2),
#                           col = "orange", main = "Bar Plot of Means with Standard Deviation",
#                           xlab = "Enhancer # in promoter hub", ylab = "Promoter hub #", las = 2)
# 添加标准差（误差条）
arrows(1:6, means - std_dev, 1:6, means + std_dev, angle = 90, code = 3, length = 0.1, col="red")

# 添加自定义的 X 轴标签
axis(1, at = 1:6, labels = labels)

# 在每个点上显示整数数值
text(1:6, means  + 10, labels = round(means, 1), cex = 0.8, pos = 3)
dev.off()

####Brain
Brain=NULL
k=0
for( TisN in c("Cerebral", "Cerebellum","Spinal"))
{
print(TisN)
Tis_b1=Pall[grepl(TisN,Pall$enh_tis),]
frequency <- table(unique(Tis_b1)$promoter)
Pname=names(which(frequency>=10))
k=k+1
Brain[[k]]=Tis_b1[match(names(which(frequency>=10)),Tis_b1$promoter),"pair"]
}

gene_overlap1=intersect(intersect(Brain[[1]]$pair,Brain[[3]]$pair),Brain[[2]]$pair)
write.table(as.data.frame(gene_overlap1), file=paste0("05Brain_10enh_genes.txt"), sep="\t", row.names=TRUE, quote=FALSE, col.names=TRUE)


####Brain
Adrenal=NULL
k=0
for( TisN in c("AdrenalCortex", "AdrenalMe"))
{
print(TisN)
Tis_b1=Pall[grepl(TisN,Pall$enh_tis),]
frequency <- table(unique(Tis_b1)$promoter)
Pname=names(which(frequency>=10))
k=k+1
Adrenal[[k]]=Tis_b1[match(names(which(frequency>=20)),Tis_b1$promoter),"pair"]
}

gene_overlap=intersect(Adrenal[[1]]$pair,Adrenal[[2]]$pair)
write.table(as.data.frame(gene_overlap), file=paste0("06Adrenal_10enh_genes.txt"), sep="\t", row.names=TRUE, quote=FALSE, col.names=TRUE)


