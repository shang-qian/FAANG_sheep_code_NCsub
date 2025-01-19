# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

#####Annotated file GFF
gff_file <- "/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff"
gff_data <- import.gff(gff_file)

chrom_lengths <- aggregate(end ~ seqnames, data = gff_data, FUN = max)
chrom_lengths = chrom_lengths[1:28,]

# 创建全基因组的 GRanges 对象
genome_gr <- GRanges(
  seqnames = chrom_lengths$seqnames,
  ranges = IRanges(start = 1, end = chrom_lengths$end)
)


# 筛选基因或转录本注释
genes <- gff_data[gff_data$type %in% c("gene")|gff_data$type %in% c("pseudogene")]
Gene_data <- data.frame(  seqnames = seqnames(genes),  start = start(genes),  end = end(genes),  strand = strand(genes),  geneID = mcols(genes)$ID  )


exon <- gff_data[gff_data$type %in% c("exon")]
Exon_data <- data.frame(  seqnames = seqnames(exon),  start = start(exon),  end = end(exon),  strand = strand(exon),  geneID = mcols(exon)$ID  )


library(GenomicRanges)
# 创建 Gene_data 和 Exon_data 的 GRanges 对象
Gene_gr<- GRanges(seqnames = Gene_data$seqnames,
                     ranges = IRanges(start = Gene_data$start, end = Gene_data$end),
                     strand = Gene_data$strand,
                     geneID = Gene_data$geneID)
Exon_gr <- GRanges(seqnames = Exon_data$seqnames,
                     ranges = IRanges(start = Exon_data$start, end = Exon_data$end),
                     strand = Exon_data$strand,
                     geneID = Exon_data$geneID)

intron_regions <- setdiff(Gene_gr, Exon_gr)




# promoter文件
promoter_df <- read.table("~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final/11Gene_promoter.info", header = F, sep = "\t")
promoter_gr <- GRanges(  seqnames = promoter_df$V1,  ranges = IRanges(start = promoter_df$V2, end = promoter_df$V3),  GeneID = promoter_df$V4)
  
# enhancer   
enhancer_df <- read.table("~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt", header = TRUE, sep = "\t")
enhancer_gr <- GRanges(  seqnames = enhancer_df$seqnames,  ranges = IRanges(start = enhancer_df$start, end = enhancer_df$end),  strand = enhancer_df$strand)
    

# 从基因组中依次扣除基因、启动子和增强子区域
# 将基因组的链信息设为未指定（保持原样）
strand(genome_gr) <- "*"


ggene_gr=Gene_gr
gpromoter_gr=promoter_gr
genhancer_gr=enhancer_gr

strand(ggene_gr) <- "*"
strand(gpromoter_gr)="*"
strand(genhancer_gr)="*"

intergenic_regions_step1 <- setdiff(genome_gr, ggene_gr)
intergenic_regions_step2 <- setdiff(intergenic_regions_step1, gpromoter_gr)

# 如果有 Enhancer_gr，扣除增强子区域
intergenic_regions <- setdiff(intergenic_regions_step2, genhancer_gr)

# 查看最终的 intergenic 区域
intergenic_regions






#SNP
SNP_df <- read.table("~/FAANG/analysis/13SNP/02New/01SNP_position_chr_all.bed", header = F, sep = "\t")
# 建立替换列表
chr_map <- c(
  "chr1" = "CM028704.1", "chr2" = "CM028705.1", "chr3" = "CM028706.1",
  "chr4" = "CM028707.1", "chr5" = "CM028708.1", "chr6" = "CM028709.1",
  "chr7" = "CM028710.1", "chr8" = "CM028711.1", "chr9" = "CM028712.1",
  "chr10" = "CM028713.1", "chr11" = "CM028714.1", "chr12" = "CM028715.1",
  "chr13" = "CM028716.1", "chr14" = "CM028717.1", "chr15" = "CM028718.1",
  "chr16" = "CM028719.1", "chr17" = "CM028720.1", "chr18" = "CM028721.1",
  "chr19" = "CM028722.1", "chr20" = "CM028723.1", "chr21" = "CM028724.1",
  "chr22" = "CM028725.1", "chr23" = "CM028726.1", "chr24" = "CM028727.1",
  "chr25" = "CM028728.1", "chr26" = "CM028729.1", "chrX" = "CM028730.1")

# 替换 V1 列中的 chr 编号
SNP_df$V1 <- unname(chr_map[SNP_df$V1])

SNP_gr= GRanges(seqnames = SNP_df$V1,  ranges = IRanges(start = SNP_df$V2, end = SNP_df$V3),  SNPname = SNP_df$V4, SNPID= SNP_df$V5)

write.table(as.data.frame(SNP_df), file = paste0("00SNP_all577834.txt"), sep = "\t", row.names = T, col.names = T, quote = FALSE)




SNP_den=function(SNP_gr,enhancer_gr)
{
SNP_enh=findOverlaps(SNP_gr,enhancer_gr)
SNP_n=length(unique(queryHits(SNP_enh)))
enh_s=enhancer_gr[unique(subjectHits(SNP_enh))]
Total_len=sum(width(enh_s))
SNP_density=SNP_n/Total_len
res_enh=c(SNP_n,length(enh_s),Total_len,SNP_density)
return(res_enh)
}

SNPden=NULL
SNPden=rbind(SNPden,SNP_den(SNP_gr,Exon_gr))
SNPden=rbind(SNPden,SNP_den(SNP_gr,intron_regions))
SNPden=rbind(SNPden,SNP_den(SNP_gr,promoter_gr))
SNPden=rbind(SNPden,SNP_den(SNP_gr,enhancer_gr))
SNPden=rbind(SNPden,SNP_den(SNP_gr,intergenic_regions))


colnames(SNPden)=c("SNP_no","feat_no","Feature_length","density")
rownames(SNPden)=c("Exon","Intron","Promoter","Enhancer","Intergenic")

write.table(as.data.frame(Exon_gr), file = paste0("01SNP_location.Exon"), sep = "\t", row.names = T, col.names = T, quote = FALSE)
write.table(as.data.frame(intron_regions), file = paste0("01SNP_location.Intron"), sep = "\t", row.names = T, col.names = T, quote = FALSE)
write.table(as.data.frame(promoter_gr), file = paste0("01SNP_location.Promoter"), sep = "\t", row.names = T, col.names = T, quote = FALSE)
write.table(as.data.frame(enhancer_gr), file = paste0("01SNP_location.Enhancer"), sep = "\t", row.names = T, col.names = T, quote = FALSE)

write.table(as.data.frame(SNPden), file = paste0("02SNP_location.density"), sep = "\t", row.names = T, col.names = T, quote = FALSE)




QTL_df <- fread("~/FAANG/analysis/13SNP/02New/03QTL_position_chr.bed",header=F)
QTL_df$V1 <- unname(chr_map[QTL_df$V1])

QTL_gr <- GRanges(seqnames = QTL_df$V1,
                     ranges = IRanges(start = QTL_df$V2, end = QTL_df$V3),
                     Anno = QTL_df$V4)



###Promoter
prom_tmp=findOverlaps(QTL_gr,promoter_gr)
QTL_data=as.data.frame(QTL_gr[queryHits(prom_tmp)])
promoter_data=as.data.frame(promoter_gr[subjectHits(prom_tmp)])

promoter_QTL=cbind(QTL_data,promoter_data)
write.table(as.data.frame(promoter_QTL), file = paste0("04Promoter_QTL.txt"), sep = "\t", row.names = T, col.names = T, quote = FALSE)


###Enhancer
enh_tmp=findOverlaps(QTL_gr,enhancer_gr)
QTL_data=as.data.frame(QTL_gr[queryHits(enh_tmp)])
enhancer_data=as.data.frame(enhancer_gr[subjectHits(enh_tmp)])

enhancer_QTL=cbind(QTL_data,enhancer_data)
write.table(as.data.frame(enhancer_QTL), file = paste0("05Enhancer_QTL.txt"), sep = "\t", row.names = T, col.names = T, quote = FALSE)


####Pairs enhancer
Pair_df <- fread("~/FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt", header = T)

regions=Pair_df$Xname
chr <- sub("_.*", "", regions)
start <- as.integer(sub(".*_(\\d+)_.*", "\\1", regions))
end <- as.integer(sub(".*_(\\d+)$", "\\1", regions))
  
pair_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end),geneID=Pair_df$promoter)

####
pair_tmp=findOverlaps(QTL_gr,pair_gr)
QTL_data=as.data.frame(QTL_gr[queryHits(pair_tmp)])
pair_data1=as.data.frame(pair_gr[subjectHits(pair_tmp)])
pair_data2=as.data.frame(Pair_df[subjectHits(pair_tmp)])

pair_QTL=cbind(QTL_data,pair_data2)
write.table(as.data.frame(pair_QTL), file = paste0("06Pair_QTL.txt"), sep = "\t", row.names = T, col.names = T, quote = FALSE)



####COMMD1
#
pdf(paste0("../07COMMD1.pdf"))
data <- data.frame(
  Tissue = c("AdrenalCortex", "AdrenalMedulla", "Bladder", "Cerebellum", "CerebralCortex", 
             "DescendingColon", "Duodenum", "Gallbladder", "HeartRightAtrium", "HeartRightVentricle", 
             "Ileum", "IleumPeyersPatch", "Lung", "LymphNodeMesenteric", "MuscleSM", "Ovary", 
             "Oviduct", "Reticulum", "RumenAtrium", "SpinalCord", "SpiralColon", "Tongue", 
             "Tonsil", "Uterus"),
  Expression = c(73.995193, 60.528011, 88.452957, 49.573303, 80.665794, 0.192342, 51.265121, 
                 0.038178, 70.277763, 75.414467, 0, 53.109875, 46.661968, 58.547806, 0, 
                 71.186226, 48.696178, 55.030457, 30.493835, 74.481728, 37.727085, 0.702589, 
                 30.855618, 0),
  Category = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0)
)

# 根据表达量对数据进行降序排序
data <- data[order(data$Expression, decreasing = TRUE), ]
# 设置颜色：类别为 1 时使用红色，否则使用蓝色
colors <- ifelse(data$Category == 1, "red", "blue")
# 绘制散点图，使用较大的点
plot(data$Expression, pch = 19, col = colors, xaxt = "n", cex = 2, # cex 控制点的大小
      xlab ="Tissue", ylab = "TPM", ylim=c(0,100),
     main = "Expression Levels in Different Tissues")

# 添加 x 轴标签
axis(1, at = 1:length(data$Tissue), labels = data$Tissue, las = 2, cex.axis = 1)

# 在每个点上显示表达值
text(1:length(data$Expression), data$Expression, 
     labels = round(data$Expression, 2), pos = 3, cex = 0.8)

# 添加图例
legend("topright", legend = c("Yes", "No"), col = c("red", "blue"), pch = 19, title = "Enhancer",cex=1.5)

dev.off()



###Density
# 示例数据
SNPden <- data.frame(
  Category = c("Intergenic", "Intron", "Exon",  "Enhancer","Promoter"),
  SNP_no = c(280977, 236253,50803, 33015,1706),
  density = c(0.000219145, 0.0002572553, 0.0005060315, 0.0016289188,0.0022746788)
)

# 1. SNP_no 的饼图280977  67190 1282150814 0.0002191450
# 定义颜色
colors <- c("skyblue", "lightgreen", "red", "orange","yellow")
# 绘制饼图
pdf(paste0("../../08Pie.pdf"))
snp_labels <- paste(SNPden$Category, "\n", SNPden$SNP_no/sum(SNPden$SNP_no), sep = "")
pie(SNPden$SNP_no, labels = snp_labels, col = colors,
    main = "Distribution of SNP_no across genomic features")
dev.off()



# 2. density 的柱状图
pdf(paste0("../../09Density.pdf"))
bar_positions <- barplot(SNPden$density, names.arg = SNPden$Category, col = colors, ylim=c(0,0.0025), width=0.3,space = 1.5,
        ylab = "SNP density (#/bp)", xlab = "Genomic feature", main = "SNP density across genomic features")
text(x = bar_positions, y = SNPden$density, labels = round(SNPden$density, 5), 
     pos = 3, cex = 0.8, col = "black")
dev.off()





