
library(data.table)
####Enh_location
file1="~/FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt"
P1=fread(file1,header=T)
file2="~/FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt"
P2=fread(file2,header=T)
Pair12= c(paste(P1$promoter,P1$enhancer,sep="_"),paste(P2$genename,P2$enh_tmp,sep="_"))
Cor_pair=sub(".*?(CM.*)", "\\1", Pair12)

EP1=P1[P1$source!="None",1:2]
EP2=P2[,1:2]
EP12=rbind(as.matrix(EP1),as.matrix(EP2))
dim(unique(EP12))
length(unique(EP12[,1]))
length(unique(EP12[,2]))

#Prom_T=unique(P1$promoter,P2$genename)
#Enh_T=unique(P1$enhancer,P2$enh_tmp)

###TPM
Gene_sample=fread("~/FAANG/analysis/10interaction/03Cor8/02Exp_gene.txt",header=T)
Gene_sample=as.data.frame(Gene_sample)
TPM=Gene_sample[,-ncol(Gene_sample)]
TPMn=paste(TPM$geneID, TPM$Pseq, TPM$Pstart, TPM$Pend,sep="_")

####Enh_location
fileenh="~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt"
Enh_all=fread(fileenh,header=T)
Enhn=paste(Enh_all$seqnames, Enh_all$start, Enh_all$end,sep="_")

EP_all=NULL
for (i in 1:nrow(P1))
{
  print(i)
  Pindex=match(P1$promoter[i],TPMn)
  TPM24=TPM[Pindex,5:28]
  #enh
  Eindex=match(P1$enhancer[i],Enhn)
  Enh24=Enh_all[Eindex,6:29]
  
  EP24=t(rbind(Enh24,TPM24))
  EP24_info=cbind(P1$promoter[i],P1$enhancer[i],EP24)
  EP_all=rbind(EP_all,EP24_info)
}

write.table(as.data.frame(EP_all), file=paste0("03EP_TPM_24tissues_all.txt"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


# 假设EP_all为数据框，将第三列和第四列提取出来
data <- as.data.frame(EP_all)
# 将第三列和第四列转换为数值型
data[, 3] <- as.numeric(data[, 3])
data[, 4] <- as.numeric(data[, 4])
# 分别提取出第三列为1和0的对应的第四列数值
group_1 <- log10(data[data[, 3] == 1, 4]+0.0001)  # 对应第3列为1的第四列数值
group_0 <- log10(data[data[, 3] == 0, 4]+0.0001)  # 对应第3列为0的第四列数值
# 进行显著性检验 (Mann-Whitney U test 非参数检验)
test_result <- wilcox.test(group_1, group_0)
# 输出检验结果
print(test_result)
# 使用 R 自带的 boxplot 函数绘图，并隐藏离群点
pdf(paste0("04Pairs_TPM_01_vioplot.pdf"))
# 绘制小提琴图并加上散点
vioplot(group_0, group_1,
        names = c("Without enhancer", "with enhancer"),
        col = c("red", "orange"),
        main = "Comparison of TPM between Group 0 and Group 1",
        ylab = "TPM",
        ylim = c(-4, 6),
        yaxt = "n"             
        )
axis(2, at = c(-4, -2, 0, 2, 4, 6), labels = c("0", expression(10^-2), "1", expression(10^2), expression(10^4), expression(10^6)))    
# 添加散点
#points(jitter(rep(1, length(group_0))), group_0, pch = 16, col = "black", cex = 0.8)
#points(jitter(rep(2, length(group_1))), group_1, pch = 16, col = "black", cex = 0.8)        
dev.off()


data <- as.data.frame(EP_all)
# 将第三列和第四列转换为数值型
data[, 3] <- as.numeric(data[, 3])
data[, 4] <- as.numeric(data[, 4])
split_result <- do.call(rbind, strsplit(as.character(data$V2), "_"))
# 将结果存储在 data 中的新列
data$Chromosome <- split_result[, 1]  # 第一列: 染色体
data$Start <- split_result[, 2]       # 第二列: 起始位置
data$End <- split_result[, 3]         # 第三列: 结束位置
data$tissue <- rownames(EP_all)

tissue_colors <- c(
  "AdrenalCortex" = "red",
  "AdrenalMedulla" = "blue",
  "Bladder" = "green",
  "Cerebellum" = "purple",
  "CerebralCortex" = "orange",
  "DescendingColon" = "pink",
  "Duodenum" = "cyan",
  "Gallbladder" = "magenta",
  "HeartRightAtrium" = "yellow",
  "HeartRightVentricle" = "darkblue",
  "Ileum" = "darkgreen",
  "IleumPeyersPatch" = "darkred",
  "Lung" = "lightblue",
  "LymphNodeMesenteric" = "lightgreen",
  "MuscleSM" = "darkorange",
  "Ovary" = "lightpink",
  "Oviduct" = "purple",
  "Reticulum" = "lightyellow",
  "RumenAtrium" = "gold",
  "SpinalCord" = "violet",
  "SpiralColon" = "salmon",
  "Tongue" = "turquoise",
  "Tonsil" = "beige",
  "Uterus" = "brown"
)

# 如果 V3 == 0 则为灰色，否则根据 tissue 映射颜色
data$color <- ifelse(data$V3 == 0, "grey", tissue_colors[data$tissue])
# 加载 ggplot2 包
library(ggplot2)

data$Start <- as.numeric(as.character(data$Start))
# 绘制散点图，按照染色体分组
pdf(paste0("05Pairs_TPM_01_genome_position20241203.pdf"))
ggplot(data, aes(x = Start / 1000000, y = log10(V4+0.0001), color = color)) +
  geom_point() +  # 绘制散点
  labs(title = "Scatter plot of TPM along Chromosome positions",
       x = "Genome position (10^6)",
       y = "TPM") +
  theme_minimal() +
  scale_color_identity() +  # 使用指定的颜色映射
  facet_wrap(~ Chromosome, scales = "free_x") + # 按染色体分面绘图
  scale_x_continuous(labels = function(x) round(x)) +  # 自定义横坐标刻度
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6), 
                     labels = c("0", expression(10^-2), "1", expression(10^2), expression(10^4), expression(10^6))) 
dev.off()  


library(ggplot2)
library(dplyr)

# 计算每个染色体的长度并生成绝对起始位置
chrom_lengths <- data %>%
  group_by(Chromosome) %>%
  summarise(max_length = max(Start))  # 计算每个染色体的长度

# 将绝对坐标加入原始数据
data_abs <- data %>%
  left_join(chrom_lengths, by = "Chromosome") %>%
  mutate(abs_position = Start)  # 这里直接保留原始坐标

# 输出到 PDF，每页一条染色体，x 轴固定范围
pdf("05Pairs_TPM_Chromosome_Scatterplots_NoGridlines20241203.pdf", width = 8, height = 6)

for (chrom in unique(data_abs$Chromosome)) {
  # 过滤当前染色体的数据
  chrom_data <- data_abs %>% filter(Chromosome == chrom)
  
  # 绘制当前染色体的图
  p <- ggplot(chrom_data, aes(x = abs_position / 1e6, y = log10(V4 + 0.0001), color = color)) +
    geom_point() +
    labs(
      title = paste("Scatter plot of TPM - Chromosome", chrom),
      x = "Genome position (Mb)",
      y = "TPM"
    ) +
    theme(
      panel.grid = element_blank(),  # 去除网格线
      panel.background = element_blank(),  # 移除背景
      axis.line = element_line(color = "black"),  # 添加坐标轴线
      axis.text = element_text(size = 20),  # 坐标轴刻度文字大小
      axis.title = element_text(size = 22)  # 坐标轴标题文字大小
    ) +
    scale_color_identity() +
    scale_x_continuous(
      limits = c(0, 300),  # 固定 x 轴范围为 0 到 300 Mb
      breaks = seq(0, 300, 50),  # 每隔 50 Mb 设置一个刻度
      labels = scales::comma  # 格式化刻度
    ) +
    scale_y_continuous(
      limits = c(-4, 4),  
      breaks = c(-4, -2, 0, 2, 4),
      labels = c("0", "0.01", "1", "1.0e2", "1.0e4")
    )

  # 输出图到 PDF 的新页面
  print(p)
}

dev.off()




# select chr1 as represent
chr1_data <- subset(data, Chromosome == "CM028704.1")

# 绘制一号染色体的散点图
pdf("05Pairs_TPM_01_chr1_genome_position.pdf")  # 修改输出文件名
ggplot(chr1_data, aes(x = Start / 1000000, y = log10(V4 + 0.0001), color = color)) +
  geom_point() +  # 绘制散点
  labs(title = "Scatter plot of TPM along Chromosome 1 positions",
       x = "Chr1 position (Mb)",
       y = "TPM") +
  theme_minimal() +
  scale_color_identity() +  # 使用指定的颜色映射
  scale_x_continuous(labels = function(x) round(x)) +  # 自定义横坐标刻度
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4, 6), 
                     labels = c("0", expression(10^-2), "1", expression(10^2), expression(10^4), expression(10^6))) 
dev.off()



# 加载 ggplot2 包
library(ggplot2)
# 为你提供的24个组织赋值颜色
tissue_colors <- c(
  "AdrenalCortex" = "red",
  "AdrenalMedulla" = "blue",
  "Bladder" = "green",
  "Cerebellum" = "purple",
  "CerebralCortex" = "orange",
  "DescendingColon" = "pink",
  "Duodenum" = "cyan",
  "Gallbladder" = "magenta",
  "HeartRightAtrium" = "yellow",
  "HeartRightVentricle" = "darkblue",
  "Ileum" = "darkgreen",
  "IleumPeyersPatch" = "darkred",
  "Lung" = "lightblue",
  "LymphNodeMesenteric" = "lightgreen",
  "MuscleSM" = "darkorange",
  "Ovary" = "lightpink",
  "Oviduct" = "darkpurple",
  "Reticulum" = "lightyellow",
  "RumenAtrium" = "gold",
  "SpinalCord" = "violet",
  "SpiralColon" = "salmon",
  "Tongue" = "turquoise",
  "Tonsil" = "beige",
  "Uterus" = "brown"
)
# 模拟数据，如果 V3 == 0 则为灰色，否则根据 tissue 映射颜色
data$color <- ifelse(data$V3 == 0, "grey", tissue_colors[data$tissue])
# 将 Start 和 End 转换为数值型
data$Start <- as.numeric(as.character(data$Start))
data$End <- as.numeric(as.character(data$End))

# 计算每个染色体的最大终点（长度）
chrom_lengths <- aggregate(End ~ Chromosome, data, max)
# 计算累积起始位置
chrom_lengths$cumulative_start <- c(0, cumsum(chrom_lengths$End[-nrow(chrom_lengths)]))
# 将每个染色体的累积起始位置加到 Start 上，形成全基因组的坐标
data <- merge(data, chrom_lengths[, c("Chromosome", "cumulative_start")], by = "Chromosome")
data$genome_start <- data$Start + data$cumulative_start
# 使用基础的 plot() 函数绘制散点图
pdf(paste0("06Pairs_TPM_01_genome_position.pdf"))
plot(
  data$genome_start, log10(data$V4+0.0001),
  col = data$color, pch = 16,
  xlab = "Chromosomes",
  ylab = "TPM",
  main = "Scatter plot of TPM along concatenated genome",
  cex = 0.5,
  xaxt = "n", yaxt = "n",  # 关闭默认的坐标轴
  ylim=c(-4,6)
)
# 在染色体分界处绘制虚线
for (i in chrom_lengths$cumulative_start[-1]) {
  abline(v = i, lty = 3, col = "black")  # 虚线
}
chromosome_labels <- c(paste0(1:26), "X")
chrom_lengths$cumulative_mid <- chrom_lengths$cumulative_start + chrom_lengths$End / 2

axis(1, at = chrom_lengths$cumulative_mid, labels = chromosome_labels, las = 1)
axis(2, at = c(-4, -2, 0, 2, 4, 6), labels = c("0", expression(10^-2), "1", expression(10^2), expression(10^4), expression(10^6)))

# 添加图例# 构建图例，包括24个组织+1个 nonexp
legend_labels <- c(names(tissue_colors), "nonexp")
legend_colors <- c(tissue_colors, "grey")

# 添加图例，显示25个颜色
legend("topright", legend = legend_labels, col = legend_colors, pch = 16, cex = 0.7)
dev.off()



