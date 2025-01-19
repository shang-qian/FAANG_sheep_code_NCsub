
library(data.table)
####Enh_location
file1="~/FAANG/analysis/10interaction/03Cor8/00oneAll/gene-ADORA1_CM028715.1_700901_700924oneAll.txt"
P1=fread(file1,header=F)
P1=as.data.frame(P1[1:24])
colnames(P1)=c("Tissue","TPM","chr12_682652_683139","chr12_737891_738081")

tissue_colors <- c(
  "AdrenalCortex" = "red",
  "AdrenalMedulla" = "blue",
  "Bladder" = "green",
  "Cerebellum" = "purple",
  "CerebralCortex" = "pink",
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

# 绘制散点图，映射颜色为组织的颜色
pdf("07gene-ADORA1_scatter_plot.pdf")
# 设置颜色：TPM为0的用灰色，TPM为1用组织的颜色
colors_chr12_682652 <- ifelse(P1$chr12_682652_683139 == 0, "grey", tissue_colors[P1$Tissue])
colors_chr12_737891 <- ifelse(P1$chr12_737891_738081 == 0, "grey", tissue_colors[P1$Tissue])
# 创建基础图（圆点表示 chr12:682652-683139）
plot(P1$chr12_682652_683139, log10(as.numeric(P1$TPM)+0.0001), 
     col = colors_chr12_682652, pch = 16, cex = 1,  # 圆点表示 chr12:682652-683139
     xlab = "chr12 positions", ylab = "TPM",
     main = "Scatter plot for both chr12 positions vs TPM",
     xlim = c(-0.5, 1.5), ylim = range(log10(as.numeric(P1$TPM)+0.0001)),
	 yaxt = "n",  
	 xaxt="n")  # 设置坐标范围
axis(2, at = c(-4, -2, 0, 2, 4, 6), labels = c("0", expression(10^-2), "1", expression(10^2), expression(10^4), expression(10^6)))    
axis(1, at = c(0, 1), labels = c("without enhancer", "with enhancer"))    
# 添加第二个 chr12 位置的散点图（三角形表示 chr12:737891-738081）
points(P1$chr12_737891_738081, log10(as.numeric(P1$TPM)+0.0001), 
       col = colors_chr12_737891, pch = 2, cex = 1.5)  # 三角形表示 chr12:737891-738081
# 为第一个 chr12 添加回归线
model1 <- lm(log10(as.numeric(P1$TPM)+0.0001) ~ P1$chr12_682652_683139)
abline(model1, col = "orange", lwd = 2)
# 在图上标出第一个 chr12 的相关系数和 R²
text(0.5,-2, paste0("chr12:682652-683139 (Cor:0.844, R²:0.698)"), col = "orange")
# 第二个 chr12 添加回归线
model2 <- lm(log10(as.numeric(P1$TPM)+0.0001) ~ P1$chr12_737891_738081)
abline(model2, col = "orange", lwd = 2)  # 第二个回归线用虚线表示
# 计算第二个 chr12 的相关系数和 R²
cor2 <- cor(log10(as.numeric(P1$TPM)+0.0001), as.numeric(P1$chr12_737891_738081))
r2_2 <- summary(model2)$r.squared
# 在图上标出第二个 chr12 的相关系数和 R²
text(0.3, 0.8 ,paste("chr12:737891-738081 (Cor:", round(cor2, 3),"R²:", round(r2_2, 3),")"), col = "orange")

# 标注组织名到 chr12_682652_683139 对应的点
text(P1$chr12_682652_683139, log10(as.numeric(P1$TPM)+0.0001), labels = paste(P1$Tissue, "TPM=", round(as.numeric(P1$TPM), 2)), pos = 2, cex = 0.8, col = colors_chr12_682652)  # pos=4 表示右侧标注
# 添加图例
legend("topright", legend = c("chr12:682652-683139", "chr12:737891-738081", "Enh = 0 (grey)"), 
        pch = c(16, 2, 16), lty=c(NA,NA,NA),lwd = 2)  # 圆点、三角形和灰色分别代表不同情况
dev.off()



