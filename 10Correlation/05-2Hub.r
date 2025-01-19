# parameters
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
####Enh_location
file1="~/FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt"
P1=fread(file1,header=T)
file2="~/FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt"
P2=fread(file2,header=T)
Pair12= c(paste(P1$promoter,P1$enhancer,sep="_"),paste(P2$genename,P2$enh_tmp,sep="_"))
Cor_pair=sub(".*?(CM.*)", "\\1", Pair12)


# 使用 lapply 函数批量处理每行数据
split_result <- lapply(Pair12, function(x) {
  # 分割字符串，按 "_CM" 切割
  parts <- strsplit(x, "_CM")[[1]]
  gene_name <- parts[1]                          # 第一部分：gene name
  promoter_info <- paste0("CM", parts[2])         # 第二部分：promoter 信息
  enhancer_info <- paste0("CM", parts[3])         # 第三部分：enhancer 信息  
  return(c(gene_name, promoter_info, enhancer_info))
})
# 将结果转为数据框
split_df <- as.data.frame(do.call(rbind, split_result), stringsAsFactors = FALSE)
colnames(split_df) <- c("GeneName", "Promoter", "Enhancer")
split_df$Pair=paste(split_df$Promoter,split_df$Enhancer,sep=":")


Final_all=split_df
length(unique(unique(Final_all)$Promoter))
length(unique(unique(Final_all)$Enhancer))


# Promoter频数数据
frequency <- table(table(unique(split_df)$Promoter))

# 生成1到30的标签
labels <- 1:20
# 计算大于20的频数总和
greater_than_30 <- sum(frequency[21:length(frequency)])
# 组合数据，将大于20的频数作为一个类别
frequency_combined <- c(frequency[1:20], greater_than_30)
frequency_combined[is.na(frequency_combined)]=0
# 设置标签，大于20的显示为 ">20"
labels_combined <- c(labels, ">20")
# 计算累积和，总和的一半和9/10的位置
cumulative_sum <- cumsum(frequency_combined)
half_total <- sum(frequency_combined) / 2
nine_tenth_total <- sum(frequency_combined) * 0.8

# 找到累积和达到总数一半和9/10的位置
half_position <- which(cumulative_sum >= half_total)[1]
nine_tenth_position <- which(cumulative_sum >= nine_tenth_total)[1]

# 创建柱状图
pdf(paste0("01Promoter_hub.pdf")) # 打开PDF设备

bar_midpoints <- barplot(frequency_combined, names.arg = labels_combined, col = "orange",
        xlab = "Enhancer # in promoter hub", ylab = "promoter hub #",
        main = "Frequency distribution of contribution of enhancer for gene",
        ylim = c(0, 2500))

# 在每个柱状图上添加对应的数目
text(bar_midpoints, frequency_combined, 
     labels = c(frequency_combined), pos = 3, cex = 0.8, col = "black")

dev.off()




#####Enhancer
library(plotrix)
# 新的频数数据
frequency <- table(table(unique(split_df)$Enhancer))
# 将大于5的频数归为一类
frequency_combined <- c(frequency[1:5], sum(frequency[6:length(frequency)]))
frequency_combined[is.na(frequency_combined)]=0
labels_combined <- c(1, 2, 3, 4, 5,">5")

# 自动计算适当的间隔
get_y_interval <- function(y_max, num_intervals) {
  rough_interval <- y_max / num_intervals
  round_interval <- round(rough_interval, -nchar(as.integer(rough_interval)) + 1)  # 取最接近的整数值
  if (round_interval == 0) round_interval <- 1  # 防止间隔为0
  return(round_interval)
}

# 计算y轴间隔和刻度
y_interval <- get_y_interval(5000, 5)
y_ticks <- seq(0, 5000, by = y_interval)

# 创建柱状图，将第一个值变为最大Y坐标的对应数
pdf(paste0("02Enhancer_hub.pdf")) # 打开PDF设备

max_frequency <- frequency_combined[1]/2
bar_midpoints <- barplot(c(max_frequency, frequency_combined[-1]), names.arg = labels_combined, col = "red",
                         xlab = "Promoter # in hub", ylab = "Enhancer hub #",
                         main = "Frequency Distribution of promoters",
                         ylim = c(0, 5000), yaxt = "n")  # 禁用默认y轴


# 添加断线和标记
axis.break(axis = 2, breakpos = 2500, style = "slash", brw = 0.02)
# 添加自定义的y轴刻度
axis(2, at = y_ticks, labels = c("0","1000","2000","6000","8000","10000"))

# 在每个柱状图上添加对应的数目
text(bar_midpoints, c(max_frequency, frequency_combined[-1]), 
     labels = c(frequency_combined), pos = 3, cex = 0.8, col = "black")
# 添加网格线
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted")

dev.off()



