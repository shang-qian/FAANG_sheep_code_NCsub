args <- commandArgs(trailingOnly = TRUE)

# 加载必要的R包
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)


print(args)

# 函数：读取peak文件并转换为GRanges对象
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

# 读取和转换peak文件
granges_list <- lapply(args[1:2], read_peaks)

# 检查是否成功读取和转换
if (any(sapply(granges_list, is.null))) {
  stop("Error in reading or converting peak files to GRanges")
}


H3K4me1 <- granges_list[[1]]
H3K27ac <- granges_list[[2]]


# 定义函数计算两两并集
compute_union <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
  gr1_unique <- gr1[-queryHits(overlaps)]
  gr2_unique <- gr2[-subjectHits(overlaps)]
  intersected <- gr1[queryHits(overlaps)]  
  union_result <- c(gr1_unique, intersected, gr2_unique)
  return(union_result)
}

# 合并所有结果
final_output <- compute_union(H3K4me1, H3K27ac)

# 去重并确保没有间隙
final_output <- reduce(final_output)

#export(final_output, "final_combined_overlap.bed")
dir.create(args[3])
write.table(as.data.frame(final_output), file=paste0(args[3],"/",args[3],"_peaks_Promoter_",args[4],".bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



