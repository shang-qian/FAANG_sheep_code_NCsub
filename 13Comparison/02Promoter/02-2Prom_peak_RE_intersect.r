args <- commandArgs(trailingOnly = TRUE)

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


Peak <- granges_list[[1]]
RE <- granges_list[[2]]

tmp=findOverlaps(RE,Peak)

Promoter=unique(RE[queryHits(tmp)])

dir.create(args[3])
write.table(as.data.frame(Promoter), file=paste0(args[3],"/01",args[3],"_promoter_",args[4],".bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



