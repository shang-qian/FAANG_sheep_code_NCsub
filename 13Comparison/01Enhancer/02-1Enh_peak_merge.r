args <- commandArgs(trailingOnly = TRUE)

library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)

print(args)

read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}


granges_list <- lapply(args[1:2], read_peaks)

H3K4me1 <- granges_list[[1]]
H3K27ac <- granges_list[[2]]



compute_union <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
  gr1_unique <- gr1[-queryHits(overlaps)]
  gr2_unique <- gr2[-subjectHits(overlaps)]
  intersected <- gr1[queryHits(overlaps)]  
  union_result <- c(gr1_unique, intersected, gr2_unique)
  return(union_result)
}


final_output <- compute_union(H3K4me1, H3K27ac)


final_output <- reduce(final_output)

#export(final_output, "final_combined_overlap.bed")
dir.create(args[3])
write.table(as.data.frame(final_output), file=paste0(args[3],"/",args[3],"_peaks_Enh_",args[4],".bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



