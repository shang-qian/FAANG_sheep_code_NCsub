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

if (any(sapply(granges_list, is.null))) {
  stop("Error in reading or converting peak files to GRanges")
}

Peak <- granges_list[[1]]
RE <- granges_list[[2]]

Enh0=intersect(Peak,RE)
Enh <- reduce(Enh0)
Enh50=Enh[which(width(Enh)>=10)]


#export(final_output, "final_combined_overlap.bed")
dir.create(args[3])
write.table(as.data.frame(Enh50), file=paste0(args[3],"/01",args[3],"_enhancer_",args[4],".bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



