# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)

# Load R packages
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)

print(args)

# function read files and transfer to GRanges
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

Tisname=args[length(args)]

granges_list <- lapply(args[1:4], read_peaks)                       
H3K4me3 <- granges_list[[1]]
H3K27ac <- granges_list[[2]]
H3K4me3P5 <-granges_list[[3]] 
CAGE <- granges_list[[4]] 
overlap_annotated=(c(H3K4me3,H3K27ac,CAGE,H3K4me3P5))


sorted_overlap_annotated <- sort(overlap_annotated)
unique_sorted_overlap_annotated <- unique(sorted_overlap_annotated)
print(unique_sorted_overlap_annotated)

write.table(as.data.frame(unique_sorted_overlap_annotated), file = paste0(Tisname,"/",Tisname,"_02final_peak_union.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




