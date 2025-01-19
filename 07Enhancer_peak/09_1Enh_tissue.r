# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)

print(args)

#library(rtracklayer)
library(data.table)
library(GenomicRanges)
library(ChIPseeker)

Tisname=args[length(args)]

# read peak files transfer to GRanges format
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

Prom_A <- fread(args[1], header = T)
Prom_A_gr <- GRanges(seqnames = Prom_A$seqnames,
                     ranges = IRanges(start = Prom_A$start, end = Prom_A$end),
                     strand = Prom_A$strand)

File_list <- lapply(args[2:4], read_peaks)    
H3K4me3=File_list[[1]]
H3K27ac=File_list[[2]]
CAGE=File_list[[3]]

# calculate overlap of Prom_A and other dataset
overlap_H3K4me3 <- findOverlaps(Prom_A_gr, H3K4me3)
overlap_H3K27ac <- findOverlaps(Prom_A_gr, H3K27ac)
overlap_CAGE <- findOverlaps(Prom_A_gr, CAGE)

Prom_A$H3K4me3_overlap <- 0
Prom_A$H3K27ac_overlap <- 0
Prom_A$CAGE_overlap <- 0

# mark 1 as the identifed peak
Prom_A$H3K4me3_overlap[queryHits(overlap_H3K4me3)] <- 1
Prom_A$H3K27ac_overlap[queryHits(overlap_H3K27ac)] <- 1
Prom_A$CAGE_overlap[queryHits(overlap_CAGE)] <- 1

# total overlap
Prom_A$total_overlap <- Prom_A$H3K4me3_overlap + Prom_A$H3K27ac_overlap + Prom_A$CAGE_overlap
Tissue_Prom_A <- Prom_A[Prom_A$total_overlap > 0, ]

write.table(as.data.frame(Tissue_Prom_A), file = paste0(Tisname,"/",Tisname,"_01All_enhancer.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

