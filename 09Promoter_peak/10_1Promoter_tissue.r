# Rscripts parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

#library(rtracklayer)
library(data.table)
library(GenomicRanges)
library(ChIPseeker)

###Tissue name
Tisname=args[length(args)]

read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

Prom_A <- fread(args[1], header = T)
Prom_A_gr <- GRanges(seqnames = Prom_A$Gseqnames,
                     ranges = IRanges(start = Prom_A$start, end = Prom_A$end),
                     strand = Prom_A$Gstrand)

File_list <- lapply(args[2:4], read_peaks)    
H3K4me3=File_list[[1]]
H3K27ac=File_list[[2]]
CAGE=File_list[[3]]

# overlap of Prom_A and other datasets
overlap_H3K4me3 <- findOverlaps(Prom_A_gr, H3K4me3)
overlap_H3K27ac <- findOverlaps(Prom_A_gr, H3K27ac)
overlap_CAGE <- findOverlaps(Prom_A_gr, CAGE)

Prom_A$H3K4me3_overlap <- 0
Prom_A$H3K27ac_overlap <- 0
Prom_A$CAGE_overlap <- 0

Prom_A$H3K4me3_overlap[queryHits(overlap_H3K4me3)] <- 1
Prom_A$H3K27ac_overlap[queryHits(overlap_H3K27ac)] <- 1
Prom_A$CAGE_overlap[queryHits(overlap_CAGE)] <- 1

Prom_A$total_overlap <- Prom_A$H3K4me3_overlap + Prom_A$H3K27ac_overlap + Prom_A$CAGE_overlap
Tissue_Prom_A <- Prom_A[Prom_A$total_overlap > 0, ]

Tissue_Prom_anno <- Tissue_Prom_A[!grepl("STRG",Tissue_Prom_A$geneID),]
Tissue_Prom_new <- Tissue_Prom_A[grepl("STRG",Tissue_Prom_A$geneID),]

write.table(as.data.frame(Tissue_Prom_A), file = paste0(Tisname,"/",Tisname,"_01All_promoter_peaks.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(Tissue_Prom_anno), file = paste0(Tisname,"/",Tisname,"_02Annotate_promoter_peaks.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(as.data.frame(Tissue_Prom_new), file = paste0(Tisname,"/",Tisname,"_03New_promoter_peaks.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
