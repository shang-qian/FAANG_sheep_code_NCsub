# Rscript paprameters
args <- commandArgs(trailingOnly = TRUE)

# load R packages
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)

print(args)

# read peak file and transfer GRanges
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

# read files
granges_list <- lapply(args[1:4], read_peaks)
Enh_A <- granges_list[[1]]                         
H3K4me1 <- granges_list[[2]]
H3K27ac <- granges_list[[3]]
CAGE <- granges_list[[4]]


O_H3K4me1=findOverlaps(H3K4me1,Enh_A)
Tissue_H3K4me1=H3K4me1[queryHits(O_H3K4me1)]

O_H3K27ac=findOverlaps(H3K27ac,Enh_A)
Tissue_H3K27ac=H3K27ac[queryHits(O_H3K27ac)]

O_CAGE=findOverlaps(CAGE,Enh_A)
Tissue_CAGE=CAGE[queryHits(O_CAGE)]

#output
Tisname=args[5]
H3K4me1_annotated <- Tissue_H3K4me1[width(Tissue_H3K4me1) >=20]
print(H3K4me1_annotated)
write.table(as.data.frame(H3K4me1_annotated), file = paste0(Tisname,"/",Tisname,"_00H3K4me1.Propeak"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

H3K27ac_annotated <- Tissue_H3K27ac[width(Tissue_H3K27ac) >=20]
print(H3K27ac_annotated)
write.table(as.data.frame(H3K27ac_annotated), file = paste0(Tisname,"/",Tisname,"_00H3K27ac.Propeak"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

CAGE_annotated <- Tissue_CAGE[width(Tissue_CAGE) >=20]
print(CAGE_annotated)
write.table(as.data.frame(CAGE_annotated), file = paste0(Tisname,"/",Tisname,"_00CAGE.Propeak"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#output
overlap_annotated=reduce(c(H3K4me1_annotated,H3K27ac_annotated,CAGE_annotated))
write.table(as.data.frame(overlap_annotated), file = paste0(Tisname,"/",Tisname,"_01final1.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




