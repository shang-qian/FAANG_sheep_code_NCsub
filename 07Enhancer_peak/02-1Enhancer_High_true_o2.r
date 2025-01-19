args <- commandArgs(trailingOnly = TRUE)

# load R package
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)


print(args)

# read peak file and tranfer to GRanges object
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

granges_list <- lapply(args[1:4], read_peaks)

# checke the files
if (any(sapply(granges_list, is.null))) {
  stop("Error in reading or converting peak files to GRanges")
}

# count peak
peak_counts <- sapply(granges_list, length)
peak_count_file <- paste0(args[5],"/",args[5], "_peak_counts.txt")
peak_count_output=data.frame(File = args[1:4], Peak_Count = peak_counts)
colnames(peak_count_output)=c("File",args[5])
write.table(peak_count_output, file = peak_count_file, sep = "\t", row.names = FALSE, quote = FALSE)

H3K4me1 <- granges_list[[3]]
H3K27ac <- granges_list[[4]]
CAGE <- granges_list[[1]]
ATAC <- granges_list[[2]]


# calculate overlap
compute_filtered_overlap <- function(gr1, gr2, lp) {
  overlap <- intersect(gr1, gr2)
  Enh_gr1=overlap[width(overlap)>=1]
  return(Enh_gr1)
}

# filter overlap
filtered_H3K4me1_H3K27ac <- compute_filtered_overlap(H3K4me1, H3K27ac, args[6])
filtered_H3K4me1_CAGE <- compute_filtered_overlap(H3K4me1, CAGE, args[6])
filtered_H3K4me1_ATAC <- compute_filtered_overlap(H3K4me1, ATAC, args[6])
filtered_H3K27ac_CAGE <- compute_filtered_overlap(H3K27ac, CAGE, args[6])

# any two overlap
compute_union <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
  gr1_unique <- gr1[-queryHits(overlaps)]
  gr2_unique <- gr2[-subjectHits(overlaps)]
  intersected <- intersect(gr1[queryHits(overlaps)], gr2[subjectHits(overlaps)])  
  union_result <- c(intersected, gr1_unique, gr2_unique)
  return(union_result)
}

# union any two overlaps
union_H3K4me1_H3K27ac_CAGE <- compute_union(filtered_H3K4me1_H3K27ac, filtered_H3K4me1_CAGE)
union_H3K4me1_H3K27ac_ATAC <- compute_union(filtered_H3K4me1_H3K27ac, filtered_H3K4me1_ATAC)
union_H3K4me1_H3K27ac_CAGE2 <- compute_union(filtered_H3K4me1_H3K27ac, filtered_H3K27ac_CAGE)
union_H3K4me1_CAGE_ATAC <- compute_union(filtered_H3K4me1_CAGE, filtered_H3K4me1_ATAC)
union_H3K4me1_CAGE_CAGE2 <- compute_union(filtered_H3K4me1_CAGE, filtered_H3K27ac_CAGE)
union_H3K27ac_CAGE_ATAC <- compute_union(filtered_H3K4me1_ATAC,filtered_H3K27ac_CAGE)

# merge all 
final_output <- c(union_H3K4me1_H3K27ac_CAGE, union_H3K4me1_H3K27ac_ATAC, union_H3K4me1_H3K27ac_CAGE2, union_H3K4me1_CAGE_ATAC, union_H3K4me1_CAGE_CAGE2, union_H3K27ac_CAGE_ATAC)

# remve duplicates
final_output <- reduce(final_output)

write.table(as.data.frame(final_output), file=paste0(args[5],"/",args[5],"_peaks_Enh_",args[6],".bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



