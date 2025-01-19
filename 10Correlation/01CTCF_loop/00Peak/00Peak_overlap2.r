args <- commandArgs(trailingOnly = TRUE)

# load R package
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)

args=c("~/FAANG/analysis/04CTCF/01Peak/02CER_F1/bowtie2/mergedLibrary/macs2/narrowPeak/CER_F1_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02CER_F2/bowtie2/mergedLibrary/macs2/narrowPeak/CER_F2_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02CER_M1/bowtie2/mergedLibrary/macs2/narrowPeak/CER_M1_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02CER_M2/bowtie2/mergedLibrary/macs2/narrowPeak/CER_M2_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02LIV_F2/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_F2_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02LIV_M1/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_M1_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02LIV_M2/bowtie2/mergedLibrary/macs2/narrowPeak/LIV_M2_CTCF_peaks.narrowPeak",
"~/FAANG/analysis/04CTCF/01Peak/02SPL_M2/bowtie2/mergedLibrary/macs2/narrowPeak/SPL_M2_CTCF_peaks.narrowPeak")
print(args)

# read files
read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

granges_list <- lapply(args[1:8], read_peaks)

# check data input
if (any(sapply(granges_list, is.null))) {
  stop("Error in reading or converting peak files to GRanges")
}

# overlap
compute_filtered_overlap <- function(gr1, gr2) {
  overlap <- intersect(gr1, gr2)
  Enh_gr1=overlap[width(overlap)>=5]
  return(Enh_gr1)
}

overlap2=list()
count=1
for (i in 1:7)
{ 
 print(i)
 for (j in (i+1):8)
 { print(paste("Comb.:",i,"-",j))
   
   overlap2[[count]]=compute_filtered_overlap(granges_list[[i]], granges_list[[j]])
   count=count+1
 
 }
}
comb_gr=do.call(c, overlap2)
merged_gr <- reduce(comb_gr)
merged_gr_filter=merged_gr[grepl("^CM",seqnames(merged_gr))]

write.table(as.data.frame(merged_gr_filter), file=paste0("01Peaks_CTCF.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


