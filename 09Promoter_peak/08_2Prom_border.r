# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)

# Load R packages
#library(rtracklayer)
###########TPM10_new TSS
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ChIPseeker)

print(args)

TTSS <- fread(args[1], header = F)
colnames(TTSS) <- c("Gseqnames","Gstart","Gend","Gwidth","Gstrand","geneID","TSS")


# peak read
peaks <- readPeakFile(args[2])
peak_counts<- makeGRangesFromDataFrame(peaks)

# Process the part with explicit "+" and "-" strands
TTSS$tss_start <- ifelse(TTSS$Gstrand == "+", pmax(as.numeric(TTSS$TSS)-3000,0), as.numeric(TTSS$TSS))  
TTSS$tss_end <- ifelse(TTSS$Gstrand == "+", as.numeric(TTSS$TSS), as.numeric(TTSS$TSS)+3000)  

tss_gr <- GRanges(
  seqnames = TTSS$Gseqnames,
  ranges = IRanges(start = TTSS$tss_start, end = TTSS$tss_end),
  strand = TTSS$Gstrand
)

data_gr=peak_counts[which(width(peak_counts)<=2000)]

overlaps <- findOverlaps(data_gr, tss_gr)# Extract the overlapping or near-matching peaks

gene_list=unique(TTSS[subjectHits(overlaps),]$geneID)

# each gene
library(parallel)
process_gene <- function(gene_id, data_gr, TTSS, overlaps) {
  overlap_TSS <- TTSS[subjectHits(overlaps),]
  gene_pos <- which(overlap_TSS$geneID == gene_id)
  tmp_TSS <- overlap_TSS[gene_pos]
  tmp_peak <- data_gr[queryHits(overlaps)[gene_pos]]  
  # 
  combined_df <- cbind(tmp_TSS, as.data.frame(tmp_peak)) 
  # calculate distance of TSS and peak，+- strand
  combined_df$distance <- ifelse(combined_df$Gstrand == "+",
                                 abs(combined_df$end - combined_df$TSS),
                                 abs(combined_df$start - combined_df$TSS))  
  # closest promoter
  selected_prom <- combined_df[which.min(combined_df$distance),]  
  return(selected_prom)
}

 
result_list <- mclapply(gene_list, function(gene_id) {process_gene(gene_id, data_gr, TTSS, overlaps)}, mc.cores = detectCores())
  
Gene_promoter=do.call(rbind, result_list)

Final_promoter=Gene_promoter[,c(1:7,15,11:13)]

write.table(as.data.frame(Final_promoter), file = paste0(args[3],"_01final_promoter.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






