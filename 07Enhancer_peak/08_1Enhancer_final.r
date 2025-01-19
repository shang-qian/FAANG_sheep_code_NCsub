# Rscript paramters
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(GenomicRanges)
#library(ChIPseeker)
#library(dplyr)

#all peak files
input_files<- list.files(path = "~/FAANG/analysis/07Enhancer_peak/06Peak_tissue_high/01All", pattern = "*02final_peak_union_enh.bed", full.names = TRUE, recursive = TRUE)
###Tissue name
tissue_names=unlist(lapply(strsplit(basename(input_files),"_"), function(x) x[1]))

#read peak files transfer to GRanges
read_peaks <- function(file) {
  prom_one <- fread(file, header = F)
  return(prom_one)
}

Prom_list <- lapply(input_files, read_peaks)    

All_prom=do.call(rbind, Prom_list)
colnames(All_prom)=c("chr","start","end","width","strand")

gr <- GRanges(seqnames = All_prom$chr,
              ranges = IRanges(start = All_prom$start, end = All_prom$end),
              strand = All_prom$strand)

reduced_gr <- reduce(gr)

#### one-enhancer
library(parallel)

## Step 1: reduced_gr to data.table
process_reduce_enh <- function(one_enh) {
  selected_start <- median(start(one_enh))
  selected_end <- median(end(one_enh))
    
  one_index=which.min(abs(start(one_enh) - selected_start) + abs(end(one_enh) - selected_end))
  Enh_r1=one_enh[one_index]
return(Enh_r1)
}
## 

process_region <- function(i, reduced_gr, gr) {
  results_list <- GRanges()
  
  current_region <- reduced_gr[i]
  matching_one <- findOverlaps(current_region, gr)
  one_enh <- gr[subjectHits(matching_one)]
  Enh_r1 <- process_reduce_enh(one_enh)
  
  results_list <- c(results_list, Enh_r1)
  
  rem <- findOverlaps(Enh_r1, one_enh)
  other_rows <- one_enh[-subjectHits(rem)]
  
  if (length(other_rows) > 0) {
    row2 <- reduce(other_rows)
    for (j in seq_along(row2)) {
      row2_overlap <- findOverlaps(row2[j], other_rows)
      one_row2 <- other_rows[subjectHits(row2_overlap)]
      Enh_r2 <- process_reduce_enh(one_row2)
      results_list <- c(results_list, Enh_r2)
    }
  }
  
  return(results_list)
}

# core number
num_cores <- detectCores() - 1  
# treat each core
results_list <- mclapply(seq_along(reduced_gr), process_region, reduced_gr = reduced_gr, gr = gr, mc.cores = num_cores)

# merge result from GRanges
final_results <- do.call(c, results_list)

write.table(as.data.frame(final_results), file = paste0("01final_enhancer_total.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

final2000=final_results[which(width(final_results)<=2000)] 
#final3001=final_results[which(width(final_results)>3000)]
write.table(as.data.frame(final3000), file = paste0("02final_enhancer_total2000.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

