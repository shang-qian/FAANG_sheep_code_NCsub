# load R packages
library(data.table)
library(parallel)
library(GenomicRanges)
library(ChIPseeker)


read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

enh_file="~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt"
enhancer_data=read_peaks(enh_file)

prom_file="~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final/06All_promoter_across_tissues.txt"
promoter_data=read_peaks(prom_file)


####Loop
# read FIMO file
fimo_data <- fread("~/FAANG/analysis/10interaction/01CTCF/01CTCF_motif/03fimo_output/fimo.tsv", header = TRUE, sep = "\t")
fimo_data=fimo_data[fimo_data$'p-value'<=1e-5]

setkey(fimo_data, sequence_name)
sequence_names <- unique(fimo_data$sequence_name)


chr_single_optimized <- function(seq_name, fimo_data) {

  df <- fimo_data[sequence_name == seq_name]
  n <- nrow(df)

  num_cores <- detectCores() - 1  

  results <- mclapply(1:(n - 1), function(i) {
  
    distances <- abs(df$start[(i + 1):n] - df$start[i])
    
    valid_indices <- which(df$strand[(i + 1):n] != df$strand[i] & distances > 1000 & distances < 2000000)
    
    if (length(valid_indices) > 0) {
      new_results <- data.table(
        sequence_name = seq_name,
        start1 = df$start[i],
        stop1 = df$stop[i],
        strand1 = df$strand[i],
        start2 = df$start[valid_indices + i],
        stop2 = df$stop[valid_indices + i],
        strand2 = df$strand[valid_indices + i],
        distance = distances[valid_indices]
      )      
      new_results[, `:=`(loop_start = pmin(start1, start2), loop_end = pmax(stop1, stop2))]     
      return(new_results)
    } else {
      return(NULL)
    }
  }, mc.cores = num_cores)
  results_final <- rbindlist(results)  
  return(results_final)
}


gr_loop_fun=function(new_loops){
gr_loops <- GRanges(
  seqnames = Rle(new_loops$sequence_name),
  ranges = IRanges(start = new_loops$loop_start, end = new_loops$loop_end),
  strand = "*", 
  loop_distance = new_loops$distance #
)
return(gr_loops)
}


library(BiocParallel)
Final_Loop_Prom_Enh <- GRanges()

process_sequence <- function(seq_name,fimo_data,promoter_data, enhancer_data) {
  print(seq_name)
  results_final <- chr_single_optimized(seq_name, fimo_data)
  
  # loop interval
  gr_loop <- unique(gr_loop_fun(results_final))
  
  # promoter
  pro_op <- findOverlaps(gr_loop, promoter_data)
  
  #### 
  number_counts <- table(queryHits(pro_op))
  unique_numbers <- as.numeric(names(number_counts[number_counts == 1]))
  loop1pro1=gr_loop[unique_numbers]
  
  # find loop1pro1 and promoter_data overlap
  loop1_pro_op <- findOverlaps(loop1pro1, promoter_data)
  
  pro1 <- unique(subjectHits(loop1_pro_op))
  
  loop_result <- GRanges()
  
  loop_result_list <- lapply(pro1, function(p) {
    tmp <- loop1pro1[which(subjectHits(loop1_pro_op) == p)]
    tmp1 <- tmp[which.min(tmp$'loop_distance')]
    mcols(tmp1) <- promoter_data[p]
    tmp1
  })
  
  loop_result <- do.call(c, loop_result_list)
  colnames(mcols(loop_result)) <- "Promoter"
  
  Final_loop_enh_overlap <- findOverlaps(loop_result, enhancer_data)
  
  Loop_pro <- loop_result[queryHits(Final_loop_enh_overlap)]
  Loop_pro$Enhancer <- enhancer_data[subjectHits(Final_loop_enh_overlap)]
  
  return(Loop_pro)
}


results_list <- list()
for (i in 1:length(sequence_names))
{
results_list[[i]]=process_sequence(sequence_names[i],fimo_data, promoter_data, enhancer_data)
}


Final_Loop_Prom_Enh <- do.call(c, results_list)

write.table(as.data.frame(Final_Loop_Prom_Enh), file=paste0("01Loop_Prom_Enh.txt"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


promoter_info <- mcols(Final_Loop_Prom_Enh)$Promoter
promoter_counts <- table(promoter_info)
promoter_counts_df <- as.data.frame(promoter_counts)
colnames(promoter_counts_df) <- c("Promoter", "Enhancer_Count")
promoter_counts_df



