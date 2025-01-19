# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)

args=c("~/FAANG/analysis/07Enhancer_peak/08Enh/02final_enhancer_total2000.txt")
library(data.table)

print(args)
Prom_A <- fread(args[1], header = T)

####
library(ggplot2)

# plot width bar
width_plot <- ggplot(Prom_A, aes(x = width)) +
  geom_histogram(binwidth = 50, fill = "orange", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Width", x = "Width", y = "Frequency") +
  theme_minimal()
# save width plot
ggsave("02width_distribution.pdf", plot = width_plot, device = "pdf", width = 8, height = 6)

library(GenomicRanges)
#library(ChIPseeker)
input_files<- list.files(path = "~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final", pattern = "*_01All_enhancer.txt", full.names = TRUE, recursive = TRUE)
###Tissue name
Tisname=sapply(strsplit(basename(input_files),"_"),function(x)x[1])

# read peak file transfer GRanges
read_peaks <- function(file) {
  one <- fread(file, header = T)
  gr <- GRanges(seqnames = one$seqnames,
                     ranges = IRanges(start = one$start, end = one$end),
                     strand = one$strand)
  return(gr)
}


File_list <- lapply(input_files, read_peaks)    
Combine_granges <- do.call(c, File_list)
Uniq_com=unique(Combine_granges)

# New matri use to save result，row is Uniq_com region，colume is one GRanges of File_list
overlap_matrix <- matrix(0, nrow = length(Uniq_com), ncol = length(File_list))

# Traversal File_list's GRanges，calculae Uniq_com overlap
for (i in seq_along(File_list)) {
  overlaps <- findOverlaps(Uniq_com, File_list[[i]])
  overlap_matrix[queryHits(overlaps), i] <- 1
}

# overlap_matrix transfer to data frame，add one column for total 
overlap_df <- as.data.frame(overlap_matrix)
colnames(overlap_df) <- Tisname
overlap_df$total_tissues <- rowSums(overlap_df)

# merge overlap_df and Uniq_com
final_df <- cbind(as.data.frame(Uniq_com), overlap_df)
write.table(final_df, file = paste0("03All_enhancer_across_tissues.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


###descript the overlap tissues peaks
common_peaks <- list()
common_peaks_lengths=common_promo_lengths <- numeric()
for(i in 1:24)
{
common_peaks[[i]]=final_df[as.numeric(final_df$total_tissues)>=i,]
common_peaks_lengths[i]=nrow(common_peaks[[i]])
common_promo_lengths[i]=nrow(final_df[as.numeric(final_df$total_tissues)==i,])
}


log_common_peaks_lengths <- common_peaks_lengths

names_arg <- 1:length(common_peaks_lengths)

log10_labels <- function(x) {
  parse(text = paste0(x,"x10^3"))
}

pdf("04Enhancer_24tissues_share.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common Enhancer", 
                         xlab = "Consistent tissue # >=", 
                         ylab = "Number of Enhancer", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 350000), 
                         yaxt = "n")

# add Y lables
axis(2, at = c(0,50000,100000,150000,200000,250000,300000,350000), labels = c("0","50000","100000","150000","200000","250000","300000","350000"))
#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, log_common_peaks_lengths + 0.1, labels = common_peaks_lengths, cex = 0.7, pos = 3)
dev.off()


names_arg <- 1:length(common_peaks_lengths)
# save pdf
pdf("05Enhancer_24tissues_share_1.pdf", width = 10, height = 6)
bar_midpoints <- barplot(common_promo_lengths, 
                         main = "Common Enhancer", 
                         xlab = "Consistent tissue #", 
                         ylab = "Number of Enhancer", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 40000), 
                         yaxt = "n")

# add Y 
axis(2, at = c(0,5000,10000,20000,30000,40000), labels = c("0","5000","10000","20000","30000","40000"))

#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, common_promo_lengths + 0.1, labels = common_promo_lengths, cex = 0.7, pos = 3)
dev.off()





