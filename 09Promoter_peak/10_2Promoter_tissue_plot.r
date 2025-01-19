# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)
#library(rtracklayer)
library(data.table)
library(GenomicRanges)
#library(ChIPseeker)

args=c("~/FAANG/analysis/08Promoter_peak/09Promoter_Final/01final_promoter_gene_TSS.txt")

print(args)
Prom_A <- fread(args[1], header = T)

####
library(ggplot2)
Prom_A$adjusted_distance <- ifelse(Prom_A$Gstrand == "+",
                                   Prom_A$end - Prom_A$TSS,
                                   Prom_A$start - Prom_A$TSS)

# plot adjusted_distance histogram
distance_plot <- ggplot(Prom_A, aes(x = adjusted_distance)) +
  geom_histogram(binwidth = 50, fill = "orange", color = "black", alpha = 0.7) +
  labs(title = "Distance to TSS", x = "Distance to TSS", y = "Frequency") +
  theme_minimal()
ggsave("04Distance_distribution.pdf", plot = distance_plot, device = "pdf", width = 8, height = 6)

# plot width histogram
width_plot <- ggplot(Prom_A, aes(x = width)) +
  geom_histogram(binwidth = 50, fill = "orange", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Width", x = "Width", y = "Frequency") +
  theme_minimal()
ggsave("05width_distribution.pdf", plot = width_plot, device = "pdf", width = 8, height = 6)


# read all peak files
input_files<- list.files(path = "~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final", pattern = "*_01All_promoter_peaks.txt", full.names = TRUE, recursive = TRUE)
###Tissue name
Tisname=sapply(strsplit(basename(input_files),"_"),function(x)x[1])

read_peaks <- function(file) {
  one <- fread(file, header = T)
  gr <- GRanges(seqnames = one$Gseqnames,
                     ranges = IRanges(start = one$start, end = one$end),
                     strand = one$Gstrand)
  return(gr)
}


File_list <- lapply(input_files, read_peaks)    
Combine_granges <- do.call(c, File_list)
Uniq_com=unique(Combine_granges)

overlap_matrix <- matrix(0, nrow = length(Uniq_com), ncol = length(File_list))

# File_list GRanges ，overlap with Uniq_com
for (i in seq_along(File_list)) {
  overlaps <- findOverlaps(Uniq_com, File_list[[i]])
  overlap_matrix[queryHits(overlaps), i] <- 1
}

# add one column: total number
overlap_df <- as.data.frame(overlap_matrix)
colnames(overlap_df) <- Tisname
overlap_df$total_tissues <- rowSums(overlap_df)

# overlap_df,Uniq_com merge
final_df <- cbind(as.data.frame(Uniq_com), overlap_df)
write.table(final_df, file = paste0("06All_promoter_across_tissues.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


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
# 
names_arg <- 1:length(common_peaks_lengths)
# labels
log10_labels <- function(x) {
  parse(text = paste0(x,"x10^3"))
}
# save PDF file
pdf("07Promoter_24tissues_share.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common Promoters", 
                         xlab = "Consistent tissue # >=", 
                         ylab = "Number of Promoters (x10^3)", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 30000), 
                         yaxt = "n")

# add Y labels
axis(2, at = c(0,5000,10000,15000,20000,25000,30000), labels = c("0","5","10","15","20","25","30"))
#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, log_common_peaks_lengths + 0.1, labels = common_peaks_lengths, cex = 0.7, pos = 3)
dev.off()

names_arg <- 1:length(common_peaks_lengths)
# save PDF file
pdf("08Promoter_24tissues_share_1.pdf", width = 10, height = 6)
bar_midpoints <- barplot(common_promo_lengths, 
                         main = "Common Promoters", 
                         xlab = "Consistent tissue #", 
                         ylab = "Number of Promoters", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 6000), 
                         yaxt = "n")

# add Y labels
axis(2, at = c(0,1000,2000,4000,6000), labels = c("0","1000","2000","4000","6000"))

#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, common_promo_lengths + 0.1, labels = common_promo_lengths, cex = 0.7, pos = 3)
dev.off()


###New promoter
input_files<- list.files(path = "~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final", pattern = "*_03New_promoter_peaks.txt", full.names = TRUE, recursive = TRUE)
###Tissue name
Tisname=sapply(strsplit(basename(input_files),"_"),function(x)x[1])

read_peaks <- function(file) {
  one <- fread(file, header = T)
  gr <- GRanges(seqnames = one$Gseqnames,
                     ranges = IRanges(start = one$start, end = one$end),
                     strand = one$Gstrand)
  return(gr)
}



File_list <- lapply(input_files, read_peaks)    
Combine_granges <- do.call(c, File_list)
Uniq_com=unique(Combine_granges)

overlap_matrix <- matrix(0, nrow = length(Uniq_com), ncol = length(File_list))

for (i in seq_along(File_list)) {
  overlaps <- findOverlaps(Uniq_com, File_list[[i]])
  overlap_matrix[queryHits(overlaps), i] <- 1
}

overlap_df <- as.data.frame(overlap_matrix)
colnames(overlap_df) <- Tisname
overlap_df$total_tissues <- rowSums(overlap_df)

# overlap_df Uniq_com merge
final_df <- cbind(as.data.frame(Uniq_com), overlap_df)



###descript the overlap tissues peaks
common_peaks <- list()
common_peaks_lengths=common_promo_lengths <- numeric()
for(i in 1:24)
{
common_peaks[[i]]=final_df[as.numeric(final_df$total_tissues)>=i,]
common_peaks_lengths[i]=nrow(common_peaks[[i]])
common_promo_lengths[i]=nrow(final_df[as.numeric(final_df$total_tissues)==i,])
}

names_arg <- 1:length(common_peaks_lengths)
# save PDF
pdf("09New_Promoter_24tissues_share.pdf", width = 10, height = 6)
bar_midpoints <- barplot(common_promo_lengths, 
                         main = "Common New promoters", 
                         xlab = "Consistent tissue #", 
                         ylab = "Number of New Promoters", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 80), 
                         yaxt = "n")

# add Y labels
axis(2, at = c(0,20,40,60,80), labels = c("0","20","40","60","80"))

#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, common_promo_lengths + 0.1, labels = common_promo_lengths, cex = 0.7, pos = 3)
dev.off()

