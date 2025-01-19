# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)
# load R package
library(GenomicRanges)
library(rtracklayer)

input_files <- list.files(path ="~/FAANG/analysis/07Enhancer_peak/02Peak_o2", pattern = "*_peaks_Enh_20.bed", full.names = TRUE, recursive = TRUE)
print(input_files)

# read BED files and transfer to GRanges
read_bed <- function(file) {
  data <- read.table(file, header = T)
  gr <- GRanges(seqnames = data[, 1],
                ranges = IRanges(start = data[, 2], end = data[, 3]))
  return(gr)
}

gr_list <- lapply(input_files, read_bed)


# extract tissue names
extract_tissue_name <- function(file_path) {
  parts <- strsplit(file_path, "/")[[1]]
  tissue_name <- parts[length(parts) - 1]
  return(tissue_name)
}
tissue_names <- sapply(input_files, extract_tissue_name)

# define GRanges object to store overlap result
all_peaks <- do.call(c, gr_list)
peak_counts <- reduce(all_peaks, with.revmap=TRUE)

##findOverlaps
presence_matrix <- matrix(0, nrow = length(peak_counts), ncol = length(gr_list))
rownames(presence_matrix) <- paste(seqnames(peak_counts), start(peak_counts), end(peak_counts), sep = "_")
colnames(presence_matrix) <- tissue_names
# matrix
for (i in seq_along(gr_list)) {
    overlaps <- findOverlaps(peak_counts, gr_list[[i]])
    presence_matrix[queryHits(overlaps), i] <- 1
}

sum_presence <- rowSums(presence_matrix)

Total_presence=cbind(presence_matrix,sum_presence)
Total_presence <- as.data.frame(Total_presence)
common_peaks <- list()
for(i in 1:length(gr_list))
{
common_peaks[[i]]=peak_counts[Total_presence$sum_presence>=i,]
new_df=Total_presence[Total_presence$sum_presence>=i,]
mcols(common_peaks[[i]]) <- cbind(mcols(common_peaks[[i]]), new_df)

}


# calculated the length
common_peaks_lengths <- lengths(common_peaks)
# calculated log10 value
log_common_peaks_lengths <- log10(common_peaks_lengths)
names_arg <- 1:length(common_peaks_lengths)

log10_labels <- function(x) {
  parse(text = paste0("10^", x))
}

# save pdf file
pdf("01Merge_common_peaks_lengths.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common Peaks", 
                         xlab = "Consistant tissue # >=", 
                         ylab = "Number of Peaks", 
                         col = "orange", 
                         names.arg = names_arg, 
                         ylim = c(0, 6), 
                         yaxt = "n")

# add  labels
axis(2, at = 0:6, labels = log10_labels(0:6))
#abline(h = 5, col = "red", lty = 2)

text(bar_midpoints, log_common_peaks_lengths + 0.01, labels = common_peaks_lengths, cex = 0.7, pos = 3)
dev.off()


#####plot cutoff
data <- common_peaks_lengths
# rate of each stats
rate_of_change <- (data[-length(data)]-data[-1])/data[-1]
# calculate the max rate
max_change_index <- which(rate_of_change<2)[1]
# obtain cutoff value
cutoff <- data[2]
# save PDF 
pdf("02common_peaks_20tissues_lengths.pdf", width = 10, height = 6)
plot(data, type = "o", col = "orange", xlab = "Tissue # ", ylab = "Number of peaks", main = "Original Data Points with Cutoff")
text(x = 1:20, y = data, labels = data, pos = 3, cex = 1, col = "black")

abline(v = 2, col = "red", lty =2)
#legend("topright", legend = paste("Cutoff:", cutoff), col = "red", lty = 2)
dev.off()




##extract the common peaks at least two samples
all_overlaps <- common_peaks[[2]]
final_intersections <- all_overlaps
df_intersections <- as.data.frame(final_intersections)
df_intersections <-  subset(df_intersections, select = -revmap)
#df_intersections$revmap <- sapply(df_intersections$revmap, function(x) paste(x, collapse = ","))
write.table(df_intersections, file = "03All27_tissues_01raw_tis2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


