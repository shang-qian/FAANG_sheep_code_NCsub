# Rscript parameters
args <- commandArgs(trailingOnly = TRUE)
# load R packages
library(GenomicRanges)
library(rtracklayer)

#Input
input_files <- list.files(path ="~/FAANG/analysis/08Promoter_peak/02Peak_o2", pattern = "*_peaks_Pro_20.bed", full.names = TRUE, recursive = TRUE)
print(input_files)

# read peak file function
read_bed <- function(file) {
  data <- read.table(file, header = T)
  gr <- GRanges(seqnames = data[, 1],
                ranges = IRanges(start = data[, 2], end = data[, 3]))
  return(gr)
}

gr_list <- lapply(input_files, read_bed)
all_peaks <- do.call(c, gr_list)
peak_counts <- reduce(all_peaks, with.revmap=TRUE)

# tissue names
extract_tissue_name <- function(file_path) {
  parts <- strsplit(basename(file_path), "_")[[1]]
  tissue_name <- parts[1]
  return(tissue_name)
}
tissue_names <- sapply(input_files, extract_tissue_name)

##findOverlaps
presence_matrix <- matrix(0, nrow = length(peak_counts), ncol = length(gr_list))
rownames(presence_matrix) <- paste(seqnames(peak_counts), start(peak_counts), end(peak_counts), sep = "_")
colnames(presence_matrix) <- tissue_names
# fill matrix with identified peak
for (i in seq_along(gr_list)) {
    overlaps <- findOverlaps(peak_counts, gr_list[[i]])
    presence_matrix[queryHits(overlaps), i] <- 1
}

sum_presence <- rowSums(presence_matrix)

Total_presence=cbind(presence_matrix,sum_presence)
Total_presence <- as.data.frame(Total_presence)

###descript the overlap tissues peaks
common_peaks <- list()
for(i in 1:20)
{
common_peaks[[i]]=peak_counts[Total_presence$sum_presence>=i,]
new_df=Total_presence[Total_presence$sum_presence>=i,]

mcols(common_peaks[[i]]) <- cbind(mcols(common_peaks[[i]]), new_df)
}




# calculate peak length
common_peaks_lengths <- lengths(common_peaks)
# calculate log10 value
log_common_peaks_lengths <- log10(common_peaks_lengths)
names_arg <- 1:length(common_peaks_lengths)
# set labels
log10_labels <- function(x) {
  parse(text = paste0("10^", x))
}
# save PDF file
pdf("01Merge_peaks_20tissues_log.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common Peaks", 
                         xlab = "Consistent tissue # >=", 
                         ylab = "Number of Peaks", 
                         col = "blue", 
                         names.arg = names_arg, 
                         ylim = c(0, log10(10^6)), 
                         yaxt = "n")

# add Y lables
axis(2, at = 0:6, labels = log10_labels(0:6))
#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, log_common_peaks_lengths + 0.1, labels = common_peaks_lengths, cex = 0.7, pos = 3)
dev.off()

#####plot cutoff
data <- common_peaks_lengths
# calculate rate change
rate_of_change <- (data[-length(data)]-data[-1])/data[-1]
max_change_index <- length(rate_of_change[rate_of_change>1])
# cutoff
cutoff <- data[2]

pdf("02Merge_peaks_20tissues_cutoff3.pdf", width = 10, height = 6)
plot(data, type = "o", col = "blue", xlab = "Tissue #", ylab = "Number of peak", main = "Original Data Points with Cutoff: 3tissues")
#text(x = 2, y = data[2], labels = data[2], pos = 3, cex = 1.8, col = "black")
text(x = 1:20, y = data, labels = data, pos = 3, cex = 1, col = "black")
abline(v = 2, col = "red", lty = 2)
#legend("topright", legend = paste("Cutoff:", cutoff), col = "red", lty = 2)
dev.off()




##extract the common peaks at least two samples
all_overlaps <- common_peaks[[2]]
final_intersections <- all_overlaps

# transfer GRanges to data.frame
df_intersections <- as.data.frame(final_intersections)
df_intersections <-  subset(df_intersections, select = -revmap)
#df_intersections$revmap <- sapply(df_intersections$revmap, function(x) paste(x, collapse = ","))
write.table(df_intersections, file = "03All20_tissues_01raw_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


