
library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

#####Annotated file read GFF
gff_file <- "/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff"
gff_data <- import.gff(gff_file)
# select gene from annotation
genes <- gff_data[gff_data$type %in% c("gene")|gff_data$type %in% c("pseudogene")]
#
Gene_data <- data.frame(
  seqnames = seqnames(genes),
  start = start(genes),
  end = end(genes),
  strand = strand(genes),
  geneID = mcols(genes)$ID  
)

# extract TSS positon function
extract_tss <- function(start, end, strand) {
  if (strand == "+") {
    return(start)
  } else if (strand == "-") {
    return(end)
  } else {
    return(NA)  # if no strand information，return NA
  }
}

Gene_data$TSS <- mapply(extract_tss, Gene_data$start, Gene_data$end, as.character(Gene_data$strand))
write.table(as.data.frame(Gene_data), file = paste0("01TSS_annotate.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# read peak
input_files<- list.files(path = "~/FAANG/analysis/08RNA_TPM/02Promoter/02TPM", pattern = "*_TPM_gene.txt", full.names = TRUE, recursive = F)

###Tissue name
tissue_names=unlist(lapply(strsplit(basename(input_files),"_"), function(x) x[1]))

AnTSS_gr <- GRanges(
  seqnames = Gene_data$seqnames,
  ranges = IRanges(start = Gene_data$start, end = Gene_data$end),
  strand = Gene_data$strand,
  geneID = Gene_data$geneID,
  TSS = Gene_data$TSS
)

# function: read bed files
read_bed <- function(file, AnTSS_gr) {
  df <- read.table(file, header = TRUE, sep = "\t")
  df$Strand[df$Strand=="."]="*"
  df_New <- df[grepl("^STRG", df$GeneID)&grepl("^CM", df$Chromosome), ]
  gr_New <- GRanges(seqnames = df_New$Chromosome,
                  ranges = IRanges(start = df_New$Start, end = df_New$End),
                  strand = df_New$Strand,
                  TPM = df_New$TPM,
                  GeneID = df_New$GeneID)
                  
  New_o=findOverlaps(gr_New, AnTSS_gr)
  all_indices <- seq_len(length(gr_New))
  # 
  non_overlapping_indices <- setdiff(all_indices, queryHits(New_o))  ###All non overlap STRG
  gr_new_intergenic <- gr_New[non_overlapping_indices]
  TPM10=sort(gr_new_intergenic[gr_new_intergenic$TPM>=10,])
  return(TPM10)
} 
 
TPM_list <- lapply(input_files, function(i){read_bed(i,AnTSS_gr)})

##CAGE
CAGE_files<- list.files(path = "~/FAANG/analysis/02CAGE/01TSS/", pattern = "*.CAGE.TSS.bed", full.names = TRUE, recursive = F)
# read BED files transfer into GRanges
read_CAGE <- function(file) {
  data <- read.table(file, header = F)
  gr <- GRanges(seqnames = data[, 1],
                ranges = IRanges(start = data[, 2], end = data[, 3]))
  return(gr)
}
#read CAGE file
CAGE_list <- lapply(CAGE_files, read_CAGE)

gr_list=list()
for (i in 1:24)
{
tmpC=findOverlaps(TPM_list[[i]],CAGE_list[[i]])
gr_list[[i]]=TPM_list[[i]][unique(queryHits(tmpC))]
}


###gr data for all 24 tissues
combined_24_TPM <- do.call(c, gr_list)


TPM24=sort(combined_24_TPM)

Re_TPM24=TPM24
strand(Re_TPM24) <- "*"
###TPM max present the multiple ranges
reduced_gr <- reduce(Re_TPM24, with.revmap = TRUE)

library(parallel)
# treat each region
process_region <- function(i,reduced_gr, TPM24) {
  original_indices <- reduced_gr$revmap[[i]]
  tmp_gr <- TPM37[original_indices]
  widths <- width(tmp_gr)
  range_diff <- max(widths) - min(widths)

  # if more than 2000，select the min one； else，select the max TPM one
  if (range_diff > 2000) {
  selected_TPM <- tmp_gr[which.min(widths)]
  Out=paste0("Wid:",min(widths))
  } else {
  selected_TPM <- tmp_gr[which.max(mcols(tmp_gr)$TPM)]
  Out=paste0("TPM:",max(mcols(tmp_gr)$TPM))
  }
  selected_TPM$GeneID <- paste0("STRG-", i, "-T", length(widths), "-",Out)
  return(selected_TPM)
}

df_final_list <- mclapply(seq_along(reduced_gr), function(i) { process_region(i, reduced_gr, TPM37)}, mc.cores = detectCores())
df_final <- do.call(c, df_final_list)



peak_counts=df_final
presence_matrix <- matrix(0, nrow = length(peak_counts), ncol = length(gr_list))
rownames(presence_matrix) <- paste(seqnames(peak_counts), start(peak_counts), end(peak_counts), sep = "_")
colnames(presence_matrix) <- tissue_names
# fill matrix
for (i in seq_along(gr_list)) {
    overlaps <- findOverlaps(peak_counts, gr_list[[i]])
    presence_matrix[queryHits(overlaps), i] <- 1
}
sum_presence <- rowSums(presence_matrix)
Total_presence=cbind(presence_matrix,sum_presence)
Total_presence <- as.data.frame(Total_presence)

###descript the overlap tissues peaks
common_peaks <- list()
for(i in seq_along(gr_list))
{
common_peaks[[i]]=peak_counts[Total_presence$sum_presence>=i,]
new_df=Total_presence[Total_presence$sum_presence>=i,]

mcols(common_peaks[[i]]) <- cbind(mcols(common_peaks[[i]]), new_df)
}



# count the length
common_peaks_lengths <- lengths(common_peaks)
#  log10 
log_common_peaks_lengths <- log10(common_peaks_lengths)

names_arg <- 1:length(common_peaks_lengths)
# labels
log10_labels <- function(x) {
  parse(text = paste0("10^", x))
}
# save PDF file
pdf("02TSS_TPM10_24tissues_log.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common new TSS", 
                         xlab = "Consistent tissue # >=", 
                         ylab = "Number of new TSS", 
                         col = "blue", 
                         names.arg = names_arg, 
                         ylim = c(0, log10(10^5)), 
                         yaxt = "n")
# add Y labels
axis(2, at = 0:5, labels = log10_labels(0:5))
#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, log_common_peaks_lengths + 0.1, labels = common_peaks_lengths, cex = 0.7, pos = 3)
dev.off()

data <- common_peaks_lengths
pdf("03TSS_TPM10_24tissues_cutoff.pdf", width = 10, height = 6)
plot(data, type = "o", col = "blue", xlab = "Tissue #", ylab = "Number of peak", main = "Original Data Points with Cutoff: 2tissues")
text(x =1:24, y = data, labels = data, pos = 3, cex = 1.8, col = "black")
abline(v = 2, col = "red", lty = 2)
#legend("topright", legend = paste("Cutoff:", cutoff), col = "red", lty = 2)
dev.off()



High1=common_peaks[[2]]
filtered_list <- lapply(gr_list, function(gr) {gr[gr$TPM >= 50]})
High2 <- do.call(c, filtered_list)

overlaps <- findOverlaps(High2, High1)
non_overlapping_indices <- setdiff(seq_along(High2), queryHits(overlaps))
uniq_High2 <- High2[non_overlapping_indices]
reduced_High2 <- reduce(uniq_High2, with.revmap=T)

High2_list <- mclapply(seq_along(reduced_High2), function(i) { process_region(i, reduced_High2, uniq_High2)}, mc.cores = detectCores())
High2_final <- do.call(c, High2_list)


presence_matrix <- matrix(0, nrow = length(High2_final), ncol = length(gr_list))
colnames(presence_matrix) <- tissue_names

for (i in seq_along(gr_list)) {
    overlaps <- findOverlaps(High2_final, gr_list[[i]])
    presence_matrix[queryHits(overlaps), i] <- 1
}
sum_presence <- rowSums(presence_matrix)
Total_presence=cbind(presence_matrix,sum_presence)
Total_presence <- as.data.frame(Total_presence)
mcols(High2_final) <- cbind(mcols(High2_final), Total_presence)

Final_TPM=sort(c(High1,High2_final))

write.table(as.data.frame(Final_TPM), file = paste0("04TPM_New331.Distribution"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



# PDF
common_TPM=numeric()
for(i in 1:24)
{ common_TPM[i]=length(which(Final_TPM$sum_presence>=i)) }
data <- common_TPM
pdf("05TSS_TPM10_24tissues_Final.pdf", width = 10, height = 6)
plot(data, type = "o", col = "blue", xlab = "Tissue #", ylab = "Number of peak", main = "Original Data Points with Cutoff: 2tissues")
text(x = 2, y = data[2], labels = data[2], pos = 3, cex = 1.8, col = "black")
abline(v = 2, col = "red", lty = 2)
#legend("topright", legend = paste("Cutoff:", cutoff), col = "red", lty = 2)
dev.off()

# log10 value
log_common_peaks_lengths <- log10(common_TPM)
names_arg <- 1:length(common_TPM)
# set labels
log10_labels <- function(x) {
  parse(text = paste0("10^", x))
}
# PDF
pdf("06TSS_TPM10_24tissues_log_final.pdf", width = 10, height = 6)
bar_midpoints <- barplot(log_common_peaks_lengths, 
                         main = "Common new TSS", 
                         xlab = "Consistent tissue # >=", 
                         ylab = "Number of new TSS", 
                         col = "blue", 
                         names.arg = names_arg, 
                         ylim = c(0, log10(10^4)), 
                         yaxt = "n")
# add Y labels
axis(2, at = 0:4, labels = log10_labels(0:4))
#abline(h = log10(10), col = "red", lty = 2)
text(bar_midpoints, log_common_peaks_lengths + 0.1, labels = common_TPM, cex = 0.7, pos = 3)
dev.off()




#########TSS
TSS10=Final_TPM[,2]

# Filter the data into two parts
TSS10_strand <- TSS10[strand(TSS10) %in% c("+", "-"), ]
TSS10_dot <- TSS10[strand(TSS10) == "*", ]

# Process the part with explicit "+" and "-" strands
TSS10_strand$TSS <- ifelse(strand(TSS10_strand) == "+", start(TSS10_strand), end(TSS10_strand))  

TSS10_dot1=TSS10_dot2=TSS10_dot
TSS10_dot1$TSS <- start(TSS10_dot1)
strand(TSS10_dot1)="+"
TSS10_dot2$TSS <- end(TSS10_dot2)
strand(TSS10_dot2)="-"

TSS_TPM10=c(TSS10_strand,TSS10_dot1,TSS10_dot2)

Ann_Fin=AnTSS_gr[grepl("^CM",seqnames(AnTSS_gr)),]
New_Fin=TSS_TPM10[grepl("^CM",seqnames(TSS_TPM10)),]

data_Fin=as.data.frame(New_Fin)
colnames(data_Fin)[6]="geneID"

Total_TSS=rbind(as.data.frame(Ann_Fin),data_Fin)

write.table(Total_TSS, file = paste0("07TSS_annotate_TPM_final.bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


 