
args <- commandArgs(trailingOnly = TRUE)

# load R packages
library(data.table)
library(GenomicRanges)


Pro_df <- fread(args[1], header = T)
TPM_df <- fread(args[2], header = T)

# 使用 merge 函数进行左连接，基于 GeneID 和 geneID 进行匹配
Anno_df <- Pro_df[!grepl("^STRG", Pro_df$geneID), ]
New_df <- Pro_df[grepl("^STRG", Pro_df$geneID), ]

Anno_TPM <- merge(Anno_df, TPM_df[, c("GeneID", "TPM")], by.x = "geneID", by.y = "GeneID", all.x = TRUE)

TPM_df$Strand[TPM_df$Strand=="."]="*"
TPM_gr <- GRanges(
  seqnames = TPM_df$Chromosome,
  ranges = IRanges(start = TPM_df$Start, end = TPM_df$End),
  strand = TPM_df$Strand,
  TPM = TPM_df$TPM,
  geneID = TPM_df$GeneID
)

New_gr=GRanges(seqnames=New_df$Gseqnames,
   ranges= IRanges(start=New_df$Gstart, end=New_df$Gend),
   strand=New_df$Gstrand)
   
New_TPM=New_df
New_TPM$TPM=NA
#  findOverlaps find overlaps
overlaps = findOverlaps(New_gr, TPM_gr)
# obtain New_gr and TPM_gr
query_hits = queryHits(overlaps)
subject_hits = subjectHits(overlaps)
# 
best_TPM <- tapply(subject_hits, query_hits, function(idx) {
  min(TPM_gr[idx]$TPM)
})
#
New_TPM$TPM[as.numeric(names(best_TPM))] = best_TPM

Pro_TPM=rbind(Anno_TPM, New_TPM)

write.table(as.data.frame(Pro_TPM), file = paste0(args[3],"/",args[3],"_01final_promoter_TPM.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
