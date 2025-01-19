args <- commandArgs(trailingOnly = TRUE)

args=c("~/FAANG/analysis/01ChIP-seq/04Filter/N_H3K27ac",
"~/FAANG/analysis/06PeakAnno",
"/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff",
"H3K27ac")

###module load R
library(ChIPseeker)
library(GenomicFeatures)
# workdirection
path=paste0(args[2],"/",args[4])
setwd(path)

####sheep
ath <- makeTxDbFromGFF(args[3],format="gff3")

base_path <- args[1]
file_suffix <- args[4]

get_file_paths <- function(base_path, pattern) {
  file_paths <- list.files(base_path, pattern = pattern, recursive = F)
  return(file_paths)
}

fileName <- get_file_paths(base_path, file_suffix)

Annofiles=vector(mode = "list", length = length(fileName))
for (i in 1:length(fileName))
{
print(i)
sample=strsplit(fileName[i], split = paste0("_",args[4],"_"))
outname=sample[[1]][1]
print(outname)
peaksfile <- fileName[i]

peak <- readPeakFile(paste0(base_path,"/",peaksfile),header=F)
#peakAnno <- annotatePeak(peak,TxDb = ath, assignGenomicAnnotation =TRUE, tssRegion=c(-2000, 2000),genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Intergenic"))
peakAnno <- annotatePeak(peak,TxDb = ath, assignGenomicAnnotation =TRUE, tssRegion=c(-2000, 2000),genomicAnnotationPriority = c("Promoter", "Exon", "Intron","Intergenic"))
       
#slotNames(peakAnno)
tmpanno=peakAnno@annoStat
newanno=tmpanno[1:4,]

newanno[,1]=c("TSS_2kb","Exon","Intron","Intergenic")
newanno[1,2]=sum(tmpanno$Frequency[grep("Promoter",tmpanno$Feature)])
newanno[2,2]=sum(tmpanno$Frequency[grep("Exon",tmpanno$Feature)])
newanno[3,2]=sum(tmpanno$Frequency[grep("Intron",tmpanno$Feature)])
newanno[4,2]=100-sum(as.numeric(newanno[1:3,2]))

rownames(newanno)=1:4

tmpanno=rbind(tmpanno,newanno)

peakAnno@annoStat=tmpanno[(nrow(tmpanno)-3):nrow(tmpanno),]
pdf(file = paste0(outname,"_",args[4],".Pie_anno.pdf"))
plotAnnoPie(peakAnno,main=paste0(outname,"_",args[4],"\nDistribution of Peaks"),line=-8)
dev.off()

Annofiles[[i]]=peakAnno
names(Annofiles)[i]=outname
}


pdf(paste0("24_TissueBar.",args[4],".pdf"))
plotAnnoBar(Annofiles,verbose=T)
dev.off()


#Annofiles[[i]]@annoStat$Frequency
#SD and mean
all_frequencies <- list()

#  Annofiles Frequency 
for (i in 1:length(Annofiles)) {
  all_frequencies[[i]] <- Annofiles[[i]]@annoStat$Frequency
}


frequency_df <- do.call(cbind, all_frequencies)
colnames(frequency_df) <- paste0("Subset_", 1:24)

# mean
mean_frequencies <- rowMeans(frequency_df)

# sd
sd_frequencies <- apply(frequency_df, 1, sd)

library(ggplot2)


summary_data <- data.frame(
  Feature = Annofiles[[1]]@annoStat$Feature,  
  Mean_Frequency = mean_frequencies,
  SD_Frequency = sd_frequencies
)

write.table(summary_data, file = paste0("02",args[4],"Mead_SD_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


feature_colors <- c(
  "TSS_2kb" = "#AFCBE3",  # Light Blue
  "Exon" = "#4682B4",    # Steel Blue
  "Intron" = "#98FB98",    # Pale Green
  "Intergenic" = "#32CD32"     # Lime Green
)


pdf(paste0("3_24TissueBar_mean_sd.",args[4],".pdf"))
ggplot(summary_data, aes(x = Feature, y = Mean_Frequency, fill = Feature)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean_Frequency - SD_Frequency, ymax = Mean_Frequency + SD_Frequency), width = 0.2, color = "black") +
  geom_text(aes(label = paste0(round(Mean_Frequency, 2), "\n±", round(SD_Frequency, 2))), 
            vjust = -1.5, size = 3.5) +
  labs(title = "Mean and Standard Deviation of Genomic Feature Across 24 tissues", 
       x = " ", y = "Frequency (%)") +
  scale_fill_manual(values = feature_colors) +  # colors
  scale_y_continuous(breaks = seq(0, 60, by = 10)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),  # 
        axis.text.y = element_text(size = 14),  # 
        axis.title.x = element_text(size = 16),  # 
        axis.title.y = element_text(size = 16))  # 
dev.off()


