args <- commandArgs(trailingOnly = TRUE)

#library(GenomicRanges)
#library(ChIPseeker)
library(rtracklayer)
print(args)

read_peaks <- function(file) {
  data <- read.table(file, header = T) 
  return(data)
}

tissue_list <- lapply(args[1:2], read_peaks)

Overlap_gene=intersect(tissue_list[[1]]$All_genes,tissue_list[[2]]$All_genes)

write.table(as.data.frame(Overlap_gene), file=paste0("01",args[3],"_All.gene"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


####Venn
library(VennDiagram)

Enhancer = tissue_list[[1]]
Promoter = tissue_list[[2]]


input <- list(  Enhancer ,  Promoter)

names(input)<-c(paste("Enhancer",length(Enhancer)),paste("Promoter",length(Promoter)))
#temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F","#4CC1EF","#FFCC4C"))
temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F"))

pdf(paste0("02",args[3],"_venn_All_genes.pdf"))
grid.draw(temp)
dev.off()


