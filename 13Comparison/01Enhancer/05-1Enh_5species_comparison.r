args <- commandArgs(trailingOnly = TRUE)

library(rtracklayer)
print(args)

# 函数：读取peak文件并转换为GRanges对象
read_peaks <- function(file) {
  data <- read.table(file, header = F) 
  return(data)
}
#读取和转换peak文件
tissue_list <- lapply(args[1:4], read_peaks)

# 函数：读取peak文件并转换为GRanges对象
read_pair <- function(file) {
  data <- read.table(file, header = T) 
  return(data)
}
#读取和转换peak文件
pair_list <- lapply(args[6:9], read_pair)

gene_list=list()
for (i in 1:4)
{
tmp=pair_list[[i]]
geneID=unique(tmp[tmp$enh_tmp %in% tissue_list[[i]]$V1,]$pair)

gene_list[[i]]=unique(sub("_CM.*", "", geneID))
}


####Venn
library(VennDiagram)
  Cerebellum = gene_list[[1]]
  Cortex = gene_list[[2]]
  Lung = gene_list[[3]]
  Muscle = gene_list[[4]]

input <- list( Cerebellum ,  Cortex ,  Lung ,  Muscle )

names(input)<-c(paste("Cerebellum",length(Cerebellum)),paste("Cortex",length(Cortex)),paste("Lung",length(Lung)),paste("Muscle",length(Muscle)))
temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F","#4CC1EF","#FFCC4C"))
#temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F"))

pdf(paste0(args[5],"/02",args[5],"_venn_5species.pdf"))
grid.draw(temp)
dev.off()


# 1. 找到Cerebellum的独特基因
Cerebellum_unique <- setdiff(Cerebellum, union(union(Cortex, Lung), Muscle))
Cerebellum_unique1=sub("gene-","",Cerebellum_unique)
write.table(as.data.frame(Cerebellum_unique1), file=paste0(args[5],"/03",args[5],"_Cerebellum_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)
# 2. 找到Cortex的独特基因
Cortex_unique <- setdiff(Cortex, union(union(Cerebellum, Lung), Muscle))
Cortex_unique1=sub("gene-","",Cortex_unique)
write.table(as.data.frame(Cortex_unique1), file=paste0(args[5],"/04",args[5],"_Cortex_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)
# 3. 找到Lung的独特基因
Lung_unique <- setdiff(Lung, union(union(Cerebellum, Cortex), Muscle))
Lung_unique1=sub("gene-","",Lung_unique)
write.table(as.data.frame(Lung_unique1), file=paste0(args[5],"/05",args[5],"_Lung_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)
# 4. 找到Muscle的独特基因
Muscle_unique <- setdiff(Muscle, union(union(Cerebellum, Cortex), Lung))
Muscle_unique1=sub("gene-","",Muscle_unique)
write.table(as.data.frame(Muscle_unique1), file=paste0(args[5],"/06",args[5],"_Muscle_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)

# 找到Cerebellum和Cortex的交集
cerebellum_cortex_common <- intersect(Cerebellum, Cortex)
# 从交集中排除在Lung和Muscle中出现的基因
unique_brain <- setdiff(cerebellum_cortex_common, union(Lung, Muscle))
brain_unique1=sub("gene-","",unique_brain)
write.table(as.data.frame(brain_unique1), file=paste0(args[5],"/07",args[5],"_Brain_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)


common_genes_all <- Reduce(intersect, list(Cerebellum, Cortex, Lung, Muscle))
common_genes_all1=sub("gene-","",common_genes_all)
write.table(as.data.frame(common_genes_all1), file=paste0(args[5],"/08",args[5],"_4tissues_common_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)


All_genes <- union(Cerebellum, union(union(Cortex, Lung), Muscle))
All_genes1=sub("gene-","",All_genes)
write.table(as.data.frame(All_genes1), file=paste0(args[5],"/09",args[5],"_4tissues_any1_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=F)




