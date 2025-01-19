
library(data.table)

args=c("~/20240401FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt",
"~/20240401FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt")

P1=fread(args[1], header = TRUE)
P2=fread(args[2], header=TRUE)

tisname=colnames(P2)[3:26]

P1$geneID=do.call(rbind,strsplit(P1$promoter,"_"))[,1]

Pall=P1
uegene=NULL
TableS9=matrix(,24,7)
for (i in 1:24)
{
TisN=tisname[i]
print(TisN)
Tis_p=Pall[Pall$"prom_tis#" %in% paste("1",TisN,sep=":"),]
Tis_pp=Tis_p
#print(nrow(Tis_pp))
p_gene=unique(Tis_pp$geneID)

TisN=tisname[i]
print(TisN)

Tis_e=Pall[Pall$"enh_tis#" %in% paste("1",TisN,sep=":"),]
Tis_ep=Tis_e[grepl(TisN,Tis_e$"prom_tis#"),]

print(nrow(Tis_ep))
e_gene=unique(Tis_ep$geneID)

p_gene_names <- sub("gene-", "", p_gene)
e_gene_names <- sub("gene-", "", e_gene)

dir.create(TisN)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/01",TisN,"_promoter_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(e_gene_names)), file = paste0(TisN,"/02",TisN,"_enhancer_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

TableS9[i,1]=TisN
TableS9[i,2]=length(p_gene_names)
TableS9[i,3]=paste(p_gene_names,collapse="|")
TableS9[i,4]=length(e_gene_names)
TableS9[i,5]=paste(e_gene_names,collapse="|")

###Venn
library(VennDiagram)
Prom_gene=p_gene_names
Pair_gene=e_gene_names
input=list(Prom_gene,Pair_gene)
names(input)<-c(paste0(TisN,"_Prom_gene",length(Prom_gene)),paste0(TisN,"_Enh_gene",length(Pair_gene)))
#temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F","#4CC1EF","#FFCC4C"))
temp=venn.diagram(input,filename = NULL,height = 3000, width = 3000,col = "transparent",fill = c("#3A5C84","#F7931F"))

pdf(paste0(TisN,"/03",TisN,"_Venn2.pdf")) 
grid.draw(temp)
dev.off()

overlapg=intersect(p_gene,e_gene)
TableS9[i,6]=length(overlapg)
TableS9[i,7]=paste(sub("gene-", "", overlapg),collapse="|")

uegene[[i]]=setdiff(e_gene,overlapg)
uegene_all=Pall[Pall$geneID %in% uegene[[i]],]

write.table(as.data.frame(uegene_all), file = paste0(TisN,"/04",TisN,"_unique_enhancer_information.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
write.table(as.data.frame(TableS9), file = paste0("TableS9_prom_enh_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



