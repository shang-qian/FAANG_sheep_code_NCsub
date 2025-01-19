
library(data.table)

args=c("~/20240401FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt",
"~/20240401FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt")

P1=fread(args[1], header = TRUE)
P2=fread(args[2], header=TRUE)

tisname=colnames(P2)[3:26]

P1$geneID=do.call(rbind,strsplit(P1$promoter,"_"))[,1]

Pall=P1
uegene=NULL
#TableS9=matrix(,24,7)
for( i in 1:24)
{
TisN=tisname[i]
#print(TisN)

Tis_e=Pall[Pall$"enh_tis#" %in% paste("1",TisN,sep=":"),]
Tis_ep=Tis_e[grepl(TisN,Tis_e$"prom_tis#"),]
#Tis_ep=Tis_e[Tis_e$"prom_tis#" %in% paste("1",TisN,sep=":"),]
print(nrow(Tis_ep))

enh_pos=do.call(rbind,strsplit(Tis_ep$enhancer,"_"))
#prom_pos=unique(enh_pos[,2:4])
enh_pos_gene=cbind(enh_pos,Tis_ep$geneID)

dir.create(TisN)
write.table(as.data.frame(enh_pos_gene), file = paste0(TisN,"/01",TisN,"_enhancer_position_gene_2.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}







