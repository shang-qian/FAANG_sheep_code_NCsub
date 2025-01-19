# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

library(data.table)

args=c("~/20240401FAANG/analysis/10interaction/04Pairs/01Exp_gene_pair_final.txt",
"~/20240401FAANG/analysis/10interaction/04Pairs/02Unexp_gene_pair_final.txt")

P1=fread(args[1], header = TRUE)
P2=fread(args[2], header=TRUE)

tisname=colnames(P2)[3:26]

P1$geneID=do.call(rbind,strsplit(P1$promoter,"_"))[,1]

Pall=P1




####Brain
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:Cerebellum",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:CerebralCortex",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:Cerebellum:CerebralCortex",]
pool=c("1:Cerebellum","1:CerebralCortex","2:Cerebellum:CerebralCortex")

p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)

p_gene=unique(c(p_gene1,p_gene2,p_gene3))

geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")

index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}

p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Brain"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


###Heart
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:HeartRightAtrium",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:HeartRightVentricle",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:HeartRightAtrium:HeartRightVentricle",]
pool=c("1:HeartRightAtrium","1:HeartRightVentricle","2:HeartRightAtrium:HeartRightVentricle")

p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)

p_gene=unique(c(p_gene1,p_gene2,p_gene3))

geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")

index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}

p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Heart"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



#Tis_p=Pall[Pall$"enh_tis#" %in% "2:AdrenalCortex:AdrenalMedulla",]
###Adrenal
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:AdrenalCortex",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:AdrenalMedulla",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:AdrenalCortex:AdrenalMedulla",]
pool=c("1:AdrenalCortex","1:AdrenalMedulla","2:AdrenalCortex:AdrenalMedulla")
p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)
p_gene=unique(c(p_gene1,p_gene2,p_gene3))
geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")
index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}
p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Adrenal"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#Pall[Pall$"enh_tis#" %in% "2:Reticulum:RumenAtrium",]
###Rumen
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:Reticulum",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:RumenAtrium",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:Reticulum:RumenAtrium",]
pool=c("1:Reticulum","1:RumenAtrium","2:Reticulum:RumenAtrium")
p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)
p_gene=unique(c(p_gene1,p_gene2,p_gene3))
geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")
index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}
p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Digestive"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




#Tis_p=Pall[Pall$"enh_tis#" %in% "2:Ileum:IleumPeyersPatch",]
#TisN="Small_intestine"
##Ileum
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:Ileum",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:IleumPeyersPatch",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:Ileum:IleumPeyersPatch",]
pool=c("1:Ileum","1:IleumPeyersPatch","2:Ileum:IleumPeyersPatch")
p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)
p_gene=unique(c(p_gene1,p_gene2,p_gene3))
geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")
index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}
p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Small_intestine"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



#Tis_p=Pall[Pall$"enh_tis#" %in% "2:DescendingColon:SpiralColon",]
#TisN="Intestinal"
###COlon
Tis_p1=Pall[Pall$"enh_tis#" %in% "1:DescendingColon",]
Tis_p2=Pall[Pall$"enh_tis#" %in% "1:SpiralColon",]
Tis_p3=Pall[Pall$"enh_tis#" %in% "2:DescendingColon:SpiralColon",]
pool=c("1:DescendingColon","1:SpiralColon","2:DescendingColon:SpiralColon")
p_gene1=unique(Tis_p1$geneID)
p_gene2=unique(Tis_p2$geneID)
p_gene3=unique(Tis_p3$geneID)
p_gene=unique(c(p_gene1,p_gene2,p_gene3))
geneT=g_list=NULL;k=0
for (i in 1:length(p_gene))
{
tmp=Pall[Pall$geneID==p_gene[i],]
enh_uniq=unique(tmp$"enh_tis#")
index=match(enh_uniq,pool)
if(length(index[!is.na(index)])>1)
{print(i)
print(p_gene[i])
k=k+1
geneT[k]=p_gene[i]
g_list[[k]]=tmp
}
}
p_gene_names <- sub("gene-", "", geneT)
Tis_ep=do.call(rbind,g_list)
TisN="Intestinal"
dir.create(TisN)
write.table(as.data.frame(Tis_ep), file = paste0(TisN,"/01",TisN,"_pairs.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(as.data.frame(sort(p_gene_names)), file = paste0(TisN,"/02",TisN,"_pairs_gene_list.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

