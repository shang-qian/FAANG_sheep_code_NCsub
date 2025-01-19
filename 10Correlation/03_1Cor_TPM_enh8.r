
args=c("~/FAANG/analysis/09Promoter_TPM/03pro_TPM")
input_files <- list.files(path =args[1], pattern = "*_01final_promoter_TPM.txt", full.names = TRUE, recursive = TRUE)
print(input_files)

tissue_names <- gsub(".*/03pro_TPM/([^/]+)/.*", "\\1", input_files)

library(data.table)
read_bed <- function(file) {
gr=fread(file, header = TRUE)
return(gr)
}

# gene_expression
GE_list <- lapply(input_files, read_bed)



# geneID
gene_names <- lapply(GE_list, function(df) paste(df$geneID,df$Gseqnames,df$start,df$end,sep="_"))
gene_names_vector <- unlist(gene_names)

###Total gene
TG=unique(gene_names_vector)
TG_split <- strsplit(TG, "_")
TG_df <- do.call(rbind, TG_split)
colnames(TG_df) <- c("geneID", "Pseq", "Pstart", "Pend")
TG_df <- as.data.frame(TG_df, stringsAsFactors = FALSE)

###
TG_df_with_TPM <- TG_df
for (i in seq_along(GE_list)) {
  print(i)
  current_GE <- GE_list[[i]][,c("geneID","TPM")]
  result <- do.call(rbind, lapply(split(current_GE, current_GE$geneID), function(x) {
  if (nrow(x) > 1) {
    x[which.max(x$TPM), ]  # retain TPM max row
  } else {
    x  # retain original
  }
  }))

  #  merge  TG_df and current_GE
  merged_data <- merge(TG_df_with_TPM, result, by = "geneID", all.x = TRUE)
  merged_data$TPM[is.na(merged_data$TPM)] <- 0 
  print(nrow(merged_data))
  colnames(merged_data)[i+4]=paste0("Tis",i)
  TG_df_with_TPM <- merged_data
}

Fs110=TG_df_with_TPM
colnames(Fs110)[5:ncol(Fs110)]=tissue_names

####expression
allexp=Fs110[,5:ncol(Fs110)]

custom_function <- function(x) {
  count_value <- length(x[x!=0])
  return(count_value)
}
allcount=apply(allexp,1,custom_function)
all_gene_exp=cbind(Fs110,allcount)

Unexp_gene=all_gene_exp[all_gene_exp[,ncol(all_gene_exp)]==0,]
write.table(Unexp_gene,"01Unexp_gene.txt",quote=F,row.names=F,sep = "\t")

Gene_sample=all_gene_exp[all_gene_exp[,ncol(all_gene_exp)]>=1,]
write.table(Gene_sample,"02Exp_gene.txt",quote=F,row.names=F,sep = "\t")

TPM=Gene_sample[,-ncol(Gene_sample)]

####Enh_location
fileenh="~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt"
Enh_all=fread(fileenh,header=T)

##CTCF_loop
Loop="~/FAANG/analysis/10interaction/01CTCF/02Loop/01Loop_Prom_Enh.txt"
Loop_all=fread(Loop,header=T)

##CAGE
CAGET="~/FAANG/analysis/10interaction/01CAGE/02Prom_enh/01CAGE_Prom_Enh_uniq.txt"
CAGE_all=fread(CAGET,header=T)
colnames(CAGE_all)=c("Pname","Pst","Pend","Pwid","Ps","Ename","Est","Eend","Ewid","Estrand")



#######
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)

Gene_interaction_fun = function(Gene_s,Enh_all,CAGE_all, Loop_all, tissue_names)
{
set.seed(2024)
Y=Gene_s
genename=paste(Y[1],Y[2],Y[3],Y[4],sep="_")
y=as.matrix(as.numeric(Y[-(1:4)]))

X_2m=Enh_all[Enh_all$seqnames==Y$Pseq & Enh_all$start>as.numeric(Y[3])-1000000 & Enh_all$end<as.numeric(Y[4])+1000000,]
tmp_enh=paste(X_2m$seqnames, X_2m$start, X_2m$end,sep="_")

loop_index=cage_index=numeric()
Loop_pro=Loop_all[Y$Pseq==Loop_all$Promoter.seqnames & as.numeric(Y[3])==Loop_all$Promoter.start-1 & as.numeric(Y[4])==Loop_all$Promoter.end,]
if(nrow(Loop_pro)>0)
 {loop_enh=paste(Loop_pro$Enhancer.seqnames, Loop_pro$Enhancer.start-1, Loop_pro$Enhancer.end, sep="_") 
  loop_index=which(tmp_enh %in% loop_enh)
 }

CAGE_pro=CAGE_all[CAGE_all$Pname==Y$Pseq & as.numeric(Y[3])==CAGE_all$Pst-1 & as.numeric(Y[4])==CAGE_all$Pend,] 
if(nrow(CAGE_pro)>0)
 {
  cage_enh=paste(CAGE_pro$Ename, CAGE_pro$Est-1, CAGE_pro$Eend, sep="_") 
  cage_index=which(tmp_enh %in% cage_enh)
 }


Tindex=unique(c(loop_index,cage_index)) 
if(length(Tindex)>0)
{
NX_2m=X_2m[-Tindex,] 

both_2m=loop_2m=cage_2m=NULL
both=intersect(loop_index,cage_index)
if(length(both)>0) 
  { both_2m=cbind(X_2m[both,],"CTCF:CAGE") 
    colnames(both_2m)[31]="Type"}
loopi=setdiff(loop_index,both)
if(length(loopi)>0) 
  { loop_2m=cbind(X_2m[loopi,],"CTCF")
    colnames(loop_2m)[31]="Type"}
cagei=setdiff(cage_index,both)
if(length(cagei)>0) 
  { cage_2m=cbind(X_2m[cagei,],"CAGE")
    colnames(cage_2m)[31]="Type"}

YX_2m=rbind(both_2m,loop_2m,cage_2m)
} else {
NX_2m=X_2m
YX_2m=NULL
}

NX_2m <- if (is.null(NX_2m)) NULL else as.matrix(NX_2m)
YX_2m <- if (is.null(YX_2m)) NULL else as.matrix(YX_2m)

#rawP=Y_X_cor(YX_2m,NX_2m,y,genename, tissue_names)
######################correlation

tmpy2=y
if(length(which(tmpy2==0))>0)
{tmpy2[tmpy2!=0]=1}

###NX_2m
if(length(NX_2m)>0)
{
NX=as.matrix(t(NX_2m)[6:(ncol(NX_2m)-1),])
vec=matrix(,ncol(NX),3)
for (i in 1:ncol(NX))
{
  if(sum(as.numeric(NX[,i]))!=24)
  {
    vectmp=cor.test(log10(y+0.0001),as.numeric(NX[,i]))
    vectmp2=cor.test(tmpy2,as.numeric(NX[,i]))
    vec[i,1]=i
    vec[i,2]=max(vectmp$estimate,vectmp2$estimate)
    vec[i,3]=min(vectmp$p.value,vectmp2$p.value)
  } 
  if(sum(as.numeric(NX[,i]))==24) 
  {
    if(length(which(tmpy2!=0))==24) 
    {
    vec[i,1]=i
    vec[i,2]=1
    vec[i,3]=0
    } 
  }
}
vec1=matrix(vec[!is.na(vec[,2]),],,3)
vecT1=matrix(vec1[vec1[,2]>=0.8&vec1[,3]<=0.01,],,3)
nx_matrix=NULL
if(nrow(vecT1)>0)
{
nx1_matrix=apply(as.matrix(NX[,vecT1[,1]]),2,as.numeric)
colnames(nx1_matrix)=paste(NX_2m[,1],NX_2m[,2],NX_2m[,3],sep="_")[vecT1[,1]]
nx_matrix=rbind(nx1_matrix,"None")
}
} else {nx_matrix=NULL}


if(length(YX_2m)>0)
{
YX=as.matrix(t(YX_2m)[6:(ncol(YX_2m)-2),])
vec=matrix(,ncol(YX),3)
for (i in 1:ncol(YX))
{
  if(sum(as.numeric(YX[,i]))!=24)
  {
    vectmp=cor.test(log10(y+0.0001),as.numeric(YX[,i]))
    vec[i,1]=i
    vec[i,2]=vectmp$estimate
    vec[i,3]=vectmp$p.value
  }
  if(sum(as.numeric(YX[,i]))==24)
  {
    if(length(which(tmpy2!=0))==24) 
    {
    vec[i,1]=i
    vec[i,2]=1
    vec[i,3]=0
    } 
  }
}
vec1=matrix(vec[!is.na(vec[,2]),],,3)
vecT2=matrix(vec1[vec1[,2]>0,],,3)
yx_matrix=NULL
if(nrow(vecT2)>0)
{
yx1_matrix=apply(as.matrix(YX[,vecT2[,1]]),2,as.numeric)
colnames(yx1_matrix)=paste(YX_2m[,1],YX_2m[,2],YX_2m[,3],sep="_")[vecT2[,1]]
source=matrix(YX_2m[,ncol(YX_2m)][vecT2[,1]],1,)
yx_matrix=rbind(yx1_matrix,source)
} 
}else {yx_matrix=NULL}


x_matrix=cbind(nx_matrix,yx_matrix)


if(length(x_matrix)>0)
{
ytpm=matrix(c(y,"TPM"),,1)
rownames(ytpm)=c(tissue_names,"Exp")
colnames(ytpm)=genename
oneAll=cbind(ytpm,x_matrix)
if (!dir.exists("00oneAll")) {
  dir.create("00oneAll")}
write.table(oneAll,paste0("00oneAll/",genename,"oneAll.txt"),quote=F)


#each interaction gene
pcoef <- matrix(,ncol(x_matrix),10)
colnames(pcoef)=c("promoter","enhancer","lm_coef","lm_pvalue","lm_R2","prom_tis#","enh_tis#","cor_coef_TPM","cor_coef_01","source")
pcoef[,1]=genename
pcoef[,2]=colnames(x_matrix)

for (i in 1:nrow(pcoef)) 
{
#  print(i)
  if(sum(as.numeric(x_matrix[1:24, i]))!=24)
  {
  lm_model <- lm(log10(as.numeric(y)+0.0001) ~ as.numeric(x_matrix[1:24, i]))  
  pcoef[i,3:4] <- summary(lm_model)$coefficients[2, c(1,4)]  # 第二行对应系数，第1/4列对应effect/p-value
  pcoef[i,5]=summary(lm_model)$adj.r.squared 
  
  cor_result <- cor.test(as.numeric(x_matrix[1:24,i]), log10(as.numeric(y)+0.0001))
  pcoef[i,8] <- cor_result$estimate

  cor_result2 <- cor.test(as.numeric(x_matrix[1:24,i]), tmpy2)
  pcoef[i,9] <- cor_result2$estimate
  } 
  
  if(sum(as.numeric(x_matrix[1:24, i]))==24)
  {
    if(length(tmpy2[tmpy2!=0])==24)
    {
     pcoef[i,3:4] <- c(1,0)
     pcoef[i,5]=1
     pcoef[i,8] <- 1
     pcoef[i,9] <- 1
    }
    if(length(tmpy2[tmpy2!=0])!=24)
    {
     pcoef[i,3:4] <- c(0,1)
     pcoef[i,5]=0.1
     pcoef[i,8] <- 0.1
     pcoef[i,9] <- 0.1
    }
  }
  
  if(as.numeric(pcoef[i,5])>=0.9)
  {
  if (!dir.exists("02Plot_R9")) {
  dir.create("02Plot_R9")} 
  pdf(paste0("02Plot_R9/",pcoef[i,1],"_",pcoef[i,2],"-TPM_scatter.pdf"))
  plot(as.numeric(x_matrix[1:24,i]), as.numeric(y), pch=16, xlab=pcoef[i,2], ylab="TPM", main = paste0(pcoef[i,1]," Scatter Plot"), xlim=c(-0.5,1.5))      # 禁止自动绘制x轴的刻度)
#  predicted_values <- predict(lm_model)
#  lines(as.numeric(x_matrix[1:24,i]), predicted_values, col = "red")
  dev.off()
  }  
  
  TPMn=tissue_names[which(y!=0)]
  pcoef[i,6]=paste0(length(TPMn), ":", paste(TPMn, collapse=":"))

  nonzerox=tissue_names[which(as.numeric(x_matrix[1:24, i])!=0)]
  pcoef[i,7] <-paste0(length(nonzerox),":",paste(nonzerox,collapse=":")) 

  pcoef[i,10] <- as.character(x_matrix[25,i])
}
if (!dir.exists("01Cor_allP")) {
  dir.create("01Cor_allP")}
write.table(pcoef,paste0("01Cor_allP/",genename,"Cor_allP.txt"),quote=F)

return(pcoef)
}
}


core.number <- detectCores()
cl <- makeCluster(getOption("cl.cores", core.number))
clusterEvalQ(cl, {library(orthopolynom)})
clusterEvalQ(cl, {library(glmnet)})

clusterExport(cl,c("Enh_all","Gene_interaction_fun","TPM","CAGE_all","Loop_all","tissue_names"),envir=environment())

####calculate each gene
#Gene_list=pbsapply(1:nrow(TPM),function(c) Gene_interaction_fun(Gene_s=TPM[c,],X=X),cl=cl)
#Gene_list=pbsapply(1:nrow(TPM),function(c) Gene_interaction_fun(Gene_s=TPM[c,],Enh_all=Enh_all,CAGE_all=CAGE_all, Loop_all=Loop_all, tissue_names),cl=cl)
Gene_list=pbsapply(1:nrow(TPM),function(c) Gene_interaction_fun(Gene_s=TPM[c,],Enh_all=Enh_all,CAGE_all=CAGE_all, Loop_all=Loop_all, tissue_names),cl=cl)

stopCluster(cl)


final_list=do.call(rbind, Gene_list)

write.table(final_list,"04TPM_interaction_raw.txt",quote=F,row.names=F)
#write.csv(final_list,"05TPM_interaction.csv",quote=F,row.names=F)

# obatin p_values
p_values <- as.numeric(final_list[, "lm_pvalue"])
# Benjamini-Hochberg adjust p values
adjusted_p_values <- p.adjust(p_values, method = "BH")

# adjust p value to final_list 
final_list <- cbind(final_list, adjusted_p_values)
# filter p < 0.05 
significant_results <- final_list[adjusted_p_values < 0.05, ]

dim(significant_results)


write.table(significant_results,"05TPM_interaction_final.txt",quote=F,row.names=F)



