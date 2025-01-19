
####Enh_location
fileenh="~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt"
Enh_all=fread(fileenh,header=T)

###calculate the effect of pairs
library(glmnet)

sigres=fread("05TPM_interaction_final.txt",header=T)

prom_ls=unique(sigres$promoter)
TPMname=paste(TPM$geneID,TPM$Pseq,TPM$Pstart,TPM$Pend,sep="_")
Enhname=paste(Enh_all$seqnames,Enh_all$start,Enh_all$end,sep="_")

f_eff=NULL
for(i in 1:length(prom_ls))
{
print(i)

TPMindex=which(TPMname==prom_ls[i])
TPMy=as.numeric(TPM[TPMindex,5:28])

enhn=sigres[sigres$promoter==prom_ls[i],]$enhancer

Enhindex=match(enhn, Enhname)
Enhx=as.matrix(Enh_all[Enhindex,6:29])
rownames(Enhx)=enhn

# y is promoter expression，x_matrix the matrix of enhancer or not
x_matrix <- t(Enhx) 
y <- log10(TPMy+0.0001)  
#y[y!=0]=1
if(ncol(x_matrix)>1)
{
result<-tryCatch({
      # 使用岭回归
      cv_ridge <- cv.glmnet(x_matrix, y, alpha = 0)
      best_lambda <- cv_ridge$lambda.min
      
      # 最佳 lambda 下拟合岭回归模型
      ridge_coefficients <- coef(glmnet(x_matrix, y, alpha = 0, lambda = best_lambda))
      
      # 提取非截距项的系数
      coef_df <- data.frame(Effect = as.numeric(ridge_coefficients[-1]), 
                            Xname = rownames(ridge_coefficients)[-1])
      return(coef_df)
    }, warning = function(w) {
    # 捕获警告并继续执行岭回归
    message("Warning during Ridge regression: ", conditionMessage(w))
    cv_ridge <- cv.glmnet(x_matrix, y, alpha = 0)
    best_lambda <- cv_ridge$lambda.min
    ridge_coefficients <- coef(glmnet(x_matrix, y, alpha = 0, lambda = best_lambda))
    coef_df <- data.frame(Effect = as.numeric(ridge_coefficients[-1]), 
                          Xname = rownames(ridge_coefficients)[-1])
    return(coef_df)   
  }, error = function(e) {
      # 如果岭回归失败，回退到逐个 enhancer 的回归分析
      message("Ridge regression failed, performing individual linear regression for each enhancer on promoter: ", prom_ls[i])
      coef_df <- data.frame(Effect = numeric(0), Xname = character(0))  # 初始化空数据框   
      # 针对每个 enhancer 进行单独回归
      for (j in 1:ncol(x_matrix)) {
        lm_model <- lm(y ~ x_matrix[, j])
        coef_df <- rbind(coef_df, data.frame(Effect = coef(lm_model)[2], Xname = colnames(x_matrix)[j]))
      }
      return(coef_df)
    })
coef_df=result    
} 

if(ncol(x_matrix)==1) 
{
lm_model <- lm(y ~ x_matrix)
coefficients <- coef(lm_model)
coef_df=data.frame(Effect=coefficients["x_matrix"],Xname=colnames(x_matrix))
}

tmp=sigres[sigres$promoter==prom_ls[i],]
ftmp=cbind(tmp,coef_df)
f_eff=rbind(f_eff,ftmp)
}


write.table(f_eff,"06TPM_interaction_Effect_final.txt",quote=F,row.names=F)



####UNexpress gene
un_gene=Unexp_gene[,1:4]
##CTCF_loop
Loop="~/FAANG/analysis/10interaction/01CTCF/02Loop/01Loop_Prom_Enh.txt"
Loop_all=fread(Loop,header=T)

##CAGE
CAGET="~/FAANG/analysis/10interaction/01CAGE/02Prom_enh/01CAGE_Prom_Enh_uniq.txt"
CAGE_all=fread(CAGET,header=T)
colnames(CAGE_all)=c("Pname","Pst","Pend","Pwid","Ps","Ename","Est","Eend","Ewid","Estrand")

pcoef=NULL
for (i in 1:nrow(un_gene))
{

Y=un_gene[i,]
genename=paste(Y[1],Y[2],Y[3],Y[4],sep="_")

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

both=intersect(loop_index,cage_index)

if(length(both)>0) 
  { 
    print(i)
    print(length(both))
    both_2m=cbind(X_2m[both,],"CTCF:CAGE") 
    print(both_2m)
    colnames(both_2m)[31]="Type"    
    enh_tmp=tmp_enh[both]
    peak=both_2m[,6:31]    
    tmp_both=cbind(genename,enh_tmp,peak)    
    pcoef=rbind(pcoef,tmp_both)  
    }
}


write.table(pcoef,"07Unexp_interaction_final.txt",quote=F,row.names=F)

enh_tis=NULL
for(i in 1: nrow(pcoef))
{
Tis_tmp=pcoef[i,3:26]

Tis_info=c(length(which(Tis_tmp==1)),colnames(Tis_tmp)[which(Tis_tmp==1)])

enh_tis[i]=paste(Tis_info,collapse=":")

}

pcoef$enh_tis=enh_tis


write.table(pcoef,"07Unexp_interaction_final.txt",quote=F,row.names=F)

