
args=c("~/FAANG/analysis/07Enhancer_peak/9Enhancer_Tissue_Final/03All_enhancer_across_tissues.txt",
"~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final/06All_promoter_across_tissues.txt",
"~/FAANG/analysis/10interaction/01CAGE/01TSS_enh_CAGE")

library(data.table)
library(parallel)
library(GenomicRanges)
library(ChIPseeker)

read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

enh_file=args[1]
enhancer_data=read_peaks(enh_file)

prom_file=args[2]
promoter_data=read_peaks(prom_file)


####CAGE
fimo_ls=list()
cagefile=list.files(path = args[3], pattern = "*CAGE_enh_pro_input.bed", full.names = TRUE)
for ( i in 1:24)
{
fimo_ls[[i]] <- fread(cagefile[i], header = F, sep = "\t")
}

fimo_data=do.call(rbind, fimo_ls)
fimo_data=unique(fimo_data)

colnames(fimo_data)=c("sequence_name","pro_st","pro_end","enh_seqname","enh_st","enh_end")
sequence_names <- unique(fimo_data$sequence_name)

gr_prom=function(new_loops){
gr_loops <- GRanges(
  seqnames = new_loops$sequence_name,
  ranges = IRanges(start = new_loops$pro_st-2000, end = new_loops$pro_end+2000),
  strand = "*" 
)
return(gr_loops)
}

gr_enhancer=function(new_loops){
gr_loops <- GRanges(
  seqnames = new_loops$enh_seqname,
  ranges = IRanges(start = new_loops$enh_st-2000, end = new_loops$enh_end+2000),
  strand = "*" 
)
return(gr_loops)
}

  gr_pro <- gr_prom(fimo_data)
  gr_enh <- gr_enhancer(fimo_data)
  pro_op <- findOverlaps(gr_pro, promoter_data)
  enh_op <- findOverlaps(gr_enh, enhancer_data)
  
  ov_p=intersect(queryHits(pro_op),queryHits(enh_op))

  result_p_e=NULL
  for (i in 1:length(ov_p))
  {
  print(paste(i,"total:",length(ov_p)))
  gr_pro[ov_p[i]]
  gr_enh[ov_p[i]]
  proID=which(queryHits(pro_op)==ov_p[i])
  promoter1=promoter_data[subjectHits(pro_op)[proID]]

  enhID=which(queryHits(enh_op)==ov_p[i])  
  enhancer1=enhancer_data[subjectHits(enh_op)[enhID]]
  for (j in 1:length(proID))
  {
  tmp_p_e=cbind(as.data.frame(promoter1[j]),as.data.frame(enhancer1))
  result_p_e=rbind(result_p_e,tmp_p_e)
  }
  }


write.table(as.data.frame(unique(result_p_e)), file=paste0("01CAGE_Prom_Enh_uniq.txt"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

