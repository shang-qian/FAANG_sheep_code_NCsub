args <- commandArgs(trailingOnly = TRUE)


library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)
print(args)
Tisname=args[2]


read_peaks <- function(file) {
  peaks <- readPeakFile(file)
  gr <- makeGRangesFromDataFrame(peaks)
  return(gr)
}

# 读取和转换peak文件
s_list <- lapply(args[1], read_peaks)
sheep <- s_list[[1]]

#cattle
c_files=list.files(path =".", pattern = "*01cattle.*liftover8.enhancer", full.names = TRUE, recursive = TRUE)
c_list <- lapply(c_files, read_peaks)
cattle=reduce(do.call(c, c_list))

#pig
p_files=list.files(path =".", pattern = "*01pig.*liftover8.enhancer", full.names = TRUE, recursive = TRUE)
p_list <- lapply(p_files, read_peaks)
pig=reduce(do.call(c, p_list))

#mouse
m_files=list.files(path =".", pattern = "*01mouse.*liftover8.enhancer", full.names = TRUE, recursive = TRUE)
m_list <- lapply(m_files, read_peaks)
mouse=reduce(do.call(c, m_list))

#human
h_files=list.files(path =".", pattern = "*01human.*liftover8.enhancer", full.names = TRUE, recursive = TRUE)
h_list <- lapply(h_files, read_peaks)
human=reduce(do.call(c, h_list))



# 定义函数计算两两并集
compute_union <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
#  gr1_unique <- gr1[-queryHits(overlaps)]
  gr2_unique <- gr2[-unique(subjectHits(overlaps))]
  intersected <- gr1[unique(queryHits(overlaps))]  
  union_result <- c(intersected, gr2_unique)
  return(union_result)
}

# 合并所有结果
cattle_output <- unique(compute_union(sheep, cattle))
pig_output <- unique(compute_union(sheep, pig))
mouse_output <- unique(compute_union(sheep, mouse))
human_output <- unique(compute_union(sheep, human))

write.table(as.data.frame(cattle_output), file=paste0(Tisname,"/02",Tisname,"cattle_peaks_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(pig_output), file=paste0(Tisname,"/02",Tisname,"pig_peaks_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(mouse_output), file=paste0(Tisname,"/02",Tisname,"mouse_peaks_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(human_output), file=paste0(Tisname,"/02",Tisname,"human_peaks_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



####bar for cattle, pig, mouse human and common 4 species
compute_overlap <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
  intersected <- gr1[queryHits(overlaps)]  
  return(intersected)
}

# 合并所有结果
cs_output <- unique(compute_overlap(sheep, cattle_output))
ps_output <- unique(compute_overlap(sheep, pig_output))
ms_output <- unique(compute_overlap(sheep, mouse_output))
hs_output <- unique(compute_overlap(sheep, human_output))

tmp0=findOverlaps(sheep,cs_output,type="equal")
SC=sheep[unique(queryHits(tmp0))]

tmp1=findOverlaps(SC,ps_output,type="equal")
CP=SC[unique(queryHits(tmp1))]

tmp2=findOverlaps(CP,ms_output,type="equal")
CPM=CP[unique(queryHits(tmp2))]

tmp3=findOverlaps(CPM,hs_output,type="equal")
CPMH=CPM[unique(queryHits(tmp3))]

com5=CPMH

any1=unique(c(cs_output,ps_output,ms_output,hs_output))


write.table(as.data.frame(cs_output), file=paste0(Tisname,"/03",Tisname,"cattle_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(ps_output), file=paste0(Tisname,"/03",Tisname,"pig_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(ms_output), file=paste0(Tisname,"/03",Tisname,"mouse_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
write.table(as.data.frame(hs_output), file=paste0(Tisname,"/03",Tisname,"human_Enh_sheep.bed1"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)



# Calculate range counts for each dataset
dataset_names <- c("sheep", "cattle", "pig", "mouse", "human", "any1", "com5")
range_counts <- c(length(sheep), length(cs_output), length(ps_output), length(ms_output), length(hs_output), length(any1), length(com5))

# Plot bar chart
pdf(paste0(Tisname,"/04",Tisname,"_barplot_5species.pdf"))
bar_positions=barplot(range_counts, names.arg = dataset_names, col = "orange",
        main = "Number of enhancers identified in Each Dataset", xlab = "Dataset", ylab = "Enhancer count", ylim=c(0,max(range_counts)+2000))
text(x = bar_positions, y = range_counts, label = range_counts, pos = 3, cex = 0.8, col = "blue")
dev.off()



###upset
granges_list <- list(
    sheep = sheep,
    cattle = cattle_output,
    pig = pig_output,
    mouse = mouse_output,
    human = human_output
)

compute_gr2 <- function(gr1, gr2) {
  overlaps <- findOverlaps(gr1, gr2)
#  gr1_unique <- gr1[-queryHits(overlaps)]
  gr2_unique <- gr2[-unique(subjectHits(overlaps))]
  return(gr2_unique)
}

cattle_gr2 <- unique(compute_gr2(sheep, cattle))
pig_gr2 <- unique(compute_gr2(sheep, pig))
mouse_gr2 <- unique(compute_gr2(sheep, mouse))
human_gr2 <- unique(compute_gr2(sheep, human))

cpmh_gr2=reduce(c(cattle_gr2,pig_gr2,mouse_gr2,human_gr2))

combined_gr <- c(granges_list[[1]],cpmh_gr2)

# 将每个物种的GRanges转换为二元列
overlap_matrix <- sapply(granges_list, function(gr) {
    as.integer(countOverlaps(combined_gr, gr) > 0)
})

# 给列加上物种名称
colnames(overlap_matrix) <- names(granges_list)

library(GenomicRanges)
library(UpSetR)
# 将二元矩阵转换为 data.frame
upset_data <- as.data.frame(overlap_matrix)

# 生成并保存 UpSet 图
pdf(paste0(Tisname,"/05",Tisname,"_species_upset_plot.pdf") ) # 保存为 PDF 文件
upset(upset_data, sets = colnames(overlap_matrix), 
      main.bar.color = "orange", sets.bar.color = "orange", 
      order.by = "freq", keep.order = TRUE,
      sets.x.label = "Number of Enhancers")
dev.off()


#####Sheep unique
us_tmp=findOverlaps(sheep,any1)
sheep_uniq=unique(sheep[-queryHits(us_tmp)])

write.table(as.data.frame(sheep_uniq), file=paste0(Tisname,"/06",Tisname,"_sheep_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


####ruminant 
csp_tmp=findOverlaps(cs_output,ps_output,type="equal")
rum1=unique(cs_output[-queryHits(csp_tmp)])
cspm_tmp=findOverlaps(rum1,ms_output,type="equal")
rum2=unique(rum1[-queryHits(cspm_tmp)])
cspmh_tmp=findOverlaps(rum2,hs_output,type="equal")
rum3=unique(rum2[-queryHits(cspmh_tmp)])

ruminant_uniq <- rum3
write.table(as.data.frame(ruminant_uniq), file=paste0(Tisname,"/07",Tisname,"_sheep_ruminant_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


####sheep_common5
common_uniq <- com5
write.table(as.data.frame(common_uniq), file=paste0(Tisname,"/08",Tisname,"_sheep_common5_uniq_enhancer.bed"), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)


