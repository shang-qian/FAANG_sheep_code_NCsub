#!/bin/bash

workdir=~/FAANG/analysis/05ChromHMM
#1Model 
cd $workdir/03Model/01Output
###M
for tissue in $(ls -d $workdir/03Model/01Output/[A-Z]*)
do
echo $tissue
TN=$(basename $tissue)
echo $TN
mkdir -p $TN/compare_model $workdir/03Model/02Evaluate/$TN/
cp $TN/[7-9]*/emissions_*.txt $TN/compare_model
cp $TN/[1][0-8]*/emissions_*.txt $TN/compare_model
java -mx4000M -jar /mnt/ceph/bmurdoch/ChromHMM/ChromHMM.jar CompareModels $workdir/03Model/01Output/$TN/18/emissions_18.txt $workdir/03Model/01Output/$TN/compare_model/ $workdir/03Model/02Evaluate/$TN/comp_18
done


cd $workdir/03Model/02Evaluate
module load R/4.2.3
R
tissue_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
all_data <- list()
for (tissue_dir in tissue_dirs) {
  file_path <- file.path(tissue_dir, "comp_18.txt")  
  # check the files
  if (file.exists(file_path)) {
    data <- as.matrix(read.table(file_path, header=TRUE, sep="\t", row.names=1))    
    all_data[[tissue_dir]] <- data
  } else {
    message(paste("File comp_18.txt not found in", tissue_dir))
  }
}


#calculate median
calculate_column_means <- function(matrix) {
  apply(matrix, 2, median)
}

column_means_list <- lapply(all_data, calculate_column_means)

column_means_matrix <- do.call(rbind, column_means_list)
rownames(column_means_matrix) <- sub("./", "", rownames(column_means_matrix))
colnames(column_means_matrix) <- paste0(7:18)

print(column_means_matrix)

library(pheatmap)
pdf("heatmap_column_means.pdf")
# Plot heatmap
pheatmap(column_means_matrix, cluster_rows = FALSE, cluster_cols = FALSE, main = "Column Means Heatmap")
dev.off()

