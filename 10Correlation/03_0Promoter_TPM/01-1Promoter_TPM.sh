#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/09Promoter_TPM/03pro_TPM
cd $workdir/09Promoter_TPM/03pro_TPM

module load R/4.2.3
while read tissue
do
echo $tissue
Promoter=~/FAANG/analysis/08Promoter_peak/10Promoter_Tissue_Final/$tissue/${tissue}_01All_promoter_peaks.txt

TPM=~/FAANG/analysis/08RNA_TPM/02Promoter/02TPM/${tissue}_TPM_gene.txt
Rscript 01-2Enhancer_tissue_TPM.r $Promoter $TPM $tissue

done < $workdir/tissue_24.txt

