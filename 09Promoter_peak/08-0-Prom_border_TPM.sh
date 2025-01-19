#!/bin/bash

##01TPM_new
workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/08Prom_border_TPM/01TPM
cd $workdir/08Promoter_peak/08Prom_border_TPM/01TPM


##02 Promoter_tissue
workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/08Prom_border_TPM/02Prom_tis
cd $workdir/08Promoter_peak/08Prom_border_TPM/02Prom_tis

while read tissue
do
echo $tissue
Pro_p1=~/FAANG/analysis/08Promoter_peak/08Prom_border_TPM/01TPM/07TSS_annotate_TPM_final.bed
Pro_p2=~/FAANG/analysis/08Promoter_peak/06Peak_tissue_high/01All/${tissue}/${tissue}_02final_peak_union.bed

ls $Pro_p1
ls $Pro_p2
#ls $Pro_p4
module load R/4.2.3
#mkdir $tissue
Rscript ~/FAANG/scripts/08Promoter_peak/08_2Prom_border.r $Pro_p1 $Pro_p2 $tissue
done < $workdir/08Promoter_peak/tissues_24.txt

