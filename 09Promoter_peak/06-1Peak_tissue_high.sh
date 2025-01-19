#!/bin/bash

###promoter
workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/06Peak_tissue_high/01All
cd $workdir/08Promoter_peak/06Peak_tissue_high/01All

while read tissue
do
echo $tissue
Pro_p1=~/FAANG/analysis/08Promoter_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00H3K4me3.Propeak
Pro_p2=~/FAANG/analysis/08Promoter_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00H3K27ac.Propeak
Pro_p3=~/FAANG/analysis/01ChIP-seq/04Filter/N_H3K4me3/${tissue}_H3K4me3.peak
Pro_p4=~/FAANG/analysis/08Promoter_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00CAGE.Propeak

ls $Pro_p1
ls $Pro_p2
ls $Pro_p3
#ls $Pro_p4
module load R/4.2.3
mkdir $tissue

Rscript ~/FAANG/scripts/08Promoter_peak/06-2Peak_tissue_high.r $Pro_p1 $Pro_p2 $Pro_p3 $Pro_p4 $tissue

done < $workdir/08Promoter_peak/tissues_24.txt

