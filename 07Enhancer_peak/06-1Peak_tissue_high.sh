#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/07Enhancer_peak/06Peak_tissue_high/01All
cd $workdir/07Enhancer_peak/06Peak_tissue_high/01All

while read tissue
do
echo $tissue
Pro_p1=~/FAANG/analysis/07Enhancer_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00H3K4me1.Propeak
Pro_p2=~/FAANG/analysis/07Enhancer_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00H3K27ac.Propeak
Pro_p3=~/FAANG/analysis/07Enhancer_peak/05H3K4me1_peak_P3/01SingleP/${tissue}.H3K4me1.P3.peak
Pro_p4=~/FAANG/analysis/07Enhancer_peak/04Mergep_intersect_tissue/${tissue}/${tissue}_00CAGE.Propeak


module load R/4.2.3
mkdir $tissue

Rscript ~/FAANG/scripts/07Enhancer_peak/06-2Peak_tissue_high.r $Pro_p1 $Pro_p2 $Pro_p3 $Pro_p4 $tissue

done < $workdir/08Promoter_peak/tissues_24.txt


