#!/bin/bash

##Enhancer_final
workdir=~/FAANG/analysis
mkdir -p $workdir/07Enhancer_peak/9Enhancer_Tissue_Final
cd $workdir/07Enhancer_peak/9Enhancer_Tissue_Final

module load R/4.2.3

while read tissue
do
echo $tissue
Pro_T=$workdir/07Enhancer_peak/08Enh/02final_enhancer_total2000.txt
Pro_H3K4me1=$workdir/01ChIP-seq/04Filter/B_H3K4me1/${tissue}_H3K4me1_peaks.broadPeak.peak
Pro_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak
Pro_CAGE=$workdir/02CAGE//${tissue}.CAGE.bed

mkdir $tissue
Rscript ~/FAANG/scripts/07Enhancer_peak/9_1Enh_tissue.r $Pro_T $Pro_H3K4me1 $Pro_H3K27ac $Pro_CAGE $tissue
done < $workdir/08Promoter_peak/tissues_24.txt

