#!/bin/bash


workdir=~/FAANG/analysis
mkdir -p $workdir/07Enhancer_peak/02Peak_o2
cd $workdir/07Enhancer_peak/02Peak_o2

###Enhancer from CAGE, ATAC, H3K4me1 and H3K27ac
while read tissue
do
echo $tissue
Enh_C=$workdir/02CAGE/${tissue}.CAGE.bed
Enh_ATAC=$workdir/03ATAC/03Filter/${tissue}.ATAC.peak
Enh_H3K4me1=$workdir/01ChIP-seq/04Filter/B_H3K4me1/${tissue}_H3K4me1_peaks.broadPeak.peak
Enh_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak

ls $Enh_C
ls $Enh_ATAC
ls $Enh_H3K4me1
ls $Enh_H3K27ac

module load R/4.2.3
mkdir $tissue

Rscript ~/FAANG/scripts/07Enhancer_peak/02-1Enhancer_High_true_o2.r $Enh_C $Enh_ATAC $Enh_H3K4me1 $Enh_H3K27ac $tissue 20

done < $workdir/07Enhancer_peak/tissues20.txt


