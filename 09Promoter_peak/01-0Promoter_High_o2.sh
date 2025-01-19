#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/02Peak_o2
cd $workdir/08Promoter_peak/02Peak_o2

module load R/4.2.3

###Promoter from CAGE, ATAC, H3K4me3
while read tissue
do
echo $tissue
Pro_C=$workdir/02CAGE/${tissue}.CAGE.TSS.bed
Pro_ATAC=$workdir/03ATAC/03Filter/${tissue}.ATAC.peak
Pro_H3K4me3=$workdir/01ChIP-seq/04Filter/N_H3K4me3/${tissue}_H3K4me3.peak
Pro_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak

ls $Pro_C
ls $Pro_ATAC
ls $Pro_H3K4me3

mkdir $tissue

Rscript 02-1Promoter_High_o2.r $Pro_C $Pro_ATAC $Pro_H3K4me3 $Pro_H3K27ac $tissue 20

done < $workdir/07Enhancer_peak/tissues20.txt


