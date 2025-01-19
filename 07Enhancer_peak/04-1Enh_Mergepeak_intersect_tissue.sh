#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/07Enhancer_peak/04Mergep_intersect_tissue
cd $workdir/07Enhancer_peak/04Mergep_intersect_tissue

while read tissue
do
echo $tissue
Enh_T=$workdir/07Enhancer_peak/03Merge_tissue/03All27_tissues_01raw_tis2.txt
Enh_H3K4me1=$workdir/01ChIP-seq/04Filter/B_H3K4me1/${tissue}_H3K4me1_peaks.broadPeak.peak
Enh_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak
Enh_CAGE=$workdir/02CAGE/${tissue}.CAGE.bed

ls $Enh_T
ls $Enh_H3K4me1
ls $Enh_H3K27ac
ls $Enh_CAGE

module load R/4.2.3
mkdir $tissue

Rscript ~/FAANG/scripts/07Enhancer_peak/04-2Enh_Mergepeak_intersect_tissue.r $Enh_T $Enh_H3K4me1 $Enh_H3K27ac $Enh_CAGE $tissue

done < $workdir/08Promoter_peak/tissues_24.txt

