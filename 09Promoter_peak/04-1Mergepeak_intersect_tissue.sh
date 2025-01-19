#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/04Mergep_intersect_tissue
cd $workdir/08Promoter_peak/04Mergep_intersect_tissue

while read tissue
do
echo $tissue
pro_T=$workdir/08Promoter_peak/03Merge_tissue/03All20_tissues_01raw_2.txt
pro_H3K4me3=$workdir/01ChIP-seq/04Filter/N_H3K4me3/${tissue}_H3K4me3.peak
pro_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak
pro_CAGE=$workdir/02CAGE/${tissue}.CAGE.TSS.bed
ls $pro_T
ls $pro_H3K4me3
ls $pro_H3K27ac
ls $pro_CAGE

module load R/4.2.3
mkdir $tissue

Rscript 04-2Mergepeak_intersect_tissue.r $pro_T $pro_H3K4me3 $pro_H3K27ac $pro_CAGE $tissue

done < $workdir/08Promoter_peak/tissues_24.txt


