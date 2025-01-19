#!/bin/bash


workdir=~/FAANG/analysis
mkdir -p $workdir/07Enhancer_peak/05H3K4me1_peak_P3
cd $workdir/07Enhancer_peak/05H3K4me1_peak_P3

mkdir 01SingleP
while read tissue
do
echo $tissue
awk '$7>=3&&$8>=3 {print $0}' $workdir/01ChIP-seq/04Filter/B_H3K4me1/${tissue}_H3K4me1_peaks.broadPeak.peak > 01SingleP/${tissue}.H3K4me1.P3.1peak
done < $workdir/08Promoter_peak/tissues_24.txt




