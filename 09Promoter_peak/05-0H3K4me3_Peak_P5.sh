#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/05H3K4me3_peak_P3
cd $workdir/08Promoter_peak/05H3K4me3_peak_P3

mkdir 01SingleP
while read tissue
do
echo $tissue
awk '$7>=3&&$8>=3 {print $0}' ~/FAANG/analysis/01ChIP-seq/04Filter/N_H3K4me3/${tissue}_H3K4me3.peak > 01SingleP/${tissue}.H3K4me3.P3.peak
done < $workdir/08Promoter_peak/tissues_24.txt


