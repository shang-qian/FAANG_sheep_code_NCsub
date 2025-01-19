#!/bin/bash

#CAGE
workdir=~/FAANG/analysis
mkdir -p $workdir/10interaction/01CAGE/TSS_enh
cd $workdir/10interaction/01CAGE/TSS_enh

while read tissue
do
echo $tissue

awk -F"," '{print $15"\t"$16"\t"$17"\t"$2"\t"$3"\t"$4}' $workdir/02CAGE/24tissues/${tissue}_links_results.csv > 01$tissue.CAGE_enh_pro.bed
awk '/^CM/ {print $0}' 01$tissue.CAGE_enh_pro.bed |sort -k1,1 -k2,2n |uniq >02$tissue.CAGE_enh_pro_input.bed 

done < $workdir/01tissue_24.txt


