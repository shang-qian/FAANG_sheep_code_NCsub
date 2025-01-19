#!/bin/bash

###
module load python
export PATH=$PATH:/mnt/ceph/bmurdoch/Homer/bin/

workdir=~/FAANG/analysis/12Comparison_sheep/04GO
cd $workdir


for i in $(ls *.csv)
do
echo $i
samplename=$(basename $i ".csv")
echo $samplename
awk -F "," -v samp=$samplename '$3<=0.05 {print samp"\t"$1"\t"$3}' $i >> 02GO.total0.05.txt 
done







