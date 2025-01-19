#!/bin/bash


workdir=~/FAANG/analysis/05ChromHMM
#1Model 
mkdir -p $workdir/03Model
cd $workdir/03Model
###M
for tissue in $(ls -d $workdir/02BinaryFile/02Output/*)
do
echo $tissue
TN=$(basename $tissue)
echo $TN
for i in {1..17}
do
echo $i
java -mx4000M -jar /mnt/ceph/bmurdoch/Kim/ChromHMM/ChromHMM.jar LearnModel -b 200 $workdir/02BinaryFile/02Output/$TN $workdir/03Model/01Output/$TN/$i $i rambouilletv2
done
done

