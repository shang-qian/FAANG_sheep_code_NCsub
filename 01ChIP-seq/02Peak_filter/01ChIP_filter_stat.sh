#!/bin/bash

#1.Chip-seq 24 tissues
workdir=~/FAANG/analysis/01ChIP-seq/03Filter
mkdir -p $workdir/N_H3K27ac $workdir/N_H3K4me3 $workdir/B_H3K27me3 $workdir/B_H3K4me1

###H3K27ac
cd $workdir/N_H3K27ac
for i in $(ls ../../02nfcore_5e3/01Narrow/NF_*/bowtie2/mergedLibrary/macs2/narrowPeak/*H3K27ac*.narrowPeak)
do
echo $i
sample=$(basename $i _peaks.narrowPeak)
echo $sample
awk '/^CM/ {print $0}' $i > ${sample}.peak
done

##H3K4me3
cd $workdir/N_H3K4me3
for i in $(ls ../../02nfcore_5e3/01Narrow/NF_*/bowtie2/mergedLibrary/macs2/narrowPeak/*H3K4me3*.narrowPeak)
do
echo $i
sample=$(basename $i _peaks.narrowPeak)
echo $sample
awk '/^CM/ {print $0}' $i > ${sample}.peak
done

##H3K4me1
cd $workdir/B_H3K4me1
for i in $(ls ../../02nfcore_5e3/02Broad/NF_*/bowtie2/mergedLibrary/macs2/broadPeak/*H3K4me1*.broadPeak)
do
echo $i
sample=$(basename $i _peaks.broadPeak)
echo $sample
awk '/^CM/ {print $0}' $i > ${sample}.peak
done

##H3K27me3
cd $workdir/B_H3K27me3
for i in $(ls ../../02nfcore_5e3/02Broad/NF_*/bowtie2/mergedLibrary/macs2/broadPeak/*H3K27me3*.broadPeak)
do
echo $i
sample=$(basename $i _peaks.broadPeak)
echo $sample
awk '/^CM/ {print $0}' $i > ${sample}.peak
done
