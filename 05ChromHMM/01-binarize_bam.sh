#!/bin/bash

workdir=~/FAANG/analysis/05ChromHMM

#1chrom size
mkdir -p $workdir/01GenomicBins
cd $workdir/01GenomicBins
cp /mnt/ceph/bmurdoch/ChromHMM/CHROMSIZES/rambouilletv2.txt RamV2.chromsize

#2 ln bam files
mkdir -p $workdir/02BinaryFile/01Input
cd $workdir/02BinaryFile/01Input

#ATAC data
for tissue in $(ls ~/FAANG/analysis/03ATAC/02nfcore/*_ATAC*/bowtie2/merged_library/*_REP1.mLb.clN.sorted.bam)
do
TN=$(basename $tissue _REP1.mLb.clN.sorted.bam)
echo $TN
mkdir $TN
ln -sf $tissue ${TN}/${TN}_ATAC_link.bam
done

#Chip-seq data H3K27ac
for mark in H3K27ac H3K4me3 Input
do
echo $mark
for tissue in $(ls ~/FAANG/analysis/01ChIP-seq/02nfcore_5e3/01Narrow/NF_*/bowtie2/mergedLibrary/*_${mark}.mLb.clN.sorted.bam)
do
TN=$(basename $tissue _${mark}.mLb.clN.sorted.bam)
echo $TN
mkdir $TN
ln -sf $tissue ${TN}/${TN}_${mark}_link.bam
done
done

for mark in H3K27me3 H3K4me1 Input
do
echo $mark
for tissue in $(ls ~/FAANG/analysis/01ChIP-seq/02nfcore_5e3/02Broad/NF_*/bowtie2/mergedLibrary/*_${mark}.mLb.clN.sorted.bam)
do
TN=$(basename $tissue _${mark}.mLb.clN.sorted.bam)
echo $TN
mkdir $TN
ln -sf $tissue ${TN}/${TN}_${mark}_link.bam
done


##3Generate the bam file table
cd $workdir/02BinaryFile/01Input
for tissue in $(ls )
do
echo $tissue
if [ $(ls $tissue | wc -l) -eq 7 ]; then
  mv $tissue 5/${tissue}
fi


##generate table
cd $workdir/02BinaryFile/01Input/5
for tissue in $(ls )
do
echo $tissue
echo $tissue |awk -v Tis=$tissue '{print Tis"\tH3K27ac\t"Tis"_H3K27ac_link.bam\t"Tis"_Input_n_link.bam\n"\
Tis"\tH3K4me3\t"Tis"_H3K4me3_link.bam\t"Tis"_Input_n_link.bam\n"\
Tis"\tH3K27me3\t"Tis"_H3K27me3_link.bam\t"Tis"_Input_link.bam\n"\
Tis"\tH3K4me1\t"Tis"_H3K4me1_link.bam\t"Tis"_Input_link.bam\n"\
Tis"\tATAC\t"Tis"_ATAC_link.bam\t"}' >$tissue/${tissue}_table.txt
done


###generate binaryfile
for tissue in $(ls -d *)
do
echo $tissue
java -mx4000M -jar /mnt/ceph/bmurdoch/ChromHMM/ChromHMM.jar BinarizeBam -b 200 $workdir/01GenomicBins/RamV2.chromsize $workdir/02BinaryFile/01Input/5/$tissue $workdir/02BinaryFile/01Input/5/$tissue/${tissue}_table.txt $workdir/02BinaryFile/02Output/$tissue
done

