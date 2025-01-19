#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/11Tissue_specific/03Homer
cd $workdir/11Tissue_specific/03Homer

module load R/4.2.3
#/mnt/ceph/sxie/20240401FAANG/scripts/11Tissue_specific/02-4Tissue_specific_homer.r


module load python
export PATH=$PATH:/mnt/ceph/bmurdoch/Homer/bin/

for i in $(ls */*_enhancer_position_gene_2.txt)
do
echo $i
#samplen=$(dirname $i _enhancer_position_gene.txt)
tissue=$(dirname $i)
sort -k1,1 -k2,2n $i >$tissue/02${tissue}_enhancer_sort_2.bed

outdir=$workdir/11Tissue_specific/03Homer/${tissue}/03homer_2
mkdir -p $outdir
peak_file=$workdir/11Tissue_specific/03Homer/$tissue/02${tissue}_enhancer_sort_2.bed
genome_sequence=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna
findMotifsGenome.pl ${peak_file} ${genome_sequence} $outdir -size 200 -p 40 -mknown /mnt/ceph/bmurdoch/Homer/data/knownTFs/vertebrates/all.motifs

done


###significant TF
workdir=~/FAANG/analysis
mkdir -p $workdir/11Tissue_specific/03Homer
cd $workdir/11Tissue_specific/03Homer

for i in $(ls */*_enhancer_position_gene.txt)
do

tissue=$(dirname $i)
echo $tissue
awk -v samp=$tissue 'BEGIN {print "Tissue\t"samp}
                     {if($3<=0.001) {print $1";"$2"\t"$3}}' $tissue/03homer/knownResults.txt |sort |uniq > $tissue/04${tissue}_noname_0.001.txt
done


cat */*_noname_0.001.txt |awk '{print $1}' |sort |uniq >Total_TF_0.txt

k=-1
for i in $(ls */*_noname_0.001.txt)
do
k=$[k+1]
m=$[k+1]
mk=$[m+1]
echo $m
sort -k1,1 -k2,2n Total_TF_${k}.txt > sort_Total_TF_${k}.txt
sort -k1,1 -k2,2n $i > ${i}_sort

join -a1 -a2 sort_Total_TF_${k}.txt ${i}_sort  \
|awk -v kk=$mk '{if($kk<0) 
                {tmp=$1; for(i=2;i<kk;i++) {tmp=tmp"\t"$i} ; print tmp"\t1" }
               if($kk>=0) {print $0} }' > Total_TF_${m}.txt
wc -l Total_TF_${m}.txt
done

awk '{tmp=$1; for (i=2;i<=NF;i++) {tmp=tmp"\t"$i}; print NF"\t"tmp}' Total_TF_24.txt >final_raw_p.txt 



module load R/4.2.3



