#!/bin/bash

##Promoter_final
workdir=~/FAANG/analysis
mkdir -p $workdir/08Promoter_peak/10Promoter_Tissue_Final
cd $workdir/08Promoter_peak/10Promoter_Tissue_Final

while read tissue
do
echo $tissue
Pro_T=$workdir/08Promoter_peak/09Promoter_Final/01final_promoter_gene_TSS.txt
Pro_H3K4me3=$workdir/01ChIP-seq/04Filter/N_H3K4me3/${tissue}_H3K4me3.peak
Pro_H3K27ac=$workdir/01ChIP-seq/04Filter/N_H3K27ac/${tissue}_H3K27ac.peak
Pro_CAGE=$workdir/02CAGE/${tissue}.CAGE.TSS.bed

ls $Pro_T
ls $Pro_H3K4me3
ls $Pro_H3K27ac
#ls $Pro_p4
module load R/4.2.3
mkdir $tissue

Rscript ~/FAANG/scripts/08Promoter_peak/10_1Promoter_tissue.r $Pro_T $Pro_H3K4me3 $Pro_H3K27ac $Pro_CAGE $tissue

done < $workdir/08Promoter_peak/tissues_24.txt




awk 'NR>1 {print $1"\t"$9"\t"$10"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7"\t"$8}' ../09Promoter_Final/01final_promoter_gene_TSS.txt > 11Gene_promoter.info
awk 'NR>1 {print $0}' 06All_promoter_across_tissues.txt >12Promoter_24tissues.txt
bedtools intersect -loj -a 12Promoter_24tissues.txt -b 11Gene_promoter.info |awk '$31!="."&&$2==$32&&$3==$33 {print $0}' >13Total_promoter25975_gene.txt

mkdir 14seq_new_prom
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5"\t"$6}' ../09Promoter_Final/01final_promoter_gene_TSS.txt |grep "STRG" > 14seq_new_prom/01New_promoter_gene.bed

bedtools getfasta -fi /mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCA_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -bed 01New_promoter_gene.bed > 02New_promoter_gene.seq
