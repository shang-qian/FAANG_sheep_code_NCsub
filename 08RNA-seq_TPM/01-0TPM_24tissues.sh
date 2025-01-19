#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/08RNA_TPM/02Promoter/01StringT
cd $workdir/08RNA_TPM/02Promoter/01StringT

#1.stringtie 
module load stringtie
while read tissue
do
echo $tissue
mkdir $tissue
cd $tissue

bamfile=/mnt/ceph/bmurdoch/Kim/Benz2616_RNA-seq_analyses/02-Mapped/merged/Benz2616_${tissue}_RNA_PE.sorted.bam
gtf=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff
stringtie -p 40 -G $gtf -B -o ${tissue}.output.gtf -A $tissue.gene_abundence.tab -C $tissue.known.cov_refs.gtf $bamfile

cd ..
done < ~/FAANG/analysis/08RNA_TPM/tissue_24.txt

###extract gene TPM
workdir=~/FAANG/analysis
mkdir -p $workdir/08RNA_TPM/02Promoter/02TPM
cd $workdir/08RNA_TPM/02Promoter/02TPM

while read tissue
do
echo $tissue
TPMfile=$workdir/08RNA_TPM/02Promoter/01StringT/${tissue}/${tissue}.gene_abundence.tab
awk 'BEGIN{FS=OFS="\t"} NR==1 {print "Chromosome", "Start", "End", "Strand", "TPM", "GeneID"} NR>1 {GeneID=($1 == ".") ? $2 : $1; print $3, $5, $6, $4, $9, GeneID}' $TPMfile \
> ${tissue}_TPM_gene.txt
done < tissue_24.txt

