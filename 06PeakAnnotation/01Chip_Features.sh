#!/bin/bash

workdir=~/FAANG/analysis
mkdir -p $workdir/06PeakAnno
cd $workdir/06PeakAnno

mkdir H3K4me1 H3K4me3 H3K27ac H3K27me3 CAGE ATAC

gff=/mnt/ceph/bmurdoch/ARS-UI_Ramb_v2.0/GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic_chrom.gff

module load R/4.2.3
R

#Rscript 